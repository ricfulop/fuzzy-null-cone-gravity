"""Bounded weak-form validation for the prepared oblate curvature blocks.

This checker is intentionally local and staged.  It validates the prepared
first-derivative potential and lapse densities against an independently derived
flat-oblate Christoffel table, linearized Einstein tensor, linearized scalar
curvature, and Green fluxes in small L=2 sectors with m=0,1,2.

It does not assemble finite radial K,E,D matrices, choose a physical endpoint
domain, infer a spectrum, or call the Paper I strong stationary exporter.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import time

import sympy as sp

import kinematic_blocks as kinematics
import weak_forms as prepared


STAGE_SCOPE = (
    "independent Christoffels",
    "derivative-slot densities at sampled rational oblate points",
    "radial endpoint flux normalization",
    "potential Green identity evaluated exactly at representative rational points",
    "lapse Green identity evaluated exactly at representative rational points",
    "continuum pure-gauge zero controls",
    "Cartesian TT sign and one-half normalization",
    "nonzero and ablation controls",
)


VALIDATION_BATCHES = 8


def emit(stage, **payload):
    print(json.dumps({"stage": stage, **payload}, sort_keys=True), flush=True)


def emit_orx(payload):
    print(f"::orx:: {json.dumps(payload, sort_keys=True)}", flush=True)


def emit_substage(batch, m, name, t0):
    now = time.monotonic()
    emit("SUBSTAGE", batch=batch, m=m, name=name, seconds=now - t0)
    return now


def normalize_rational(value):
    return sp.cancel(value)


def chart_conjugate(value):
    return value.subs(sp.I, -sp.I)


def matrix_chart_conjugate(matrix):
    return matrix.applyfunc(chart_conjugate)


def derivative(value, coord, data):
    return field_derivative(value, coord, data, data["m"])


def field_derivative(value, coord, data, azimuthal_m):
    if coord == 0:
        return sp.diff(value, data["r"])
    if coord == 1:
        return sp.diff(value, data["x"])
    return sp.I * azimuthal_m * value


def background_derivative(value, coord, data):
    return field_derivative(value, coord, data, 0)


def exact_zero(expression, x):
    """Chart exact zero on -1<x<1."""
    expression = sp.sympify(expression).replace(
        lambda term: (term.is_Pow and term.exp.is_Rational
                      and term.exp.q == 2 and term.base.is_polynomial(x)),
        lambda term: sp.Pow(sp.expand(term.base), term.exp),
    )
    d = 1 - x**2
    for weight in (sp.S.One, sp.sqrt(d), 1 / sp.sqrt(d)):
        value = sp.factor(sp.cancel(sp.powsimp(expression * weight, force=False)))
        if value == 0 or value.is_zero is True:
            return True
    return False


def assert_zero(expression, x, label):
    if not exact_zero(expression, x):
        raise AssertionError(f"{label}: {expression}")


def assert_zero_matrix(matrix, x, label):
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            assert_zero(matrix[i, j], x, f"{label} entry {(i, j)}")


def independent_connection(data):
    metric = data["metric"]
    inverse = metric.inv()
    gamma = [[[sp.S.Zero for _ in range(3)] for _ in range(3)] for _ in range(3)]
    for a in range(3):
        for b in range(3):
            for c in range(3):
                gamma[a][b][c] = sp.simplify(sum(
                    inverse[a, d] * (
                        background_derivative(metric[d, c], b, data)
                        + background_derivative(metric[d, b], c, data)
                        - background_derivative(metric[b, c], d, data)
                    )
                    for d in range(3)
                ) / 2)
    return gamma


def tensor_jet(tensor, data, value_slot, derivative_slot, gamma):
    h = value_slot * tensor
    nabla = []
    for a in range(3):
        partial = tensor.diff(data["r"]) if a == 0 else (
            tensor.diff(data["x"]) if a == 1 else sp.I * data["m"] * tensor
        )
        first = value_slot * partial
        if a == 0:
            first += derivative_slot * tensor
        first -= sp.Matrix(3, 3, lambda i, j: sum(
            gamma[b][a][i] * h[b, j] + gamma[b][a][j] * h[i, b]
            for b in range(3)
        ))
        nabla.append(first)
    inverse = data["inverse_diagonal"]
    trace = sum(inverse[i] * h[i, i] for i in range(3))
    trace_gradient = [sum(inverse[i] * nabla[a][i, i] for i in range(3))
                      for a in range(3)]
    divergence_lower = [sum(inverse[i] * nabla[i][i, j] for i in range(3))
                        for j in range(3)]
    return {
        "h": h,
        "nabla": nabla,
        "trace": trace,
        "trace_gradient": trace_gradient,
        "divergence_lower": divergence_lower,
    }


def linearized_christoffel(jet, data):
    inverse = data["inverse_diagonal"]
    return sp.MutableDenseNDimArray([
        normalize_rational(inverse[i] * (
            jet["nabla"][j][i, k]
            + jet["nabla"][k][i, j]
            - jet["nabla"][i][j, k]
        ) / 2)
        for i in range(3) for j in range(3) for k in range(3)
    ], (3, 3, 3))


def covariant_derivative_delta_gamma(delta_gamma, first_index, data, gamma):
    return sp.MutableDenseNDimArray([
        normalize_rational(
            derivative(delta_gamma[a, b, c], first_index, data)
            + sum(gamma[a][first_index][d] * delta_gamma[d, b, c]
                  - gamma[d][first_index][b] * delta_gamma[a, d, c]
                  - gamma[d][first_index][c] * delta_gamma[a, b, d]
                  for d in range(3))
        )
        for a in range(3) for b in range(3) for c in range(3)
    ], (3, 3, 3))


def varied_ricci_lower(jet, data, gamma):
    delta_gamma = linearized_christoffel(jet, data)
    covariant = [
        covariant_derivative_delta_gamma(delta_gamma, i, data, gamma)
        for i in range(3)
    ]
    return sp.Matrix(3, 3, lambda i, j: normalize_rational(
        sum(covariant[a][a, i, j] - covariant[j][a, i, a]
            for a in range(3))
    ))


def scalar_and_einstein_from_ricci(ricci, data):
    scalar = normalize_rational(
        sum(data["inverse_diagonal"][i] * ricci[i, i] for i in range(3))
    )
    einstein = sp.Matrix(3, 3, lambda i, j: normalize_rational(
        ricci[i, j] - data["metric"][i, j] * scalar / 2
    ))
    return scalar, einstein


def linearized_scalar_curvature(jet, data, gamma):
    ricci = varied_ricci_lower(jet, data, gamma)
    scalar, _ = scalar_and_einstein_from_ricci(ricci, data)
    return scalar


def linearized_einstein_lower(jet, data, gamma):
    ricci = varied_ricci_lower(jet, data, gamma)
    _, einstein = scalar_and_einstein_from_ricci(ricci, data)
    return einstein


def lower_two_covariant_derivative(jet, data, gamma):
    second = []
    for a in range(3):
        row = []
        for b in range(3):
            row.append(
                sp.Matrix(3, 3, lambda i, j:
                          derivative(jet["nabla"][b][i, j], a, data)
                          - sum(gamma[c][a][b] * jet["nabla"][c][i, j]
                                for c in range(3))
                          - sum(gamma[c][a][i] * jet["nabla"][b][c, j]
                                + gamma[c][a][j] * jet["nabla"][b][i, c]
                                for c in range(3)))
            )
        second.append(row)
    return second


def scalar_gradient(value, data):
    return [derivative(value, i, data) for i in range(3)]


def scalar_hessian(value, data, gamma):
    grad = scalar_gradient(value, data)
    return [
        [
            derivative(grad[j], i, data)
            - sum(gamma[c][i][j] * grad[c] for c in range(3))
            for j in range(3)
        ]
        for i in range(3)
    ]


def covector_divergence(covector, data, gamma):
    inverse = data["inverse_diagonal"]
    return sum(
        inverse[i] * (derivative(covector[i], i, data)
                      - sum(gamma[c][i][i] * covector[c] for c in range(3)))
        for i in range(3)
    )


def formula_scalar_curvature(jet, data, gamma):
    return covector_divergence(jet["divergence_lower"], data, gamma) - sum(
        data["inverse_diagonal"][i] * scalar_hessian(jet["trace"], data, gamma)[i][i]
        for i in range(3)
    )


def formula_einstein_lower(jet, data, gamma):
    inverse = data["inverse_diagonal"]
    second = lower_two_covariant_derivative(jet, data, gamma)
    hessian_trace = scalar_hessian(jet["trace"], data, gamma)
    laplace_trace = sum(inverse[a] * hessian_trace[a][a] for a in range(3))
    div_div = covector_divergence(jet["divergence_lower"], data, gamma)
    result = sp.zeros(3)
    for i in range(3):
        for j in range(3):
            laplace_h = sum(inverse[a] * second[a][a][i, j]
                            for a in range(3))
            grad_div_ij = (
                derivative(jet["divergence_lower"][j], i, data)
                - sum(gamma[c][i][j] * jet["divergence_lower"][c]
                      for c in range(3))
            )
            grad_div_ji = (
                derivative(jet["divergence_lower"][i], j, data)
                - sum(gamma[c][j][i] * jet["divergence_lower"][c]
                      for c in range(3))
            )
            result[i, j] = (
                -laplace_h
                - hessian_trace[i][j]
                + grad_div_ij
                + grad_div_ji
                - data["metric"][i, j] * div_div
                + data["metric"][i, j] * laplace_trace
            ) / 2
    return result


def strong_potential_density(left, right, data, gamma):
    einstein = linearized_einstein_lower(right, data, gamma)
    inverse = data["inverse_diagonal"]
    return data["volume_density"] * sum(
        inverse[i] * inverse[j] * chart_conjugate(left["h"][i, j]) * einstein[i, j]
        for i in range(3) for j in range(3)
    )


def strong_lapse_density(y, f_slot, right, data, gamma):
    value = f_slot * y
    return data["volume_density"] * chart_conjugate(value) * (
        linearized_scalar_curvature(right, data, gamma) / 2
    )


def independent_lapse_density(y, f_slot, f_prime_slot, right, data):
    gradient_f = [
        f_prime_slot * y,
        f_slot * sp.diff(y, data["x"]),
        sp.I * data["m"] * f_slot * y,
    ]
    return -data["volume_density"] * sum(
        data["inverse_diagonal"][j] * chart_conjugate(gradient_f[j])
        * (right["divergence_lower"][j] - right["trace_gradient"][j])
        for j in range(3)
    ) / 2


def potential_flux_vector(left, right, data):
    inverse = data["inverse_diagonal"]
    h = matrix_chart_conjugate(left["h"])
    trace = chart_conjugate(left["trace"])
    flux = []
    for a in range(3):
        flux_a = -inverse[a] * sum(
            inverse[i] * inverse[j] * h[i, j] * right["nabla"][a][i, j]
            for i in range(3) for j in range(3)
        ) / 2
        flux_a += sum(
            inverse[a] * inverse[j] * h[a, j]
            * (right["divergence_lower"][j] - right["trace_gradient"][j] / 2)
            for j in range(3)
        )
        flux_a += trace * inverse[a] * (
            -right["divergence_lower"][a] + right["trace_gradient"][a]
        ) / 2
        flux.append(flux_a)
    return flux


def lapse_flux_vector(y, f_slot, right, data):
    inverse = data["inverse_diagonal"]
    return [
        chart_conjugate(f_slot * y) * inverse[a]
        * (right["divergence_lower"][a] - right["trace_gradient"][a]) / 2
        for a in range(3)
    ]


def volume_divergence(vector, data, azimuthal_m=0):
    density = data["volume_density"]
    return sum(field_derivative(density * vector[a], a, data, azimuthal_m)
               for a in range(3))


def direct_lie(metric, vector, data):
    return sp.Matrix(3, 3, lambda i, j:
                     sum(vector[a] * background_derivative(metric[i, j], a, data)
                         for a in range(3))
                     + sum(metric[k, j] * derivative(vector[k], i, data)
                           + metric[i, k] * derivative(vector[k], j, data)
                           for k in range(3)))


def sample_labels(data):
    """Small representative set: scalar, vector, tensor when present."""
    labels = data["spatial_labels"]
    wanted = {"lambda", "X", "V_polar", "V_axial", "C_polar", "C_axial"}
    chosen = []
    seen = set()
    for index, label in enumerate(labels):
        if label[0] in wanted and label[0] not in seen:
            chosen.append(index)
            seen.add(label[0])
    return chosen


def sample_point(data, m):
    return {
        data["r"]: sp.Rational(3, 2) + sp.Rational(m, 7),
        data["x"]: sp.Rational(1, 5),
        data["epsilon"]: sp.Rational(1, 4),
    }


def polynomial_profiles(data, label):
    r = data["r"]
    profiles = {
        "left": r**2 - 2 * r + sp.Rational(5, 4),
        "right_a": r**3 - r + sp.Rational(3, 2),
        "right_b": r**2 + r / 3 + sp.Rational(2, 5),
    }
    profile = profiles[label]
    return profile, sp.diff(profile, r)


def write_checkpoint(output_dir, stage, payload, started):
    if output_dir is None:
        return
    record = {
        "stage": stage,
        "status": "passed",
        "seconds": time.monotonic() - started,
        "payload": payload,
    }
    with (output_dir / "weak-form-progress.jsonl").open("a") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def emit_validation_progress(completed_batches):
    emit_orx({
        "phase": "validation",
        "step": completed_batches,
        "of": VALIDATION_BATCHES,
        "unit": "batches",
        "m": {"completed_batches": completed_batches},
    })


def verify_derivative_slots(output_dir=None, started=None, completed_batches=0):
    counts = {}
    for m in range(3):
        emit("STAGE", name="weak_density_sector", m=m)
        data = kinematics.sector(2, m)
        gamma = independent_connection(data)
        prepared_gamma = prepared.connection(data)
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    assert_zero(gamma[a][b][c] - prepared_gamma[a][b][c],
                                data["x"], f"Christoffel m={m} {(a, b, c)}")

        tensors = data["spatial_tangents"]
        for index in sample_labels(data):
            tensor = tensors[index]
            for value_slot, derivative_slot in ((1, 0), (0, 1), (2, -3)):
                expected = tensor_jet(tensor, data, value_slot, derivative_slot, gamma)
                actual = prepared.tensor_jet(
                    tensor, data, value_slot, derivative_slot, prepared_gamma
                )
                assert_zero_matrix(actual["h"] - expected["h"], data["x"],
                                   f"jet h m={m} index={index}")
                for a in range(3):
                    assert_zero_matrix(actual["nabla"][a] - expected["nabla"][a],
                                       data["x"],
                                       f"jet nabla m={m} index={index} a={a}")

        point = sample_point(data, m)
        left_profile = polynomial_profiles(data, "left")
        right_profile = polynomial_profiles(data, "right_a")
        left = [prepared.tensor_jet(tensors[i], data, *left_profile, prepared_gamma)
                for i in sample_labels(data)]
        right = [prepared.tensor_jet(tensors[i], data, *right_profile, prepared_gamma)
                 for i in sample_labels(data)]
        right_cache = []
        for j, right_jet in enumerate(right):
            cached = {
                "einstein": linearized_einstein_lower(right_jet, data, gamma),
                "einstein_formula": formula_einstein_lower(right_jet, data, gamma),
                "scalar": linearized_scalar_curvature(right_jet, data, gamma),
                "scalar_formula": formula_scalar_curvature(right_jet, data, gamma),
            }
            assert_zero_matrix(
                (cached["einstein"] - cached["einstein_formula"]).subs(point),
                data["x"], f"Ricci-variation Einstein m={m} {j}",
            )
            assert_zero(
                (cached["scalar"] - cached["scalar_formula"]).subs(point),
                data["x"], f"Ricci-variation scalar m={m} {j}",
            )
            right_cache.append(cached)
        for i, left_jet in enumerate(left):
            for j, right_jet in enumerate(right):
                expected = data["volume_density"] * sum(
                    data["inverse_diagonal"][a] * data["inverse_diagonal"][b]
                    * chart_conjugate(left_jet["h"][a, b])
                    * right_cache[j]["einstein"][a, b]
                    for a in range(3) for b in range(3)
                )
                weak = prepared.potential_density(left_jet, right_jet, data)
                flux = volume_divergence(potential_flux_vector(left_jet, right_jet, data), data, 0)
                assert_zero((expected - weak - flux).subs(point), data["x"],
                            f"potential Green density m={m} {(i, j)}")
                assert_zero(
                    (prepared.potential_boundary_flux(left_jet, right_jet, data)
                     - data["volume_density"]
                     * potential_flux_vector(left_jet, right_jet, data)[0]).subs(point),
                    data["x"], f"radial potential flux m={m} {(i, j)}"
                )

        for ell_index, (_, ell, _) in enumerate(data["lapse_labels"]):
            y = kinematics.harmonic_amplitudes(ell, data["m"], data["x"])[0]
            for j, right_jet in enumerate(right):
                for f_slot, f_prime_slot in ((5, 0), (0, -7)):
                    expected_lapse = independent_lapse_density(
                        y, f_slot, f_prime_slot, right_jet, data
                    )
                    actual_lapse = prepared.lapse_density(
                        y, f_slot, f_prime_slot, right_jet, data
                    )
                    assert_zero(
                        (actual_lapse - expected_lapse).subs(point),
                        data["x"],
                        f"lapse slot density m={m} {(ell_index, j)}"
                        f" {(f_slot, f_prime_slot)}",
                    )
                f, f_prime = polynomial_profiles(data, "right_b")
                strong = (
                    data["volume_density"] * chart_conjugate(f * y)
                    * right_cache[j]["scalar"] / 2
                )
                weak = prepared.lapse_density(y, f, f_prime, right_jet, data)
                flux = volume_divergence(lapse_flux_vector(y, f, right_jet, data), data, 0)
                assert_zero((strong - weak - flux).subs(point), data["x"],
                            f"lapse Green density m={m} {(ell_index, j)}")
                assert_zero(
                    (prepared.lapse_boundary_flux(y, f, right_jet, data)
                     - data["volume_density"]
                     * lapse_flux_vector(y, f, right_jet, data)[0]).subs(point),
                    data["x"], f"radial lapse flux m={m} {(ell_index, j)}"
                )

        counts[f"m={m}"] = {
            "sampled_spatial_tensors": len(left),
            "potential_pairs": len(left) * len(right),
            "lapse_pairs": len(data["lapse_labels"]) * len(right),
            "selected_spatial_labels": [data["spatial_labels"][i] for i in sample_labels(data)],
            "rational_point": {str(key): str(value) for key, value in point.items()},
            "left_profile": str(left_profile[0]),
            "right_profile": str(right_profile[0]),
            "lapse_profile": str(polynomial_profiles(data, "right_b")[0]),
            "green_evidence": "exact at the recorded point, not a global symbolic identity certificate",
        }
        print(f"PASS weak density and Green identities L=2 m={m}", flush=True)
        write_checkpoint(output_dir, f"derivative_slots_m{m}", counts[f"m={m}"], started)
        completed_batches += 1
        emit_validation_progress(completed_batches)
    return counts, completed_batches


def verify_gauge_controls(output_dir=None, started=None, completed_batches=0):
    results = []
    for m in range(3):
        emit("STAGE", name="pure_gauge_control", m=m)
        batch_started = time.monotonic()
        data = kinematics.sector(2, m)
        gamma = independent_connection(data)
        batch_started = emit_substage(
            "pure_gauge_control", m, "connection", batch_started
        )
        r, x = data["r"], data["x"]
        y = kinematics.harmonic_amplitudes(max(m, 1), m, x)[0]
        vector = sp.Matrix([r * (2 - r) * y, x * y, sp.I * y])
        gauge = direct_lie(data["metric"], vector, data)
        batch_started = emit_substage(
            "pure_gauge_control", m, "coordinate_lie", batch_started
        )
        jet = tensor_jet(gauge, data, 1, 0, gamma)
        batch_started = emit_substage(
            "pure_gauge_control", m, "tensor_jet", batch_started
        )
        ricci = varied_ricci_lower(jet, data, gamma)
        batch_started = emit_substage(
            "pure_gauge_control", m, "varied_ricci", batch_started
        )
        scalar, einstein = scalar_and_einstein_from_ricci(ricci, data)
        batch_started = emit_substage(
            "pure_gauge_control", m, "scalar_einstein", batch_started
        )
        for i in range(3):
            for j in range(3):
                assert_zero(einstein[i, j], x, f"pure gauge Einstein m={m} {(i, j)}")
        assert_zero(scalar, x, f"pure gauge scalar curvature m={m}")
        emit_substage("pure_gauge_control", m, "exact_zero", batch_started)
        results.append(f"m={m}")
        print(f"PASS continuum pure-gauge zero m={m}", flush=True)
        write_checkpoint(output_dir, f"pure_gauge_m{m}", {"m": m}, started)
        completed_batches += 1
        emit_validation_progress(completed_batches)
    return results, completed_batches


def verify_cartesian_tt(output_dir=None, started=None, completed_batches=0):
    emit("STAGE", name="cartesian_tt")
    z = sp.Symbol("z", real=True)
    f = z**3 - 2 * z
    fp = sp.diff(f, z)
    fpp = sp.diff(f, z, 2)
    y = sp.Symbol("y", real=True)
    data = {
        "r": z,
        "x": y,
        "m": 0,
        "metric": sp.eye(3),
        "inverse_diagonal": [sp.S.One, sp.S.One, sp.S.One],
        "volume_density": sp.S.One,
    }
    gamma = independent_connection(data)
    h = sp.diag(0, 1, -1)
    left = tensor_jet(h, data, f, fp, gamma)
    right = tensor_jet(h, data, f, fp, gamma)
    strong = strong_potential_density(left, right, data, gamma)
    weak = prepared.potential_density(left, right, data)
    flux = volume_divergence(potential_flux_vector(left, right, data), data, 0)
    assert sp.expand(strong - weak - flux) == 0
    assert sp.expand(weak - fp**2) == 0
    assert sp.expand(strong + f*fpp) == 0

    bad_delta_g = fpp * h / 2
    bad_strong = sum(f * h[i, j] * bad_delta_g[i, j]
                     for i in range(3) for j in range(3))
    assert sp.expand(bad_strong - weak - flux) != 0
    print("PASS Cartesian TT sign and one-half normalization; wrong-sign control fails", flush=True)
    result = {"profile": "z**3 - 2*z", "wrong_sign_rejected": True}
    write_checkpoint(output_dir, "cartesian_tt", result, started)
    completed_batches += 1
    emit_validation_progress(completed_batches)
    return result, completed_batches


def verify_hermitian_ablation(output_dir=None, started=None, completed_batches=0):
    emit("STAGE", name="ablation_controls")
    data = kinematics.sector(2, 0)
    gamma = prepared.connection(data)
    tensors = data["spatial_tangents"]
    left = prepared.tensor_jet(tensors[0], data, *polynomial_profiles(data, "left"), gamma)
    right = prepared.tensor_jet(tensors[1], data, *polynomial_profiles(data, "right_b"), gamma)
    correct = prepared.potential_density(left, right, data)

    inverse = data["inverse_diagonal"]
    derivative_pair = sum(
        inverse[a] * inverse[i] * inverse[j]
        * chart_conjugate(left["nabla"][a][i, j]) * right["nabla"][a][i, j]
        for a in range(3) for i in range(3) for j in range(3)
    )
    divergence_pair = sum(inverse[j] * chart_conjugate(left["divergence_lower"][j])
                          * right["divergence_lower"][j] for j in range(3))
    broken = data["volume_density"] * (derivative_pair / 2 - divergence_pair)
    assert not exact_zero((correct - broken).subs(sample_point(data, 1)), data["x"])
    print("PASS nonzero ablation control rejects omitted cross/trace terms", flush=True)
    result = {"omitted_cross_trace_terms_rejected": True}
    write_checkpoint(output_dir, "hermitian_ablation", result, started)
    completed_batches += 1
    emit_validation_progress(completed_batches)
    return result, completed_batches


def verify(output_dir=None):
    started = time.monotonic()
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "weak-form-progress.jsonl").write_text("")
    completed_batches = 0
    emit_validation_progress(completed_batches)
    here = Path(__file__).resolve().parent
    derivative_slot_counts, completed_batches = verify_derivative_slots(
        output_dir, started, completed_batches
    )
    gauge_controls, completed_batches = verify_gauge_controls(
        output_dir, started, completed_batches
    )
    cartesian_tt, completed_batches = verify_cartesian_tt(
        output_dir, started, completed_batches
    )
    ablation_controls, completed_batches = verify_hermitian_ablation(
        output_dir, started, completed_batches
    )
    result = {
        "status": "passed",
        "scope": STAGE_SCOPE,
        "physics_source_commit": "157ec0f04866ea36509b325ec65f0a68890f9b17",
        "radial_domain": "none; local derivative-slot and endpoint-flux checks only",
        "finite_assembly_performed": False,
        "operator_spectrum_performed": False,
        "nonlinear_theorem_claimed": False,
        "derivative_slot_counts": derivative_slot_counts,
        "gauge_controls": gauge_controls,
        "cartesian_tt": cartesian_tt,
        "ablation_controls": ablation_controls,
        "python": platform.python_version(),
        "dependencies": {
            name: importlib.metadata.version(name) for name in ("sympy", "numpy")
        },
        "source_sha256": {
            str(path.relative_to(here)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                Path(__file__).resolve(),
                here / "weak_forms.py",
                here / "potential_hessian.py",
                here / "kinematic_blocks.py",
            )
        },
        "seconds": time.monotonic() - started,
    }
    print("PASS bounded weak-form validation", flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    emit_orx({"phase": "setup", "step": 1, "of": 1, "unit": "environment"})
    checks = verify(args.output_dir)
    emit_orx({
        "phase": "export",
        "step": 0,
        "of": 1,
        "unit": "checker result",
        "m": {"elapsed_seconds": checks["seconds"]},
    })
    (args.output_dir / "weak-form-checks.json").write_text(
        json.dumps(checks, indent=2) + "\n"
    )
    emit_orx({
        "phase": "export",
        "step": 1,
        "of": 1,
        "unit": "checker result",
        "m": {"validation_passed": 1, "elapsed_seconds": checks["seconds"]},
    })
    print(json.dumps(checks, indent=2))
