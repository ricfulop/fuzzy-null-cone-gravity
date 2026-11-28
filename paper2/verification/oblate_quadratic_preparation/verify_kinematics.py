"""Prepared local operator check; not executed during preparation.

This checks tensor maps and canonical normalization, not a gravity spectrum.
Requires separate user approval before running. The pinned reference is read
only through its field-construction function; its verifiers are not called.
"""

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import time

import sympy as sp

import kinematic_blocks as prepared


REFERENCE_HASH = "55ba3d4ba6c44104668362258aa96465bd485fef8f407f65283e12aa0c643ebe"


def exact_zero(expression, x):
    """Cancel on the chart -1<x<1; no tolerances or numerical probes."""
    # Expand only the polynomial base, retaining the original half-integer
    # exponent. Equal radicands then have the same representation without
    # splitting radicals or assuming sqrt(z**2) == z.
    expression = sp.sympify(expression).replace(
        lambda term: (term.is_Pow and term.exp.is_Rational
                      and term.exp.q == 2 and term.base.is_polynomial(x)),
        lambda term: sp.Pow(sp.expand(term.base), term.exp),
    )
    d = 1 - x**2
    for weight in (sp.S.One, sp.sqrt(d), 1 / sp.sqrt(d)):
        value = sp.cancel(sp.powsimp(expression * weight, force=False))
        value = sp.simplify(value)
        if value == 0 or value.is_zero is True:
            return True
    return False


def verify_exact_zero_controls():
    """Replay the archived residual and reject genuinely nonzero identities."""
    x = sp.Symbol("x", real=True)
    factored_root = sp.Pow(-(x - 1) * (x + 1), sp.S.Half, evaluate=False)
    expanded_root = sp.Pow(1 - x**2, sp.S.Half, evaluate=False)
    residual = sp.Add(
        sp.Mul(3, sp.I, factored_root, x - 1, x + 1, evaluate=False),
        sp.Mul(-3, sp.I, expanded_root, x**2 - 1, evaluate=False),
        evaluate=False,
    )
    zero_controls = {
        "archived_m1_reference_residual": residual,
        "unevaluated_zero_times_I": sp.Mul(0, sp.I, evaluate=False),
        "reciprocal_chart_roots": 1 / factored_root - 1 / expanded_root,
    }
    nonzero_controls = {
        "unit_offset": residual + 1,
        "exact_small_offset": residual + sp.Rational(1, 10**30),
        "chart_polynomial_offset": residual + 1 - x**2,
        "wrong_relative_radical_sign": factored_root + expanded_root,
        "square_root_branch": sp.Pow(x**2, sp.S.Half, evaluate=False) - x,
    }
    for label, expression in zero_controls.items():
        assert exact_zero(expression, x), f"Exact-zero regression failed: {label}"
    for label, expression in nonzero_controls.items():
        assert not exact_zero(expression, x), f"Nonzero control accepted: {label}"
    print("PASS exact-zero controls: archived residual and structural zeros; nonzero and branch controls rejected", flush=True)
    return {"zero_controls": list(zero_controls), "nonzero_controls": list(nonzero_controls)}


def assert_zero_matrix(matrix, x, label):
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            if not exact_zero(matrix[i, j], x):
                raise AssertionError(f"{label} entry {(i, j)}: {matrix[i, j]}")


def direct_lie(metric, vector, r, x, m):
    """Coordinate Lie derivative from its definition, including radial jets."""
    def derivative_amplitude(value, index):
        if index == 0:
            return sp.diff(value, r)
        if index == 1:
            return sp.diff(value, x)
        return sp.I * m * value

    result = sp.zeros(3)
    for i in range(3):
        for j in range(3):
            transport = vector[0] * sp.diff(metric[i, j], r) + vector[1] * sp.diff(metric[i, j], x)
            result[i, j] = transport + sum(
                metric[k, j] * derivative_amplitude(vector[k], i)
                + metric[i, k] * derivative_amplitude(vector[k], j)
                for k in range(3)
            )
    return result


def direct_metric_tangent(data, label):
    """Differentiate the actual nonlinear spatial coefficient map."""
    r, x = data["r"], data["x"]
    kind, ell, m = label
    y, polar, axial, electric, magnetic = prepared.harmonic_amplitudes(ell, m, x)
    parameter = sp.Symbol("variation", real=True)
    delta_lambda = y if kind == "lambda" else 0
    delta_gamma = sp.zeros(2)
    delta_vector = sp.zeros(2, 1)
    if kind == "X":
        delta_gamma = sp.diag(y / (1 - x**2), (1 - x**2) * y)
    elif kind.startswith("V_"):
        delta_vector = polar if kind == "V_polar" else axial
    elif kind.startswith("C_"):
        delta_gamma = electric if kind == "C_polar" else magnetic
    gamma = data["metric"][1:, 1:] + parameter * delta_gamma
    vector = parameter * delta_vector
    metric = sp.zeros(3)
    metric[0, 0] = data["metric"][0, 0] + parameter * delta_lambda + (vector.T * gamma * vector)[0]
    mixed = gamma * vector
    metric[0, 1] = metric[1, 0] = mixed[0]
    metric[0, 2] = metric[2, 0] = mixed[1]
    metric[1:, 1:] = gamma
    return metric.diff(parameter).subs(parameter, 0)


def verify():
    started = time.monotonic()
    zero_controls = verify_exact_zero_controls()
    here = Path(__file__).resolve().parent
    reference = here / "reference/all_modes_157ec0f.py"
    if hashlib.sha256(reference.read_bytes()).hexdigest() != REFERENCE_HASH:
        raise ValueError("Pinned Paper I reference hash does not match")
    spec = importlib.util.spec_from_file_location("paper1_pinned_fields", reference)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    counts = {"complex_spatial_tangents": 0, "complex_shift_columns": 0, "r1_reference_blocks": 0}

    for m in range(3):
        data = prepared.sector(2, m)
        r, x, epsilon = data["r"], data["x"], data["epsilon"]
        T = data["spatial_tangents"]
        U0 = data["shift_columns_without_profile_derivatives"]
        J1 = data["shift_derivative_injection"]
        for label, tensor in zip(data["spatial_labels"], T):
            assert_zero_matrix(tensor - direct_metric_tangent(data, label), x, f"tangent {label}")
            counts["complex_spatial_tangents"] += 1

        profile = sp.Function("profile")(r)
        for j, (kind, ell, order) in enumerate(data["shift_labels"]):
            y, polar, axial, _, _ = prepared.harmonic_amplitudes(ell, order, x)
            if kind == "b":
                vector = sp.Matrix([profile * y / data["metric"][0, 0], 0, 0])
            else:
                angular = polar if kind == "W_polar" else axial
                vector = sp.Matrix([0, profile * angular[0], profile * angular[1]])
            direct = direct_lie(data["metric"], vector, r, x, m)
            derivative_tensor = sum((T[i] * J1[i, j] for i in range(len(T))), sp.zeros(3))
            candidate = profile * U0[j] + sp.diff(profile, r) * derivative_tensor
            assert_zero_matrix(direct - candidate, x, f"Lie column {(kind, ell, order)}")
            counts["complex_shift_columns"] += 1

        old_T, old_labels, old_U, old_shift_labels, _ = module.fields(2, m)
        assert [label[:2] for label in data["spatial_labels"]] == [tuple(label) for label in old_labels]
        assert [label[:2] for label in data["shift_labels"]] == [tuple(label) for label in old_shift_labels]
        substitute = {module.x: x, module.eps: epsilon}
        for tensors, originals in ((T, old_T), (U0, old_U)):
            for tensor, original in zip(tensors, originals):
                assert_zero_matrix(tensor.subs(r, 1) - original.subs(substitute), x, f"r=1, m={m}")
        counts["r1_reference_blocks"] += 1
        print(f"PASS geometry and pinned-reference sector L=2 m={m}", flush=True)

    # Exercise the actual density path. At a deformed rational point compare
    # with an independent full-matrix contraction, including its determinant.
    data = prepared.sector(2, 2)
    r, x, epsilon = data["r"], data["x"], data["epsilon"]
    T = data["spatial_tangents"]
    U = data["shift_columns_without_profile_derivatives"]
    densities = prepared.angular_density_blocks(data)
    point = {r: sp.Integer(2), x: sp.Rational(1, 3), epsilon: sp.Rational(1, 5)}
    metric = data["metric"].subs(point)
    inverse = metric.inv()
    volume = sp.sqrt(metric.det())
    for key, left, right in (("S_density", T, T), ("C0_density", T, U), ("U0_norm_density", U, U)):
        direct = sp.Matrix(len(left), len(right), lambda i, j:
            volume * (sp.trace(inverse * left[i].subs(point).conjugate() * inverse * right[j].subs(point))
                      - sp.trace(inverse * left[i].subs(point).conjugate()) * sp.trace(inverse * right[j].subs(point))))
        assert_zero_matrix(densities[key].subs(point) - direct, x, f"deformed point {key}")

    def exact_angular_mean(expression):
        polynomial = sp.Poly(sp.cancel(sp.powsimp(expression, force=True)), x)
        return sum(coefficient / (powers[0] + 1) for powers, coefficient in polynomial.terms() if powers[0] % 2 == 0)

    round_subs = {r: sp.Integer(2), epsilon: sp.S.Zero}
    integrated = {key: matrix.subs(round_subs).applyfunc(exact_angular_mean)
                  for key, matrix in densities.items()}
    # Map r=2 tensors to the pinned r=1 frame, retaining volume r^2.
    frame = sp.diag(1, sp.Rational(1, 2), sp.Rational(1, 2))
    for key, left, right in (("S_density", T, T), ("C0_density", T, U), ("U0_norm_density", U, U)):
        reference_gram = sp.Matrix(len(left), len(right), lambda i, j: 4 * module.round_pair(
            (frame * left[i].subs(round_subs) * frame).subs(x, module.x),
            (frame * right[j].subs(round_subs) * frame).subs(x, module.x)))
        assert_zero_matrix(integrated[key] - reference_gram, x, f"round integrated {key}")

    # Pass derivative columns through the same density/integration path before
    # projection; the extended normal Hessian must have zero derivative legs.
    J1 = data["shift_derivative_injection"]
    derivative_columns = [sum((T[i] * J1[i, j] for i in range(len(T))), sp.zeros(3)) for j in range(J1.cols)]
    extended = U + derivative_columns
    Cext = sp.Matrix([[prepared.dewitt_density(a, b, data) for b in extended] for a in T]).subs(round_subs).applyfunc(exact_angular_mean)
    Hext = sp.Matrix([[prepared.dewitt_density(a, b, data) for b in extended] for a in extended]).subs(round_subs).applyfunc(exact_angular_mean)
    Sround = integrated["S_density"]
    assert_zero_matrix(Cext[:, len(U):] - Sround*J1, x, "integrated derivative coupling")
    Aext = (Cext.conjugate().T * Sround.inv() * Cext - Hext) / 2
    assert_zero_matrix(Aext[:, len(U):], x, "normal-Hessian derivative columns")
    assert_zero_matrix(Aext[len(U):, :], x, "normal-Hessian derivative rows")
    print("PASS DeWitt densities at a deformed point and integrated round sector; derivative Schur legs vanish", flush=True)

    # Canonical normalization control with an indefinite invertible spatial Gram.
    S = sp.Matrix([[0, -2], [-2, -2]])
    C = sp.Matrix([[1, 2], [3, -1]])
    H = sp.Matrix([[4, 1], [1, 5]])
    blocks = prepared.canonical_angular_blocks(S, C, H, sp.zeros(2))
    p = sp.Matrix(sp.symbols("p0:2", real=True))
    shift = sp.Matrix(sp.symbols("s0:2", real=True))
    velocity = S.inv() * (2 * p + C * shift)
    lagrangian = (velocity.T * S * velocity)[0] / 4 - (velocity.T * C * shift)[0] / 2 + (shift.T * H * shift)[0] / 4
    direct_hamiltonian = (p.T * velocity)[0] - lagrangian
    assembled = (p.T * blocks["K_pp"] * p)[0] / 2 + (p.T * blocks["J0"] * shift)[0] + (shift.T * blocks["A"] * shift)[0] / 2
    assert sp.expand(direct_hamiltonian - assembled) == 0
    omitted_shift_norm = assembled + (shift.T * H * shift)[0] / 4
    assert sp.expand(direct_hamiltonian - omitted_shift_norm) != 0
    print("PASS indefinite Legendre normalization; omitted shift-norm control fails as intended", flush=True)

    phi = sp.Symbol("phi", real=True)
    real_fields = sp.Matrix([sp.sqrt(2) * sp.cos(phi), -sp.sqrt(2) * sp.sin(phi), sp.sqrt(2) * sp.sin(phi), sp.sqrt(2) * sp.cos(phi)])
    direct_real_gram = (real_fields * real_fields.T).applyfunc(lambda value: sp.integrate(value, (phi, 0, 2 * sp.pi)) / (2 * sp.pi))
    complex_gram = sp.Matrix([[1, sp.I], [-sp.I, 1]])
    assert prepared.realify_block(complex_gram, 1) == direct_real_gram
    print("PASS positive-sine realification including imaginary cross block", flush=True)
    return {
        "status": "passed", "counts": counts,
        "exact_zero_controls": zero_controls,
        "legendre_control": "passed; manufactured indefinite Gram, not gravitational data",
        "realification_control": "passed; synthetic imaginary cross block",
        "density_controls": "deformed rational point, exact integrated L=2,m=2 round r=2 sector, and derivative-column Schur cancellation",
        "physics_source_commit": "157ec0f04866ea36509b325ec65f0a68890f9b17",
        "python": platform.python_version(),
        "dependencies": {name: importlib.metadata.version(name) for name in ("sympy", "numpy")},
        "source_sha256": {str(path.relative_to(here)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (Path(__file__).resolve(), here / "kinematic_blocks.py", reference)},
        "seconds": time.monotonic() - started,
        "scope": "Kinematic tensor maps, selected density/integration controls, normalization, and r=1 reference; no potential Hessian, radial spectrum, nonlinear consistency, or physical mode count",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "kinematic-checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
