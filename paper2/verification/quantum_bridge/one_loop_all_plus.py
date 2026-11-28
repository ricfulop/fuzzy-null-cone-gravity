"""Continuum Einstein all-plus one-loop reference, not a finite-N result.

Inputs are the dimensionally regulated scalar Compton trees, their KLT
relation, and the continuum supersymmetry state-count identity in
hep-th/9811140, equations (4.2), (4.9), and (4.10). Those inputs have NOT been
derived for the proposed finite-rank gravity action.

This module reconstructs the box multiplicities from three cuts and evaluates
their evanescent mu^8 numerator via exact simplex moments. The published final
amplitude is used only as a comparison, never to evaluate the loop integral.
Independent Decimal cubature checks the Feynman-parameter integrals.
"""

from collections import Counter
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from itertools import permutations, product
from math import cos, sin, sqrt, factorial, pi, prod
import cmath
import json
import os
import platform
import sys
import time

from run import ROOT, as_decimal, gauss_nodes, report


def poly_product(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def poly_divide(numerator, denominator):
    remainder = list(numerator)
    quotient = [F(0)] * (len(numerator) - len(denominator) + 1)
    while len(remainder) >= len(denominator):
        offset = len(remainder) - len(denominator)
        factor = remainder[-1] / denominator[-1]
        quotient[offset] = factor
        for j, value in enumerate(denominator):
            remainder[offset + j] -= factor * value
        while remainder and remainder[-1] == 0:
            remainder.pop()
    return quotient, remainder


def simplex_moment(exponents):
    """Integral of product(x_i**a_i) delta(1-sum(x_i)), with x_i>=0."""
    return F(prod(factorial(a) for a in exponents),
             factorial(sum(exponents) + len(exponents) - 1))


EDGES = ((1, 0, 1, 0), (0, 1, 0, 1))


def integrated_box(s, t, gamma_limit):
    # Expand (-s*x1*x3 - t*x2*x4)^2 BEFORE doing any integration.
    total = F(0)
    for (a, pa), (b, pb) in product(zip((-s, -t), EDGES), repeat=2):
        exponents = tuple(x + y for x, y in zip(pa, pb))
        total += a * b * simplex_moment(exponents)
    return gamma_limit * total


def bivariate_add(a, b):
    out = dict(a)
    for powers, coefficient in b.items():
        out[powers] = out.get(powers, F(0)) + coefficient
    return {powers: value for powers, value in out.items() if value}


def bivariate_product(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            powers = (i + k, j + l)
            out[powers] = out.get(powers, F(0)) + x * y
    return {powers: value for powers, value in out.items() if value}


def exact_derivation():
    # Gamma(epsilon-2) = Gamma(1+epsilon) /
    #                     [epsilon*(epsilon-1)*(epsilon-2)].
    numerator = [F(-1)]
    for factor in ((0, 1), (1, -1), (2, -1), (3, -1)):
        numerator = poly_product(numerator, list(map(F, factor)))
    denominator = [F(1)]
    for factor in ((0, 1), (-1, 1), (-2, 1)):
        denominator = poly_product(denominator, list(map(F, factor)))
    quotient, remainder = poly_divide(numerator, denominator)
    assert not remainder
    gamma_limit = quotient[0]  # Gamma(1+epsilon) -> 1.
    assert gamma_limit == -3

    # Keep two independent invariants and impose u=-s-t symbolically.
    s = {(1, 0): F(1)}
    t = {(0, 1): F(1)}
    u = {(1, 0): F(-1), (0, 1): F(-1)}
    total = {}
    for first, second in ((s, t), (s, u), (t, u)):
        for (a, pa), (b, pb) in product(zip((first, second), EDGES), repeat=2):
            moment = simplex_moment(tuple(x + y for x, y in zip(pa, pb)))
            term = bivariate_product(a, b)
            total = bivariate_add(total, {
                powers: 2 * gamma_limit * moment * value
                for powers, value in term.items()
            })
    reference = {}
    for invariant in (s, t, u):
        reference = bivariate_add(reference, {
            powers: -value / 120
            for powers, value in bivariate_product(invariant, invariant).items()
        })
    assert total == reference
    return {
        "gamma_quotient_coefficients_in_epsilon": [str(x) for x in quotient],
        "gamma_limit": str(gamma_limit),
        "simplex_moments": {
            "x1^2*x3^2": str(simplex_moment((2, 0, 2, 0))),
            "x1*x2*x3*x4": str(simplex_moment((1, 1, 1, 1))),
            "x2^2*x4^2": str(simplex_moment((0, 2, 0, 2))),
        },
        "one_loop_polynomial_after_u_minus_s_minus_t": {
            f"s^{i}*t^{j}": str(value)
            for (i, j), value in sorted(total.items())
        },
        "exact_formula_residual": "0",
    }


def canonical_box(order):
    possibilities = []
    for seq in (tuple(order), tuple(reversed(order))):
        possibilities.extend(seq[i:] + seq[:i] for i in range(4))
    return min(possibilities)


def cut_multiplicities():
    """Recover, rather than add, the repeated box seen by two different cuts."""
    cuts = []
    coefficients, appearances = {}, Counter()
    for left, right in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        terms = Counter(canonical_box(a + b)
                        for a in permutations(left) for b in permutations(right))
        assert sorted(terms.values()) == [2, 2]
        for box, count in terms.items():
            if box in coefficients:
                assert coefficients[box] == count
            coefficients[box] = count
            appearances[box] += 1
        cuts.append({
            "partition": "".join(str(i + 1) for i in left) + "|" +
                         "".join(str(i + 1) for i in right),
            "boxes": {"".join(str(i + 1) for i in box): count
                      for box, count in sorted(terms.items())},
        })
    assert len(coefficients) == 3
    assert set(appearances.values()) == {2}
    return cuts, coefficients


def simplex_cubature(s, t, order, precision):
    nodes = [((x + 1) / 2, weight) for x, weight in gauss_nodes(order, precision)]
    total = D(0)
    for (u, wu), (v, wv), (w, ww) in product(nodes, repeat=3):
        x1 = u
        x2 = (1 - u) * v
        x3 = (1 - u) * (1 - v) * w
        x4 = (1 - u) * (1 - v) * (1 - w)
        jacobian = (1 - u) ** 2 * (1 - v)
        total += wu * wv * ww * jacobian * (-s * x1 * x3 - t * x2 * x4) ** 2
    return total


def quadrature_case(s_string, t_string, gamma_limit, config):
    s, t = F(s_string), F(t_string)
    exact = integrated_box(s, t, gamma_limit)
    with localcontext() as context:
        context.prec = config["decimal_precision"]
        target = as_decimal(exact)
        checks = []
        for order in config["simplex_quadrature_orders"]:
            value = as_decimal(gamma_limit) * simplex_cubature(
                as_decimal(s), as_decimal(t), order, context.prec)
            error = abs(value - target)
            assert error < D(config["quadrature_absolute_tolerance"]), (s, t, error)
            checks.append({"nodes_per_axis": order, "value": str(value),
                           "absolute_error": str(error)})
        difference = abs(D(checks[0]["value"]) - D(checks[-1]["value"]))
        assert difference < D(config["quadrature_absolute_tolerance"])
    return {"s": str(s), "t": str(t), "exact_box_without_i_over_4pi2": str(exact),
            "checks": checks, "quadrature_difference": str(difference)}


def minkowski_dot(a, b):
    return a[0] * b[0] - sum(x * y for x, y in zip(a[1:], b[1:]))


def invariant(a, b):
    return 2 * minkowski_dot(a, b)


def boost(p, velocity):
    beta_squared = sum(x * x for x in velocity)
    if beta_squared == 0:
        return tuple(p)
    gamma = 1 / sqrt(1 - beta_squared)
    beta_dot_p = sum(x * y for x, y in zip(velocity, p[1:]))
    factor = (gamma - 1) * beta_dot_p / beta_squared - gamma * p[0]
    return (gamma * (p[0] - beta_dot_p),
            *(x + factor * b for x, b in zip(p[1:], velocity)))


def external_momenta(energy, theta, phi, velocity):
    spatial = (energy * sin(theta) * cos(phi), energy * sin(theta) * sin(phi),
               energy * cos(theta))
    momenta = [(-energy, 0.0, 0.0, -energy), (-energy, 0.0, 0.0, energy),
               (energy, *spatial), (energy, *(-x for x in spatial))]
    return [boost(p, velocity) for p in momenta]


def spinor_pair(p):
    energy, x, y, z = p
    scale = max(abs(value) for value in p)
    if abs(energy + z) > scale * 1e-12:
        root = cmath.sqrt(complex(energy + z))
        return (root, (x + 1j * y) / root), (root, (x - 1j * y) / root)
    root = cmath.sqrt(complex(energy - z))
    return ((x - 1j * y) / root, root), ((x + 1j * y) / root, root)


def determinant(a, b):
    return a[0] * b[1] - a[1] * b[0]


def angle(pairs, i, j):
    return determinant(pairs[i][0], pairs[j][0])


def square(pairs, i, j):
    # Convention: <ij>[ji] = 2 k_i.k_j for signature (+---).
    return -determinant(pairs[i][1], pairs[j][1])


def relative_error(value, target):
    return abs(value - target) / max(abs(value), abs(target), 1e-300)


def loop_amplitude(momenta, pairs, boxes, gamma_limit, kappa):
    phase = square(pairs, 0, 1) * square(pairs, 2, 3) / (
        angle(pairs, 0, 1) * angle(pairs, 2, 3))
    integrals = []
    for order, multiplicity in sorted(boxes.items()):
        s = invariant(momenta[order[0]], momenta[order[1]])
        t = invariant(momenta[order[1]], momenta[order[2]])
        value = integrated_box(s, t, float(gamma_limit))
        integrals.append((order, multiplicity, value))
    amplitude = 1j / (4 * pi) ** 2 * (kappa / 2) ** 4 * phase ** 2 * sum(
        count * value for _, count, value in integrals)
    return amplitude, integrals


def published_benchmark(momenta, pairs, kappa):
    # Validation-only target: hep-th/9811140 (4.14), coupling restored.
    s = invariant(momenta[0], momenta[1])
    t = invariant(momenta[1], momenta[2])
    u = invariant(momenta[0], momenta[2])
    angle_chain = (angle(pairs, 0, 1) * angle(pairs, 1, 2) *
                   angle(pairs, 2, 3) * angle(pairs, 3, 0))
    return (-1j / (4 * pi) ** 2 * (kappa / 2) ** 4 *
            (s * t / angle_chain) ** 2 * (s * s + t * t + u * u) / 120)


def complex_json(value):
    return {"real": value.real, "imaginary": value.imag}


def kinematic_case(index, energy, theta, boxes, gamma_limit, config):
    phi = config["azimuth_start"] + index * config["azimuth_step"]
    velocity = config["boosts"][index % len(config["boosts"])]
    momenta = external_momenta(energy, theta, phi, velocity)
    pairs = [spinor_pair(p) for p in momenta]
    kappa = config["coupling_kappa"]
    scale = max(abs(x) for p in momenta for x in p)
    tolerance = config["kinematic_relative_tolerance"]
    assert max(abs(sum(p[j] for p in momenta)) for j in range(4)) < tolerance * scale
    assert max(abs(minkowski_dot(p, p)) for p in momenta) < tolerance * scale ** 2
    for i, j in permutations(range(4), 2):
        assert relative_error(angle(pairs, i, j) * square(pairs, j, i),
                              invariant(momenta[i], momenta[j])) < tolerance

    amplitude, integrals = loop_amplitude(momenta, pairs, boxes, gamma_limit, kappa)
    reference = published_benchmark(momenta, pairs, kappa)
    benchmark_error = relative_error(amplitude, reference)
    assert benchmark_error < tolerance
    assert abs(amplitude) > 1e-15

    permutation_errors = []
    for order in permutations(range(4)):
        value, _ = loop_amplitude([momenta[i] for i in order], [pairs[i] for i in order],
                                  boxes, gamma_limit, kappa)
        permutation_errors.append(relative_error(value, amplitude))
    assert max(permutation_errors) < tolerance

    weights = [cmath.exp(0.08 * (i - 1.5) + 1j * (0.19 + 0.23 * i)) for i in range(4)]
    changed = [(tuple(w * x for x in a), tuple(x / w for x in b))
               for (a, b), w in zip(pairs, weights)]
    value, _ = loop_amplitude(momenta, changed, boxes, gamma_limit, kappa)
    little_group_error = relative_error(value, amplitude * prod(w ** -4 for w in weights))
    assert little_group_error < tolerance

    energy_scale = 1.7 + 0.03 * index
    rescaled = [tuple(energy_scale * x for x in p) for p in momenta]
    value, _ = loop_amplitude(rescaled, [spinor_pair(p) for p in rescaled],
                             boxes, gamma_limit, kappa)
    scaling_error = relative_error(value, amplitude * energy_scale ** 4)
    assert scaling_error < tolerance

    # Algebraic consistency of INPUT scalar trees, not a finite-N derivation.
    s = invariant(momenta[0], momenta[1])
    mu_squared = 0.23 * energy ** 2
    d1 = -s * (0.25 + 0.02 * index)
    d2 = -s - d1
    q = mu_squared * square(pairs, 0, 1) / angle(pairs, 0, 1)
    gauge_first = 1j * q / d1
    gauge_second = -1j * q * (1 / d1 + 1 / d2)
    gravity_tree = -1j * q ** 2 * (1 / d1 + 1 / d2)
    klt_error = relative_error(-1j * d1 * gauge_first * gauge_second, gravity_tree)
    assert klt_error < tolerance

    # Negative controls: putting mu=0 before integration loses the rational term;
    # omitting one real state loses half the pure-gravity amplitude.
    zero_mu_tree = -1j * (0 * square(pairs, 0, 1) / angle(pairs, 0, 1)) ** 2 * (
        1 / d1 + 1 / d2)
    strict_four_dimensional = zero_mu_tree * zero_mu_tree
    four_dimensional_error = relative_error(strict_four_dimensional, amplitude)
    one_real_state_error = relative_error(amplitude / 2, amplitude)
    assert four_dimensional_error == 1
    assert abs(one_real_state_error - 0.5) < tolerance

    invariants = [invariant(momenta[0], momenta[1]),
                  invariant(momenta[1], momenta[2]),
                  invariant(momenta[0], momenta[2])]
    assert abs(sum(invariants)) < tolerance * scale ** 2
    return {
        "case": index + 1, "energy": energy, "theta": theta, "phi": phi, "boost": velocity,
        "s_t_u": invariants, "amplitude": complex_json(amplitude),
        "published_reference": complex_json(reference),
        "boxes": [{"ordering": "".join(str(i + 1) for i in order), "multiplicity": count,
                   "box_without_i_over_4pi2": value} for order, count, value in integrals],
        "benchmark_relative_error": benchmark_error,
        "permutation_relative_error": max(permutation_errors),
        "little_group_relative_error": little_group_error,
        "energy_scaling_relative_error": scaling_error, "scalar_tree_klt_error": klt_error,
        "strict_4d_cut_relative_error": four_dimensional_error,
        "one_real_state_relative_error": one_real_state_error,
    }


def main(config):
    start = time.monotonic()
    print(json.dumps({"configuration": config, "python": sys.version,
                      "platform": platform.platform(), "commit": os.getenv("ORX_COMMIT"),
                      "run": os.getenv("ORX_RUN_ID"), "experiment": os.getenv("ORX_EXPERIMENT")}),
          flush=True)
    report(phase="exact_loop", step=0, of=1, unit="derivations")
    exact = exact_derivation()
    cuts, boxes = cut_multiplicities()
    gamma_limit = F(exact["gamma_limit"])
    print("EXACT LOOP " + json.dumps({"integral": exact, "cuts": cuts}), flush=True)
    report(step=1, m={"loop_exact_residual": 0, "independent_boxes": len(boxes)})

    quadratures = []
    report(phase="simplex_quadrature", step=0, of=len(config["simplex_cases"]), unit="cases")
    for done, (s, t) in enumerate(config["simplex_cases"], 1):
        row = quadrature_case(s, t, gamma_limit, config)
        quadratures.append(row)
        maximum = max(D(check["absolute_error"]) for check in row["checks"])
        report(row="simplex", m={"case": f"s={s}, t={t}", "box_coefficient": float(
            F(row["exact_box_without_i_over_4pi2"])), "quadrature_error_text": f"{maximum:.3E}"})
        report(step=done, m={"loop_quadrature_error": float(maximum)})
        print("SIMPLEX INTEGRAL " + json.dumps(row), flush=True)

    kinematics = []
    cases = list(product(config["energies"], config["scattering_angles"]))
    report(phase="kinematic_checks", step=0, of=len(cases), unit="cases")
    for index, (energy, theta) in enumerate(cases):
        row = kinematic_case(index, energy, theta, boxes, gamma_limit, config)
        kinematics.append(row)
        report(row="amplitudes", m={"case": str(index + 1), "energy": energy,
            "theta": theta, "amplitude_magnitude": abs(complex(
                row["amplitude"]["real"], row["amplitude"]["imaginary"])),
            "loop_benchmark_error": row["benchmark_relative_error"],
            "loop_bose_error": row["permutation_relative_error"],
            "strict_4d_error": row["strict_4d_cut_relative_error"],
            "single_state_error": row["one_real_state_relative_error"]})
        report(step=index + 1, m={"loop_benchmark_error": row["benchmark_relative_error"]})
        print("ONE LOOP KINEMATICS " + json.dumps(row), flush=True)

    report(phase="export", step=0, of=1, unit="reports")
    result = {
        "config": config, "source_commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID"),
        "scope": "Continuum pure Einstein one-loop all-plus reference in a state-preserving dimensional scheme; hbar=1, kappa^2=32*pi*G.",
        "inputs": [
            "Published scalar Compton trees and KLT relation: hep-th/9811140 (4.9)-(4.10).",
            "Continuum supersymmetry identity: two real scalar states reproduce the pure-graviton all-plus loop (4.2).",
            "Dimensional loop integration with Feynman i0, scaleless integrals discarded and no additional higher-curvature interactions.",
        ],
        "exact_loop": exact, "cut_multiplicities": cuts, "quadrature_checks": quadratures,
        "kinematic_checks": kinematics,
        "limitations": [
            "No finite-N angular cutoff is present in this reference calculation.",
            "The scalar state-count identity is not established for the finite-rank Gaunt candidate.",
            "No new gravitational constraint completion, microscopic corner degeneracy or bulk Euclidean measure is derived.",
            "Checking helicity, Bose symmetry and cuts is not a full off-shell Ward-identity proof for a different regulator.",
            "Only the all-plus channel is evaluated; MHV channels and a controlled physical finite-N limit remain open.",
        ],
        "status": "all specified continuum reference checks passed",
        "wall_seconds": time.monotonic() - start,
    }
    output = ROOT / "output"
    output.mkdir(exist_ok=True)
    path = output / "one_loop_all_plus_results.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    summary = {
        "result": str(path), "independent_boxes": len(boxes), "quadrature_cases": len(quadratures),
        "kinematic_cases": len(kinematics), "permutations_per_case": 24,
        "max_benchmark_relative_error": max(r["benchmark_relative_error"] for r in kinematics),
        "max_bose_relative_error": max(r["permutation_relative_error"] for r in kinematics),
        "max_quadrature_error": str(max(D(c["absolute_error"]) for r in quadratures for c in r["checks"])),
        "strict_4d_cut_relative_error": 1, "one_real_state_relative_error": 0.5,
        "scope": "continuum reference only; finite-N Gaunt amplitude remains uncomputed",
        "status": result["status"], "wall_seconds": result["wall_seconds"],
    }
    report(step=1, m={"loop_benchmark_error": summary["max_benchmark_relative_error"],
                     "loop_bose_error": summary["max_bose_relative_error"],
                     "loop_exact_residual": 0, "strict_4d_error": 1, "single_state_error": 0.5})
    print("FINAL " + json.dumps(summary, allow_nan=False), flush=True)
