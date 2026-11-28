"""Independent polynomial-oscillator controls for phase_space_loop.

These are finite-mode combinatorial, quadrature and canonical-coordinate
bookkeeping checks. They are NOT a finite-N graviton amplitude, a check of the
gravitational constraint reduction, or a proof of gravitational consistency.
No files are written and no calculation is started on import.
"""

from fractions import Fraction as F
from itertools import combinations, permutations
from math import cos, factorial, fsum, isfinite, sin

try:
    from .phase_space_loop import integrate_vertex, ordered_partitions
except ImportError:
    from phase_space_loop import integrate_vertex, ordered_partitions


def _json_complex(value):
    return {"real": float(value.real), "imaginary": float(value.imag)}


def _comparison(value, reference, relative_tolerance, absolute_tolerance):
    error = abs(value - reference)
    return {
        "value": _json_complex(value), "reference": _json_complex(complex(reference)),
        "absolute_error": error,
        "relative_error": error / abs(reference) if reference != 0 else None,
        "passed": error <= absolute_tolerance + relative_tolerance * abs(reference),
    }


def _subsets(n):
    return (frozenset(subset) for size in range(1, n + 1)
            for subset in combinations(range(n), size))


def _scalar_blocks(n, couplings):
    return {subset: ((complex(couplings.get(len(subset) + 2, 0)), 0j), (0j, 0j))
            for subset in _subsets(n)}


def _static_reference(omega, couplings, n):
    """Exact rational series of (1/2)*sqrt(V''(q)), independent of diagrams.

    V(q)=omega^2*q^2/2 + sum_m g_m*q^m/m!, for m=3,...,6.
    If sqrt(V'')=sum a_j*q^j, solve its square coefficient by coefficient.
    Return the n-th derivative, not the ordinary Taylor coefficient.
    """
    series = [F(omega)]
    for degree in range(1, n + 1):
        hessian = F(couplings.get(degree + 2, 0)) / factorial(degree)
        convolution = sum((series[j] * series[degree - j] for j in range(1, degree)), F(0))
        series.append((hessian - convolution) / (2 * series[0]))
    return F(factorial(n), 2) * series[n]


def _bubble_integral(a, b, omega):
    """Integral dnu/(2pi) / [(nu^2+a^2)((nu+Omega)^2+b^2)]."""
    return (a + b) / (2 * a * b * ((a + b) ** 2 + omega ** 2))


def _quartic_blocks(n, direction, coupling):
    # H4=g*(u.Q)^4/4!, external Q vectors chosen so u.z_r=1.
    size = 2 * len(direction)
    block = [[0j] * size for _ in range(size)]
    for i, left in enumerate(direction):
        for j, right in enumerate(direction):
            block[2 * i][2 * j] = coupling * left * right
    return {frozenset(subset): block for subset in combinations(range(n), 2)}


def _quartic_reference(frequencies, external, direction, coupling):
    # Six ordered (2,2) partitions group into the three exchange channels.
    return -coupling ** 2 / 2 * fsum(
        direction[i] ** 2 * direction[j] ** 2
        * _bubble_integral(a, b, external[0] + external[partner])
        for partner in (1, 2, 3)
        for i, a in enumerate(frequencies) for j, b in enumerate(frequencies))


def _two_point_case(frequencies, external_frequency, order):
    """Pure-Q cubic loop plus a quartic Q/P tadpole, evaluated independently.

    H3=.4*(v.Q)^3/3!, H4=.7*(r.z)^4/4!, with distinct external directions.
    This exercises exactly the two-point interface used by the reduced model.
    """
    legs = ((1.0, 0.2, -0.3, 0.1), (0.4, 0.7, 0.1, -0.2))
    v, r = (0.8, 0.4), (-0.2, 0.6, 0.7, -0.4)
    blocks = {}
    for label, leg in enumerate(legs):
        coefficient = 0.4 * sum(v[i] * leg[2 * i] for i in range(2))
        matrix = [[0j] * 4 for _ in range(4)]
        for i in range(2):
            for j in range(2):
                matrix[2 * i][2 * j] = coefficient * (v[i] * v[j])
        blocks[frozenset((label,))] = matrix
    coefficient = 0.7
    for leg in legs:
        coefficient *= sum(a * b for a, b in zip(r, leg))
    tadpole_block = [[coefficient * (a * b) for b in r] for a in r]
    blocks[frozenset((0, 1))] = tadpole_block
    tadpole = sum(tadpole_block[2 * i][2 * i] / (4 * omega)
                  + omega * tadpole_block[2 * i + 1][2 * i + 1] / 4
                  for i, omega in enumerate(frequencies))
    b0, b1 = blocks[frozenset((0,))], blocks[frozenset((1,))]
    bubble = -sum(b0[2 * i][2 * j] * b1[2 * j][2 * i]
                  * _bubble_integral(a, b, external_frequency)
                  for i, a in enumerate(frequencies) for j, b in enumerate(frequencies)) / 2
    value = integrate_vertex(frequencies, (external_frequency, -external_frequency), blocks, order)
    return value, tadpole + bubble, tadpole, bubble


def _mixed_blocks(n):
    """H_m=g_m*(v_m.z)^m/m! gives integrable, mixed Q/P Hessian blocks.

    The real polynomial Hamiltonian has two oscillators. Distinct complex
    external directions exercise complex multilinearity; each degree uses a
    different real direction so the loop factors need not commute.
    """
    directions = {
        3: (0.8, -0.3, 0.4, 0.2), 4: (-0.2, 0.6, 0.7, -0.4),
        5: (0.5, 0.2, -0.6, 0.3), 6: (0.3, -0.5, 0.2, 0.8),
    }
    couplings = {3: 0.4, 4: 0.7, 5: -0.3, 6: 0.9}
    legs = ((1, 0.2j, -0.3, 0.1), (0.4, 0.7, 0.1j, -0.2),
            (-0.2, 0.1, 0.8, 0.3j), (0.3j, -0.4, 0.2, 0.9))[:n]
    result = {}
    for subset in _subsets(n):
        degree = len(subset) + 2
        vector = directions[degree]
        coefficient = couplings[degree]
        for label in sorted(subset):
            coefficient *= sum(a * b for a, b in zip(vector, legs[label]))
        # Form each symmetric entry once, to avoid arithmetic-order asymmetry.
        matrix = [[0j] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i, 4):
                matrix[i][j] = matrix[j][i] = coefficient * (vector[i] * vector[j])
        result[subset] = matrix
    return result


def _canonical_rotation(frequencies, angles):
    # Old coordinates z=T*z_new. T preserves both the symplectic form and H2.
    size = 2 * len(frequencies)
    matrix = [[0.0] * size for _ in range(size)]
    for i, (omega, angle) in enumerate(zip(frequencies, angles)):
        q, p = 2 * i, 2 * i + 1
        matrix[q][q] = matrix[p][p] = cos(angle)
        matrix[q][p] = sin(angle) / omega
        matrix[p][q] = -omega * sin(angle)
    return matrix


def _congruence(matrix, transform):
    size = len(transform)
    out = [[0j] * size for _ in range(size)]
    for i in range(size):
        for j in range(i, size):
            out[i][j] = out[j][i] = sum(
                transform[a][i] * matrix[a][b] * transform[b][j]
                for a in range(size) for b in range(size))
    return out


def _rotation_residuals(frequencies, transform):
    size = len(transform)
    symplectic = [[0.0] * size for _ in range(size)]
    free = [[0.0] * size for _ in range(size)]
    for i, omega in enumerate(frequencies):
        q, p = 2 * i, 2 * i + 1
        symplectic[q][p], symplectic[p][q] = 1, -1
        free[q][q], free[p][p] = omega ** 2, 1

    def residual(matrix):
        return max(abs(sum(transform[a][i] * matrix[a][b] * transform[b][j]
                           for a in range(size) for b in range(size)) - matrix[i][j])
                   for i in range(size) for j in range(size))

    return {"symplectic_residual": residual(symplectic),
            "free_hamiltonian_residual": residual(free)}


def _permuted_inputs(external, blocks, permutation):
    # permutation[new_label] is the corresponding old label.
    old_to_new = {old: new for new, old in enumerate(permutation)}
    return (tuple(external[old] for old in permutation),
            {frozenset(old_to_new[old] for old in subset): block
             for subset, block in blocks.items()})


def _checks_at_order(order, relative_tolerance, absolute_tolerance):
    def compare(value, reference):
        return _comparison(value, reference, relative_tolerance, absolute_tolerance)

    couplings = {3: F(2, 5), 4: F(7, 10), 5: F(-3, 10), 6: F(9, 10)}
    cases = (
        ("all_H3_through_H6", couplings),
        ("cubic_only_box", {3: couplings[3]}),
        ("quartic_only_bubble", {4: couplings[4]}),
        ("cubic_quartic_triangle_plus_box_and_bubble", {3: couplings[3], 4: couplings[4]}),
        ("cubic_quintic_pair_plus_box", {3: couplings[3], 5: couplings[5]}),
        ("sextic_only_tadpole", {6: couplings[6]}),
    )
    omega = F(5, 4)
    scalar = []
    for name, interactions in cases:
        for n in (2, 3, 4):
            exact = _static_reference(omega, interactions, n)
            blocks = _scalar_blocks(n, interactions)
            value = integrate_vertex([omega], [0] * n, blocks, order)
            scalar.append({"case": name, "n": n, "exact_derivative": str(exact),
                           **compare(value, float(exact))})

    external_cases = ((0.2, -0.45, 0.8, -0.55), (0.35, 0.65, -0.25, -0.75))
    bubble = []
    for frequencies, direction in (((1.25,), (1.0,)), ((0.8, 1.4), (0.6, 0.8))):
        for external in external_cases:
            blocks = _quartic_blocks(4, direction, 0.7)
            value = integrate_vertex(frequencies, external, blocks, order)
            reference = _quartic_reference(frequencies, external, direction, 0.7)
            bubble.append({"frequencies": list(frequencies), "direction": list(direction),
                           "external_frequencies": list(external), "frequency_sum": fsum(external),
                           **compare(value, reference)})

    frequencies = (0.8, 1.4)
    two_point = []
    for external_frequency in (0.0, 0.6):
        value, reference, tadpole, pair = _two_point_case(frequencies, external_frequency, order)
        two_point.append({"frequencies": list(frequencies),
                          "external_frequencies": [external_frequency, -external_frequency],
                          "analytic_tadpole": _json_complex(complex(tadpole)),
                          "analytic_bubble": _json_complex(complex(pair)),
                          **compare(value, reference)})
    external = external_cases[0]
    blocks = _mixed_blocks(4)
    baseline = integrate_vertex(frequencies, external, blocks, order)
    bose = []
    for permutation in permutations(range(4)):
        permuted_external, permuted_blocks = _permuted_inputs(external, blocks, permutation)
        value = integrate_vertex(frequencies, permuted_external, permuted_blocks, order)
        bose.append({"permutation": list(permutation), **compare(value, baseline)})

    rotations = []
    for angles in ((0.37, -0.61), (1.5707963267948966, 0.29)):
        transform = _canonical_rotation(frequencies, angles)
        for n, ext in ((2, (0.3, -0.3)), (3, (0.2, 0.35, -0.55)), (4, external)):
            original = _mixed_blocks(n)
            transformed = {subset: _congruence(block, transform) for subset, block in original.items()}
            reference = integrate_vertex(frequencies, ext, original, order)
            value = integrate_vertex(frequencies, ext, transformed, order)
            rotations.append({"n": n, "angles": list(angles), "external_frequencies": list(ext),
                              **_rotation_residuals(frequencies, transform), **compare(value, reference)})

    # Also compare a Q/P-mixed oscillator directly to the independent static
    # square-root answer; invariance alone would not fix the normalization.
    scalar_rotations = []
    transform = _canonical_rotation((float(omega),), (0.71,))
    for n in (2, 3, 4):
        transformed = {subset: _congruence(block, transform)
                       for subset, block in _scalar_blocks(n, couplings).items()}
        value = integrate_vertex((omega,), [0] * n, transformed, order)
        scalar_rotations.append({"n": n, **compare(value, float(_static_reference(omega, couplings, n)))})

    # The exact-zero optimization must agree with omission, not a threshold.
    sparse = {subset: block for subset, block in _scalar_blocks(4, {4: F(7, 10)}).items()
              if len(subset) == 2}
    with_zeros = _scalar_blocks(4, {4: F(7, 10)})
    zeros = compare(integrate_vertex((omega,), (0, 0, 0, 0), with_zeros, order),
                    integrate_vertex((omega,), (0, 0, 0, 0), sparse, order))
    rows = scalar + bubble + two_point + bose + rotations + scalar_rotations + [zeros]
    return {
        "quadrature_order": order, "scalar_static": scalar, "quartic_bubbles": bubble,
        "two_point_multimode": two_point,
        "four_leg_permutations": bose, "canonical_rotations": rotations,
        "rotated_scalar_static": scalar_rotations, "exact_zero_block_omission": zeros,
        "passed": all(row["passed"] for row in rows),
        "max_absolute_error": max(row["absolute_error"] for row in rows),
        "max_relative_error_nonzero_references": max(
            row["relative_error"] for row in rows if row["relative_error"] is not None),
        "mixed_four_point": _json_complex(baseline),
    }


def run_controls(config):
    """Return JSON-serializable controls; no output files or launches.

    Config keys: loop_quadrature_orders (default [32,64]),
    loop_relative_tolerance (default 1e-9), loop_absolute_tolerance (1e-12).
    Orders are sorted and deduplicated. All orders are reported; final success
    requires the highest-order controls, and (when available) agreement of
    independent test values at the two highest orders. Lower-order failures
    remain visible and are not hidden by the final status.
    """
    orders = config.get("loop_quadrature_orders", [32, 64])
    if not orders or any(isinstance(order, bool) or not isinstance(order, int) or order < 2 for order in orders):
        raise ValueError("loop_quadrature_orders must contain integers at least 2")
    orders = sorted(set(orders))
    relative_tolerance = float(config.get("loop_relative_tolerance", 1e-9))
    absolute_tolerance = float(config.get("loop_absolute_tolerance", 1e-12))
    if any(not isfinite(value) or value < 0 for value in (relative_tolerance, absolute_tolerance)):
        raise ValueError("loop tolerances must be finite and nonnegative")
    counts = {str(n): [sum(len(partition) == k for partition in ordered_partitions(n))
                       for k in range(1, n + 1)] for n in (2, 3, 4)}
    counts_passed = counts == {"2": [1, 2], "3": [1, 6, 6], "4": [1, 14, 36, 24]}
    results = [_checks_at_order(order, relative_tolerance, absolute_tolerance) for order in orders]

    convergence = []
    if len(results) > 1:
        low, high = results[-2:]
        for group in ("scalar_static", "quartic_bubbles", "two_point_multimode", "canonical_rotations", "rotated_scalar_static"):
            for index, (left, right) in enumerate(zip(low[group], high[group])):
                a = complex(left["value"]["real"], left["value"]["imaginary"])
                b = complex(right["value"]["real"], right["value"]["imaginary"])
                convergence.append({"group": group, "index": index,
                                    "lower_order": orders[-2], "higher_order": orders[-1],
                                    **_comparison(a, b, relative_tolerance, absolute_tolerance)})
        left, right = low["mixed_four_point"], high["mixed_four_point"]
        convergence.append({"group": "mixed_four_point", "index": 0,
                            "lower_order": orders[-2], "higher_order": orders[-1],
                            **_comparison(complex(left["real"], left["imaginary"]),
                                          complex(right["real"], right["imaginary"]),
                                          relative_tolerance, absolute_tolerance)})
    passed = counts_passed and results[-1]["passed"] and all(row["passed"] for row in convergence)
    return {
        "scope": "Finite-mode oscillator bookkeeping controls only; NOT a finite-N graviton amplitude.",
        "hbar": 1, "canonical_order": "Q0,P0,Q1,P1,...",
        "quadrature_orders": orders, "relative_tolerance": relative_tolerance,
        "absolute_tolerance": absolute_tolerance,
        "ordered_partition_counts": counts, "ordered_partition_counts_passed": counts_passed,
        "quadrature_results": results, "convergence": convergence,
        "convergence_checked": len(results) > 1,
        "passed": passed,
        "status": "all specified bookkeeping controls passed" if passed else "bookkeeping control failure",
        "assumptions": [
            "Midpoint/Weyl phase-space perturbation theory, positive free frequencies, and no added counterterms.",
            "H_m=g_m*q^m/m!; B_S are derivatives, so sqrt(V'') Taylor coefficients are multiplied by n!.",
            "Static references are exact rational local Taylor expansions, not claims of global polynomial stability.",
            "Bubble references use the analytic whole-line two-propagator integral, including unequal internal frequencies.",
            "Canonical rotations preserve the free Hamiltonian and symplectic form; B transforms as T^transpose B T with external directions transformed consistently.",
            "Permutation checks relabel both external frequencies and every Hessian derivative block.",
            "Two quadrature orders provide a convergence diagnostic, not a certified integration error bound.",
        ],
        "limitations": [
            "No gravitational H3-H6 coefficients are generated or validated by these oscillator controls.",
            "No exchange graphs, external-leg normalization, Lorentzian continuation or physical scattering limit are supplied.",
            "No claim of full gravitational consistency or interacting reflection positivity follows.",
        ],
    }
