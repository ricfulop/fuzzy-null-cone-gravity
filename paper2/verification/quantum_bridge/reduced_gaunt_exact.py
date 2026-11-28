"""Exact polynomial algebra for a declared reduction-before-cutoff candidate.

Ten smooth l=2 TT profiles (five electric and five magnetic) are a selected
finite radial/angular basis, not all radial degrees of freedom at N=3.
Polynomials multiply exp(-r^2/2); their products integrate with exp(-r^2).
All angular products are complete pointwise/Gaunt contractions.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import factorial, prod

from run import inner, legendre

ONE = {(0, 0, 0): F(1)}


def add(a, b, factor=F(1)):
    result = dict(a)
    for powers, value in b.items():
        result[powers] = result.get(powers, F(0)) + factor * value
    return {powers: value for powers, value in result.items() if value}


def scale(a, factor):
    return {powers: factor * value for powers, value in a.items() if factor * value}


def multiply(a, b):
    result = {}
    for pa, va in a.items():
        for pb, vb in b.items():
            powers = tuple(x + y for x, y in zip(pa, pb))
            result[powers] = result.get(powers, F(0)) + va * vb
    return {powers: value for powers, value in result.items() if value}


def derivative(poly, axis):
    result = {}
    for powers, value in poly.items():
        if powers[axis]:
            new = list(powers)
            new[axis] -= 1
            result[tuple(new)] = powers[axis] * value
    return result


def gaussian_derivative(poly, axis, rate=F(1)):
    powers = [0, 0, 0]
    powers[axis] = 1
    return add(derivative(poly, axis), multiply({tuple(powers): -rate}, poly))


@lru_cache(None)
def gaussian_multi(powers, rate=F(1)):
    if sum(powers) == 0:
        return ONE
    axis = next(i for i, power in enumerate(powers) if power)
    prior = list(powers)
    prior[axis] -= 1
    return gaussian_derivative(gaussian_multi(tuple(prior), rate), axis, rate)


def multi(indices):
    return tuple(indices.count(i) for i in range(3))


def tt_tensor(polarization):
    """(P_i^k P_j^l - P_ij P^kl/2) E_kl exp(-r^2/2)."""
    laplacian_squared = {}
    laplacian_hessian = {}
    for i in range(3):
        for j in range(3):
            term = {}
            for a in range(3):
                term = add(term, gaussian_multi(multi([i, j, a, a])))
            laplacian_hessian[i, j] = term
    for a in range(3):
        laplacian_squared = add(laplacian_squared, laplacian_hessian[a, a])
    laplacian_contraction = {}
    for k in range(3):
        for l in range(3):
            laplacian_contraction = add(laplacian_contraction,
                laplacian_hessian[k, l], polarization[k][l])
    result = []
    for i in range(3):
        row = []
        for j in range(3):
            value = scale(laplacian_squared, polarization[i][j])
            for k in range(3):
                value = add(value, laplacian_hessian[i, k], -polarization[k][j])
                value = add(value, laplacian_hessian[j, k], -polarization[k][i])
                for l in range(3):
                    value = add(value, gaussian_multi(multi([i, j, k, l])),
                                polarization[k][l] / 2)
            if i == j:
                value = add(value, laplacian_contraction, F(1, 2))
            row.append(value)
        result.append(row)
    return result


def levi_civita(i, j, k):
    if len({i, j, k}) < 3:
        return 0
    return 1 if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1


def symmetric_curl(tensor):
    result = []
    for i in range(3):
        row = []
        for j in range(3):
            value = {}
            for a in range(3):
                for k in range(3):
                    value = add(value, gaussian_derivative(tensor[k][j], a),
                                F(levi_civita(i, a, k), 2))
                    value = add(value, gaussian_derivative(tensor[i][k], a),
                                F(levi_civita(j, a, k), 2))
            row.append(value)
        result.append(row)
    return result


def assert_tt(tensor):
    trace = {}
    for i in range(3):
        trace = add(trace, tensor[i][i])
        divergence = {}
        for j in range(3):
            assert tensor[i][j] == tensor[j][i]
            divergence = add(divergence, gaussian_derivative(tensor[i][j], j))
        assert not divergence, (i, divergence)
    assert not trace


def odd_double_factorial(n):
    return prod(range(1, n + 1, 2))


@lru_cache(None)
def gaussian_moment(powers):
    """Integral of x**powers exp(-r^2) divided by pi**(3/2)."""
    if any(power % 2 for power in powers):
        return F(0)
    return F(prod(odd_double_factorial(power - 1) for power in powers),
             2 ** (sum(powers) // 2))


def gaussian_inner(a, b):
    total = F(0)
    for pa, va in a.items():
        for pb, vb in b.items():
            powers = tuple(x + y for x, y in zip(pa, pb))
            total += va * vb * gaussian_moment(powers)
    return total


def tensor_inner(a, b):
    return sum((gaussian_inner(a[i][j], b[i][j])
                for i in range(3) for j in range(3)), F(0))


def tensor_gradient(tensor, axis):
    return [[gaussian_derivative(tensor[i][j], axis) for j in range(3)]
            for i in range(3)]


def five_polarizations():
    entries = [
        [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
        [[F(-1, 2), 0, 0], [0, F(-1, 2), 0], [0, 0, 1]],
        [[0, 1, 0], [1, 0, 0], [0, 0, 0]],
        [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
        [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
    ]
    return [[[F(x) for x in row] for row in matrix] for matrix in entries]


def exact_mode_basis(progress=None):
    electric = [tt_tensor(e) for e in five_polarizations()]
    modes = electric + [symmetric_curl(t) for t in electric]
    gradients = []
    for tensor in modes:
        assert_tt(tensor)
        gradients.append([tensor_gradient(tensor, axis) for axis in range(3)])
    gram, stiffness = [], []
    for i, a in enumerate(modes):
        gram.append([tensor_inner(a, b) for b in modes])
        stiffness.append([
            sum((tensor_inner(gradients[i][axis], gradients[j][axis])
                 for axis in range(3)), F(0)) for j in range(len(modes))
        ])
        if progress:
            progress(i + 1)
    for i in range(len(modes)):
        assert gram[i][i] > 0
        expected = F(11, 2) if i < 5 else F(13, 2)
        assert stiffness[i][i] / gram[i][i] == expected
        for j in range(len(modes)):
            if i != j:
                assert gram[i][j] == stiffness[i][j] == 0
    assert all(gram[i + 5][i + 5] == F(11, 2) * gram[i][i] for i in range(5))
    assert gram[1][1] == F(567, 16)
    # Exact dual-Gram canonical pairing E_alpha . sigma_beta = delta_ab.
    pairing = [[gram[i][j] / gram[j][j] for j in range(len(modes))]
               for i in range(len(modes))]
    assert all(value == int(i == j)
               for i, row in enumerate(pairing) for j, value in enumerate(row))
    return modes, {
        "angular_rank_N": 3, "maximum_free_total_angular_momentum": 2,
        "electric_profiles": 5, "magnetic_profiles": 5,
        "radial_scope": "one selected smooth profile in each parity sector; not a complete radial basis",
        "gram_divided_by_pi_to_3_over_2": [[str(x) for x in row] for row in gram],
        "stiffness_divided_by_pi_to_3_over_2": [[str(x) for x in row] for row in stiffness],
        "frequency_squared": [str(stiffness[i][i] / gram[i][i]) for i in range(len(modes))],
        "trace_divergence_and_dual_pairing_residual": "0",
        "units": "packet width a=1 and 16*pi*G=1",
    }


def tensor_square(tensor):
    result = {}
    for row in tensor:
        for component in row:
            result = add(result, multiply(component, component))
    return result


def axisymmetric_source(tensor):
    source = tensor_square(tensor)
    rotation = add(multiply({(0, 1, 0): F(-1)}, derivative(source, 0)),
                   multiply({(1, 0, 0): F(1)}, derivative(source, 1)))
    assert not rotation
    radial_angular = {}
    for (a, b, c), value in source.items():
        assert a % 2 == b % 2 == c % 2 == 0
        if b:
            continue  # Evaluate the invariant scalar on the y=0 meridian.
        radial_power = a + c
        polynomial = radial_angular.setdefault(radial_power, [F(0)] * 9)
        for j in range(a // 2 + 1):
            coefficient = F(factorial(a // 2),
                            factorial(j) * factorial(a // 2 - j))
            polynomial[2 * j + c] += value * (-1) ** j * coefficient
    gaunt = {}
    for ell in range(9):
        radial = {power: (2 * ell + 1) * inner(tuple(poly), legendre(ell))
                  for power, poly in radial_angular.items()}
        radial = {power: value for power, value in radial.items() if value}
        if radial:
            gaunt[ell] = radial
    expected = {
        0: {0: F(54), 2: F(-72), 4: F(48), 6: F(-48, 5), 8: F(3, 5)},
        2: {2: F(72), 4: F(-444, 7), 6: F(96, 7), 8: F(-6, 7)},
        4: {4: F(108, 7), 6: F(-144, 35), 8: F(9, 35)},
    }
    assert gaunt == expected
    # Positive radial majorant B bounds sigma^2 and its Newton potential.
    potential_bound = sum((abs(value) * F(factorial(power // 2), 16)
                          for radial in gaunt.values() for power, value in radial.items()), F(0))
    assert potential_bound == F(2301, 56)
    norm = sum(value * gaussian_moment(powers) for powers, value in source.items())
    assert norm == F(567, 16)
    central_potential = sum((value * F(factorial(power // 2), 16)
                             for power, value in gaunt[0].items()), F(0))
    assert central_potential == F(87, 40)

    # Fourier transform: F[sigma^2] = pi^(3/2) exp(-k^2/4) Q(k).
    fourier = {}
    for powers, value in source.items():
        fourier = add(fourier, gaussian_multi(powers, F(1, 2)),
                      value * (-1) ** (sum(powers) // 2))
    assert fourier[(0, 0, 0)] == norm
    coulomb = F(0)
    for powers, value in multiply(fourier, fourier).items():
        if any(power % 2 for power in powers):
            continue
        n = sum(powers) // 2
        coulomb += value * F(prod(odd_double_factorial(power - 1) for power in powers),
                             16 * (2 * n + 1))
    assert coulomb > 0
    return gaunt, potential_bound, {
        "source_harmonics": {str(ell): {str(power): str(value)
                            for power, value in sorted(radial.items())}
                            for ell, radial in sorted(gaunt.items())},
        "source_norm_divided_by_pi_to_3_over_2": str(norm),
        "u2_at_origin": str(central_potential),
        "uniform_u2_upper_bound": str(potential_bound),
        "fourier_coulomb_integral_divided_by_pi_to_3_over_2_over_sqrt2": str(coulomb),
        "energy_quartic_coefficient_divided_by_pi_to_3_over_2_over_sqrt2": str(-7 * coulomb),
        "source_field_units": "sigma=epsilon*T; T is the electric axisymmetric l=2 profile",
    }
