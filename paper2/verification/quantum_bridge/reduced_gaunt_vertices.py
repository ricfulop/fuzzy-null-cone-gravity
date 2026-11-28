"""Exact H2, H3 and H4 of the declared reduction-first canonical candidate.

Q and S are flat TT tensors, each a polynomial times exp(-r^2/2).
S = sum_alpha (G0^-1 p)^alpha e_alpha, not the canonical coefficient p.
Products and Gaussian moments retain every generated harmonic. No binary
angular truncation is used. See reports/quantum-bridge/reduced-vertices.md.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import prod

from reduced_gaunt_exact import (
    ONE, add, scale, multiply, gaussian_derivative, gaussian_multi,
    odd_double_factorial, tt_tensor, symmetric_curl, assert_tt,
)


def psum(terms):
    out = {}
    for term in terms:
        out = add(out, term)
    return out


def zero_matrix():
    return [[{} for _ in range(3)] for _ in range(3)]


def matrix_add(a, b, factor=F(1)):
    return [[add(a[i][j], b[i][j], factor) for j in range(3)] for i in range(3)]


def matrix_scale(a, factor):
    return [[scale(x, factor) for x in row] for row in a]


def matrix_multiply(a, b):
    return [[psum(multiply(a[i][k], b[k][j]) for k in range(3))
             for j in range(3)] for i in range(3)]


def trace(a):
    return psum(a[i][i] for i in range(3))


def contract(a, b):
    return psum(multiply(a[i][j], b[i][j]) for i in range(3) for j in range(3))


@lru_cache(None)
def moment(powers, order):
    """Normalized moment with weight exp(-order*r^2/2)."""
    if any(n % 2 for n in powers):
        return F(0)
    return F(prod(odd_double_factorial(n - 1) for n in powers),
             order ** (sum(powers) // 2))


def integrate(poly, order):
    """Integral divided by (2*pi/order)^(3/2)."""
    return sum((v * moment(p, order) for p, v in poly.items()), F(0))


def pg(poly, axis, order):
    return gaussian_derivative(poly, axis, F(order))


def tensor_from_labels(electric, magnetic):
    result = matrix_add(tt_tensor(electric), symmetric_curl(tt_tensor(magnetic)))
    assert_tt(result)
    return result


def christoffel_series(q):
    """Coefficients Gamma[1:4] of the ordinary metric exp(lambda Q)."""
    q2 = matrix_multiply(q, q)
    q3 = matrix_multiply(q2, q)
    eye = [[ONE if i == j else {} for j in range(3)] for i in range(3)]
    metric = [eye, q, matrix_scale(q2, F(1, 2)), matrix_scale(q3, F(1, 6))]
    inverse = [matrix_scale(g, F((-1) ** n)) for n, g in enumerate(metric)]
    first_kind = [None]
    for n in range(1, 4):
        first_kind.append([[[scale(psum((pg(metric[n][j][l], i, n),
            pg(metric[n][i][l], j, n), scale(pg(metric[n][i][j], l, n), -1))), F(1, 2))
            for j in range(3)] for i in range(3)] for l in range(3)])
    gamma = [None]
    for n in range(1, 4):
        gamma.append([[[psum(multiply(inverse[a][k][l], first_kind[n-a][l][i][j])
                            for a in range(n) for l in range(3))
                       for j in range(3)] for i in range(3)] for k in range(3)])
        for i in range(3):
            assert not psum(gamma[n][k][i][k] for k in range(3)), (n, i)
    return metric, gamma


def gamma_pair(a, b):
    return [[psum(multiply(a[k][i][l], b[l][j][k])
                  for k in range(3) for l in range(3))
             for j in range(3)] for i in range(3)]


def curvature_data(q):
    metric, gamma = christoffel_series(q)
    p11, p12 = gamma_pair(gamma[1], gamma[1]), gamma_pair(gamma[1], gamma[2])
    l2 = trace(p11)
    l3 = add(scale(trace(p12), 2), contract(q, p11), -1)
    l4 = psum((scale(trace(gamma_pair(gamma[1], gamma[3])), 2),
               trace(gamma_pair(gamma[2], gamma[2])),
               scale(contract(q, p12), -2), contract(metric[2], p11)))
    # Pointwise R2, including its total derivatives. It is NOT replaceable
    # by the integrated Gamma-Gamma density in the conformal source.
    r2 = psum(pg(gamma[2][k][i][i], k, 2) for k in range(3) for i in range(3))
    for i in range(3):
        for j in range(3):
            r2 = add(r2, multiply(q[i][j], psum(
                pg(gamma[1][k][i][j], k, 1) for k in range(3))), -1)
    r2 = add(r2, l2, -1)
    grad = psum(multiply(pg(q[i][j], a, 1), pg(q[i][j], a, 1))
                for a in range(3) for i in range(3) for j in range(3))
    assert integrate(r2, 2) == integrate(l2, 2) == -integrate(grad, 2) / 4
    return {"gamma": gamma, "metric": metric,
            "r2": r2, "l2": l2, "l3": l3, "l4": l4}


def york_source_and_raw_gram(q, s, geometry):
    """Compare direct covariant divergence with its reduced TT expression."""
    qs, sq = matrix_multiply(q, s), matrix_multiply(s, q)
    e1 = matrix_scale(matrix_add(qs, sq), F(1, 2))
    gamma1 = geometry["gamma"][1]
    direct, compact = [], []
    for i in range(3):
        value = psum(pg(e1[i][j], j, 2) for j in range(3))
        value = add(value, psum(multiply(q[j][k], pg(s[i][j], k, 1))
                                for j in range(3) for k in range(3)), -1)
        value = add(value, psum(multiply(gamma1[l][j][i], s[l][j])
                                for l in range(3) for j in range(3)), -1)
        value = add(value, psum(multiply(gamma1[l][j][j], s[i][l])
                                for l in range(3) for j in range(3)), -1)
        direct.append(value)
        compact.append(scale(psum(
            add(add(multiply(s[l][j], pg(q[j][i], l, 1)),
                    multiply(s[j][l], pg(q[j][l], i, 1)), -1),
                    multiply(q[j][l], pg(s[i][j], l, 1)), -1)
            for l in range(3) for j in range(3)), F(1, 2)))
        assert direct[-1] == compact[-1]

    comm = matrix_add(qs, sq, -1)
    comm_squared_trace = trace(matrix_multiply(comm, comm))
    q2 = matrix_scale(geometry["metric"][2], 2)
    e2 = matrix_scale(matrix_add(matrix_add(matrix_multiply(q2, s),
        matrix_multiply(qs, q)), matrix_multiply(s, q2)), F(1, 6))
    a1 = matrix_add(e1, qs, -1)
    a2 = matrix_add(matrix_add(e2, matrix_multiply(q, e1), -1),
                    matrix_scale(matrix_multiply(q2, s), F(1, 2)))
    raw1 = scale(trace(matrix_multiply(s, a1)), 2)
    raw2 = add(scale(trace(matrix_multiply(s, a2)), 2), trace(matrix_multiply(a1, a1)))
    assert not raw1
    assert raw2 == scale(comm_squared_trace, F(-1, 12))
    return compact, comm_squared_trace, raw2


def fourier(poly):
    """Real and imaginary Q in F[P exp(-r^2)]=pi^(3/2) exp(-k^2/4) Q.

    Convention exp(-ik.x): multiplying by x is i*d/dk. Odd polynomials
    therefore require the imaginary component; dropping it changes mixed
    electric/magnetic couplings.
    """
    real, imag = {}, {}
    for powers, value in poly.items():
        n = sum(powers)
        term = gaussian_multi(powers, F(1, 2))
        if n % 2:
            imag = add(imag, term, value * (-1) ** ((n - 1) // 2))
        else:
            real = add(real, term, value * (-1) ** (n // 2))
    return real, imag


@lru_cache(None)
def inverse_laplacian_moment(powers, inverse_power):
    if any(n % 2 for n in powers):
        return F(0)
    n = sum(powers) // 2
    numerator = prod(odd_double_factorial(k - 1) for k in powers)
    if inverse_power == 0:
        return F(numerator)
    if inverse_power == 1:
        return F(numerator, 2*n+1)
    if inverse_power == 2 and n >= 1:
        return F(numerator, (2*n+1)*(2*n-1))
    raise ValueError("Nonintegrable monomial or unsupported inverse power")


def fourier_inner(a, b, inverse_power=1):
    value = F(0)
    for ac, bc in zip(a, b):
        for powers, coefficient in multiply(ac, bc).items():
            value += coefficient * inverse_laplacian_moment(powers, inverse_power)
    return value


def coulomb(a, b):
    """<a,C b> / (pi^(3/2)/sqrt(2)), C=(-Delta)^-1/8."""
    return fourier_inner(fourier(a), fourier(b)) / 16


def york_norm(vector):
    """||L W||^2 / (pi^(3/2)/sqrt(2)); Delta_L W=vector."""
    vf = [fourier(component) for component in vector]
    k = [{tuple(int(i == a) for i in range(3)): F(1)} for a in range(3)]
    k2 = psum(multiply(x, x) for x in k)
    k4 = multiply(k2, k2)
    div = tuple(psum(multiply(k[a], vf[a][part]) for a in range(3)) for part in range(2))
    # Clear denominators and verify the longitudinal inverse on the entire
    # Fourier polynomial, not only at selected momenta.
    for part in range(2):
        wnum = [add(scale(multiply(k2, vf[i][part]), -1),
                    scale(multiply(k[i], div[part]), F(1, 4))) for i in range(3)]
        kd = psum(multiply(k[i], wnum[i]) for i in range(3))
        for i in range(3):
            residual = psum((scale(multiply(k2, wnum[i]), -1),
                             scale(multiply(k[i], kd), F(-1, 3)),
                             scale(multiply(k4, vf[i][part]), -1)))
            assert not residual
    result = sum((fourier_inner(x, x) for x in vf), F(0)) - fourier_inner(div, div, 2) / 4
    assert result >= 0
    return result


def hamiltonian_coefficients(q, s, progress=None):
    geometry = curvature_data(q)
    if progress:
        progress("curvature")
    vector, comm, raw2 = york_source_and_raw_gram(q, s, geometry)
    if progress:
        progress("Gram_and_York_source")
    source = contract(s, s)
    r2 = geometry["r2"]
    c_ss, c_sr, c_rr = coulomb(source, source), coulomb(source, r2), coulomb(r2, r2)
    if progress:
        progress("scalar_Green_contractions")
    longitudinal = york_norm(vector)
    if progress:
        progress("longitudinal_Green_contraction")
    # Local fourth-order integrals have unit (pi/2)^(3/2), half the
    # common quartic unit pi^(3/2)/sqrt(2).
    local_curvature = -integrate(geometry["l4"], 4) / 2
    commutator = integrate(comm, 4) / 24
    h4q = local_curvature + c_rr
    h4mixed = commutator + longitudinal + 6*c_sr
    h4p = -7*c_ss
    coefficients = {
        "H2_configuration": -integrate(geometry["l2"], 2),
        "H2_momentum": integrate(source, 2),
        "H3_configuration": -integrate(geometry["l3"], 3),
        "H3_mixed": F(0),
        "H4_configuration": h4q,
        "H4_mixed": h4mixed,
        "H4_momentum": h4p,
        "quartic_local_curvature": local_curvature,
        "quartic_curvature_coulomb": c_rr,
        "quartic_raw_commutator": commutator,
        "quartic_longitudinal": longitudinal,
        "quartic_scalar_mixed": 6*c_sr,
    }
    if q == s:
        expected = F(3, 64) * integrate(multiply(source, source), 4)
        assert longitudinal == expected
        assert not comm
        assert all(vector[i] == scale(pg(source, i, 2), F(-1, 4)) for i in range(3))
    # Parseval provides a kernel-free normalization check for both Fourier
    # parities independently of the Coulomb denominator formula.
    assert fourier_inner(fourier(source), fourier(source), 0) == integrate(multiply(source, source), 4)
    return coefficients, {**geometry, "source": source, "vector": vector, "raw2": raw2}


def label_matrix(values):
    return [[F(x) for x in row] for row in values]


def label_add(a, b, factor=F(1)):
    return [[a[i][j]+factor*b[i][j] for j in range(3)] for i in range(3)]


def label_multiply(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(3)), F(0))
             for j in range(3)] for i in range(3)]


def label_trace(a):
    return sum((a[i][i] for i in range(3)), F(0))


def invariants(e, b):
    ee, bb, eb = label_multiply(e, e), label_multiply(b, b), label_multiply(e, b)
    te, tb, teb = label_trace(ee), label_trace(bb), label_trace(eb)
    return {
        "cubic": [label_trace(label_multiply(ee, e)), label_trace(label_multiply(e, bb))],
        "quartic": [te*te, tb*tb, te*tb, teb*teb, label_trace(label_multiply(ee, bb))],
    }


def solve_rational(matrix, rhs):
    n = len(rhs)
    a = [[F(x) for x in row] + [F(y)] for row, y in zip(matrix, rhs)]
    assert all(len(row) == n+1 for row in a)
    for i in range(n):
        pivot = next(k for k in range(i, n) if a[k][i])
        a[i], a[pivot] = a[pivot], a[i]
        divisor = a[i][i]
        a[i] = [x/divisor for x in a[i]]
        for k in range(n):
            if k != i:
                factor = a[k][i]
                a[k] = [x-factor*y for x, y in zip(a[k], a[i])]
    return [row[-1] for row in a]


def test_cases():
    z = label_matrix([[0]*3 for _ in range(3)])
    a = label_matrix([[F(-1, 2), 0, 0], [0, F(-1, 2), 0], [0, 0, 1]])
    d = label_matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]])
    xy = label_matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    xz = label_matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    yz = label_matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    generic = (label_add(d, xy, F(1, 2)), label_add(a, xz, F(1, 3)),
               label_add(a, d, F(1, 4)), [[x/2 for x in row] for row in yz])
    rot = label_matrix([[F(1, 3), F(-2, 3), F(2, 3)],
                        [F(2, 3), F(2, 3), F(1, 3)],
                        [F(-2, 3), F(1, 3), F(2, 3)]])
    rt = [list(row) for row in zip(*rot)]
    assert label_multiply(rot, rt) == label_matrix([[int(i == j) for j in range(3)] for i in range(3)])
    rotated = tuple(label_multiply(label_multiply(rot, x), rt) for x in generic)
    return [
        ("electric_axis", (a, z, a, z)),
        ("magnetic_axis", (z, a, z, a)),
        ("parallel_parities", (a, a, a, [[-x for x in row] for row in a])),
        ("orthogonal_labels", (a, d, d, z)),
        ("noncommuting_labels", (a, xz, xz, a)),
        ("held_out_mixed", generic),
        ("held_out_rotated", rotated),
    ]
