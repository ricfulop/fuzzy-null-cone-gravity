"""Exact local Gaussian-polynomial ingredients for gravitational H5/H6.

Inputs q and s are 3-by-3 symmetric flat-TT tensors of sparse polynomial
dicts with rational coefficients. The physical tensors are q*exp(-r^2/2)
and s*exp(-r^2/2); s is the dual-Gram momentum field, not canonical p.
No spatial projection or approximate coefficient pruning is performed.

All returned polynomials omit their Gaussian envelope:
  metric[n], inverse[n], gamma[n], R[n], L[n]: exp(-n*r^2/2);
  E[n], b[n]: exp(-(n+1)*r^2/2), with n counting q powers at fixed s.
These are ordinary homogeneous coefficients, not derivative tensors with
extra factorials. Gamma indexing is gamma[n][upper][lower][lower].
R is pointwise scalar curvature. L is a Gamma--Gamma *density* whose
spatial integral equals that of R; L cannot replace R in scalar sources.

Public entry points:
  curvature_data(q, local_R_order=4, density_order=6, progress=None)
  york_data(q, s, max_order=3, geometry=None, progress=None)
  local_geometry(q, s=None, local_R_order=4, density_order=6,
                 york_order=3, progress=None)
  compare_legacy(q, s=None, geometry=None)
The combined entry point shares the connection. Standalone tangent_series
does not construct curvature or Christoffels. Requested orders are respected:
curvature_data(q, 2, 2) constructs only gamma[1:3], R1/R2, and L2.
Set either curvature order to zero to omit that family. No function runs on
import; comparisons must be invoked by the tracked calculation runner.
Progress callbacks receive one dict with stage, order, and maximum_order.

This module does not solve York/scalar inverses or evaluate the full reduced
Hamiltonian. Archived third/fourth-order routines are left unchanged.
"""

from fractions import Fraction as F
from math import factorial

from reduced_gaunt_exact import ONE, add, multiply, scale, assert_tt
from reduced_gaunt_vertices import (
    psum, zero_matrix, matrix_multiply, matrix_scale, contract, pg,
    gamma_pair,
)


def _order(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(name + " must be a nonnegative integer")
    return value


def _check_tensor(tensor, name):
    if len(tensor) != 3 or any(len(row) != 3 for row in tensor):
        raise ValueError(name + " must be a 3-by-3 polynomial tensor")
    for row in tensor:
        for polynomial in row:
            if not isinstance(polynomial, dict):
                raise TypeError(name + " entries must be sparse polynomial dicts")
            for powers, value in polynomial.items():
                if (not isinstance(powers, tuple) or len(powers) != 3
                        or any(type(p) is not int or p < 0 for p in powers)):
                    raise ValueError("polynomial powers must be triples of nonnegative integers")
                if not isinstance(value, (int, F)):
                    raise TypeError("exact polynomial coefficients must be int or Fraction")
    assert_tt(tensor)


def _progress(callback, stage, order, maximum):
    if callback is not None:
        callback({"stage": stage, "order": order, "maximum_order": maximum})


def _identity():
    return [[dict(ONE) if i == j else {} for j in range(3)] for i in range(3)]


def _connection(q, max_order, progress):
    metric = [_identity()]
    inverse = [_identity()]
    first_kind, gamma = [None], [None]
    for n in range(1, max_order + 1):
        # g_n = g_(n-1) Q/n, with its full order-n Gaussian envelope.
        metric.append(matrix_scale(matrix_multiply(metric[-1], q), F(1, n)))
        inverse.append(matrix_scale(metric[n], F((-1) ** n)))
        derivatives = [[[pg(metric[n][i][j], axis, n) for j in range(3)]
                        for i in range(3)] for axis in range(3)]
        lowered = [zero_matrix() for _ in range(3)]
        for ell in range(3):
            for i in range(3):
                for j in range(i, 3):
                    value = scale(psum((derivatives[i][j][ell],
                                        derivatives[j][i][ell],
                                        scale(derivatives[ell][i][j], -1))), F(1, 2))
                    lowered[ell][i][j] = lowered[ell][j][i] = value
        first_kind.append(lowered)
        raised = [zero_matrix() for _ in range(3)]
        for k in range(3):
            for i in range(3):
                for j in range(i, 3):
                    value = psum(multiply(inverse[a][k][ell], first_kind[n - a][ell][i][j])
                                 for a in range(n) for ell in range(3))
                    raised[k][i][j] = raised[k][j][i] = value
        gamma.append(raised)
        # det(exp(Q))=1 implies Gamma^k_(ik)=0 pointwise, order by order.
        for i in range(3):
            if psum(raised[k][i][k] for k in range(3)):
                raise AssertionError(("nonzero Christoffel trace", n, i))
        _progress(progress, "christoffel", n, max_order)
    return {
        "q": q, "metric": metric, "inverse": inverse, "gamma": gamma,
        "connection_order": max_order,
        "checks": {"christoffel_trace_zero": True},
        "gaussian_rates": {
            "metric": {n: n for n in range(max_order + 1)},
            "inverse": {n: n for n in range(max_order + 1)},
            "gamma": {n: n for n in range(1, max_order + 1)},
        },
    }


def christoffel_series(q, max_order=5, progress=None):
    """Return g_n, g^-1_n and Gamma_n through max_order, with exact trace checks.

    Unlike the archived function of the same name, this returns a dictionary.
    The zero connection is represented by gamma[0]=None; metric[0]=I.
    Treat returned coefficient dictionaries as immutable (symmetric entries
    may share storage). Input q must be flat TT, as required by this model.
    """
    _order(max_order, "max_order")
    _check_tensor(q, "q")
    return _connection(q, max_order, progress)


def _curvature(geometry, local_R_order, density_order, progress):
    inverse, gamma = geometry["inverse"], geometry["gamma"]
    pair_cache, density_cache = {}, {}

    def density(n):
        if n not in density_cache:
            terms = []
            for a in range(n - 1):
                for b in range(1, (n - a) // 2 + 1):
                    c = n - a - b
                    if (b, c) not in pair_cache:
                        pair_cache[b, c] = gamma_pair(gamma[b], gamma[c])
                    # Pair(c,b) is Pair(b,c)^T. The inverse metric is
                    # symmetric, so their contracted densities are equal.
                    factor = 1 if b == c else 2
                    terms.append(scale(contract(inverse[a], pair_cache[b, c]), factor))
            density_cache[n] = psum(terms)
        return density_cache[n]

    ricci_derivative = {}
    pointwise = {}
    for n in range(1, local_R_order + 1):
        ricci_derivative[n] = [[psum(pg(gamma[n][k][i][j], k, n) for k in range(3))
                                for j in range(3)] for i in range(3)]
        # Unimodular R = g^ij (partial_k Gamma^k_ij - Gamma^k_il Gamma^l_jk).
        derivative_term = psum(contract(inverse[a], ricci_derivative[n - a])
                               for a in range(n))
        pointwise[n] = add(derivative_term, density(n), -1)
        if n == 1 and pointwise[n]:
            raise AssertionError("pointwise R1 must vanish for a flat TT q")
        _progress(progress, "pointwise_curvature", n, local_R_order)
    densities = {}
    for n in range(2, density_order + 1):
        densities[n] = density(n)
        _progress(progress, "integrated_curvature_density", n, density_order)
    geometry.update({"R": pointwise, "L": densities,
                     "local_R_order": local_R_order, "density_order": density_order})
    geometry["checks"]["pointwise_R1_zero"] = not pointwise[1] if local_R_order else None
    geometry["gaussian_rates"].update({"R": {n: n for n in pointwise},
                                       "L": {n: n for n in densities}})
    # Low-order aliases make comparison with archived results straightforward.
    geometry.update({"r" + str(n): value for n, value in pointwise.items()})
    geometry.update({"l" + str(n): value for n, value in densities.items()})
    return geometry


def curvature_data(q, local_R_order=4, density_order=6, progress=None):
    """Exact pointwise R1..R_local_R_order and integrated densities L2..L_density_order.

    Density order d needs Gamma only through d-1; pointwise order r needs
    Gamma through r. In particular defaults never construct pointwise R5/R6.
    Integral normalization: integrate(L[n], n) from reduced_gaunt_vertices
    returns the spatial integral divided by (2*pi/n)^(3/2), exactly as before.
    """
    _order(local_R_order, "local_R_order")
    _order(density_order, "density_order")
    _check_tensor(q, "q")
    connection_order = max(local_R_order, density_order - 1, 0)
    geometry = _connection(q, connection_order, progress)
    return _curvature(geometry, local_R_order, density_order, progress)


def _compatible(q, geometry, max_order):
    if geometry.get("q") != q:
        raise ValueError("supplied geometry belongs to a different q")
    if geometry.get("connection_order", -1) < max_order:
        raise ValueError("supplied geometry has insufficient connection order")


def _tangent(q, s, max_order, metric=None):
    if metric is not None:
        powers = [matrix_scale(metric[n], F(factorial(n))) for n in range(max_order + 1)]
    else:
        powers = [_identity()]
        for _ in range(max_order):
            powers.append(matrix_multiply(powers[-1], q))
    # Cache Q^a S; each ordered placement in the Frechet derivative remains.
    left = [matrix_multiply(power, s) for power in powers]
    result = []
    for n in range(max_order + 1):
        terms = [matrix_multiply(left[a], powers[n - a]) for a in range(n + 1)]
        coefficient = F(1, factorial(n + 1))
        result.append([[scale(psum(term[i][j] for term in terms), coefficient)
                        for j in range(3)] for i in range(3)])
    return result


def tangent_series(q, s, max_order=3):
    """Return E[n]=sum_(a=0)^n Q^a S Q^(n-a)/(n+1)! without computing geometry.

    E[n] has Gaussian rate n+1, although its index counts only powers of q.
    """
    _order(max_order, "max_order")
    _check_tensor(q, "q")
    _check_tensor(s, "s")
    return _tangent(q, s, max_order)


def _york(s, geometry, max_order, progress):
    inverse, gamma = geometry["inverse"], geometry["gamma"]
    tangent = _tangent(geometry["q"], s, max_order, geometry["metric"])
    sources = []
    covariant_derivative = []
    for n in range(max_order + 1):
        # coefficient n of nabla_k E_ij, with n powers of q and one s.
        derivative = []
        for k in range(3):
            matrix = zero_matrix()
            for i in range(3):
                for j in range(i, 3):
                    connection_term = psum(
                        add(multiply(gamma[b][ell][k][i], tangent[n - b][ell][j]),
                            multiply(gamma[b][ell][k][j], tangent[n - b][i][ell]))
                        for b in range(1, n + 1) for ell in range(3))
                    value = add(pg(tangent[n][i][j], k, n + 1), connection_term, -1)
                    matrix[i][j] = matrix[j][i] = value
            derivative.append(matrix)
        covariant_derivative.append(derivative)
        sources.append([psum(multiply(inverse[a][j][k], covariant_derivative[n - a][k][i][j])
                             for a in range(n + 1) for j in range(3) for k in range(3))
                        for i in range(3)])
        if n == 0 and any(sources[n]):
            raise AssertionError("b0 must vanish for a flat TT s")
        _progress(progress, "york_source", n, max_order)
    return {
        "E": tangent, "b": sources, "york_order": max_order,
        "gaussian_rates": {"E": {n: n + 1 for n in range(max_order + 1)},
                           "b": {n: n + 1 for n in range(max_order + 1)}},
        "checks": {"flat_divergence_b0_zero": True},
    }


def york_data(q, s, max_order=3, geometry=None, progress=None):
    """Direct curved divergence b_n = [g^jk nabla_k E_ij]_n, through q^n s.

    b[n] is a three-component covector, not a vector with a silently raised
    index. It has Gaussian rate n+1. Supply a connection/curvature dictionary
    from this module to share the expensive geometry; it is not modified.
    No York inverse or b-independent scalar curvature is constructed here.
    """
    _order(max_order, "max_order")
    _check_tensor(q, "q")
    _check_tensor(s, "s")
    if geometry is None:
        geometry = _connection(q, max_order, progress)
    else:
        _compatible(q, geometry, max_order)
    return _york(s, geometry, max_order, progress)


def local_geometry(q, s=None, local_R_order=4, density_order=6,
                   york_order=3, progress=None):
    """Shared-connection wrapper returning curvature and, if supplied, s-dependent data.

    R/L dictionaries and E/b lists are keyed/indexed by homogeneous order.
    The gaussian_rates dictionary records the envelope rate for every family.
    """
    for name, value in (("local_R_order", local_R_order), ("density_order", density_order),
                        ("york_order", york_order)):
        _order(value, name)
    _check_tensor(q, "q")
    if s is not None:
        _check_tensor(s, "s")
    maximum = max(local_R_order, density_order - 1, york_order if s is not None else 0, 0)
    geometry = _connection(q, maximum, progress)
    _curvature(geometry, local_R_order, density_order, progress)
    if s is not None:
        york = _york(s, geometry, york_order, progress)
        geometry.update({key: york[key] for key in ("E", "b", "york_order")})
        geometry["gaussian_rates"].update(york["gaussian_rates"])
        geometry["checks"].update(york["checks"])
    return geometry


def compare_legacy(q, s=None, geometry=None):
    """Explicit tracked-run comparison against archived r2, l2/l3/l4 and optional b1.

    With no supplied result, request ONLY R through 2, L through 4, and b
    through 1. No fifth/sixth-order calculation is needed for this check.
    Returns JSON-serializable exact equality checks and residual term counts;
    this comparison is never performed automatically at import or construction.
    """
    from reduced_gaunt_vertices import curvature_data as archived_curvature
    from reduced_gaunt_vertices import york_source_and_raw_gram

    if geometry is None:
        geometry = local_geometry(q, s, local_R_order=2, density_order=4, york_order=1)
    _compatible(q, geometry, 3)
    if geometry.get("local_R_order", 0) < 2 or geometry.get("density_order", 0) < 4:
        raise ValueError("legacy comparison needs R2 and L2 through L4")
    archived = archived_curvature(q)
    residuals = {name: len(add(geometry[name], archived[name], -1))
                 for name in ("r2", "l2", "l3", "l4")}
    if s is not None:
        # Recompute b1 from the supplied q geometry and this s, so a reused
        # combined dictionary cannot silently compare a different momentum.
        current = york_data(q, s, max_order=1, geometry=geometry)["b"][1]
        previous, _, _ = york_source_and_raw_gram(q, s, archived)
        residuals.update({"b1_" + str(i): len(add(current[i], previous[i], -1))
                          for i in range(3)})
    checks = {name: count == 0 for name, count in residuals.items()}
    return {"passed": all(checks.values()), "checks": checks,
            "residual_monomials": residuals}
