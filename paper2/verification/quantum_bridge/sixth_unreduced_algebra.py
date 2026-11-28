"""Pointwise original recursion and energy for the independent H6 reference.

All inputs are physical Cartesian field values, including Gaussian envelopes.
Q and S are raw flat-TT tensors; S is the dual-Gram momentum field. Spatial
jets use derivative-first indices. Metric inverse coefficients are ordinary
Taylor coefficients: inverse1=-Q and inverse2=Q^2/2, with their physical
spatial derivatives. No factorial or canonical-coordinate factor is added.

This module performs no inversion, angular projection, quadrature, or work
on import. It uses neither the production reduced sources nor its HEAT,
Fourier, or pairing implementations. The caller supplies independently solved
W1, W2, u, u3, u4 and the physical canonical correction Ec. These functions
implement the unreduced reference report; numerical validation is pending.

Energy output has nine spatial terms. The separate global g^T G0^-1 g term
must be added once after integration. There is no separate c-vector term in
this representation. Full finite-radius exterior completion remains with
the caller, including the inhomogeneous exterior source in the u4 solve.
"""

from collections.abc import Mapping
from math import fsum, isfinite


SPATIAL_FAMILIES = (
    "raw_gram", "york_second_response", "york_third_source", "york_A2_strong",
    "curvature", "scalar_B4_u", "scalar_B3_u3", "scalar_B2_u4", "scalar_s2_u2",
)


def _real(value, name):
    if isinstance(value, bool):
        raise TypeError(name+" must be a real finite scalar")
    result = float(value)
    if not isfinite(result):
        raise ValueError(name+" must be finite")
    return result


def _array(value, shape, name):
    if not shape:
        return _real(value, name)
    if not isinstance(value, (tuple, list)) or len(value) != shape[0]:
        raise ValueError(name+" has an incorrect Cartesian shape")
    return tuple(_array(item, shape[1:], name) for item in value)


def _matrix(value, name):
    return _array(value, (3, 3), name)


def _vector(value, name):
    return _array(value, (3,), name)


def _multiply(a, b):
    return tuple(tuple(fsum(a[i][k]*b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def _combine(*terms):
    return tuple(tuple(fsum(factor*matrix[i][j] for factor, matrix in terms)
                       for j in range(3)) for i in range(3))


def _trace(a):
    return fsum(a[i][i] for i in range(3))


def _dot(a, b):
    return fsum(a[i]*b[i] for i in range(3))


def _squared_norm(a):
    return fsum(a[i][j]**2 for i in range(3) for j in range(3))


def _checked(record):
    """Reject arithmetic overflow before a caller records a finite density."""
    for name, value in record.items():
        if isinstance(value, (tuple, list)):
            for item in value:
                _checked({name: item})
        else:
            _real(value, name)
    return record


def scalar_metric_operator(inverse_coefficient, derivative, scalar_jet):
    """Return div(A grad f), with derivative[a][i][j]=partial_a A^ij.

    The derivative term is retained explicitly; no TT simplification is used.
    scalar_jet must contain physical gradient and Hessian arrays.
    """
    a = _matrix(inverse_coefficient, "inverse_coefficient")
    da = _array(derivative, (3, 3, 3), "derivative")
    gradient = _vector(scalar_jet["gradient"], "gradient")
    hessian = _matrix(scalar_jet["hessian"], "hessian")
    return _real(fsum(a[i][j]*hessian[i][j]+da[i][i][j]*gradient[j]
                      for i in range(3) for j in range(3)), "metric operator")


def third_sources(*, q, s, t0, r2, r3, b2, a1w1,
                  u_jet, inverse1, dinverse1):
    """Build s3, F3 and the actual W2 source from raw tangent products.

    T0=L0 W1 and a1w1=A1 W1. Returns r=b2-A1W1 and F3 for
    W2=Y r, u3=C F3; B2/B3 are also supplied for the energy ledger.
    """
    q, s, t0 = (_matrix(value, name) for value, name in
                 ((q, "q"), (s, "s"), (t0, "t0")))
    r2, r3 = _real(r2, "r2"), _real(r3, "r3")
    b2, a1w1 = _vector(b2, "b2"), _vector(a1w1, "a1w1")
    qs, sq = _multiply(q, s), _multiply(s, q)
    e1 = _combine((.5, qs), (.5, sq))
    z1 = _combine((1, e1), (-1, t0))
    u1 = _combine((1, z1), (-1, qs))
    s2 = _squared_norm(s)
    s3 = 2*_trace(_multiply(s, u1))
    d1u = scalar_metric_operator(inverse1, dinverse1, u_jet)
    return _checked({"s2": s2, "s3": s3, "B2": r2+7*s2,
                     "B3": r3+7*s3, "D1u": d1u,
                     "F3": s3-r3+8*d1u,
                     "r": tuple(b2[i]-a1w1[i] for i in range(3)),
                     "E1": e1, "Z1": z1, "U1": u1})


def fourth_sources(*, q, s, t0, t1, l0w2, ec, r2, r4,
                   u_jet, u3_jet, inverse1, dinverse1, inverse2, dinverse2):
    """Build the original s4 and F4, without adjoint elimination of u4.

    Ec=sum e_alpha (G0^-1 g)^alpha has its defining plus sign; Z2
    subtracts it. T1=L1 W1 is distinct from L0 W2. tr(U1^2) below
    is a matrix trace, not a squared Frobenius norm.
    """
    q, s, t0, t1, l0w2, ec = (_matrix(value, name) for value, name in
        ((q, "q"), (s, "s"), (t0, "t0"), (t1, "t1"),
         (l0w2, "l0w2"), (ec, "ec")))
    r2, r4 = _real(r2, "r2"), _real(r4, "r4")
    q2 = _multiply(q, q)
    qs, sq = _multiply(q, s), _multiply(s, q)
    e1 = _combine((.5, qs), (.5, sq))
    e2 = _combine((1/6, _multiply(q2, s)), (1/6, _multiply(qs, q)),
                  (1/6, _multiply(s, q2)))
    z1 = _combine((1, e1), (-1, t0))
    z2 = _combine((1, e2), (-1, l0w2), (-1, t1), (-1, ec))
    u1 = _combine((1, z1), (-1, qs))
    u2 = _combine((1, z2), (-1, _multiply(q, z1)), (.5, _multiply(q2, s)))
    s2 = _squared_norm(s)
    s4 = _trace(_multiply(u1, u1))+2*_trace(_multiply(s, u2))
    d1u3 = scalar_metric_operator(inverse1, dinverse1, u3_jet)
    d2u = scalar_metric_operator(inverse2, dinverse2, u_jet)
    b2 = r2+7*s2
    f4 = s4-r4+8*d1u3+8*d2u-b2*_real(u_jet["value"], "u")
    return _checked({"s4": s4, "B4": r4+7*s4, "F4": f4,
                     "D1u3": d1u3, "D2u": d2u,
                     "E2": e2, "Z2": z2, "U2": u2})


def energy_densities(*, q, s, sources, w1, w2, b3, a2w1, l6, u, u3, u4):
    """Nine signed spatial densities; add global g^T G0^-1 g separately.

    sources is the union of third_sources and fourth_sources. The returned
    terms include the original -B2*u4 and -B3*u3 without reduced bilinears.
    """
    q, s = _matrix(q, "q"), _matrix(s, "s")
    w1, w2, b3, a2w1, r = (_vector(value, name) for value, name in
        ((w1, "w1"), (w2, "w2"), (b3, "b3"), (a2w1, "a2w1"), (sources["r"], "r")))
    u, u3, u4, l6 = (_real(value, name) for value, name in
                     ((u, "u"), (u3, "u3"), (u4, "u4"), (l6, "l6")))
    commutator = _combine((1, _multiply(q, s)), (-1, _multiply(s, q)))
    double = _combine((1, _multiply(q, commutator)), (-1, _multiply(commutator, q)))
    return _checked({
        "raw_gram": -_trace(_multiply(double, double))/360,
        "york_second_response": -2*_dot(r, w2),
        "york_third_source": -4*_dot(b3, w1),
        "york_A2_strong": 2*_dot(w1, a2w1),
        "curvature": -l6,
        "scalar_B4_u": -_real(sources["B4"], "B4")*u,
        "scalar_B3_u3": -_real(sources["B3"], "B3")*u3,
        "scalar_B2_u4": -_real(sources["B2"], "B2")*u4,
        "scalar_s2_u2": 28*_real(sources["s2"], "s2")*u*u,
    })


def york_norm_flux(w, longitudinal, normal):
    """Surface density +2 W dot (T n) in the finite-ball norm identity."""
    w, t, n = _vector(w, "w"), _matrix(longitudinal, "longitudinal"), _vector(normal, "normal")
    return _real(2*fsum(w[i]*t[i][j]*n[j] for i in range(3) for j in range(3)), "York flux")


def assemble_integrated_energy(spatial_integrals, gram_square):
    """Add the separately supplied global Gram square exactly once.

    The caller must first integrate and exterior-complete every spatial
    family in physical units. This function supplies no integration or
    convergence verdict; an extra canonical/c-vector key is an error.
    """
    if not isinstance(spatial_integrals, Mapping) or set(spatial_integrals) != set(SPATIAL_FAMILIES):
        raise ValueError("Exactly the nine unreduced spatial families are required")
    ledger = {name: _real(spatial_integrals[name], name) for name in SPATIAL_FAMILIES}
    ledger["gram_square"] = _real(gram_square, "gram_square")
    return {"ledger": ledger, "total": _real(fsum(ledger.values()), "total"),
            "absolute_term_sum": _real(fsum(abs(value) for value in ledger.values()), "absolute_term_sum"),
            "global_gram_square_count": 1, "integration_and_exterior_completion_checked": False}


def a2_flux(*, q, w1, t0, t1, t2, normal):
    """Surface density for integral(2 W1 A2W1 - five weak densities).

    Multiply by the physical surface measure. No pointwise strong/weak
    equality is implied or tested by this flux function.
    """
    q, t0, t1, t2 = (_matrix(value, name) for value, name in
                     ((q, "q"), (t0, "t0"), (t1, "t1"), (t2, "t2")))
    w, n = _vector(w1, "w1"), _vector(normal, "normal")
    current = _combine((1, t2), (-1, _multiply(t1, q)),
                       (.5, _multiply(t0, _multiply(q, q))))
    return _real(2*fsum(w[i]*current[i][j]*n[j] for i in range(3) for j in range(3)), "A2 flux")


def scalar_adjoint_flux(v_jet, u4_jet, normal):
    """Surface density 8(v dn u4-u4 dn v); no radius/solid angle included."""
    n = _vector(normal, "normal")
    v, u4 = _real(v_jet["value"], "v"), _real(u4_jet["value"], "u4")
    dv, du4 = _vector(v_jet["gradient"], "grad v"), _vector(u4_jet["gradient"], "grad u4")
    return _real(8*(v*_dot(du4, n)-u4*_dot(dv, n)), "scalar adjoint flux")
