"""Exact localized scalar/York bilinear sources for the reduced H6 formula.

Use prepare_sixth_sources(data), where data is an existing dictionary from
reduced_fifth_order.prepare_fields with pointwise R3 and York b2 available.
No geometry, first inverse, nested inverse, energy or integral is recomputed.
Preparation of spatial derivative jets is separate from source assembly:
    jets = prepare_source_jets(data)
    result = assemble_sixth_sources(jets)
The wrapper performs these two steps only when explicitly called. There is
no field/source construction, numerical evaluation or checking on import.

Returned direct keys F,K are GaussianHeatSource scalars. r,t,J_L,J_R,A1W
are three-component tuples of covector sources; J_R aliases r. T1 is a
symmetric 3-by-3 tuple of GaussianHeatSource covariant tensor components.
statistics contains only structural, JSON-serializable counts and rates;
the source objects themselves retain exact Fraction coefficients.

Definitions (physical fields, with all Gaussian envelopes included):
    W=W1, T=L0 W, T1=L1 W, scalar w=7*u+v,
    F=-2*S:T-R3-8*Q:Hess(u),
    K=R3-14*S:T-8*Q:Hess(v),
    r_i=b2_i-partial_j T1_ij+Q_jk partial_k T_ij
                          +(partial_i Q_jk)*T_jk/2,
    t_i=S_ij partial_j w,  J_L=r+2*t,  J_R=r.
The compact A1 formula used in r requires flat TT Q and symmetric T.
S is the flat TT dual-Gram momentum field, not the canonical coefficient p.
W remains CONTRAVARIANT and fixed in L1/A1; all vector-source outputs are
covariant. No assumption of axial or single-parity data is made.

The external Gaussian exponents beta in GaussianHeatSource are 1/2 for
Q,S and their physical spatial derivatives, and 3/2 for R3,b2. They are
not the derivative-helper rates, which are twice beta. All supplied first
HeatFields have inner rate a=1; that rate stays inside I_m and is never
added to beta. Source differentiation includes both envelope and heat terms.

These sources retain the momentum-dependent H6 bilinear sectors. They are
not the full H6 energy: local terms, canonical Gram corrections, and the
no-Gaussian triple-Heat tail term remain separate inputs to that calculation.
"""

from collections.abc import Mapping
from fractions import Fraction as F

from gaussian_green_fields import HeatField
from gaussian_heat_sources import GaussianHeatSource, from_heat, from_polynomial
from reduced_gaunt_vertices import pg


def _shape(tree, shape, name):
    if not shape:
        return
    if not isinstance(tree, (tuple, list)) or len(tree) != shape[0]:
        raise ValueError(name + " must have shape " + str(shape))
    for i, item in enumerate(tree):
        _shape(item, shape[1:], name + "[" + str(i) + "]")


def _heat_tree(tree, shape, name):
    _shape(tree, shape, name)
    if shape:
        for i, item in enumerate(tree):
            _heat_tree(item, shape[1:], name + "[" + str(i) + "]")
    elif not isinstance(tree, HeatField) or tree.a != F(1):
        raise ValueError(name + " must be a supplied first HeatField with exact inner rate a=1")


def _symmetric_hessian(gradient):
    result = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(i, 3):
            result[i][j] = result[j][i] = gradient[i].derivative(j)
    return tuple(tuple(row) for row in result)


def prepare_source_jets(data):
    """Differentiate supplied PG/Heat fields only, without rebuilding geometry.

    Returns the raw q,s,R3,b2 plus W,dW,T,div_W,dq,dT,hess_u,hess_v,grad_w.
    dq[k][i][j] uses pg(q_ij,k,1), with exp(-r^2/2) omitted from that
    polynomial. dT[k][i][j] is an analytic HeatField derivative including
    its complete heat dependence. grad_w means the gradient of scalar 7u+v,
    distinct from the contravariant vector W and its gradient dW.
    """
    if not isinstance(data, Mapping):
        raise TypeError("data must be an existing prepared H5 mapping")
    required = ("q", "s", "geometry", "york", "w", "dw", "t0", "divergence_w", "grad_u", "grad_v")
    missing = [name for name in required if name not in data]
    if missing:
        raise ValueError("prepared H5 data is missing: " + ", ".join(missing))
    geometry, york = data["geometry"], data["york"]
    if geometry.get("q") != data["q"]:
        raise ValueError("supplied geometry and q do not agree")
    if 3 not in geometry.get("R", {}) or len(york.get("b", ())) < 3:
        raise ValueError("supplied data must already contain pointwise R3 and York b2")
    q, s, b2 = data["q"], data["s"], york["b"][2]
    _shape(q, (3, 3), "q")
    _shape(s, (3, 3), "s")
    _shape(b2, (3,), "b2")
    for name, shape in (("w", (3,)), ("dw", (3, 3)), ("t0", (3, 3)),
                        ("divergence_w", ()), ("grad_u", (3,)), ("grad_v", (3,))):
        _heat_tree(data[name], shape, name)

    dq = tuple(tuple(tuple(pg(q[i][j], k, 1) for j in range(3)) for i in range(3)) for k in range(3))
    T = data["t0"]
    dT = tuple(tuple(tuple(T[i][j].derivative(k) for j in range(3)) for i in range(3)) for k in range(3))
    return {
        "q": q, "s": s, "R3": geometry["R"][3], "b2": b2,
        "W": data["w"], "dW": data["dw"], "T": T, "div_W": data["divergence_w"],
        "dq": dq, "dT": dT,
        "hess_u": _symmetric_hessian(data["grad_u"]),
        "hess_v": _symmetric_hessian(data["grad_v"]),
        "grad_w": tuple(data["grad_u"][j].scale(7).add(data["grad_v"][j]) for j in range(3)),
    }


def _sum(weighted_sources):
    """Merge exact source keys once per contraction, avoiding successive copies."""
    terms = {}
    for factor, source in weighted_sources:
        for key, coefficient in source.terms.items():
            coefficient = terms.get(key, F(0)) + factor * coefficient
            if coefficient:
                terms[key] = coefficient
            else:
                terms.pop(key, None)
    return GaussianHeatSource(terms)


def _pg_heat(poly, heat):
    return from_heat(heat, poly, F(1, 2))


def assemble_sixth_sources(jets):
    """Assemble exact sources from prepare_source_jets output, without inverses.

    The divergence of T1 is obtained by differentiating GaussianHeatSource
    objects. This differentiates their external Gaussian as well as their
    heat moments, independently of any numerical finite-metric A1 evaluator.
    """
    q, s, T = jets["q"], jets["s"], jets["T"]
    W, dW, dq = jets["W"], jets["dW"], jets["dq"]
    T1 = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(i, 3):
            terms = [(F(-2, 3), _pg_heat(q[i][j], jets["div_W"]))]
            for k in range(3):
                terms.extend(((1, _pg_heat(dq[k][i][j], W[k])),
                              (1, _pg_heat(q[i][k], dW[j][k])),
                              (1, _pg_heat(q[j][k], dW[i][k]))))
            T1[i][j] = T1[j][i] = _sum(terms)
    T1 = tuple(tuple(row) for row in T1)

    st = _sum((1, _pg_heat(s[i][j], T[i][j])) for i in range(3) for j in range(3))
    qhu = _sum((1, _pg_heat(q[i][j], jets["hess_u"][i][j])) for i in range(3) for j in range(3))
    qhv = _sum((1, _pg_heat(q[i][j], jets["hess_v"][i][j])) for i in range(3) for j in range(3))
    R3 = from_polynomial(jets["R3"], F(3, 2))
    scalar_F = _sum(((-2, st), (-1, R3), (-8, qhu)))
    scalar_K = _sum(((1, R3), (-14, st), (-8, qhv)))

    A1W, r, t = [], [], []
    for i in range(3):
        terms = [(1, T1[i][j].derivative(j)) for j in range(3)]
        for j in range(3):
            for k in range(3):
                terms.append((-1, _pg_heat(q[j][k], jets["dT"][k][i][j])))
                terms.append((F(-1, 2), _pg_heat(dq[i][j][k], T[j][k])))
        A1W.append(_sum(terms))
        r.append(_sum(((1, from_polynomial(jets["b2"][i], F(3, 2))), (-1, A1W[-1]))))
        t.append(_sum((1, _pg_heat(s[i][j], jets["grad_w"][j])) for j in range(3)))
    A1W, r, t = tuple(A1W), tuple(r), tuple(t)
    J_L = tuple(r[i].add(t[i], 2) for i in range(3))
    statistics = {"F": scalar_F.statistics(), "K": scalar_K.statistics(),
                  "T1": [[entry.statistics() for entry in row] for row in T1]}
    for name, vector in (("r", r), ("t", t), ("J_L", J_L), ("J_R", r), ("A1W", A1W)):
        statistics[name] = [source.statistics() for source in vector]
    return {
        "F": scalar_F, "K": scalar_K, "r": r, "t": t,
        "J_L": J_L, "J_R": r, "T1": T1, "A1W": A1W,
        "statistics": statistics,
        "conventions": {
            "gaussian_exponents_beta": {"Q_S_and_derivatives": "1/2", "R3_b2": "3/2"},
            "inner_heat_rate_a": "1",
            "vector_type": "W contravariant; r,t,J_L,J_R,A1W covariant",
            "scalar_pairing": "-<K,C F>, C=(-Delta)^-1/8; not evaluated",
            "york_pairing": "-2<J_L,A0^-1 J_R>, A0=Delta_L; not evaluated",
            "scope": "Exact H6 bilinear source algebra only; no H6 energy or loop computed",
        },
    }


def prepare_sixth_sources(data):
    """Convenience wrapper: supplied-field derivative preparation, then exact assembly."""
    return assemble_sixth_sources(prepare_source_jets(data))
