"""Exact assembly of the positive-Gaussian local families of reduced H6.

Implements equations (5)/(6) of the sixth-order-bilinear-and-tail-reduction
report. This is assembly only: no numerical energy, quadrature, first or
nested Green solve, or loop amplitude is computed. Nothing is constructed
on import. All prepared geometry and first HeatFields are supplied by the
caller; geometry/source construction is never invoked here.

API:
  prepare_local_fields(data, sources) -> prepared local algebra
  local_families(data, sources, max_pair_products=20_000_000,
                 progress=None, prepared=None) -> {name: PairingAccumulator}
  merge_families(families, max_pair_products=20_000_000) -> PairingAccumulator
  triple_tail(data, max_pair_products=20_000_000, progress=None)
      -> (TripleHeatTail, tuple of per-component structural statistics)

data has the H5 prepare_fields layout, upgraded by the caller to pointwise
R4, integrated density L6, and York b3. sources is the corresponding output
of reduced_sixth_sources.prepare_sixth_sources, including its exact T1.
Both must describe the same raw flat-TT Q,S direction. S is dual-Gram
momentum, not the canonical coefficient p. All first HeatFields have a=1.

The returned local/merged accumulators EXCLUDE:
  g^T G0^-1 g + 2 c^T G0^-1 g, with c_alpha=int (7u+v) S:e_alpha;
  -2 <J_L,A0^-1 J_R> and -<K,C F>;
  -int (7u+v) T:T, which is assembled separately by triple_tail.
No E_c term is added to the local densities. The canonical correction must
be assembled once, outside this module, and must not be added again through
a separately constructed s4. A merged local result is NOT the full H6.

Physical envelopes are explicit: T=L0 W has none, while
T1=exp(-r^2/2)*HT1 and T2=exp(-r^2)*HT2. Stripping these known envelopes
does not change the internal rate a of a HeatField. Every PairingAccumulator
here has beta>0 and at most two HeatFields; the zero-beta triple term never
enters that backend. Q/S parity or axial symmetry is not used to remove terms.
"""

from collections.abc import Mapping
from fractions import Fraction as F
from numbers import Rational

from gaussian_green_fields import HeatField
from gaussian_heat_sources import GaussianHeatSource
from gaussian_heat_pairings import PairingAccumulator
from reduced_gaunt_exact import ONE, add, assert_tt, multiply
from reduced_gaunt_vertices import matrix_add, matrix_multiply, matrix_scale, pg, trace


LOCAL_FAMILIES = (
    "curvature", "raw_gram_fourth", "york_source_b3",
    "york_A2_T1_square", "york_A2_T_T2", "york_A2_Q_T_T1",
    "york_A2_Q2_T_square", "york_A2_QTQT", "scalar_curvature",
    "weighted_SQT", "weighted_S_T1", "weighted_commutator",
    "scalar_gradient", "scalar_B2", "scalar_s2_u2",
)


def _tensor(tensor, name):
    if not isinstance(tensor, (list, tuple)) or len(tensor) != 3:
        raise ValueError(name + " must be a 3-by-3 polynomial tensor")
    for row in tensor:
        if not isinstance(row, (list, tuple)) or len(row) != 3:
            raise ValueError(name + " must be a 3-by-3 polynomial tensor")
        for poly in row:
            if not isinstance(poly, Mapping):
                raise TypeError(name + " components must be polynomial mappings")
            for powers, coefficient in poly.items():
                if (not isinstance(powers, tuple) or len(powers) != 3
                        or any(type(p) is not int or p < 0 for p in powers)):
                    raise ValueError(name + " has invalid Cartesian powers")
                if isinstance(coefficient, bool) or not isinstance(coefficient, Rational):
                    raise TypeError(name + " must have exact rational coefficients")
    assert_tt(tensor)


def _heat_tree(tree, shape, name):
    if shape:
        if not isinstance(tree, (tuple, list)) or len(tree) != shape[0]:
            raise ValueError(name + " must have shape " + str(shape))
        for i, child in enumerate(tree):
            _heat_tree(child, shape[1:], name + "[" + str(i) + "]")
    elif not isinstance(tree, HeatField) or tree.a != 1:
        raise ValueError(name + " must be a first HeatField with exact inner rate one")


def _validate(data, sources):
    if not isinstance(data, Mapping) or not isinstance(sources, Mapping):
        raise TypeError("data and sources must be supplied mappings")
    _tensor(data["q"], "q")
    _tensor(data["s"], "s")
    geometry, york = data["geometry"], data["york"]
    if geometry.get("q") != data["q"]:
        raise ValueError("supplied geometry does not belong to q")
    if 4 not in geometry.get("R", {}) or 6 not in geometry.get("L", {}):
        raise ValueError("local H6 assembly requires supplied R4 and L6")
    if len(geometry.get("metric", ())) < 3 or len(york.get("b", ())) < 4:
        raise ValueError("local H6 assembly requires supplied g2 and York b3")
    for name, shape in (("u", ()), ("v", ()), ("w", (3,)), ("dw", (3, 3)),
                        ("divergence_w", ()), ("t0", (3, 3)), ("grad_u", (3,)), ("grad_v", (3,))):
        _heat_tree(data[name], shape, name)
    T1 = sources.get("T1")
    if not isinstance(T1, (tuple, list)) or len(T1) != 3:
        raise ValueError("sources must contain its prepared 3-by-3 T1")
    for row in T1:
        if not isinstance(row, (tuple, list)) or len(row) != 3:
            raise ValueError("source T1 must be 3-by-3")
        for source in row:
            if not isinstance(source, GaussianHeatSource):
                raise TypeError("source T1 components must be GaussianHeatSource")
    for i, row in enumerate(T1):
        for j, source in enumerate(row):
            if source.terms != T1[j][i].terms:
                raise ValueError("source T1 must be symmetric")


def _accumulate(terms, key, coefficient):
    coefficient = terms.get(key, F(0)) + coefficient
    if coefficient:
        terms[key] = coefficient
    else:
        terms.pop(key, None)


def _strip_T1(source):
    """Remove exactly beta=1/2; reject pure PG or any unexpected heat group."""
    terms = {}
    for (beta, a, m), polynomial in source.grouped_terms().items():
        if beta != F(1, 2) or a != F(1) or type(m) is not int or m < 0:
            raise ValueError("T1 compression requires only beta=1/2,a=1 heat groups")
        for powers, coefficient in polynomial.items():
            _accumulate(terms, (powers, m), coefficient)
    return HeatField(F(1), terms)


def _sum_poly_heat(weighted_products):
    """Exact P*HeatField sum with its common external Gaussian already removed."""
    terms = {}
    for factor, poly, field in weighted_products:
        for p, c in poly.items():
            for (q, m), d in field.terms.items():
                powers = tuple(left + right for left, right in zip(p, q))
                _accumulate(terms, (powers, m), factor * c * d)
    return HeatField(F(1), terms)


def prepare_local_fields(data, sources):
    """Prepare HT1, HT2 and local polynomials, without rebuilding source fields.

    This only strips known envelopes, multiplies supplied exact polynomials,
    and differentiates the supplied g2 polynomial with Gaussian rate two.
    It does not call prepare_fields, local_geometry, york_data, or an inverse.
    The caller must ensure sources and all supplied first fields correspond
    to the same q,s; their construction history is not encoded in source keys.
    """
    _validate(data, sources)
    q, s, geometry = data["q"], data["s"], data["geometry"]
    g2 = geometry["metric"][2]
    q2 = matrix_multiply(q, q)
    if g2 != matrix_scale(q2, F(1, 2)):
        raise ValueError("supplied g2 is not Q^2/2 for this q")
    HT1 = tuple(tuple(_strip_T1(source) for source in row) for row in sources["T1"])
    HT2 = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(i, 3):
            products = [(F(-2, 3), g2[i][j], data["divergence_w"])]
            for k in range(3):
                products.extend(((1, pg(g2[i][j], k, 2), data["w"][k]),
                                 (1, g2[i][k], data["dw"][j][k]),
                                 (1, g2[j][k], data["dw"][i][k])))
            HT2[i][j] = HT2[j][i] = _sum_poly_heat(products)
    comm = matrix_add(matrix_multiply(q, s), matrix_multiply(s, q), -1)
    double_comm = matrix_add(matrix_multiply(q, comm), matrix_multiply(comm, q), -1)
    return {
        "q": q, "s": s, "source_T1": sources["T1"], "Q2": q2,
        "HT1": HT1, "HT2": tuple(tuple(row) for row in HT2),
        "scalar_w": data["u"].scale(7).add(data["v"]),
        "commutator_square_trace": trace(matrix_multiply(comm, comm)),
        "double_commutator_square_trace": trace(matrix_multiply(double_comm, double_comm)),
        "SQ": matrix_multiply(s, q),
        "B2": add(geometry["R"][2], data["s2"], 7),
    }


def _add(accumulator, polynomial, beta, fields=(), factor=F(1)):
    if polynomial and factor and all(field.terms for field in fields):
        accumulator.add(polynomial, beta, fields=fields, factor=F(factor))


def local_families(data, sources, max_pair_products=20_000_000, progress=None, prepared=None):
    """Assemble 15 signed positive-Gaussian PairingAccumulators, not a full H6.

    Five separate families sum to 2<W,A2W>, permitting termwise weak-form
    audits. progress receives stage/family/completed/of after each family.
    Supply prepared=prepare_local_fields(data,sources) to reuse HT1/HT2.
    """
    if prepared is None:
        prepared = prepare_local_fields(data, sources)
    elif (prepared["q"] != data["q"] or prepared["s"] != data["s"]
          or prepared["source_T1"] != sources["T1"]):
        raise ValueError("prepared local fields do not match the supplied q,s,T1")
    out = {name: PairingAccumulator(max_pair_products=max_pair_products) for name in LOCAL_FAMILIES}
    def done(name):
        if progress is not None:
            progress({"stage": "sixth_local_assembly", "family": name,
                      "completed": LOCAL_FAMILIES.index(name) + 1, "of": len(LOCAL_FAMILIES)})
    q, s, T = data["q"], data["s"], data["t0"]
    q2, HT1, HT2, w = prepared["Q2"], prepared["HT1"], prepared["HT2"], prepared["scalar_w"]
    geometry = data["geometry"]
    _add(out["curvature"], geometry["L"][6], F(3), factor=-1)
    done("curvature")
    _add(out["raw_gram_fourth"], prepared["double_commutator_square_trace"], F(3), factor=F(-1, 360))
    done("raw_gram_fourth")
    for i in range(3):
        _add(out["york_source_b3"], data["york"]["b"][3][i], F(2), (data["w"][i],), -4)
    done("york_source_b3")
    for i in range(3):
        for j in range(3):
            _add(out["york_A2_T1_square"], ONE, F(1), (HT1[i][j], HT1[i][j]), -1)
    done("york_A2_T1_square")
    for i in range(3):
        for j in range(3):
            _add(out["york_A2_T_T2"], ONE, F(1), (T[i][j], HT2[i][j]), -2)
    done("york_A2_T_T2")
    for i in range(3):
        for j in range(3):
            for k in range(3):
                _add(out["york_A2_Q_T_T1"], q[i][j], F(1), (T[j][k], HT1[k][i]), 4)
    done("york_A2_Q_T_T1")
    for i in range(3):
        for j in range(3):
            for k in range(3):
                _add(out["york_A2_Q2_T_square"], q2[i][j], F(1), (T[j][k], T[k][i]), -1)
    done("york_A2_Q2_T_square")
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for ell in range(3):
                    _add(out["york_A2_QTQT"], multiply(q[i][j], q[k][ell]), F(1),
                         (T[j][k], T[ell][i]), -1)
    done("york_A2_QTQT")
    _add(out["scalar_curvature"], geometry["R"][4], F(2), (data["v"].add(data["u"], -1),))
    done("scalar_curvature")
    for i in range(3):
        for j in range(3):
            _add(out["weighted_SQT"], prepared["SQ"][i][j], F(1), (w, T[j][i]), -2)
    done("weighted_SQT")
    for i in range(3):
        for j in range(3):
            _add(out["weighted_S_T1"], s[i][j], F(1), (w, HT1[i][j]), 2)
    done("weighted_S_T1")
    _add(out["weighted_commutator"], prepared["commutator_square_trace"], F(2), (w,), F(1, 12))
    done("weighted_commutator")
    for i in range(3):
        for j in range(3):
            _add(out["scalar_gradient"], q2[i][j], F(1), (data["grad_v"][i], data["grad_u"][j]), 4)
    done("scalar_gradient")
    _add(out["scalar_B2"], prepared["B2"], F(1), (data["v"], data["u"]))
    done("scalar_B2")
    _add(out["scalar_s2_u2"], data["s2"], F(1), (data["u"], data["u"]), 28)
    done("scalar_s2_u2")
    return out


def merge_families(families, max_pair_products=20_000_000):
    """Merge exact local bins without repeating angular contractions or adding exclusions."""
    total = PairingAccumulator(max_pair_products=max_pair_products)
    for accumulator in families.values():
        total.merge(accumulator)
    return total


def triple_tail(data, max_pair_products=20_000_000, progress=None):
    """Assemble ONLY -int (7u+v) T:T with its separate zero-Gaussian backend.

    Return (summed TripleHeatTail, per-component work records). Off-diagonal
    tensor entries have multiplicity two, so six i<=j contractions suffice.
    No radial/exterior evaluation, tail bound, or falloff test is run here.
    """
    from algebraic_heat_tail import TripleHeatTail, collapse_triple_heat

    _heat_tree(data["u"], (), "u")
    _heat_tree(data["v"], (), "v")
    _heat_tree(data["t0"], (3, 3), "t0")
    T = data["t0"]
    if any(T[i][j].terms != T[j][i].terms for i in range(3) for j in range(3)):
        raise ValueError("T must be symmetric to use the off-diagonal multiplicities")
    w = data["u"].scale(7).add(data["v"])
    total, rows = TripleHeatTail({}), []
    for i in range(3):
        for j in range(i, 3):
            factor = -1 if i == j else -2
            part, statistics = collapse_triple_heat((w, T[i][j], T[i][j]), factor=F(factor),
                                                    max_pair_products=max_pair_products)
            total = total.add(part)
            rows.append({"i": i, "j": j, "factor": factor, "statistics": statistics})
            if progress is not None:
                progress({"stage": "sixth_triple_tail_assembly", "component": [i, j],
                          "completed": len(rows), "of": 6})
    return total, tuple(rows)
