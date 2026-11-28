"""Exact finite-vector canonical correction in the reduced sixth-order energy.

Public API (nothing is constructed or evaluated on import):
  canonical_terms(data, basis, labels=None, *, gram_normalized=None,
                  precision=96, progress=None, ...compiler budgets...)
      -> JSON-serializable exact factors, normalized strings, physical floats
  canonical_controls(data, basis, result, labels, *, quartic_coefficients,
                     check_reciprocity=True, progress=None, ...budgets...)
      -> JSON-serializable exact component and contraction checks

data is supplied H5 data: q,s,geometry with pointwise R2, s2, york.b[1].
No geometry or first HeatField is rebuilt. basis is the supplied raw ten-mode
polynomial list, ordered E0..E4,B0..B4 as in exact_mode_basis. Each tensor
omits exp(-r^2/2); S is the dual-Gram momentum field, not canonical p.
Optional labels=(Eq,Bq,Es,Bs) are the four STF label matrices in that same
raw basis. Their reconstructed Q/S are checked without constructing modes.
If labels are omitted, coordinates follow from exact Gaussian projection.
gram_normalized, when supplied, is G0/pi^(3/2) and is checked against exact
moments of the supplied basis. Its entries may be Fraction/int or strings
from exact_mode_basis metadata. No diagonal-Gram assumption is made.

The physical vectors are
  g_i=-int Tr([Q,e_i][Q,S])/12 + 2<b1(Q,e_i),Y b1(Q,S)>,
  c_i=<S:e_i,C(14*s2-6*R2)>,  C=(-Delta)^-1/8, Y=Delta_L^-1.
Both Fourier endpoints have external beta=1 and no heat parameter. For a
compiled term c*A^i*B^j/(A+B)^k/sqrt(A+B), A=B=1 exactly. Its rational
coefficient in U4=pi^(3/2)/sqrt(2) is c/2^k, so evaluation uses NO quadrature.
The raw commutator has beta=2: its integral divided by U4 is integrate(P,4)/2.

Return conventions:
  g_bar,c_bar = exact rational strings for g/U4,c/U4;
  Ec_bar = Gbar^-1*g_bar, with Gbar=G0/pi^(3/2);
  Ec_value = physical coefficients Ec_bar/sqrt(2) in E_c=sum_i e_i Ec_i.
E_c has this positive sign; the momentum expansion SUBTRACTS E_c.
  gram_square_exact = g_bar^T Gbar^-1 g_bar;
  canonical_cross_exact = 2*c_bar^T Gbar^-1 g_bar;
  total_exact = their sum.
These last three rational factors multiply U6=pi^(3/2)/2. Corresponding
*_normalized_decimal strings use U6; *_value are physical binary64 values.
The parallel *_without_pi_decimal strings factor out only pi^(3/2).

This is ONLY g^T G0^-1 g + 2*c^T G0^-1 g. It excludes every local density,
the triple-Heat tail, and both nested scalar/York bilinears. Add it once to
the other H6 families; do not add a second weighted S:E_c contribution.
All phase-space sectors are retained. For pure S (Q=0), g=0 and the energy
correction vanishes, while c need not vanish. No amplitude is computed.

Controls use already supplied complete_quartic.calculate_coefficients
output, not a new reconstruction. The mixed H4 polynomial is quadratic in
S, permitting its exact two-point derivative identity. The mixed+pure-P
quartic derivative uses the exact five-point formula, not a two-point
approximation. Reciprocity shares the new compiler; the component targets
use the archived scalar Fourier and independent quartic representation.
Passing these checks does not independently validate a Green solution or
the complete mixed H6 energy.
"""

from collections.abc import Mapping
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import isfinite, pi, sqrt
from numbers import Rational

from gaussian_fourier_bilinears import compile_scalar_bilinear, compile_york_bilinear
from gaussian_heat_sources import from_polynomial
from reduced_gaunt_exact import add, five_polarizations, multiply, scale, tensor_inner
from reduced_gaunt_vertices import (
    contract, coulomb, integrate, invariants, label_add, matrix_add,
    matrix_multiply, pg, psum, trace, york_norm,
)


def _fraction(value, name):
    if isinstance(value, bool) or not isinstance(value, (Rational, str)):
        raise TypeError(name + " must be exact rational (or a rational string)")
    return F(value)


def _poly(value, name):
    if not isinstance(value, Mapping):
        raise TypeError(name + " must be a polynomial mapping")
    for powers, coefficient in value.items():
        if (not isinstance(powers, tuple) or len(powers) != 3
                or any(type(p) is not int or p < 0 for p in powers)):
            raise ValueError(name + " has invalid Cartesian powers")
        if isinstance(coefficient, bool) or not isinstance(coefficient, Rational):
            raise TypeError(name + " polynomial coefficients must be exact rationals")


def _tensor(value, name):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(name + " must be a 3-by-3 polynomial tensor")
    if any(not isinstance(row, (list, tuple)) or len(row) != 3 for row in value):
        raise ValueError(name + " must be a 3-by-3 polynomial tensor")
    for i in range(3):
        for j in range(3):
            _poly(value[i][j], name)
            if value[i][j] != value[j][i]:
                raise ValueError(name + " is not symmetric")
        if psum(pg(value[i][j], j, 1) for j in range(3)):
            raise ValueError(name + " is not flat transverse")
    if trace(value):
        raise ValueError(name + " is not flat traceless")


def _labels(labels):
    if not isinstance(labels, (list, tuple)) or len(labels) != 4:
        raise ValueError("labels must be (Eq,Bq,Es,Bs)")
    result = []
    for matrix in labels:
        if (not isinstance(matrix, (list, tuple)) or len(matrix) != 3
                or any(not isinstance(row, (list, tuple)) or len(row) != 3 for row in matrix)):
            raise ValueError("Each label must be a 3-by-3 STF matrix")
        a = tuple(tuple(_fraction(x, "label") for x in row) for row in matrix)
        if any(a[i][j] != a[j][i] for i in range(3) for j in range(3)) or sum(a[i][i] for i in range(3)):
            raise ValueError("Each label must be symmetric and traceless")
        result.append(a)
    return tuple(result)


def _label_coordinates(matrix):
    # Polarizations are diag(1,-1,0), diag(-1/2,-1/2,1), xy, xz, yz.
    return ((matrix[0][0] - matrix[1][1]) / 2, matrix[2][2],
            matrix[0][1], matrix[0][2], matrix[1][2])


def _combination(basis, coefficients):
    return [[psum(scale(mode[i][j], c) for mode, c in zip(basis, coefficients) if c)
             for j in range(3)] for i in range(3)]


def _same_tensor(left, right):
    return all(left[i][j] == right[i][j] for i in range(3) for j in range(3))


def _zero_tensor(tensor):
    return not any(coefficient for row in tensor for poly in row for coefficient in poly.values())


def _ldl(matrix):
    """Exact positive-definite factorization; keep a general full Gram matrix."""
    size = len(matrix)
    lower, diagonal = [[F(int(i == j)) for j in range(size)] for i in range(size)], []
    for i in range(size):
        d = matrix[i][i] - sum((lower[i][k] ** 2 * diagonal[k] for k in range(i)), F(0))
        if d <= 0:
            raise ValueError("The supplied basis Gram matrix is not positive definite")
        diagonal.append(d)
        for j in range(i + 1, size):
            lower[j][i] = (matrix[j][i] - sum(
                (lower[j][k] * lower[i][k] * diagonal[k] for k in range(i)), F(0))) / d
    return lower, diagonal


def _solve(factorization, rhs):
    lower, diagonal = factorization
    y = []
    for i, b in enumerate(rhs):
        y.append(b - sum((lower[i][j] * y[j] for j in range(i)), F(0)))
    x = [F(0)] * len(rhs)
    for i in reversed(range(len(rhs))):
        x[i] = y[i] / diagonal[i] - sum(
            (lower[j][i] * x[j] for j in range(i + 1, len(rhs))), F(0))
    return tuple(x)


def _first_york_source(q, s, dq):
    """Compact flat-TT b1; physical derivative helper rate is one (beta=1/2)."""
    ds = [[[pg(s[i][j], k, 1) for j in range(3)] for i in range(3)] for k in range(3)]
    return tuple(scale(psum(add(add(multiply(s[l][j], dq[l][j][i]),
        multiply(s[j][l], dq[i][j][l]), -1), multiply(q[j][l], ds[l][i][j]), -1)
        for l in range(3) for j in range(3)), F(1, 2)) for i in range(3))


def _commutator(left, right):
    return matrix_add(matrix_multiply(left, right), matrix_multiply(right, left), -1)


def _trace_product(left, right):
    return psum(multiply(left[i][j], right[j][i]) for i in range(3) for j in range(3))


def _prepare(data, basis, labels, gram_normalized, progress):
    if not isinstance(data, Mapping):
        raise TypeError("data must be supplied H5 data")
    q, s = data["q"], data["s"]
    _tensor(q, "q")
    _tensor(s, "s")
    basis = tuple(basis)
    if len(basis) != 10:
        raise ValueError("The raw ten-mode basis must contain exactly ten tensors")
    for i, mode in enumerate(basis):
        _tensor(mode, "basis[" + str(i) + "]")
    geometry, york = data["geometry"], data["york"]
    if geometry.get("q") != q or 2 not in geometry.get("R", {}):
        raise ValueError("data must supply pointwise R2 for this q")
    r2 = geometry["R"][2]
    _poly(r2, "R2")
    s2 = contract(s, s)
    if data.get("s2") != s2:
        raise ValueError("supplied s2 differs from the full symmetric contraction S:S")
    if len(york.get("b", ())) < 2 or len(york["b"][1]) != 3:
        raise ValueError("data must supply the three polynomial b1 components")
    b1 = tuple(york["b"][1])
    for component in b1:
        _poly(component, "b1")
    dq = [[[pg(q[i][j], k, 1) for j in range(3)] for i in range(3)] for k in range(3)]
    if b1 != _first_york_source(q, s, dq):
        raise ValueError("supplied b1 differs from the compact flat-TT source for q,s")
    gram = [[F(0)] * 10 for _ in range(10)]
    for i in range(10):
        for j in range(i, 10):
            gram[i][j] = gram[j][i] = tensor_inner(basis[i], basis[j])
        if progress:
            progress({"stage": "canonical_gram_row", "completed": i + 1, "of": 10})
    if gram_normalized is not None:
        if (len(gram_normalized) != 10 or any(len(row) != 10 for row in gram_normalized)):
            raise ValueError("gram_normalized must be 10-by-10")
        supplied = [[_fraction(x, "Gram entry") for x in row] for row in gram_normalized]
        if supplied != gram:
            raise ValueError("supplied G0/pi^(3/2) differs from exact basis Gaussian moments")
    factorization = _ldl(gram)
    labels = None if labels is None else _labels(labels)
    if labels is None:
        q_coeff = _solve(factorization, [tensor_inner(e, q) for e in basis])
        s_coeff = _solve(factorization, [tensor_inner(e, s) for e in basis])
    else:
        q_coeff = _label_coordinates(labels[0]) + _label_coordinates(labels[1])
        s_coeff = _label_coordinates(labels[2]) + _label_coordinates(labels[3])
    if not _same_tensor(_combination(basis, q_coeff), q):
        raise ValueError("q is not the stated combination of the supplied raw basis")
    if not _same_tensor(_combination(basis, s_coeff), s):
        raise ValueError("s is not the stated combination of the supplied raw basis")
    return {"q": q, "s": s, "R2": r2, "s2": s2, "b1": b1, "dq": dq,
            "basis": basis, "gram": gram, "factorization": factorization,
            "q_coefficients": q_coeff, "s_coefficients": s_coeff, "labels": labels,
            "commutator": _commutator(q, s), "w_source": add(scale(s2, 14), r2, -6)}


def _options(max_pair_products, max_fourier_terms, max_terms, max_degree, progress):
    options = dict(max_pair_products=max_pair_products, max_fourier_terms=max_fourier_terms,
                   max_terms=max_terms, max_degree=max_degree)
    if any(type(value) is not int or value < 1 for value in options.values()):
        raise ValueError("Every compiler budget must be a positive integer")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable or None")
    return options


def _progress(progress, index, family):
    if progress is None:
        return None
    return lambda event: progress({**event, "canonical_index": index, "canonical_family": family})


def _unit_rate_coefficient(kernel):
    """Exact pure-PG zero-dimensional compiler value divided by U4."""
    value = F(0)
    for block, terms in kernel.blocks.items():
        if block != ((F(1), None), (F(1), None)):
            raise ValueError("Canonical bilinears require two pure-PG beta=1 endpoints")
        for (m, n, _i, _j, k), coefficient in terms.items():
            if m != -1 or n != -1:
                raise ValueError("Canonical bilinears must have no heat integration variables")
            value += coefficient * F(2) ** (-k)
    return value


def _decimal(value, precision):
    with localcontext() as ctx:
        ctx.prec = precision
        return str(D(value.numerator) / D(value.denominator))


def _physical(value, unit):
    try:
        result = float(value) * unit
    except OverflowError as exc:
        raise ArithmeticError("Canonical physical float exceeds binary64 range") from exc
    if not isfinite(result) or (value and result == 0):
        raise ArithmeticError("Canonical physical value is nonfinite or underflows binary64")
    return result


def canonical_terms(data, basis, labels=None, *, gram_normalized=None, precision=96,
                    progress=None, max_pair_products=20_000_000,
                    max_fourier_terms=200_000, max_terms=200_000, max_degree=64):
    """Evaluate just the two canonical H6 families using exact rational algebra.

Physical assembly keys are gram_square_value, canonical_cross_value,
total_value. Fractions are serialized as strings; no tolerance is needed
for the algebra. Decimal strings round exact rational normalized factors
to precision significant digits and do not approximate pi. No scientific
calls occur unless this function is explicitly invoked by the caller.
"""
    options = _options(max_pair_products, max_fourier_terms, max_terms, max_degree, progress)
    if type(precision) is not int or precision < 16:
        raise ValueError("precision must be an integer >=16")
    prepared = _prepare(data, basis, labels, gram_normalized, progress)
    base_b = tuple(from_polynomial(poly, F(1)) for poly in prepared["b1"])
    w_source = from_polynomial(prepared["w_source"], F(1))
    g, c, raw, york, statistics = [], [], [], [], []
    for index, mode in enumerate(prepared["basis"]):
        if progress:
            progress({"stage": "canonical_component_start", "index": index, "of": 10})
        b = _first_york_source(prepared["q"], mode, prepared["dq"])
        vector = tuple(from_polynomial(poly, F(1)) for poly in b)
        y_kernel = compile_york_bilinear(vector, base_b, **options,
            progress=_progress(progress, index, "g_york"))
        endpoint = contract(prepared["s"], mode)
        c_kernel = compile_scalar_bilinear(from_polynomial(endpoint, F(1)), w_source,
            **options, progress=_progress(progress, index, "c"))
        raw.append(-integrate(_trace_product(_commutator(prepared["q"], mode),
                                             prepared["commutator"]), 4) / 24)
        york.append(2 * _unit_rate_coefficient(y_kernel))
        g.append(raw[-1] + york[-1])
        c.append(_unit_rate_coefficient(c_kernel))
        statistics.append({"index": index, "york": y_kernel.statistics(),
                           "scalar": c_kernel.statistics()})
        if progress:
            progress({"stage": "canonical_component_complete", "completed": index + 1, "of": 10,
                      "g_bar": str(g[-1]), "c_bar": str(c[-1])})
    ec = _solve(prepared["factorization"], g)
    square = sum((x * y for x, y in zip(g, ec)), F(0))
    cross = 2 * sum((x * y for x, y in zip(c, ec)), F(0))
    result = {
        "g_bar": [str(x) for x in g], "c_bar": [str(x) for x in c],
        "g_raw_bar": [str(x) for x in raw], "g_york_bar": [str(x) for x in york],
        "Ec_bar": [str(x) for x in ec],
        "g_normalized_decimal": [_decimal(x, precision) for x in g],
        "c_normalized_decimal": [_decimal(x, precision) for x in c],
        "Ec_normalized_decimal": [_decimal(x, precision) for x in ec],
        "g_value": [_physical(x, pi ** 1.5 / sqrt(2)) for x in g],
        "c_value": [_physical(x, pi ** 1.5 / sqrt(2)) for x in c],
        "Ec_value": [_physical(x, 1 / sqrt(2)) for x in ec],
        "gram_normalized": [[str(x) for x in row] for row in prepared["gram"]],
        "q_coefficients": [str(x) for x in prepared["q_coefficients"]],
        "s_coefficients": [str(x) for x in prepared["s_coefficients"]],
        "labels_Eq_Bq_Es_Bs": None if prepared["labels"] is None else
            [[[str(x) for x in row] for row in matrix] for matrix in prepared["labels"]],
        "normalization": {
            "g_bar_c_bar_g_raw_bar_g_york_bar": "multiply by pi^(3/2)/sqrt(2)",
            "Ec_bar": "multiply by 1/sqrt(2); + sign in E_c; the momentum correction subtracts E_c",
            "gram_normalized": "multiply by pi^(3/2)",
            "family_exact_and_normalized_decimal": "multiply by pi^(3/2)/2",
            "family_without_pi_decimal": "multiply by pi^(3/2)",
            "value_fields": "physical binary64 values, including all pi factors",
        },
        "statistics": {"basis_size": 10, "retained_internal_modes": 10,
                       "heat_integration_dimension": 0, "precision": precision,
                       "compiler_options": options, "components": statistics,
                       "supplied_gram_checked": gram_normalized is not None},
        "scope": "Only the two finite-vector canonical H6 terms; not the full energy or a loop amplitude",
    }
    for name, value in (("gram_square", square), ("canonical_cross", cross), ("total", square + cross)):
        result[name + "_exact"] = str(value)
        result[name + "_normalized_decimal"] = _decimal(value, precision)
        result[name + "_without_pi_decimal"] = _decimal(value / 2, precision)
        result[name + "_value"] = _physical(value, pi ** 1.5 / 2)
    return result


def _result_vector(result, name):
    values = result[name]
    if not isinstance(values, (list, tuple)) or len(values) != 10:
        raise ValueError("result." + name + " must have ten exact entries")
    return tuple(_fraction(value, name) for value in values)


def canonical_controls(data, basis, result, labels, *, quartic_coefficients,
                       check_reciprocity=True, progress=None, max_pair_products=20_000_000,
                       max_fourier_terms=200_000, max_terms=200_000, max_degree=64):
    """Exact gates against supplied quartic coefficients and archived Fourier.

quartic_coefficients is the full calculate_coefficients result with keys
mixed and momentum, not a path and not a command to reconstruct it. For
each raw basis component the independent targets in U4 units are
  gamma_i=6*C(S:e_i,R2)-[H4mixed(S+e_i)-H4mixed(S-e_i)]/4,
  kappa_i=14*C(S:e_i,s2)-6*C(S:e_i,R2),
  gamma_i+kappa_i=-D_i(H4mixed+H4pureP)/2.
The derivative of the degree-four polynomial is evaluated exactly using
[H(-2)-8H(-1)+8H(1)-H(2)]/12 about S. It is not a numerical finite
difference approximation. Pure-sector checks apply only to the supplied
direction; this function never constructs extra data cases or Green fields.
"""
    from complete_quartic import evaluate_mixed

    options = _options(max_pair_products, max_fourier_terms, max_terms, max_degree, progress)
    if type(check_reciprocity) is not bool:
        raise TypeError("check_reciprocity must be bool")
    if not isinstance(quartic_coefficients, Mapping):
        raise TypeError("quartic_coefficients must be the supplied reconstructed coefficient mapping")
    mixed = quartic_coefficients["mixed"]
    # Reject missing/short invariant vectors instead of silently truncating zip.
    if set(mixed["paired"]) != {"00", "01", "10", "11"}:
        raise ValueError("Quartic mixed paired coefficients must contain all four species")
    mixed = {"paired": {key: tuple(_fraction(x, "quartic coefficient") for x in values)
                        for key, values in mixed["paired"].items()},
             "cross": tuple(_fraction(x, "quartic coefficient") for x in mixed["cross"])}
    momentum = tuple(_fraction(x, "quartic coefficient") for x in quartic_coefficients["momentum"])
    if (any(len(x) != 3 for x in mixed["paired"].values()) or len(mixed["cross"]) != 5
            or len(momentum) != 5):
        raise ValueError("Incomplete quartic mixed/momentum coefficient vectors")
    labels = _labels(labels)
    prepared = _prepare(data, basis, labels, result["gram_normalized"], progress)
    g, c = _result_vector(result, "g_bar"), _result_vector(result, "c_bar")
    raw, york = _result_vector(result, "g_raw_bar"), _result_vector(result, "g_york_bar")
    ec = _result_vector(result, "Ec_bar")
    rows = []

    def check(name, actual, reference, index=None):
        residual = actual - reference
        row = {"name": name, "actual": str(actual), "reference": str(reference),
               "residual": str(residual), "passed": residual == 0}
        if index is not None:
            row["index"] = index
        rows.append(row)

    def total_quartic(current_labels):
        return evaluate_mixed(mixed, current_labels) + sum((coefficient * invariant
            for coefficient, invariant in zip(momentum,
                invariants(current_labels[2], current_labels[3])["quartic"])), F(0))

    polarizations = five_polarizations()
    base_b = tuple(from_polynomial(poly, F(1)) for poly in prepared["b1"]) if check_reciprocity else None
    w_source = from_polynomial(prepared["w_source"], F(1)) if check_reciprocity else None
    for index, mode in enumerate(prepared["basis"]):
        endpoint = contract(prepared["s"], mode)
        c_r2, c_s2 = coulomb(endpoint, prepared["R2"]), coulomb(endpoint, prepared["s2"])
        shifted = {}
        for step in (-2, -1, 1, 2):
            current = list(labels)
            position = 2 + index // 5
            current[position] = label_add(current[position], polarizations[index % 5], F(step))
            shifted[step] = tuple(current)
        mixed_difference = (evaluate_mixed(mixed, shifted[1]) - evaluate_mixed(mixed, shifted[-1])) / 4
        derivative = (total_quartic(shifted[-2]) - 8 * total_quartic(shifted[-1])
                      + 8 * total_quartic(shifted[1]) - total_quartic(shifted[2])) / 12
        check("g_component_quartic", g[index], 6 * c_r2 - mixed_difference, index)
        check("c_component_archived_scalar", c[index], 14 * c_s2 - 6 * c_r2, index)
        check("g_plus_c_quartic_derivative", g[index] + c[index], -derivative / 2, index)
        check("g_raw_plus_york", g[index], raw[index] + york[index], index)
        check("Ec_Gram_equation", sum((prepared["gram"][index][j] * ec[j] for j in range(10)), F(0)), g[index], index)
        if check_reciprocity:
            b = _first_york_source(prepared["q"], mode, prepared["dq"])
            y_reverse = compile_york_bilinear(base_b, tuple(from_polynomial(poly, F(1)) for poly in b),
                **options, progress=_progress(progress, index, "york_reciprocity"))
            c_reverse = compile_scalar_bilinear(w_source, from_polynomial(endpoint, F(1)),
                **options, progress=_progress(progress, index, "scalar_reciprocity"))
            check("york_reciprocity", york[index], 2 * _unit_rate_coefficient(y_reverse), index)
            check("scalar_reciprocity", c[index], _unit_rate_coefficient(c_reverse), index)
        if progress:
            progress({"stage": "canonical_control_component_complete", "completed": index + 1,
                      "of": 10, "failure_count_so_far": sum(not row["passed"] for row in rows)})
    a = prepared["s_coefficients"]
    a_g = sum((x * y for x, y in zip(a, g)), F(0))
    a_c = sum((x * y for x, y in zip(a, c)), F(0))
    c_sr, c_ss = coulomb(prepared["s2"], prepared["R2"]), coulomb(prepared["s2"], prepared["s2"])
    check("g_directional_quartic", a_g, 6 * c_sr - evaluate_mixed(mixed, labels))
    comm_square = _trace_product(prepared["commutator"], prepared["commutator"])
    check("g_directional_archived_York", a_g, -integrate(comm_square, 4) / 24 - york_norm(prepared["b1"]))
    check("c_directional", a_c, 14 * c_ss - 6 * c_sr)
    square = sum((x * y for x, y in zip(g, ec)), F(0))
    cross = 2 * sum((x * y for x, y in zip(c, ec)), F(0))
    check("gram_square_contraction", _fraction(result["gram_square_exact"], "gram square"), square)
    check("canonical_cross_contraction", _fraction(result["canonical_cross_exact"], "canonical cross"), cross)
    check("total_contraction", _fraction(result["total_exact"], "canonical total"), square + cross)
    if _zero_tensor(prepared["q"]):
        for index in range(10):
            check("pure_momentum_g_zero", g[index], F(0), index)
        check("pure_momentum_canonical_total_zero", square + cross, F(0))
    if _zero_tensor(prepared["s"]):
        for index in range(10):
            check("pure_configuration_g_zero", g[index], F(0), index)
            check("pure_configuration_c_zero", c[index], F(0), index)
        check("pure_configuration_canonical_total_zero", square + cross, F(0))
    failures = sum(not row["passed"] for row in rows)
    return {"rows": rows, "passed": failures == 0, "check_count": len(rows), "failure_count": failures,
            "maximum_absolute_exact_residual": str(max((abs(F(row["residual"])) for row in rows), default=F(0))),
            "normalization": "Vector checks in U4=pi^(3/2)/sqrt(2); final contractions in U6=pi^(3/2)/2",
            "component_targets": "Supplied independent quartic coefficients plus archived scalar Fourier",
            "reciprocity_checked": check_reciprocity,
            "limitations": ["Reciprocity shares the new compiler; it is not an independent Green solve",
                            "Pure-sector gates cover only a supplied pure-sector case",
                            "No complete H6 or finite-mode loop validation is implied"]}
