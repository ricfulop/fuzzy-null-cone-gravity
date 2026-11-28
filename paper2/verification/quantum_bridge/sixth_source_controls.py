"""Pointwise checks of exact H6 sources against physical metric/HeatField jets.

Public API:
  evaluate_sources(flat_sources, point) -> tuple of finite binary64 values
  run_source_checks(data, sources, points, scaled_tolerance=1e-8)
      -> JSON-serializable {rows, passed, maximum, ...}

The source evaluator reads exact GaussianHeatSource terms directly. It does
not call source-construction code, use Gaussian-mixture quadrature, or omit
external exp(-beta*r^2) factors. All sources at a point share moment tables
for each inner rate a and Cartesian power tables.

The comparison independently evaluates supplied metric/inverse/connection
polynomials and physical spatial jets, then calls york_series(max_order=1).
W, dW, ddW and scalar gradients/Hessians are obtained from the ORIGINAL
data['w'], data['u'], data['v'] HeatFields. The source-preparation derivative
jets and stored t0 are not used to manufacture the reference. All symmetric
matrix contractions explicitly sum all nine i,j entries.

This is an algebra/envelope/indexing control, not an independent Green solve:
both paths share the analytic HeatField derivative algebra and heat_integrals.
The pointwise check does not validate the source data's geometric derivation,
nested bilinear integrals, spatial tails, full H6, or a loop amplitude.
Numerical work and derivative-field construction occur only when explicitly
called by the tracked runner; nothing is constructed/evaluated on import.
"""

from collections.abc import Mapping
from math import exp, fsum, isfinite

from gaussian_green_fields import HeatField, evaluate_fields, heat_integrals
from gaussian_heat_sources import GaussianHeatSource
from reduced_gaunt_vertices import pg
from york_operator_fields import york_series


def _finite(value, name):
    value = float(value)
    if not isfinite(value):
        raise ArithmeticError(name + " is nonfinite")
    return value


def _point(point):
    point = tuple(_finite(value, "Cartesian coordinate") for value in point)
    if len(point) != 3:
        raise ValueError("point must contain three Cartesian coordinates")
    return point


def _power_tables(point, maxima):
    tables = []
    for value, maximum in zip(point, maxima):
        row = [1.0]
        for _ in range(maximum):
            row.append(_finite(row[-1] * value, "Cartesian monomial power"))
        tables.append(row)
    return tables


def evaluate_sources(flat_sources, point):
    """Evaluate source entries literally, with shared moment/envelope tables.

    A term (beta,None,None,powers) has no heat factor. A heat term has
    exp(-beta*r^2)*I_m(a*r^2); a is not added to beta. Empty sources return
    zero. An empty source list returns an empty tuple after point validation.
    Invalid/nonfinite inputs, intermediates, and results raise exceptions.
    """
    point = _point(point)
    sources = tuple(flat_sources)
    if any(not isinstance(source, GaussianHeatSource) for source in sources):
        raise TypeError("flat_sources must contain GaussianHeatSource objects")
    maxima = [0, 0, 0]
    heat_orders, betas = {}, set()
    for source in sources:
        for beta, a, m, powers in source.terms:
            betas.add(beta)
            for i in range(3):
                maxima[i] = max(maxima[i], powers[i])
            if a is not None:
                heat_orders[a] = max(heat_orders.get(a, 0), m)
    powers_at = _power_tables(point, maxima)
    r2 = _finite(fsum(x * x for x in point), "radius squared")
    envelopes = {beta: exp(-_finite(float(beta) * r2, "Gaussian exponent")) for beta in betas}
    moments = {a: heat_integrals(_finite(float(a) * r2, "heat argument"), maximum)
               for a, maximum in heat_orders.items()}
    output = []
    for source in sources:
        summands = []
        for (beta, a, m, powers), coefficient in source.terms.items():
            value = _finite(coefficient, "source coefficient")
            for i in range(3):
                value *= powers_at[i][powers[i]]
            value *= envelopes[beta]
            if a is not None:
                value *= moments[a][m]
            summands.append(_finite(value, "source term"))
        output.append(_finite(fsum(summands), "source value"))
    return tuple(output)


def _polynomial_value(poly, point, rate):
    """Physical P*exp(-rate*r^2/2), with rate zero allowed for g0/inverse0."""
    radius_squared = _finite(fsum(x * x for x in point), "radius squared")
    envelope = exp(-rate * radius_squared / 2)
    terms = []
    for powers, coefficient in poly.items():
        value = _finite(coefficient, "polynomial coefficient")
        for coordinate, power in zip(point, powers):
            value *= coordinate ** power
        terms.append(_finite(value * envelope, "physical polynomial term"))
    return _finite(fsum(terms), "physical polynomial value")


def _poly_tree(tree, depth, point, rate):
    if depth == 0:
        return _polynomial_value(tree, point, rate)
    return [_poly_tree(child, depth - 1, point, rate) for child in tree]


def _source_layout(sources):
    """Validate all required components and remember their evaluation indices."""
    flat = []
    def visit(tree, shape, name):
        if shape:
            if not isinstance(tree, (tuple, list)) or len(tree) != shape[0]:
                raise ValueError(name + " must have shape " + str(shape))
            return [visit(child, shape[1:], name + "[" + str(i) + "]") for i, child in enumerate(tree)]
        if not isinstance(tree, GaussianHeatSource):
            raise TypeError(name + " must be a GaussianHeatSource")
        flat.append(tree)
        return len(flat) - 1
    layout = {name: visit(sources[name], shape, name) for name, shape in
              (("F", ()), ("K", ()), ("A1W", (3,)), ("T1", (3, 3)),
               ("r", (3,)), ("t", (3,)), ("J_L", (3,)), ("J_R", (3,)))}
    return flat, layout


def _restore(layout, values):
    if isinstance(layout, int):
        return values[layout]
    return [_restore(child, values) for child in layout]


def _original_heat_jets(data):
    """Re-differentiate original fields; do not reuse the assembly's cached jets."""
    W, u, v = data["w"], data["u"], data["v"]
    if not isinstance(W, (tuple, list)) or len(W) != 3:
        raise ValueError("data.w must contain three contravariant vector HeatFields")
    if any(not isinstance(field, HeatField) for field in (*W, u, v)):
        raise TypeError("original W,u,v must be HeatFields")
    dW = [[W[k].derivative(a) for k in range(3)] for a in range(3)]
    ddW = [[[W[k].derivative(a).derivative(b) for k in range(3)] for b in range(3)] for a in range(3)]
    gu, gv = u.gradient(), v.gradient()
    # Evaluate all nine entries independently, without copying the upper half.
    hu = [[u.derivative(i).derivative(j) for j in range(3)] for i in range(3)]
    hv = [[v.derivative(i).derivative(j) for j in range(3)] for i in range(3)]
    fields = []
    def store(tree, depth):
        if depth:
            return [store(child, depth - 1) for child in tree]
        fields.append(tree)
        return len(fields) - 1
    refs = {name: store(tree, depth) for name, tree, depth in
            (("W", W, 1), ("dW", dW, 2), ("ddW", ddW, 3),
             ("grad_u", gu, 1), ("grad_v", gv, 1), ("hess_u", hu, 2), ("hess_v", hv, 2))}
    return fields, refs


def _metric_jets(geometry):
    """Spatial derivatives of supplied series coefficients only; no Taylor geometry."""
    first, second = [], []
    for n in (0, 1):
        metric = geometry["metric"][n]
        derivative = [[[pg(metric[i][j], a, n) for j in range(3)] for i in range(3)] for a in range(3)]
        hessian = [[[[pg(derivative[a][i][j], b, n) for j in range(3)] for i in range(3)]
                    for b in range(3)] for a in range(3)]
        first.append(derivative)
        second.append(hessian)
    return first, second


def _comparison(actual, reference, tolerance):
    actual, reference = _finite(actual, "comparison actual"), _finite(reference, "comparison reference")
    residual = _finite(actual - reference, "comparison residual")
    scaled = abs(residual) / max(1.0, abs(reference))
    return {"actual": actual, "reference": reference, "residual": residual,
            "scaled_error": scaled, "passed": scaled <= tolerance}


def run_source_checks(data, sources, points, scaled_tolerance=1e-8):
    """Compare every supplied source component with a separate pointwise reference.

    One output row per point contains actual/reference values, component checks,
    full-index scalar contractions, passed and maximum. Every comparison enters
    row and global pass/failure aggregation; no failed component is ignored.
    Errors use abs(actual-reference)/max(1,abs(reference)). Empty point lists
    are rejected to prevent vacuous success. Nonfinite arithmetic raises.
    """
    tolerance = _finite(scaled_tolerance, "scaled_tolerance")
    if tolerance < 0:
        raise ValueError("scaled_tolerance must be nonnegative")
    points = tuple(_point(point) for point in points)
    if not points:
        raise ValueError("at least one source-check point is required")
    if not isinstance(data, Mapping) or not isinstance(sources, Mapping):
        raise TypeError("data and sources must be supplied mappings")
    geometry = data["geometry"]
    if geometry.get("q") != data["q"]:
        raise ValueError("supplied geometry and q do not agree")
    if (3 not in geometry.get("R", {}) or len(data["york"].get("b", ())) < 3
            or any(len(geometry.get(name, ())) < 2 for name in ("metric", "inverse", "gamma"))):
        raise ValueError("source checks require supplied R3,b2 and metric/inverse/Gamma through order one")
    flat_sources, source_layout = _source_layout(sources)
    heat_fields, heat_refs = _original_heat_jets(data)
    metric_first, metric_second = _metric_jets(geometry)
    rows = []
    for point in points:
        source_values = evaluate_sources(flat_sources, point)
        actual = {name: _restore(ref, source_values) for name, ref in source_layout.items()}
        # This separate evaluator call shares moments across all reference jets.
        heat_values = evaluate_fields(heat_fields, point)
        physical = {name: _restore(ref, heat_values) for name, ref in heat_refs.items()}
        metric = [_poly_tree(geometry["metric"][n], 2, point, n) for n in (0, 1)]
        inverse = [_poly_tree(geometry["inverse"][n], 2, point, n) for n in (0, 1)]
        dmetric = [_poly_tree(metric_first[n], 3, point, n) for n in (0, 1)]
        ddmetric = [_poly_tree(metric_second[n], 4, point, n) for n in (0, 1)]
        gamma = [None, _poly_tree(geometry["gamma"][1], 3, point, 1)]
        york = york_series(metric, dmetric, ddmetric, inverse, gamma,
                           physical["W"], physical["dW"], physical["ddW"], max_order=1)
        q, s = _poly_tree(data["q"], 2, point, 1), _poly_tree(data["s"], 2, point, 1)
        R3 = _polynomial_value(geometry["R"][3], point, 3)
        b2 = _poly_tree(data["york"]["b"][2], 1, point, 3)
        T, T1 = york["L"]
        st = _finite(fsum(s[i][j] * T[i][j] for i in range(3) for j in range(3)), "S:T")
        qhu = _finite(fsum(q[i][j] * physical["hess_u"][i][j]
                          for i in range(3) for j in range(3)), "Q:Hess u")
        qhv = _finite(fsum(q[i][j] * physical["hess_v"][i][j]
                          for i in range(3) for j in range(3)), "Q:Hess v")
        grad_w = [7 * physical["grad_u"][j] + physical["grad_v"][j] for j in range(3)]
        r = [b2[i] - york["A"][1][i] for i in range(3)]
        t = [fsum(s[i][j] * grad_w[j] for j in range(3)) for i in range(3)]
        reference = {"F": -2 * st - R3 - 8 * qhu, "K": R3 - 14 * st - 8 * qhv,
                     "A1W": york["A"][1], "T1": T1, "r": r, "t": t,
                     "J_L": [r[i] + 2 * t[i] for i in range(3)], "J_R": list(r)}
        checks = {}
        for name in ("F", "K"):
            checks[name] = _comparison(actual[name], reference[name], tolerance)
        for name in ("A1W", "r", "t", "J_L", "J_R"):
            for i in range(3):
                checks[name + "[" + str(i) + "]"] = _comparison(actual[name][i], reference[name][i], tolerance)
        for i in range(3):
            for j in range(3):
                checks["T1[" + str(i) + "][" + str(j) + "]"] = _comparison(actual["T1"][i][j], T1[i][j], tolerance)
        qt = _finite(fsum(q[i][j] * T[j][i] for i in range(3) for j in range(3)), "Tr(QT)")
        checks["source_T1_trace"] = _comparison(fsum(actual["T1"][i][i] for i in range(3)), qt, tolerance)
        checks["physical_T1_trace"] = _comparison(fsum(T1[i][i] for i in range(3)), qt, tolerance)
        for i in range(3):
            for j in range(i + 1, 3):
                checks["source_T1_symmetry_" + str(i) + str(j)] = _comparison(actual["T1"][i][j], actual["T1"][j][i], tolerance)
                checks["physical_T1_symmetry_" + str(i) + str(j)] = _comparison(T1[i][j], T1[j][i], tolerance)
        maximum = max(check["scaled_error"] for check in checks.values())
        failures = sum(not check["passed"] for check in checks.values())
        rows.append({"point": list(point), "actual": actual, "reference": reference, "checks": checks,
                     "contractions": {"S_dot_T": st, "Q_hess_u": qhu, "Q_hess_v": qhv,
                                      "Tr_QT": qt, "R3": R3},
                     "passed": failures == 0, "maximum": maximum, "failure_count": failures})
    return {
        "rows": rows, "passed": all(row["passed"] for row in rows),
        "maximum": max(row["maximum"] for row in rows), "scaled_tolerance": tolerance,
        "check_count": sum(len(row["checks"]) for row in rows),
        "failure_count": sum(row["failure_count"] for row in rows),
        "scale_definition": "abs(actual-reference)/max(1,abs(reference))",
        "scope": "Pointwise sixth-order source algebra, Gaussian envelopes and York indexing only",
        "limitations": [
            "Both routes share the original first HeatFields, their analytic derivative algebra and heat_integrals.",
            "The reference independently assembles physical jets and the full two-connection York operator; it is not an independent Green solution.",
            "Supplied geometry and b2 are not recomputed or independently derived here.",
            "Passing sampled points does not prove global equality or validate nested integrals, tails, full H6, or scattering.",
            "Binary64 evaluation can lose accuracy for strongly cancelling coefficients; no interval certificate is supplied.",
        ],
    }
