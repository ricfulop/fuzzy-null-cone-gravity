"""Pointwise checks of the stripped first/second York tensors used in H6.

run_local_checks(data, sources, prepared_local, points, scaled_tolerance=1e-8)
returns JSON-serializable rows/passed/maximum plus strict check/failure counts.
prepared_local is supplied reduced_sixth_local.prepare_local_fields output.
The actual tensors are exp(-r^2/2)*HT1 and exp(-r^2)*HT2, including their
external envelopes. The reference independently evaluates metric, inverse,
Gamma and freshly differentiated metric coefficients through order two,
then calls york_series(max_order=2). W,dW,ddW are re-derived from ORIGINAL
data['w']; prepared_local HT1/HT2, data.dw, data.t0 and cached divergence
never manufacture the reference. No geometry or Green inverse is recomputed.

The trace check is Tr(T2)=Tr(Q*T1)-Tr(Q^2*T0)/2. Every matrix contraction
sums all Cartesian indices, including both symmetric off-diagonal entries.
The source T1 is also evaluated literally to cross-check its compression.

Each row separately records physical weak/strong A2 densities:
  weak = -T1:T1-2*T0:T2+4*Tr(Q*T0*T1)
         -Tr(Q^2*T0^2)-Tr(Q*T0*Q*T0),
  strong = 2*W^i*(A2 W)_i.
These are NOT pointwise equal; their difference includes a divergence.
Neither their difference nor an implied integrated equality enters pass/fail.
Integrated integration-by-parts and boundary-flux checks remain pending.

This is an algebra/envelope/indexing control with shared HeatField derivative
algebra and moment evaluation, not an independent Green solution or full-H6
validation. All derivative construction and numerical work occur only inside
the requested function call. Importing this module performs no calculation.
"""

from collections.abc import Mapping
from math import exp, fsum

from gaussian_green_fields import HeatField, evaluate_fields
from gaussian_heat_sources import GaussianHeatSource
from reduced_gaunt_vertices import pg
from sixth_source_controls import (
    _comparison, _finite, _point, _poly_tree, _restore, evaluate_sources,
)
from york_operator_fields import york_series


def _tree(tree, shape, name, leaf_type):
    if shape:
        if not isinstance(tree, (tuple, list)) or len(tree) != shape[0]:
            raise ValueError(name + " must have shape " + str(shape))
        for i, child in enumerate(tree):
            _tree(child, shape[1:], name + "[" + str(i) + "]", leaf_type)
    elif not isinstance(tree, leaf_type):
        raise TypeError(name + " has an invalid leaf type")
    elif leaf_type is HeatField and tree.a != 1:
        raise ValueError(name + " must have first-Heat inner rate a=1")


def _validate(data, sources, prepared):
    if any(not isinstance(value, Mapping) for value in (data, sources, prepared)):
        raise TypeError("data, sources and prepared_local must be supplied mappings")
    geometry = data["geometry"]
    if geometry.get("q") != data["q"]:
        raise ValueError("supplied geometry does not describe data.q")
    if any(len(geometry.get(name, ())) < 3 for name in ("metric", "inverse", "gamma")):
        raise ValueError("metric, inverse and Gamma coefficients through order two are required")
    if prepared.get("q") != data["q"] or prepared.get("s") != data["s"]:
        raise ValueError("prepared local fields must describe the supplied q,s")
    _tree(data["w"], (3,), "original W", HeatField)
    _tree(sources["T1"], (3, 3), "sources.T1", GaussianHeatSource)
    _tree(prepared["source_T1"], (3, 3), "prepared.source_T1", GaussianHeatSource)
    for name in ("HT1", "HT2"):
        _tree(prepared[name], (3, 3), "prepared." + name, HeatField)
    if any(prepared["source_T1"][i][j].terms != sources["T1"][i][j].terms
           for i in range(3) for j in range(3)):
        raise ValueError("prepared source_T1 does not match the supplied sources")


def _metric_jets(geometry):
    """Differentiate supplied g_n, n=0,1,2; pg's rate is n, not beta=n/2."""
    first, second = [], []
    for n in range(3):
        metric = geometry["metric"][n]
        derivative = [[[pg(metric[i][j], a, n) for j in range(3)]
                       for i in range(3)] for a in range(3)]
        hessian = [[[[pg(derivative[a][i][j], b, n) for j in range(3)]
                     for i in range(3)] for b in range(3)] for a in range(3)]
        first.append(derivative)
        second.append(hessian)
    return first, second


def _heat_layout(data, prepared):
    """Share evaluator moments, but keep actual fields separate from original jets."""
    W = data["w"]
    dW = [[W[k].derivative(a) for k in range(3)] for a in range(3)]
    ddW = [[[W[k].derivative(a).derivative(b) for k in range(3)]
            for b in range(3)] for a in range(3)]
    fields = []

    def store(tree, depth):
        if depth:
            return [store(child, depth - 1) for child in tree]
        fields.append(tree)
        return len(fields) - 1

    refs = {name: store(tree, depth) for name, tree, depth in
            (("W", W, 1), ("dW", dW, 2), ("ddW", ddW, 3),
             ("HT1", prepared["HT1"], 2), ("HT2", prepared["HT2"], 2))}
    return fields, refs


def _finite_tree(tree, depth, name):
    if depth:
        return [_finite_tree(child, depth - 1, name) for child in tree]
    return _finite(tree, name)


def _matrix_product(left, right):
    return [[_finite(fsum(left[i][k] * right[k][j] for k in range(3)), "matrix product")
             for j in range(3)] for i in range(3)]


def _trace_product(left, right):
    return _finite(fsum(left[i][j] * right[j][i]
                        for i in range(3) for j in range(3)), "trace product")


def _densities(q, q2, T, T1, T2, W, A2):
    """Physical reference tensors only; no integration or equality assertion."""
    qt = _matrix_product(q, T)
    families = {
        "york_A2_T1_square": -fsum(T1[i][j] * T1[i][j] for i in range(3) for j in range(3)),
        "york_A2_T_T2": -2 * fsum(T[i][j] * T2[i][j] for i in range(3) for j in range(3)),
        "york_A2_Q_T_T1": 4 * _trace_product(qt, T1),
        "york_A2_Q2_T_square": -_trace_product(q2, _matrix_product(T, T)),
        "york_A2_QTQT": -_trace_product(qt, qt),
    }
    families = {name: _finite(value, name) for name, value in families.items()}
    weak = _finite(fsum(families.values()), "weak A2 density")
    strong = _finite(2 * fsum(W[i] * A2[i] for i in range(3)), "strong A2 density")
    return {"weak_A2_families": families, "weak_A2": weak, "strong_A2": strong,
            "strong_minus_weak": _finite(strong - weak, "strong/weak density difference"),
            "strong_weak_equality_checked": False,
            "interpretation": "The difference includes a divergence; integrated IBP/flux remains unchecked"}


def run_local_checks(data, sources, prepared_local, points, scaled_tolerance=1e-8):
    """Check all T1/T2 components, T2 symmetry and trace; record A2 densities.

Returns one row per point with actual/reference tensors, component checks,
contractions, densities, passed/maximum/failure_count. Every component
comparison enters row/global aggregation. Scaled residuals are
abs(actual-reference)/max(1,abs(reference)); nonfinite values and empty
point lists are rejected. Strong-minus-weak density is diagnostic only.
"""
    tolerance = _finite(scaled_tolerance, "scaled_tolerance")
    if tolerance < 0:
        raise ValueError("scaled_tolerance must be nonnegative")
    points = tuple(_point(point) for point in points)
    if not points:
        raise ValueError("at least one local-check point is required")
    _validate(data, sources, prepared_local)
    geometry = data["geometry"]
    metric_first, metric_second = _metric_jets(geometry)
    heat_fields, heat_refs = _heat_layout(data, prepared_local)
    source_fields = tuple(sources["T1"][i][j] for i in range(3) for j in range(3))
    rows = []
    for point in points:
        values = tuple(_finite(value, "HeatField value") for value in evaluate_fields(heat_fields, point))
        physical = {name: _restore(ref, values) for name, ref in heat_refs.items()}
        source_values = evaluate_sources(source_fields, point)
        r2 = _finite(fsum(x * x for x in point), "radius squared")
        actual = {"T1": [[_finite(value * exp(-r2 / 2), "restored T1") for value in row]
                         for row in physical["HT1"]],
                  "T2": [[_finite(value * exp(-r2), "restored T2") for value in row]
                         for row in physical["HT2"]],
                  "source_T1": [list(source_values[3 * i:3 * i + 3]) for i in range(3)]}
        metric = [_poly_tree(geometry["metric"][n], 2, point, n) for n in range(3)]
        inverse = [_poly_tree(geometry["inverse"][n], 2, point, n) for n in range(3)]
        dmetric = [_poly_tree(metric_first[n], 3, point, n) for n in range(3)]
        ddmetric = [_poly_tree(metric_second[n], 4, point, n) for n in range(3)]
        gamma = [None if n == 0 and geometry["gamma"][n] is None else
                 _poly_tree(geometry["gamma"][n], 3, point, n) for n in range(3)]
        york = york_series(metric, dmetric, ddmetric, inverse, gamma,
                           physical["W"], physical["dW"], physical["ddW"], max_order=2)
        T, T1, T2 = _finite_tree(york["L"], 3, "physical York tensor")
        A2 = _finite_tree(york["A"][2], 1, "physical A2 W")
        q = _poly_tree(data["q"], 2, point, 1)
        q2 = _matrix_product(q, q)
        reference = {"T0": T, "T1": T1, "T2": T2, "A2W": A2, "W": physical["W"]}
        checks = {}
        for name, target in (("T1", T1), ("T2", T2), ("source_T1", T1)):
            for i in range(3):
                for j in range(3):
                    checks[name + "[" + str(i) + "][" + str(j) + "]"] = _comparison(
                        actual[name][i][j], target[i][j], tolerance)
        qt1, q2t = _trace_product(q, T1), _trace_product(q2, T)
        trace_target = _finite(qt1 - q2t / 2, "T2 trace target")
        checks["restored_T2_trace"] = _comparison(fsum(actual["T2"][i][i] for i in range(3)), trace_target, tolerance)
        checks["physical_T2_trace"] = _comparison(fsum(T2[i][i] for i in range(3)), trace_target, tolerance)
        for i in range(3):
            for j in range(i + 1, 3):
                checks["restored_T2_symmetry_" + str(i) + str(j)] = _comparison(
                    actual["T2"][i][j], actual["T2"][j][i], tolerance)
                checks["physical_T2_symmetry_" + str(i) + str(j)] = _comparison(T2[i][j], T2[j][i], tolerance)
        maximum = max(check["scaled_error"] for check in checks.values())
        failures = sum(not check["passed"] for check in checks.values())
        rows.append({"point": list(point), "actual": actual, "reference": reference, "checks": checks,
                     "contractions": {"Tr_Q_T1": qt1, "Tr_Q2_T0": q2t, "T2_trace_target": trace_target},
                     "densities": _densities(q, q2, T, T1, T2, physical["W"], A2),
                     "passed": failures == 0, "maximum": maximum, "failure_count": failures})
    return {"rows": rows, "passed": all(row["passed"] for row in rows),
            "maximum": max(row["maximum"] for row in rows), "scaled_tolerance": tolerance,
            "check_count": sum(len(row["checks"]) for row in rows),
            "failure_count": sum(row["failure_count"] for row in rows),
            "scale_definition": "abs(actual-reference)/max(1,abs(reference))",
            "full_functional_validation": False, "integrated_IBP_flux_checked": False,
            "scope": "Pointwise stripped York tensors, Gaussian envelopes, trace and symmetry only",
            "limitations": [
                "Reference jets use original W but share HeatField derivative algebra and heat_integrals.",
                "Supplied metric/inverse/Gamma coefficients are evaluated, not independently reconstructed.",
                "Strong/weak densities need not agree pointwise and their difference is not gated.",
                "No independent Green solution, integrated boundary flux, complete H6 or loop amplitude is validated.",
                "Sampled binary64 checks provide no global identity or interval certificate.",
            ]}
