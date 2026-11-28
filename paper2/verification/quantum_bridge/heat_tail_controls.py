"""Bounded controls for exact triple-HEAT collapse and its algebraic exterior.

run_tail_controls(config) -> JSON-compatible synthetic control results.
check_gravitational_tail(data, name, config, progress=None) -> JSON results
    for supplied first-response fields; no geometry or inverses are constructed.

All settings use the tail_* prefix; defaults are declared in _settings below.
The source monopoles are rational integrals with pi**(3/2) removed. The
pointwise sphere comparison includes 4*pi*r**2 and the off-diagonal tensor
multiplicity two. It is never compared with an integrated remainder bound.

The sphere route shares the HEAT-moment evaluator with the radial route, so
it checks angular contraction and assembly, not independent HEAT evaluation.
Exterior bounds evaluated in binary64 are diagnostics, not certificates.
No control, field construction, geometry, quadrature or evaluation runs on import.
"""

from collections.abc import Mapping
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from itertools import permutations
from math import erfc, exp, factorial, fsum, isfinite, pi

from algebraic_heat_tail import (
    TripleHeatTail, collapse_triple_heat, moment_exterior_coefficients,
)
from gaussian_green_fields import HeatField
from h5_spatial_reference import _CompiledFields
from reduced_gaunt_exact import ONE, add, gaussian_inner, scale
from reduced_sixth_local import triple_tail


def _positive_float(value, name):
    value = float(value)
    if not isfinite(value) or value <= 0:
        raise ValueError(name + " must be finite and strictly positive")
    return value


def _positive_int(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name + " must be a positive integer")
    return value


def _radii(values, name):
    result = tuple(_positive_float(value, name) for value in values)
    if not result or len(result) > 64:
        raise ValueError(name + " must contain between one and 64 radii")
    if any(not isfinite(3*r*r) or r*r == 0 for r in result):
        raise ValueError(name + " has a squared radius outside binary64 range")
    return result


def _settings(config):
    if not isinstance(config, Mapping):
        raise TypeError("config must be a mapping")
    settings = {
        "tail_control_tolerance": 1e-11,
        "tail_high_m_relative_tolerance": 1e-12,
        "tail_decimal_precision": 80,
        "tail_identity_radii": (0.5, 1, 2, 4),
        "tail_spatial_tolerance": 1e-8,
        "tail_angular_tolerance": 1e-8,
        "tail_radii": (0.5, 1, 2, 4, 8, 16),
        "tail_angular_orders": ((8, 16), (12, 32)),
        "tail_exterior_radii": (8, 12, 16),
        "tail_max_pair_products": 20_000_000,
        "tail_max_exterior_products": 1_000_000,
        "tail_max_angular_points": 100_000,
    }
    settings.update({key: config[key] for key in settings if key in config})
    for key in ("tail_control_tolerance", "tail_high_m_relative_tolerance",
                "tail_spatial_tolerance", "tail_angular_tolerance"):
        settings[key] = _positive_float(settings[key], key)
    for key in ("tail_max_pair_products", "tail_max_exterior_products", "tail_max_angular_points"):
        settings[key] = _positive_int(settings[key], key)
    precision = settings["tail_decimal_precision"]
    if type(precision) is not int or not 32 <= precision <= 256:
        raise ValueError("tail_decimal_precision must be an integer between 32 and 256")
    for key in ("tail_identity_radii", "tail_radii", "tail_exterior_radii"):
        settings[key] = list(_radii(settings[key], key))
    grids = tuple(tuple(pair) for pair in settings["tail_angular_orders"])
    if len(grids) < 2 or len(grids) > 16:
        raise ValueError("tail_angular_orders needs between two and 16 refinements")
    if any(len(pair) != 2 or any(type(order) is not int for order in pair)
           or pair[0] < 7 or pair[1] < 13 for pair in grids):
        raise ValueError("Angular j<=12 requires n_mu>=7 and n_phi>=13")
    if any(hi[0] < lo[0] or hi[1] < lo[1] or hi == lo for lo, hi in zip(grids, grids[1:])):
        raise ValueError("Angular grids must refine without decreasing either order")
    points = len(settings["tail_radii"])*sum(n_mu*n_phi for n_mu, n_phi in grids)
    if points > settings["tail_max_angular_points"]:
        raise ValueError("Planned full-sphere controls exceed tail_max_angular_points")
    settings["tail_angular_orders"] = [list(pair) for pair in grids]
    settings["planned_gravitational_angular_points"] = points
    settings["angular_bandlimit_assumed"] = 12
    settings["scaled_error_convention"] = "abs(a-b)/max(1,abs(a),abs(b))"
    settings["monopole_normalization"] = "gaussian_inner removes pi^(3/2)"
    return settings


def _bins(coefficients):
    return [{"key": list(key), "coefficient": str(value)}
            for key, value in sorted(coefficients.items())]


class _Checks:
    def __init__(self):
        self.rows = []

    def record(self, name, passed, **details):
        row = {"name": name, "passed": bool(passed), **details}
        self.rows.append(row)
        return row

    def exact(self, name, actual, expected):
        return self.record(name, actual == expected,
                           actual=_bins(actual), expected=_bins(expected))

    def near(self, name, actual, expected, tolerance, **details):
        if not isfinite(actual) or not isfinite(expected):
            return self.record(name, False, actual=str(actual), expected=str(expected),
                               error="Nonfinite comparison", **details)
        error = abs(actual-expected)
        scaled = error/max(1.0, abs(actual), abs(expected))
        return self.record(name, scaled <= tolerance, actual=actual, expected=expected,
                           absolute_error=error, scaled_error=scaled,
                           tolerance=tolerance, **details)

    def raises(self, name, function, exception, required_text=None):
        try:
            function()
        except exception as exc:
            return self.record(name, required_text is None or required_text.lower() in str(exc).lower(),
                               exception=type(exc).__name__, message=str(exc))
        except Exception as exc:
            return self.record(name, False, exception=type(exc).__name__, message=str(exc))
        return self.record(name, False, error="Expected exception was not raised")

    def result(self, settings, **details):
        return {
            "passed": all(row["passed"] for row in self.rows),
            "check_count": len(self.rows),
            "failed_checks": [row["name"] for row in self.rows if not row["passed"]],
            "maximum_scaled_error": max((row.get("scaled_error", 0.0) for row in self.rows), default=0.0),
            "checks": self.rows, "settings": settings,
            "accuracy_certified": False, **details,
        }


def _decimal_pi(precision):
    """Independent Machin formula; no tail recurrence or binary64 pi input."""
    with localcontext() as ctx:
        ctx.prec = precision+12
        threshold = D(10)**(-precision-6)

        def inverse_atan(denominator):
            x = D(1)/D(denominator)
            x2 = x*x
            power = x
            value = x
            for index in range(1, 4*precision+64):
                power *= -x2
                term = power/D(2*index+1)
                value += term
                if abs(term) < threshold:
                    return value
            raise ArithmeticError("Independent Decimal arctangent did not converge")

        value = 16*inverse_atan(5)-4*inverse_atan(239)
        ctx.prec = precision
        return +value


def run_tail_controls(config):
    """Synthetic exact identities and numerical controls, executed only by caller."""
    settings = _settings(config)
    checks = _Checks()
    tolerance = settings["tail_control_tolerance"]
    product_cap = settings["tail_max_pair_products"]
    exterior_cap = settings["tail_max_exterior_products"]
    known_moments = (
        {(1, -1): F(1, 2)},
        {(1, -3): F(1, 4), (0, -2): -F(1, 2)},
        {(1, -5): F(3, 8), (0, -4): -F(3, 4), (0, -2): -F(1, 2)},
    )
    for moment, expected in enumerate(known_moments):
        checks.exact("moment_recurrence_m"+str(moment),
                     dict(moment_exterior_coefficients(moment)), expected)

    simple = TripleHeatTail({(0, 0, 0, 1): F(1)})
    cancellation = TripleHeatTail({(2, 0, 0, 1): F(1), (0, 0, 0, 0): -F(1, 2)})
    derivative_tail = TripleHeatTail({(0, 0, 0, 1): F(3),
                                     (2, 0, 1, 1): F(-4), (2, 0, 0, 2): F(-2)})
    slow = TripleHeatTail({(0, 0, 0, 0): F(1)})
    checks.exact("simple_exterior_coefficients", dict(simple.exterior_coefficients(exterior_cap)),
                 {(3, -3): F(1, 16), (2, -2): -F(1, 8)})
    checks.exact("cancelled_exterior_coefficients", dict(cancellation.exterior_coefficients(exterior_cap)),
                 {(2, 0): -F(1, 8)})
    checks.exact("derivative_tail_exterior_coefficients", dict(derivative_tail.exterior_coefficients(exterior_cap)),
                 {(3, -3): -F(1, 8), (2, -2): F(1, 2), (2, 0): F(1, 4), (1, -1): -F(1, 2)})
    checks.raises("uncancelled_I0_cubed_falloff", lambda: slow.exterior_integral(
        2, max_products=exterior_cap), ValueError, "falloff")
    exterior_controls = []
    for radius in (1.0, 2.0, 4.0):
        simple_summary = simple.exterior_integral(radius, max_products=exterior_cap)
        checks.near("simple_algebraic_integral_R"+str(radius),
                    simple_summary["algebraic_integral"], pi**2.5/(8*radius*radius), tolerance)
        summary = cancellation.exterior_integral(radius, max_products=exterior_cap)
        complement = erfc(radius)
        full_exterior = -pi**2.5/12*complement*(3-3*complement+complement*complement)
        difference = abs(full_exterior-summary["algebraic_integral"])
        bound = summary["absolute_remainder_bound"]
        checks.record("cancelled_full_integral_bound_R"+str(radius),
                      summary["algebraic_integral"] == 0 and bound > 0
                      and difference <= bound*(1+tolerance),
                      exact_full_exterior_float=full_exterior, algebraic_integral=summary["algebraic_integral"],
                      absolute_difference=difference, absolute_remainder_bound=bound,
                      bound_ratio=difference/bound if bound else None,
                      relative_rounding_slack=tolerance)
        # This density is 4*pi*d[r**3*I0(r**2)**2*I1(r**2)]/dr.
        # Its full exterior has a nonzero A[3,-3], testing the erf**3
        # replacement in the remainder bound without quadrature or subtraction
        # of two nearly equal algebraic/full tails.
        derivative_summary = derivative_tail.exterior_integral(radius, max_products=exterior_cap)
        derivative_algebraic = -pi**2.5/(4*radius*radius)
        exponential = exp(-radius*radius)
        erf_value = 1-complement
        stable_remainder = (pi**2.5/(4*radius*radius)*complement
                            *(3-3*complement+complement*complement)
                            + pi*pi/(2*radius)*exponential*erf_value*erf_value)
        i0_closed = pi**0.5*erf_value/(2*radius)
        i1_closed = pi**0.5*erf_value/(4*radius**3)-exponential/(2*radius*radius)
        derivative_full = -4*pi*radius**3*i0_closed*i0_closed*i1_closed
        checks.near("derivative_algebraic_integral_R"+str(radius),
                    derivative_summary["algebraic_integral"], derivative_algebraic, tolerance)
        checks.near("derivative_closed_full_integral_R"+str(radius),
                    derivative_full, derivative_algebraic+stable_remainder, tolerance)
        derivative_bound = derivative_summary["absolute_remainder_bound"]
        checks.record("derivative_nonzero_A3_remainder_bound_R"+str(radius),
                      stable_remainder > 0 and derivative_bound > 0
                      and stable_remainder <= derivative_bound*(1+tolerance),
                      stable_positive_remainder=stable_remainder,
                      absolute_remainder_bound=derivative_bound,
                      bound_ratio=stable_remainder/derivative_bound if derivative_bound else None,
                      relative_rounding_slack=tolerance)
        exterior_controls.append({"radius": radius, "simple": simple_summary,
                                  "cancellation": summary, "closed_full_exterior": full_exterior,
                                  "derivative_tail": derivative_summary,
                                  "derivative_closed_full_exterior": derivative_full,
                                  "derivative_stable_positive_remainder": stable_remainder})

    for name, tail in (("simple", simple), ("cancelled", cancellation),
                       ("derivative", derivative_tail), ("slow", slow)):
        for radius in settings["tail_identity_radii"]:
            checks.near(name+"_radial_exterior_R"+str(radius), tail.evaluate_radial(radius),
                        tail.evaluate_exterior(radius, exterior_cap), tolerance)
    for radius in (1e200, 1e-200):
        for method_name in ("evaluate_radial", "evaluate_exterior", "exterior_integral"):
            checks.raises(method_name+"_radius_"+str(radius),
                          lambda method_name=method_name, radius=radius: getattr(simple, method_name)(radius),
                          ValueError, "radius" if method_name == "evaluate_radial" else None)

    x0 = HeatField(F(1), {((1, 0, 0), 0): F(1)})
    y1 = HeatField(F(1), {((0, 1, 0), 1): F(1)})
    xy2 = HeatField(F(1), {((1, 1, 0), 2): F(1)})
    one1 = HeatField(F(1), {((0, 0, 0), 1): F(1)})
    collapse_statistics = []
    for label, fields, expected in (
        ("distinct", (x0, y1, xy2), {(4, 0, 1, 2): F(1, 15)}),
        ("identical", (x0, x0, one1), {(2, 0, 0, 1): F(1, 3)}),
    ):
        for index, order in enumerate(permutations(range(3))):
            tail, statistics = collapse_triple_heat(tuple(fields[i] for i in order),
                                                   max_pair_products=product_cap)
            checks.exact(label+"_permutation_"+str(index), dict(tail.coefficients), expected)
            collapse_statistics.append({"case": label, "permutation": list(order), **statistics})
        checks.exact(label+"_signed_merge_zero", dict(tail.add(tail, -1).coefficients), {})
        checks.exact(label+"_signed_scale_merge", dict(tail.scale(-3).add(tail, 2).coefficients),
                     {key: -value for key, value in expected.items()})
    # Unlike one-term fields, this forces the identical-pair cross term's
    # factor two to survive the angular average: H**3=(I0+x**2*I1)**3.
    mixed_identical = HeatField(F(1), {((0, 0, 0), 0): F(1), ((2, 0, 0), 1): F(1)})
    mixed_tail, mixed_statistics = collapse_triple_heat((mixed_identical,)*3,
                                                       max_pair_products=product_cap)
    checks.exact("identical_pair_surviving_cross_terms", dict(mixed_tail.coefficients),
                 {(0, 0, 0, 0): F(1), (2, 0, 0, 1): F(1),
                  (4, 0, 1, 1): F(3, 5), (6, 1, 1, 1): F(1, 7)})
    checks.record("identical_pair_optimization_exercised", mixed_statistics.get("identical_pair") is True,
                  statistics=mixed_statistics)
    collapse_statistics.append({"case": "identical_pair_surviving_cross_terms", **mixed_statistics})
    checks.raises("triple_product_cap", lambda: collapse_triple_heat((x0, y1, xy2),
                  max_pair_products=1), RuntimeError, "budget")
    checks.raises("exterior_product_cap", lambda: simple.exterior_coefficients(max_products=1),
                  RuntimeError, "budget")
    other_rate = HeatField(F(2), {((0, 0, 0), 0): F(1)})
    checks.raises("nonunit_heat_rate", lambda: collapse_triple_heat((x0, y1, other_rate),
                  max_pair_products=product_cap), ValueError, "rates")

    # (179)!! = 180!/(2**90*90!), independent of the moment recurrence.
    # A[3,-181] = (179)!!/2**93, and integrating r**(-181) gives R**(-180)/180.
    moment = 90
    high = TripleHeatTail({(0, 0, 0, moment): F(1)})
    expected_coefficient = F(factorial(180), 2**183*factorial(90))
    actual_coefficient = high.exterior_coefficients(exterior_cap).get((3, -181), F(0))
    checks.record("high_m_exact_factorial_coefficient", actual_coefficient == expected_coefficient,
                  actual=str(actual_coefficient), expected=str(expected_coefficient))
    summary = high.exterior_integral(100, max_products=exterior_cap)
    normalized = F(factorial(180), 2**181*factorial(90)*180)*F(100)**(-180)
    precision = settings["tail_decimal_precision"]
    with localcontext() as ctx:
        ctx.prec = precision
        decimal_pi = _decimal_pi(precision)
        target = D(normalized.numerator)/D(normalized.denominator)*decimal_pi**2*decimal_pi.sqrt()
        actual = D.from_float(summary["algebraic_integral"])
        relative_error = abs(actual-target)/abs(target)
        high_tolerance = D(str(settings["tail_high_m_relative_tolerance"]))
        checks.record("high_m_nonzero_algebraic_integral", actual > 0 and target > 0
                      and relative_error <= high_tolerance,
                      actual_decimal=str(actual), target_decimal=str(target),
                      normalized_factorial_target=str(normalized), relative_error_decimal=str(relative_error),
                      tolerance_decimal=str(high_tolerance), target_method="Closed factorial Fraction and Decimal Machin pi")
    return checks.result(settings, exterior_controls=exterior_controls,
                         collapse_statistics=collapse_statistics, high_m_exterior=summary,
                         limitations=["Numerical bound evaluation is not interval certified",
                                      "Radial/exterior point identity is distinct from an integrated exterior error",
                                      "High-m target uses a closed factorial formula independent of tail recurrence"])


def check_gravitational_tail(data, name, config, progress=None):
    """Check supplied -w*T:T, w=7u+v, without constructing geometry or sources."""
    settings = _settings(config)
    if not isinstance(data, Mapping):
        raise TypeError("data must contain supplied first-response fields")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable or None")
    name = str(name)
    checks = _Checks()

    def notify(record):
        if progress is not None:
            progress({"name": name, **record})

    notify({"stage": "tail_assembly_start"})
    tail, component_statistics = triple_tail(data, max_pair_products=settings["tail_max_pair_products"],
                                             progress=notify)
    exterior = tail.exterior_coefficients(settings["tail_max_exterior_products"])
    monopoles = tuple(gaussian_inner(poly, ONE) for poly in data["york"]["b"][1])
    if len(monopoles) != 3:
        raise ValueError("The supplied b1 must have three components")
    scalar_source = add(scale(data["s2"], 14), data["geometry"]["R"][2], -6)
    scalar_monopole = gaussian_inner(scalar_source, ONE)
    predicted = -F(15, 4096)*scalar_monopole*sum((value*value for value in monopoles), F(0))
    observed = exterior.get((3, -3), F(0))
    checks.record("monopole_leading_coefficient", observed == predicted,
                  actual=str(observed), expected=str(predicted),
                  formula="-15*Mhbar*sum(Mbar_i**2)/4096")
    bad = {key: value for key, value in exterior.items() if key[0] == 3 and key[1] > -3}
    checks.exact("exact_algebraic_falloff", bad, {})
    if name == "monopole_mixed":
        expected_vector = (F(0), F(2079, 32), F(0))
        checks.record("named_monopole_vector", monopoles == expected_vector,
                      actual=[str(value) for value in monopoles], expected=[str(value) for value in expected_vector])
        checks.record("named_scalar_monopole", scalar_monopole == F(251559, 64),
                      actual=str(scalar_monopole), expected=str(F(251559, 64)))
        checks.record("named_nonzero_leading_tail", observed != 0, coefficient=str(observed))
    notify({"stage": "tail_exact_monopoles", "input_bins": len(tail.coefficients),
            "exterior_bins": len(exterior), "leading_coefficient": str(observed)})

    compiler = _CompiledFields()
    scalar_w = data["u"].scale(7).add(data["v"])
    w_ref = compiler.heat(scalar_w)
    tensor_refs = []
    for i in range(3):
        for j in range(i, 3):
            tensor_refs.append((compiler.heat(data["t0"][i][j]), 1 if i == j else 2))
            notify({"stage": "tail_sphere_compile", "component": [i, j]})
    angular_grids = [compiler.angular_grid(n_mu, n_phi) for n_mu, n_phi in settings["tail_angular_orders"]]
    radial_rows = []
    for radius in settings["tail_radii"]:
        collapsed = tail.evaluate_radial(radius)
        compiled = compiler.at_radius(radius)
        levels = []
        for index, (orders, grid) in enumerate(zip(settings["tail_angular_orders"], angular_grids)):
            contributions = []
            for weight, monomials in grid:
                values = compiler.evaluate(compiled, monomials)
                density = -values[w_ref]*fsum(multiplicity*values[reference]**2
                                            for reference, multiplicity in tensor_refs)
                contributions.append(radius*radius*weight*density)
            sphere = fsum(contributions)
            absolute = fsum(abs(value) for value in contributions)
            comparison = checks.near("sphere_vs_collapse_R"+str(radius)+"_level"+str(index+1),
                                     sphere, collapsed, settings["tail_spatial_tolerance"],
                                     radius=radius, angular_orders=orders)
            refinement = None
            if levels:
                refinement = checks.near("angular_refinement_R"+str(radius)+"_level"+str(index+1),
                                         sphere, levels[-1]["radial_density"], settings["tail_angular_tolerance"],
                                         radius=radius, previous_orders=levels[-1]["angular_orders"], angular_orders=orders)
            ratio = absolute/abs(sphere) if sphere else None
            levels.append({"angular_orders": orders, "nodes": len(grid), "radial_density": sphere,
                           "absolute_node_sum": absolute,
                           "cancellation_ratio": ratio if ratio is None or isfinite(ratio) else None,
                           "cancellation_ratio_overflow": ratio is not None and not isfinite(ratio),
                           "zero_sum_with_nonzero_terms": sphere == 0 and absolute != 0,
                           "collapse_comparison": comparison, "angular_refinement": refinement})
            notify({"stage": "tail_sphere_radius", "radius": radius, "level": index+1,
                    "angular_orders": orders, "radial_density": sphere, "scaled_error": comparison.get("scaled_error")})
        radial_rows.append({"radius": radius, "collapsed_radial_density": collapsed, "levels": levels})

    summaries = []
    if not bad:
        for radius in settings["tail_exterior_radii"]:
            summary = tail.exterior_integral(radius, required_decay_power=-3,
                                             max_products=settings["tail_max_exterior_products"])
            summaries.append(summary)
            checks.record("finite_exterior_summary_R"+str(radius),
                          isfinite(summary["algebraic_integral"])
                          and isfinite(summary["absolute_remainder_bound"])
                          and summary["absolute_remainder_bound"] >= 0,
                          radius=radius)
            if name == "monopole_mixed":
                checks.record("named_nonzero_exterior_R"+str(radius),
                              summary["algebraic_integral"] != 0, radius=radius,
                              algebraic_integral=summary["algebraic_integral"])
            notify({"stage": "tail_exterior_summary", "radius": radius,
                    "algebraic_integral": summary["algebraic_integral"],
                    "absolute_remainder_bound": summary["absolute_remainder_bound"]})
    return checks.result(settings, name=name,
                         monopoles_normalized=[str(value) for value in monopoles],
                         scalar_monopole_normalized=str(scalar_monopole),
                         leading_coefficient=str(observed), predicted_leading_coefficient=str(predicted),
                         input_bins=len(tail.coefficients), exterior_coefficients=_bins(exterior),
                         component_statistics=list(component_statistics),
                         sphere_statistics=compiler.statistics(), radial_checks=radial_rows,
                         exterior_summaries=summaries, exterior_skipped_for_bad_falloff=bool(bad),
                         limitations=["Supplied data must be the intended j<=12 total-density case",
                                      "Sphere and collapsed radial evaluations share the HEAT-moment evaluator",
                                      "No independent geometry or Green fields were constructed",
                                      "Pointwise density discrepancies are not compared with integrated tail bounds",
                                      "Exterior remainder inequality is evaluated in binary64, not interval certified"])
