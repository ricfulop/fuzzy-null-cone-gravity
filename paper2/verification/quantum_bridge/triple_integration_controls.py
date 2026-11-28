"""Tracked controls for the full Decimal triple-HEAT integrator.

run_triple_integration_controls(config, progress=None) -> JSON-compatible result.
All scientific work is inside that call or its helpers; import is inert.
No gravitational geometry, sources, contractions, files or runs are created.

Full-integral controls use two small exact bin dictionaries, two independently
chosen exterior radii, and the SAME tolerances at both radii. The normalized
quantity is physical/(4*pi); the cancellation target is -pi**(3/2)/48 and
the derivative target is exactly zero. Physical targets are also checked.
Their pi is computed with Decimal AGM, independently of production Machin pi.

Moment controls call the production private API
  _heat_moments(Decimal_z, maximum_moment, precision,
                _Work(max_heat_steps, max_series_terms))
and compare with green_field_check._decimal_reference's independent positive
Gaussian integration. Identical exact binary64 arguments are represented by
Decimal.from_float on the production side. Working/check precisions are fixed
at 64/96; the reference uses 96. Quadrature convergence and omitted reference
tails are explicitly gated. This is numerical validation, not certification.
"""

from collections.abc import Mapping
from decimal import Decimal as D, InvalidOperation, localcontext
from fractions import Fraction as F

from algebraic_heat_tail import TripleHeatTail
from green_field_check import _decimal_reference
from triple_heat_integration import _Work, _heat_moments, integrate_triple_heat


BASE_PRECISION = 64
CHECK_PRECISION = 96
TARGET_PRECISION = 112


def _number(value, name, allow_zero=False):
    if isinstance(value, bool):
        raise TypeError(name+" must be a finite decimal number")
    try:
        result = value if isinstance(value, D) else D(str(value))
    except (InvalidOperation, TypeError, ValueError) as exc:
        raise ValueError(name+" must be a finite decimal number") from exc
    if not result.is_finite() or result < 0 or (result == 0 and not allow_zero):
        raise ValueError(name+" must be finite and "+("nonnegative" if allow_zero else "positive"))
    return result


def _settings(config):
    if not isinstance(config, Mapping):
        raise TypeError("config must be a mapping")
    defaults = {
        "triple_control_radii": (10, 12),
        "triple_control_panel_width": 2,
        "triple_control_orders": (16, 24, 32),
        "triple_control_absolute_tolerance": "1e-10",
        "triple_control_relative_tolerance": "1e-9",
        "triple_control_reference_orders": (24, 48),
        "triple_control_moment_relative_tolerance": "1e-24",
        "triple_control_reference_relative_tolerance": "1e-24",
        "triple_control_reference_tail_relative_tolerance": "1e-24",
        "triple_control_max_term_evaluations": 20_000_000,
        "triple_control_max_heat_steps": 5_000_000,
        "triple_control_max_series_terms": 4096,
        "triple_control_max_exterior_products": 1_000_000,
    }
    settings = {key: config.get(key, value) for key, value in defaults.items()}
    for key in ("triple_control_panel_width", "triple_control_absolute_tolerance",
                "triple_control_moment_relative_tolerance", "triple_control_reference_relative_tolerance",
                "triple_control_reference_tail_relative_tolerance"):
        settings[key] = str(_number(settings[key], key))
    key = "triple_control_relative_tolerance"
    settings[key] = str(_number(settings[key], key, allow_zero=True))
    radii = tuple(_number(radius, "triple_control_radii") for radius in settings["triple_control_radii"])
    if len(radii) != 2 or radii[0] >= radii[1]:
        raise ValueError("triple_control_radii must contain exactly two increasing positive radii")
    settings["triple_control_radii"] = [str(radius) for radius in radii]
    for key, maximum_length, maximum_order in (("triple_control_orders", 4, 128),
                                               ("triple_control_reference_orders", 3, 96)):
        orders = tuple(settings[key])
        if (not 2 <= len(orders) <= maximum_length
                or any(type(order) is not int or not 8 <= order <= maximum_order for order in orders)
                or any(high <= low for low, high in zip(orders, orders[1:]))):
            raise ValueError(key+" must be a bounded strictly increasing order list")
        settings[key] = list(orders)
    for key in ("triple_control_max_term_evaluations", "triple_control_max_heat_steps",
                "triple_control_max_series_terms", "triple_control_max_exterior_products"):
        if type(settings[key]) is not int or settings[key] < 1:
            raise ValueError(key+" must be a positive integer")
    settings.update({"precision": BASE_PRECISION, "check_precision": CHECK_PRECISION,
                     "reference_precision": CHECK_PRECISION, "target_precision": TARGET_PRECISION,
                     "ABSOLUTE_TOLERANCE_UNITS": "FULL physical integral divided by 4*pi",
                     "target_pi_method": "Independent Decimal arithmetic-geometric mean",
                     "radius_refinement_rule": "Same atol/rtol and target-based budget at both radii; no relaxation"})
    return settings


def _agm_pi(precision):
    """Gauss-Legendre AGM pi, separate from the integrator's Machin series."""
    with localcontext() as ctx:
        ctx.prec = precision+16
        a, b, t, multiplier = D(1), D(1)/D(2).sqrt(), D(1)/4, D(1)
        tolerance = D(10)**(-precision-8)
        for _ in range(20):
            average = (a+b)/2
            geometric = (a*b).sqrt()
            t -= multiplier*(a-average)**2
            a, b, multiplier = average, geometric, 2*multiplier
            if abs(a-b) < tolerance:
                value = (a+b)**2/(4*t)
                ctx.prec = precision
                return +value
        raise ArithmeticError("Independent AGM pi did not converge")


def _origin_control(settings):
    rows = []
    for precision in (BASE_PRECISION, CHECK_PRECISION):
        work = _Work(settings["triple_control_max_heat_steps"], settings["triple_control_max_series_terms"])
        values, remainder, seed_terms = _heat_moments(D(0), 90, precision, work)
        with localcontext() as ctx:
            # The helper documents twelve guard digits for moments. Compare
            # its analytic origin branch at that exact arithmetic precision.
            ctx.prec = precision+12
            differences = tuple(abs(value-D(1)/(2*m+1)) for m, value in enumerate(values))
            rows.append({"precision": precision, "moment_range": [0, 90],
                         "all_origin_values_exact_at_declared_precision": not any(differences),
                         "maximum_absolute_difference_decimal": str(max(differences)),
                         "seed_relative_tail_decimal": str(remainder), "seed_terms": seed_terms,
                         "work": work.statistics(),
                         "passed": not any(differences) and remainder == 0 and seed_terms == 0})
    return {"name": "origin_closed_moments", "passed": all(row["passed"] for row in rows), "levels": rows}


def _moment_control(z_float, maximum, selected, settings):
    estimates = []
    for precision in (BASE_PRECISION, CHECK_PRECISION):
        work = _Work(settings["triple_control_max_heat_steps"], settings["triple_control_max_series_terms"])
        moments, tail, terms = _heat_moments(D.from_float(z_float), maximum, precision, work)
        estimates.append((moments, tail, terms, work.statistics()))
    reference_orders = settings["triple_control_reference_orders"]
    references = [_decimal_reference(z_float, maximum, order, CHECK_PRECISION) for order in reference_orders]
    comparisons = []
    with localcontext() as ctx:
        ctx.prec = TARGET_PRECISION
        moment_tolerance = D(settings["triple_control_moment_relative_tolerance"])
        refinement_tolerance = D(settings["triple_control_reference_relative_tolerance"])
        tail_tolerance = D(settings["triple_control_reference_tail_relative_tolerance"])
        for m in selected:
            reference = references[-1][0][m]
            reference_tail = references[-1][1][m]
            if not reference.is_finite() or reference <= 0:
                raise ArithmeticError("Independent moment reference is not finite and positive")
            denominator = abs(reference)
            changes = [abs(new[0][m]-old[0][m])/denominator for old, new in zip(references, references[1:])]
            errors = [(abs(level[0][m]-reference)+reference_tail)/denominator for level in estimates]
            precision_change = abs(estimates[1][0][m]-estimates[0][0][m])/denominator
            relative_tail = reference_tail/denominator
            seed_bounds_ok = all(level[1] <= D(10)**(-precision-6)
                                 for level, precision in zip(estimates, (BASE_PRECISION, CHECK_PRECISION)))
            passed = (all(error <= moment_tolerance for error in errors)
                      and precision_change <= moment_tolerance and changes[-1] <= refinement_tolerance
                      and relative_tail <= tail_tolerance and seed_bounds_ok)
            comparisons.append({
                "m": m, "passed": passed,
                "working_value_decimal": str(estimates[0][0][m]),
                "checked_value_decimal": str(estimates[1][0][m]),
                "reference_values_decimal": [str(level[0][m]) for level in references],
                "reference_absolute_tail_bounds_decimal": [str(level[1][m]) for level in references],
                "comparison_relative_indicators_decimal": [str(error) for error in errors],
                "precision_relative_change_decimal": str(precision_change),
                "reference_relative_changes_decimal": [str(change) for change in changes],
                "reference_relative_tail_bound_decimal": str(relative_tail),
                "seed_relative_bounds_obey_stop_criterion": seed_bounds_ok,
            })
    return {"name": "moments_z="+str(z_float)+"_mmax="+str(maximum),
            "z_exact_decimal": str(D.from_float(z_float)), "maximum_moment": maximum,
            "selected_moments": list(selected), "passed": all(row["passed"] for row in comparisons),
            "reference_orders": reference_orders, "reference_panels": [row[2] for row in references],
            "seed_relative_tail_decimals": [str(row[1]) for row in estimates],
            "seed_term_counts": [row[2] for row in estimates], "work": [row[3] for row in estimates],
            "comparisons": comparisons}


def _series_cap_control(settings):
    work = _Work(settings["triple_control_max_heat_steps"], 1)
    try:
        _heat_moments(D(2), 3, BASE_PRECISION, work)
    except ValueError as exc:
        return {"name": "kummer_series_cap_rejection", "passed": "max_series_terms" in str(exc),
                "exception": type(exc).__name__, "message": str(exc), "work": work.statistics()}
    return {"name": "kummer_series_cap_rejection", "passed": False,
            "message": "Expected max_series_terms rejection did not occur", "work": work.statistics()}


def _integral_control(name, bins, normalized_target, physical_target, target_pi, settings, progress):
    tail = TripleHeatTail(bins)
    rows = []
    with localcontext() as ctx:
        ctx.prec = TARGET_PRECISION
        fixed_budget = D(settings["triple_control_absolute_tolerance"])+D(
            settings["triple_control_relative_tolerance"])*abs(normalized_target)
        physical_budget = 4*target_pi*fixed_budget
    for radius in settings["triple_control_radii"]:
        def notify(record):
            if progress is not None:
                progress({"control": name, "control_radius": radius, **record})

        try:
            result = integrate_triple_heat(
                tail, radius=radius, panel_width=settings["triple_control_panel_width"],
                orders=settings["triple_control_orders"], precision=BASE_PRECISION, check_precision=CHECK_PRECISION,
                absolute_tolerance=settings["triple_control_absolute_tolerance"],
                relative_tolerance=settings["triple_control_relative_tolerance"],
                max_term_evaluations=settings["triple_control_max_term_evaluations"],
                max_heat_steps=settings["triple_control_max_heat_steps"],
                max_series_terms=settings["triple_control_max_series_terms"],
                max_exterior_products=settings["triple_control_max_exterior_products"], progress=notify)
        except (ArithmeticError, ValueError, TypeError, RuntimeError, AssertionError) as exc:
            rows.append({"radius": radius, "passed": False,
                         "exception": type(exc).__name__, "message": str(exc)})
            notify({"stage": "triple_control_radius_complete", "radius": radius,
                    "passed": rows[-1]["passed"], "radius_result": rows[-1]})
            continue
        with localcontext() as ctx:
            ctx.prec = TARGET_PRECISION
            value = D(result["normalized_value_decimal"])
            physical = D(result["physical_value_decimal"])
            error, physical_error = abs(value-normalized_target), abs(physical-physical_target)
            indicator = D(result["combined_indicator_normalized_decimal"])
            rows.append({"radius": radius,
                         "passed": result["converged"] and error <= fixed_budget and physical_error <= physical_budget,
                         "absolute_normalized_error_decimal": str(error),
                         "absolute_physical_error_decimal": str(physical_error),
                         "fixed_normalized_budget_decimal": str(fixed_budget),
                         "fixed_physical_budget_decimal": str(physical_budget),
                         "observed_error_within_numerical_indicator": error <= indicator,
                         "result": result})
        notify({"stage": "triple_control_radius_complete", "radius": radius,
                "passed": rows[-1]["passed"], "radius_result": rows[-1]})
    if all("result" in row for row in rows):
        with localcontext() as ctx:
            ctx.prec = TARGET_PRECISION
            delta = abs(D(rows[1]["result"]["normalized_value_decimal"])-D(rows[0]["result"]["normalized_value_decimal"]))
            physical_delta = abs(D(rows[1]["result"]["physical_value_decimal"])-D(rows[0]["result"]["physical_value_decimal"]))
            radius_passed = delta <= fixed_budget and physical_delta <= physical_budget
        radius_refinement = {"passed": radius_passed,
                             "absolute_normalized_change_decimal": str(delta),
                             "absolute_physical_change_decimal": str(physical_delta),
                             "same_fixed_budget_at_both_radii": True}
    else:
        radius_passed = False
        radius_refinement = {"passed": False, "skipped": True,
                             "reason": "At least one radius evaluation failed; other radius evidence retained"}
    return {"name": name, "passed": all(row["passed"] for row in rows) and radius_passed,
            "input_bins": [{"bin": list(key), "coefficient": str(coefficient)} for key, coefficient in sorted(bins.items())],
            "normalized_target_decimal": str(normalized_target), "physical_target_decimal": str(physical_target),
            "fixed_normalized_budget_decimal": str(fixed_budget), "fixed_physical_budget_decimal": str(physical_budget),
            "radii": rows, "radius_refinement": radius_refinement}


def run_triple_integration_controls(config, progress=None):
    """Execute small full-integral and independent-moment controls when called.

    Invalid configuration raises. Case exceptions are retained as failed rows
    so one failure does not erase the other control evidence. No run is launched.
    Completed control callbacks include control_result; each completed radius
    callback includes radius_result, including a caught exception row. Both
    payloads are emitted immediately so the driver can checkpoint partial work.
    """
    settings = _settings(config)
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable or None")
    failures = []

    def perform(stage, name, function):
        if progress is not None:
            progress({"stage": stage+"_start", "control": name})
        try:
            row = function()
        except (ArithmeticError, ValueError, TypeError, RuntimeError, AssertionError) as exc:
            row = {"name": name, "passed": False, "exception": type(exc).__name__, "message": str(exc)}
        if not row["passed"]:
            failures.append({"stage": stage, "name": name,
                             "exception": row.get("exception"), "message": row.get("message")})
        if progress is not None:
            progress({"stage": stage+"_complete", "control": name,
                      "passed": row["passed"], "control_result": row})
        return row

    origin = perform("triple_origin_control", "origin_closed_moments", lambda: _origin_control(settings))
    specifications = ((0.0, 3, (0, 1, 3)), (0.125, 3, (0, 1, 3)),
                      (2.0, 17, (0, 2, 7, 17)), (144.0, 90, (0, 1, 17, 90)),
                      (800.0, 90, (0, 1, 30, 90)))
    moments = [perform("triple_moment_control", "z="+str(z)+"_mmax="+str(maximum),
                       lambda z=z, maximum=maximum, selected=selected: _moment_control(z, maximum, selected, settings))
               for z, maximum, selected in specifications]
    series_cap = perform("triple_series_cap", "kummer_series_cap_rejection", lambda: _series_cap_control(settings))
    target_pi = _agm_pi(TARGET_PRECISION)
    with localcontext() as ctx:
        ctx.prec = TARGET_PRECISION
        normalized_target = -target_pi*target_pi.sqrt()/48
        physical_target = -target_pi*target_pi*target_pi.sqrt()/12
    cancellation = {(2, 0, 0, 1): F(1), (0, 0, 0, 0): -F(1, 2)}
    derivative = {(0, 0, 0, 1): F(3), (2, 0, 1, 1): F(-4), (2, 0, 0, 2): F(-2)}
    integrals = [
        perform("triple_full_integral", "closed_negative_integral", lambda: _integral_control(
            "closed_negative_integral", cancellation, normalized_target, physical_target, target_pi, settings, progress)),
        perform("triple_full_integral", "exact_zero_derivative_integral", lambda: _integral_control(
            "exact_zero_derivative_integral", derivative, D(0), D(0), target_pi, settings, progress)),
    ]
    return {"passed": not failures, "failed_controls": failures,
            "control_count": 2+len(moments)+len(integrals),
            "settings": settings, "origin": origin, "moment_controls": moments,
            "series_cap_control": series_cap, "full_integral_controls": integrals,
            "independent_agm_pi_decimal": str(target_pi), "accuracy_certified": False,
            "limitations": ["Full integral errors and radius changes use fixed physical/(4*pi) budgets",
                            "Observed numerical-indicator coverage is diagnostic, not a certified error bound",
                            "Independent reference moments use Gaussian quadrature rather than the Kummer recurrence",
                            "Reference quadrature refinement is diagnostic; its omitted exterior tail has a separate bound",
                            "No gravitational fields, geometry or complete H6 energy are tested here"]}
