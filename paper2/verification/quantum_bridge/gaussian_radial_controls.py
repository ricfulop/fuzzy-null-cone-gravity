"""Tracked public-API controls for independent Gaussian radial companions.

run_gaussian_radial_controls(config, progress=None) -> JSON-compatible result.
No FD solve is imported or called. The one wrapper attachment check constructs
ScalarSourceSpec only. Import performs no source construction or evaluation.

Moment references are independently integrated at 96 decimal digits with
composite Gauss--Legendre orders24/48 and a positive Gaussian exterior bound.
run.gauss_nodes returns HALF the usual [-1,1] weights, so a panel uses its
FULL width. Refinement is an error diagnostic, not a rigorous quadrature
remainder theorem. Provider Decimal interval endpoints, rather than binary64
display values, enter the 1e-22 relative gates. Source corrections use their
separate 2e-12 binary64 gate. No total accuracy certification is claimed.

Complete control_result rows are emitted at gaussian_radial_control_completed
before advancing, or at gaussian_radial_control_exception before rethrowing.
Completed moment/radius/channel rows are also emitted immediately as
level_result at gaussian_radial_control_level_completed; they do not add
control groups. An exception control_result retains its completed levels.
"""

from collections.abc import Mapping
from decimal import Decimal as D, localcontext, ROUND_CEILING
from fractions import Fraction as F
from math import factorial, isfinite

from run import gauss_nodes
from gaussian_radial_reference import GaussianRadialSource
from spherical_green_reference import ExteriorSource, OriginSeries
from spherical_york_reference import ScalarSourceSpec


REFERENCE_PRECISION = 96
MOMENT_POWERS = (0, 1, 2, 3, 8, 20, 40, 68)
MOMENT_CASES = (("1", "0"), ("1", "1/8"), ("3/2", "1"),
                ("1", "4"), ("1", "8"), ("1", "12"))


def _decimal(value):
    value = F(value)
    return D(value.numerator)/D(value.denominator)


def _power(value, exponent):
    return D(1) if exponent == 0 else value**exponent


def _settings(config):
    if not isinstance(config, Mapping):
        raise TypeError("config must be a mapping")
    defaults = {
        "gaussian_radial_control_provider_precision": 64,
        "gaussian_radial_control_orders": (24, 48),
        "gaussian_radial_control_panel_width": "0.5",
        "gaussian_radial_control_moment_relative_tolerance": "1e-22",
        "gaussian_radial_control_refinement_relative_tolerance": "1e-22",
        "gaussian_radial_control_tail_relative_tolerance": "1e-22",
        "gaussian_radial_control_correction_tolerance": 2e-12,
        "gaussian_radial_control_max_reference_terms": 500_000,
        "gaussian_radial_control_max_reference_panels": 128,
        "gaussian_radial_control_max_upper_steps": 32,
        "gaussian_radial_control_max_source_objects": 32,
        "gaussian_radial_control_max_moment_steps": 100_000,
        "gaussian_radial_control_max_series_terms": 2048,
        "gaussian_radial_control_max_point_evaluations": 20_000,
        "gaussian_radial_control_max_products": 100_000,
        "gaussian_radial_control_max_terms": 10_000,
    }
    out = {key: config.get(key, value) for key, value in defaults.items()}
    for key in ("gaussian_radial_control_panel_width", "gaussian_radial_control_moment_relative_tolerance",
                "gaussian_radial_control_refinement_relative_tolerance", "gaussian_radial_control_tail_relative_tolerance"):
        value = D(str(out[key]))
        if not value.is_finite() or value <= 0:
            raise ValueError(key+" must be finite and positive")
        out[key] = str(value)
    value = out["gaussian_radial_control_correction_tolerance"]
    if isinstance(value, bool) or not isfinite(float(value)) or float(value) <= 0:
        raise ValueError("correction tolerance must be finite and positive")
    out["gaussian_radial_control_correction_tolerance"] = float(value)
    for key in defaults:
        if "max_" in key or key.endswith("precision"):
            if type(out[key]) is not int or out[key] < 1:
                raise ValueError(key+" must be a positive integer")
    if not 32 <= out["gaussian_radial_control_provider_precision"] <= 128:
        raise ValueError("provider precision must lie in32..128")
    orders = tuple(out["gaussian_radial_control_orders"])
    if (len(orders) != 2 or any(type(order) is not int or not 16 <= order <= 96 for order in orders)
            or orders[0] >= orders[1]):
        raise ValueError("reference orders must be two increasing integers in16..96")
    out["gaussian_radial_control_orders"] = list(orders)
    out.update({"reference_precision": REFERENCE_PRECISION, "cases": [list(row) for row in MOMENT_CASES],
                "moment_powers": list(MOMENT_POWERS), "accuracy_certified": False,
                "moment_error_units": "relative to positive independent moment; Decimal enclosure used",
                "correction_error_units": "abs(actual-target)/max(1,abs(target))",
                "reference_tail_formula": "U^(n+1)*exp(-beta*U^2)/(2*beta*U^2-n), denominator positive",
                "reference_tail_arithmetic_safety_factor": "1+1e-80",
                "cutoff_rule": "increase U by1/sqrt(beta) until every analytic tail <= tolerance*positive lower bound/16"})
    return out


def _agm_pi():
    """Independent AGM pi; the provider uses Machin arctangents."""
    a, b, t, factor = D(1), D(1)/D(2).sqrt(), D(1)/4, D(1)
    for _ in range(20):
        new_a, new_b = (a+b)/2, (a*b).sqrt()
        t -= factor*(a-new_a)**2
        a, b, factor = new_a, new_b, 2*factor
        if abs(a-b) < D("1e-88"):
            return (a+b)**2/(4*t)
    raise ArithmeticError("Independent AGM pi did not converge")


class _Reference:
    def __init__(self, settings):
        self.settings = settings
        self.terms = self.nodes = self.tail_evaluations = 0
        self.cache = {}

    def _reserve(self, terms):
        if self.terms+terms > self.settings["gaussian_radial_control_max_reference_terms"]:
            raise ValueError("Independent reference exceeds max_reference_terms")
        self.terms += terms

    def _tail(self, n, beta, upper):
        self.tail_evaluations += 1
        denominator = 2*beta*upper*upper-n
        if upper <= 0 or denominator <= 0:
            raise ValueError("Positive Gaussian tail inequality is outside its domain")
        return upper**(n+1)*(-beta*upper*upper).exp()/denominator*(1+D("1e-80"))

    def _cutoff(self, beta, lower, powers):
        # f_n(r)=r^n exp(-beta*r²) is log-concave on r>0. Its minimum on
        # [a,b] is at an endpoint, hence this is a positive integral lower
        # bound on a subinterval of [lower,infinity). The factor1/2 leaves
        # ample numerical margin at96-digit arithmetic; no provider used.
        positive_lower = {}
        for n in powers:
            a = max(lower, (D(n+1)/(2*beta)).sqrt())
            delta = 1/(1+2*beta*a+D(n)/a)
            b = a+delta
            positive_lower[n] = delta*min(_power(a, n)*(-beta*a*a).exp(),
                                           _power(b, n)*(-beta*b*b).exp())/2
            if positive_lower[n] <= 0:
                raise ArithmeticError("Independent moment lower bound is not positive")
        upper = max(lower+1, D(8)/beta.sqrt())
        target = D(self.settings["gaussian_radial_control_tail_relative_tolerance"])/16
        for step in range(self.settings["gaussian_radial_control_max_upper_steps"]):
            tails = {n: self._tail(n, beta, upper) for n in powers}
            if all(tails[n] <= target*positive_lower[n] for n in powers):
                return upper, tails, step+1, positive_lower
            upper += 1/beta.sqrt()
        raise ValueError("Independent tail selection exceeds max_upper_steps")

    def integrate(self, beta, lower, upper, powers, order):
        if upper < lower:
            raise ValueError("Reference interval endpoints are reversed")
        if upper == lower:
            return {n: D(0) for n in powers}, 0
        width_limit = D(self.settings["gaussian_radial_control_panel_width"])
        panels = int(((upper-lower)/width_limit).to_integral_value(rounding=ROUND_CEILING))
        if panels > self.settings["gaussian_radial_control_max_reference_panels"]:
            raise ValueError("Independent reference exceeds max_reference_panels")
        self._reserve(panels*order*len(powers))
        nodes = gauss_nodes(order, REFERENCE_PRECISION)
        width = (upper-lower)/panels
        totals = {n: D(0) for n in powers}
        for panel in range(panels):
            midpoint = lower+(D(panel)+D("0.5"))*width
            for node, weight in nodes:
                radius = midpoint+width*node/2
                common = width*weight*(-beta*radius*radius).exp()
                for n in powers:
                    totals[n] += common*_power(radius, n)
        self.nodes += panels*order
        return totals, panels

    def moments(self, beta_fraction, radius_fraction, powers):
        powers = tuple(sorted(set(powers)))
        key = F(beta_fraction), F(radius_fraction), powers
        if key in self.cache:
            return self.cache[key]
        beta, radius = _decimal(beta_fraction), _decimal(radius_fraction)
        upper, tails, steps, lower_bounds = self._cutoff(beta, radius, powers)
        levels, panels = [], []
        for order in self.settings["gaussian_radial_control_orders"]:
            values, panel_count = self.integrate(beta, radius, upper, powers, order)
            levels.append(values)
            panels.append(panel_count)
        result = {"values": levels[-1], "coarse": levels[0], "tails": tails,
                  "differences": {n: abs(levels[1][n]-levels[0][n]) for n in powers},
                  "upper": upper, "panels": panels, "cutoff_steps": steps, "positive_lower_bounds": lower_bounds}
        self.cache[key] = result
        return result

    def finite(self, beta, lower, upper, powers):
        levels = [self.integrate(_decimal(beta), _decimal(lower), _decimal(upper), powers, order)[0]
                  for order in self.settings["gaussian_radial_control_orders"]]
        return {"values": levels[1], "coarse": levels[0], "tails": {n: D(0) for n in powers},
                "differences": {n: abs(levels[1][n]-levels[0][n]) for n in powers}}

    def statistics(self):
        return {"reference_term_evaluations": self.terms,
                "max_reference_terms": self.settings["gaussian_radial_control_max_reference_terms"],
                "reference_nodes": self.nodes, "tail_bound_evaluations": self.tail_evaluations,
                "cached_integral_sets": len(self.cache)}


class _Harness:
    def __init__(self, settings, progress):
        self.settings, self.progress = settings, progress
        self.reference = _Reference(settings)
        self.rows, self.sources = [], []
        self.source_attempts = 0
        self.exception_emitted = False
        self.current_control, self.level_rows = None, []

    def source(self, polynomial, beta=F(1), ell_max=0, **overrides):
        if self.source_attempts >= self.settings["gaussian_radial_control_max_source_objects"]:
            raise ValueError("Gaussian controls exceed max_source_objects")
        self.source_attempts += 1
        kwargs = {"precision": self.settings["gaussian_radial_control_provider_precision"],
                  "max_products": self.settings["gaussian_radial_control_max_products"],
                  "max_terms": self.settings["gaussian_radial_control_max_terms"],
                  "max_moment_steps": self.settings["gaussian_radial_control_max_moment_steps"],
                  "max_series_terms": self.settings["gaussian_radial_control_max_series_terms"],
                  "max_point_evaluations": self.settings["gaussian_radial_control_max_point_evaluations"]}
        kwargs.update(overrides)
        source = GaussianRadialSource(polynomial, beta, ell_max=ell_max, **kwargs)
        self.sources.append(source)
        return source

    def record(self, collection, row):
        collection.append(row)
        self.level_rows.append(row)
        if self.progress is not None:
            self.progress({"stage": "gaussian_radial_control_level_completed",
                           "control": self.current_control, "level_result": row})

    def run(self, name, function):
        self.current_control, self.level_rows = name, []
        try:
            with localcontext() as context:
                context.prec = REFERENCE_PRECISION
                row = {"name": name, **function()}
            if type(row.get("passed")) is not bool:
                raise TypeError("Every control needs a boolean passed")
        except Exception as exc:
            row = {"name": name, "passed": False, "exception": type(exc).__name__,
                   "message": str(exc), "reference_work": self.reference.statistics(),
                   "completed_controls": [row["name"] for row in self.rows],
                   "completed_levels": self.level_rows}
            self.rows.append(row)
            self.exception_emitted = True
            if self.progress is not None:
                self.progress({"stage": "gaussian_radial_control_exception", "control_result": row})
            raise
        self.rows.append(row)
        if self.progress is not None:
            self.progress({"stage": "gaussian_radial_control_completed", "control_result": row})
        return row


def _closed_moment(n, beta, radius, reference_pi):
    if n % 2:
        k, z = (n-1)//2, beta*radius*radius
        total = sum((_power(z, j)/factorial(j) for j in range(k+1)), D(0))
        return D(factorial(k))*(-z).exp()*total/(2*beta**(k+1))
    if radius == 0:
        k = n//2
        return (reference_pi/beta).sqrt()*D(factorial(2*k))/(D(2)**(2*k+1)*factorial(k)*beta**k)
    return None


def _moment_control(beta, radius, harness):
    source = harness.source({(0, 0, 0): F(1)}, beta)
    measured = source.moment_enclosures(radius, MOMENT_POWERS)
    reference = harness.reference.moments(beta, radius, MOMENT_POWERS)
    pi_value = _agm_pi()
    settings = harness.settings
    tolerance = D(settings["gaussian_radial_control_moment_relative_tolerance"])
    refinement_tolerance = D(settings["gaussian_radial_control_refinement_relative_tolerance"])
    tail_tolerance = D(settings["gaussian_radial_control_tail_relative_tolerance"])
    rows = []
    for record in measured["moments"]:
        n = record["power"]
        lo, hi = D(record["lower_decimal"]), D(record["upper_decimal"])
        target = reference["values"][n]
        if (not all(value.is_finite() for value in (target, lo, hi))
                or not 0 < target or not 0 <= lo <= hi):
            raise ArithmeticError("Moment comparison requires positive finite enclosures/reference")
        midpoint, half_width = (lo+hi)/2, (hi-lo)/2
        tail, change = reference["tails"][n], reference["differences"][n]
        indicator = (abs(midpoint-target)+half_width+tail)/target
        closed = _closed_moment(n, _decimal(beta), _decimal(radius), pi_value)
        closed_error = None if closed is None else (abs(midpoint-closed)+half_width)/closed
        compatible = lo <= target+tail+change and hi >= target-change
        display = D.from_float(record["value"])
        declared = D.from_float(record["absolute_error_bound"])
        display_covers_interval = declared >= max(D(0), display-lo, hi-display)
        passed = (indicator <= tolerance and change/target <= refinement_tolerance
                  and tail/target <= tail_tolerance and compatible and display_covers_interval
                  and (closed_error is None or closed_error <= tolerance))
        harness.record(rows, {"power": n, "passed": passed, "provider": record,
                     "reference_decimal": str(target), "reference_coarse_decimal": str(reference["coarse"][n]),
                     "relative_comparison_indicator": str(indicator), "relative_refinement_change": str(change/target),
                     "reference_tail_bound_decimal": str(tail), "reference_relative_tail_bound": str(tail/target),
                     "provider_relative_enclosure_width": str((hi-lo)/target),
                     "independent_reference_compatible_with_enclosure": compatible,
                     "binary64_display_bound_covers_enclosure": display_covers_interval,
                     "closed_target_decimal": None if closed is None else str(closed),
                     "relative_closed_target_error": None if closed_error is None else str(closed_error)})
    return {"passed": all(row["passed"] for row in rows), "beta": str(beta), "radius": str(radius),
            "moments": rows, "reference_upper": str(reference["upper"]), "reference_panels": reference["panels"],
            "cutoff_steps": reference["cutoff_steps"], "provider_statistics": source.statistics(),
            "reference_work": harness.reference.statistics(), "accuracy_certified": False}


def _linear(reference, coefficients):
    value = sum((_decimal(c)*reference["values"][n] for n, c in coefficients.items()), D(0))
    uncertainty = sum((abs(_decimal(c))*(reference["tails"][n]+reference["differences"][n])
                       for n, c in coefficients.items()), D(0))
    return value, uncertainty


def _comparison(actual, target, tolerance, uncertainty=D(0), declared_error=None):
    actual = float(actual)
    if not isfinite(actual) or not target.is_finite() or not uncertainty.is_finite():
        raise ArithmeticError("Nonfinite Gaussian correction comparison")
    difference = abs(D.from_float(actual)-target)
    scale = max(D(1), abs(target))
    indicator = (difference+uncertainty)/scale
    covered = None if declared_error is None else difference <= D.from_float(float(declared_error))+uncertainty+D("1e-80")*scale
    return {"value": actual, "target_decimal": str(target), "reference_uncertainty_decimal": str(uncertainty),
            "scaled_indicator_decimal": str(indicator), "declared_arithmetic_bound": declared_error,
            "arithmetic_bound_consistent_with_reference": covered,
            "passed": indicator <= D(str(tolerance)) and covered is not False}


def _boundary_control(harness):
    # P=xz*(1-2r²/3); in H_21=-3xz its exact radial polynomial is
    # c0=-1/3, c1=2/9. The source is signed and contains a radial trace.
    polynomial = {(1, 0, 1): F(1), (3, 0, 1): F(-2, 3),
                  (1, 2, 1): F(-2, 3), (1, 0, 3): F(-2, 3)}
    source = harness.source(polynomial, ell_max=2)
    provider = source.exterior(F(2))
    expected_coefficients = {(2, 1, 0): F(-1, 3), (2, 1, 1): F(2, 9)}
    powers = (1, 3, 4, 6, 8, 10)
    at_R = harness.reference.moments(F(1), F(2), powers)
    at_r = harness.reference.moments(F(1), F(3), powers)
    finite = harness.reference.finite(F(1), F(2), F(3), (6, 8))
    B_R, eB_R = _linear(at_R, {1: F(-1, 3), 3: F(2, 9)})
    T_R, eT_R = _linear(at_R, {6: F(-1, 3), 8: F(2, 9)})
    B_r, eB_r = _linear(at_r, {1: F(-1, 3), 3: F(2, 9)})
    A_r, eA_r = _linear(finite, {6: F(-1, 3), 8: F(2, 9)})
    tolerance = harness.settings["gaussian_radial_control_correction_tolerance"]
    rows = []
    for alpha in (F(1), F(1, 8)):
        ai, d, lam = _decimal(alpha), D(5), D("2.5")
        for radius, side in ((F(0), None), (F(1), None), (F(2), "left"), (F(2), "right"), (F(3), None)):
            s = _decimal(radius*radius)
            if radius < 2 or side == "left":
                targets, errors = [ai*B_R/d, D(0), D(0), D(0)], [ai*eB_R/d, D(0), D(0), D(0)]
            else:
                A, eA = (D(0), D(0)) if radius == 2 else (A_r, eA_r)
                B, eB = (B_R, eB_R) if radius == 2 else (B_r, eB_r)
                power = _decimal(radius)**(-5)
                F0 = (-s).exp()*(-D(1)/3+D(2)*s/9)
                F1 = (-s).exp()*(D(5)/9-D(2)*s/9)
                targets = [ai/d*(power*A+B), -ai*power*A/(2*s),
                           ai*(lam+1)*power*A/(2*s*s)-ai*F0/(4*s),
                           -ai*(lam+1)*(lam+2)*power*A/(2*s**3)+ai*(lam+2)*F0/(4*s*s)-ai*F1/(4*s)]
                errors = [ai/d*(power*eA+eB), ai*power*eA/(2*s),
                          ai*(lam+1)*power*eA/(2*s*s), ai*(lam+1)*(lam+2)*power*eA/(2*s**3)]
            values = provider.outside_response_regular_jet(2, 1, radius, alpha, order=3, side=side)
            bounds = provider.response_error_bounds(2, 1, radius, alpha, order=3, side=side)
            comparisons = [_comparison(value, target, tolerance, error, bound)
                           for value, target, error, bound in zip(values, targets, errors, bounds)]
            harness.record(rows, {"alpha": str(alpha), "radius": str(radius), "side": side,
                         "jets": comparisons, "passed": all(row["passed"] for row in comparisons)})
        shift = _comparison(provider.origin_shift(2, 1, alpha), ai*B_R/d, tolerance, ai*eB_R/d,
                            provider.origin_shift_error_bound(2, 1, alpha))
        term = provider.algebraic_response_terms(2, 1, alpha)
        coefficient = _comparison(term[(3, 0)], ai*T_R/d, tolerance, ai*eT_R/d,
                                  provider.algebraic_response_error_bounds(2, 1, alpha)[(3, 0)])
        # Full outside response minus its Gaussian remainder must contain the
        # asymptotic homogeneous term once, including all regular derivatives.
        full = provider.outside_response_regular_jet(2, 1, F(3), alpha, order=3)
        remainder = provider.remainder_response_regular_jet(2, 1, F(3), alpha, order=3)
        decomposition = []
        multiplier = D(1)
        for j in range(4):
            if j:
                multiplier *= -(lam+j-1)
            target = ai*T_R/d*multiplier*D(3)**(-5-2*j)
            uncertainty = ai*eT_R/d*abs(multiplier)*D(3)**(-5-2*j)
            decomposition.append(_comparison(full[j]-remainder[j], target, tolerance, uncertainty))
        harness.record(rows, {"alpha": str(alpha), "kind": "origin_shift_and_algebraic_split", "origin_shift": shift,
                     "algebraic_coefficient": coefficient, "decomposition_jets": decomposition,
                     "passed": shift["passed"] and coefficient["passed"] and all(row["passed"] for row in decomposition)})
    return {"passed": dict(source.coefficients) == expected_coefficients and all(row["passed"] for row in rows),
            "exact_source_coefficients_match": dict(source.coefficients) == expected_coefficients,
            "rows": rows, "provider_provenance": provider.provenance(),
            "scope": "full outside-source correction and one-sided source jump; no FD solve"}


def _derivative_channels(harness):
    # Parent H_1,-1=-y. For axis y, d_y H=-H00 and
    # Hplus_y=(x²-2y²+z²)/3=H20/3+H22/6, independently written.
    source = harness.source({(0, 1, 0): F(-1)}, ell_max=1)
    parent = source.exterior(F(2))
    tolerance = harness.settings["gaussian_radial_control_correction_tolerance"]
    expected = {(0, 0): ("lower", F(-1)), (2, 0): ("upper", F(1, 3)), (2, 2): ("upper", F(1, 6))}
    rows = []
    for parent_alpha in (F(1), F(1, 8)):
        provider = parent.response_remainder_derivative(1, alpha=parent_alpha)
        mode_match = provider.mode_keys == tuple(sorted(expected))
        at_R = harness.reference.moments(F(1), F(2), (1, 3, 4, 6, 8, 10))
        for mode, (branch, angular) in expected.items():
            L, m = mode
            factor = parent_alpha*angular
            if branch == "lower":
                B, eB = _linear(at_R, {3: F(1, 6), 1: F(-4, 6)})
                M, eM = _linear(at_R, {4: F(1, 9), 1: F(-8, 9)})
            else:
                B, eB = _linear(at_R, {4: F(1, 24), 1: F(-1, 3)})
                M, eM = _linear(at_R, {6: F(1, 2), 4: F(-2)})
            factor_d = _decimal(factor)
            shift = _comparison(provider.origin_shift(L, m, 1), factor_d*B/(2*L+1), tolerance,
                                abs(factor_d)*eB/(2*L+1), provider.origin_shift_error_bound(L, m, 1))
            term = provider.algebraic_response_terms(L, m, 1)
            moment = _comparison(term[(L+1, 0)], factor_d*M/(2*L+1), tolerance,
                                 abs(factor_d)*eM/(2*L+1), provider.algebraic_response_error_bounds(L, m, 1)[(L+1, 0)])
            samples = []
            for radius in (F(2), F(3)):
                reference = harness.reference.moments(F(1), radius, (1, 3, 4, 6, 8, 10))
                r, s = _decimal(radius), _decimal(radius*radius)
                envelope = (-s).exp()
                if branch == "lower":
                    targets = [reference["values"][1]/3, -envelope/6, envelope/6]
                    errors = [(reference["tails"][1]+reference["differences"][1])/3, D(0), D(0)]
                else:
                    T = reference["values"][4]
                    eT = reference["tails"][4]+reference["differences"][4]
                    targets = [r**(-5)*T, -D("2.5")*r**(-7)*T-envelope/(2*s),
                               D("8.75")*r**(-9)*T+D("1.75")*envelope/(s*s)+envelope/(2*s)]
                    errors = [r**(-5)*eT, D("2.5")*r**(-7)*eT, D("8.75")*r**(-9)*eT]
                values = provider.source_regular_jet(L, m, radius, order=2)
                bounds = provider.source_error_bounds(L, m, radius, order=2)
                comparisons = [_comparison(value, factor_d*target, tolerance, abs(factor_d)*error, bound)
                               for value, target, error, bound in zip(values, targets, errors, bounds)]
                harness.record(samples, {"parent_alpha": str(parent_alpha), "mode": list(mode),
                                "radius": str(radius), "jets": comparisons,
                                "passed": all(row["passed"] for row in comparisons)})
            harness.record(rows, {"parent_alpha": str(parent_alpha), "mode": list(mode), "branch": branch,
                         "exact_angular_factor": str(angular), "complete_modes_match": mode_match,
                         "chi_origin_shift_from_Bh": shift, "chi_algebraic_coefficient_from_Mh": moment,
                         "source_samples": samples,
                         "passed": mode_match and shift["passed"] and moment["passed"] and all(row["passed"] for row in samples)})
    return {"passed": all(row["passed"] for row in rows), "channels": rows,
            "reference": "Nonaxial H_1,-1 derivative, explicit lower/upper harmonics and independently integrated Gaussian moments",
            "source_statistics": source.statistics()}


def _trace_factors(harness):
    # P=1-3xz+2y*r² = H00+H21-2r² H_1,-1.
    poly = {(0, 0, 0): F(1), (1, 0, 1): F(-3), (2, 1, 0): F(2),
            (0, 3, 0): F(2), (0, 1, 2): F(2)}
    rows = []
    for beta in (F(1), F(3, 2)):
        source = harness.source(poly, beta, ell_max=2)
        actual = source.asymptotic_moments()
        expected = {(0, 0, 0): F(1, 4)/beta, (0, 0, 1): F(3, 8)/beta**2,
                    (2, 1, 0): F(15, 16)/beta**3, (2, 1, 1): F(105, 32)/beta**4,
                    (1, -1, 0): F(-15, 8)/beta**3, (1, -1, 1): F(-105, 16)/beta**4}
        finite = source.finite_moments(F(0))
        prefactor = (_agm_pi()/_decimal(beta)).sqrt()
        comparisons = []
        for moment in finite["moments"]:
            l, m = moment["mode"]
            key = l, m, moment["trace"]
            target = _decimal(expected[key])*prefactor
            comparison = _comparison(moment["all_space"]["value"], target,
                         harness.settings["gaussian_radial_control_correction_tolerance"],
                         declared_error=moment["all_space"]["absolute_error_bound"])
            comparisons.append({"lmt": list(key), **comparison})
        exact = dict(actual) == expected
        harness.record(rows, {"beta": str(beta), "exact_symbolic_coefficients_match": exact,
                     "coefficients": [{"lmt": list(key), "coefficient": str(value)} for key, value in sorted(actual.items())],
                     "prefactor": "sqrt(pi/beta)", "numeric_checks": comparisons,
                     "passed": exact and all(row["passed"] for row in comparisons)})
    return {"passed": all(row["passed"] for row in rows), "rates": rows,
            "trace_levels_checked": [0, 1], "normalization": "Unnormalized H radial moments; no extra norm or4pi"}


def _reject(function, fragments):
    try:
        function()
    except (ValueError, TypeError, RuntimeError, ArithmeticError, OverflowError) as exc:
        return {"passed": any(fragment in str(exc) for fragment in fragments),
                "exception": type(exc).__name__, "message": str(exc)}
    return {"passed": False, "message": "Expected rejection did not occur"}


def _guards(harness):
    source = harness.source({(0, 0, 0): F(1)})
    bad_anchor = source.exterior(F(1, 10))
    grid_radius = float(F(1, 10))
    good_anchor = source.exterior(grid_radius)
    origin = source.origin_data(0, 0, 4, F(1, 20), derivative_order=4)
    origin_series = {(0, 0): OriginSeries(origin["patch_radius_float"], origin["coefficients"], origin["error_bounds_float_s"])}
    samples = source.sample((0.0, 0.05, grid_radius))
    def attach(provider):
        exterior = ExteriorSource(grid_radius, {}, provider.weighted_outer_bounds())
        return ScalarSourceSpec(samples, origin_series, exterior,
                                residual_policy="exact_correction", residual_provider=provider)
    rejection = _reject(lambda: attach(bad_anchor), ("exact radius", "anchor"))
    attach(good_anchor)
    exact_anchor_ok = (bad_anchor.radius == good_anchor.radius and bad_anchor.radius_exact == F(1, 10)
                       and good_anchor.radius_exact == F.from_float(grid_radius)
                       and bad_anchor.radius_exact != good_anchor.radius_exact)
    rows = []
    harness.record(rows, {"name": "exact_anchor_mismatch", "anchor_identity_checked": exact_anchor_ok,
                         **rejection, "passed": rejection["passed"] and exact_anchor_ok})
    harness.record(rows, {"name": "missing_boundary_side", **_reject(
        lambda: good_anchor.outside_response_regular_jet(0, 0, grid_radius, 1, order=3), ("side",))})
    harness.record(rows, {"name": "remainder_origin_domain", **_reject(
        lambda: good_anchor.remainder_response_regular_jet(0, 0, 0, 1), ("exterior", "origin"))})
    derivative = good_anchor.response_remainder_derivative(0)
    harness.record(rows, {"name": "deeper_derivative_rejection", **_reject(
        lambda: derivative.response_remainder_derivative(0), ("unsupported",))})
    harness.record(rows, {"name": "argument_range_cap", **_reject(
        lambda: source.moment_enclosures(1000, (0,)), ("max_z",))})
    harness.record(rows, {"name": "exact_work_cap", **_reject(
        lambda: harness.source({(0, 0, 0): F(1)}, max_products=1), ("max_products",))})
    limited = harness.source({(0, 0, 0): F(1)}, max_moment_steps=1)
    harness.record(rows, {"name": "moment_step_cap", **_reject(
        lambda: limited.moment_enclosures(0, (0,)), ("max_moment_steps",))})
    series_limited = harness.source({(0, 0, 0): F(1)}, max_series_terms=1)
    harness.record(rows, {"name": "series_term_cap", **_reject(
        lambda: series_limited.moment_enclosures(0, (0,)), ("max_series_terms",))})
    point_limited = harness.source({(0, 0, 0): F(1)}, max_point_evaluations=1)
    point_limited.moment_enclosures(0, (0,))
    harness.record(rows, {"name": "point_evaluation_cap", **_reject(
        lambda: point_limited.moment_enclosures(1, (0,)), ("max_point_evaluations",))})
    return {"passed": all(row["passed"] for row in rows), "checks": rows,
            "exact_anchor_policy": "Provider from grid.radius attaches; incompatible exact rational anchor rejects",
            "fd_solves_performed": 0}


def run_gaussian_radial_controls(config, progress=None):
    harness = None
    try:
        settings = _settings(config)
        harness = _Harness(settings, progress)
        for beta, radius in MOMENT_CASES:
            harness.run("moments_beta="+beta+"_R="+radius,
                        lambda beta=F(beta), radius=F(radius): _moment_control(beta, radius, harness))
        harness.run("outside_response_boundary_jets", lambda: _boundary_control(harness))
        harness.run("both_divergence_derivative_channels", lambda: _derivative_channels(harness))
        harness.run("exact_trace0_trace1_moments", lambda: _trace_factors(harness))
        harness.run("anchor_domain_and_work_guards", lambda: _guards(harness))
        return {"passed": all(row["passed"] for row in harness.rows), "settings": settings,
                "controls": harness.rows, "reference_work": harness.reference.statistics(),
                "source_attempts": harness.source_attempts,
                "source_statistics": [source.statistics() for source in harness.sources],
                "accuracy_certified": False, "fd_solves_performed": 0,
                "scope": "PG moments and exterior-response companions only; actual gravitational source/FD checks are separate"}
    except Exception as exc:
        if harness is None or not harness.exception_emitted:
            row = {"name": "settings_or_preparation", "passed": False, "exception": type(exc).__name__,
                   "message": str(exc), "reference_work": None if harness is None else harness.reference.statistics()}
            if progress is not None:
                progress({"stage": "gaussian_radial_control_exception", "control_result": row})
        raise
