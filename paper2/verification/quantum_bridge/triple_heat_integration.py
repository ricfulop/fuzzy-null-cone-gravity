"""Full integral of a supplied rate-one TripleHeatTail, in physical/(4*pi) units.

integrate_triple_heat(tail, *, radius=12, panel_width=2, orders=(16,24,32),
    precision=64, check_precision=96, absolute_tolerance=1e-10,
    relative_tolerance=1e-9, ..., progress=None) -> JSON-compatible diagnostics.

For tail.coefficients[(D,m,n,p)]=c, the normalized inner integrand is
    sum c*r**(D+2)*I_m(r*r)*I_n(r*r)*I_p(r*r).
There is NO additional radial measure or 4*pi in this inner expression.
The normalized algebraic exterior is pi**(3/2)*sum A[3,nu]*R**(nu+1)/(-nu-1).
All absolute tolerances and acceptance arithmetic use these SAME /4*pi units.
Physical Decimal and binary64 display values are returned separately.

The exact exterior coefficients come from the unchanged TripleHeatTail API.
Algebraic powers and their signed sum remain Fraction until conversion to
Decimal. The same Gaussian remainder inequality is evaluated in Decimal,
with independent Machin pi. No binary64 heat_integrals or erf is used.

I_m(z) is evaluated by the positive Kummer series at the largest needed m,
then positive downward recurrence. A geometric series-tail estimate is retained.
Independent radial-order and Decimal-precision changes are summed termwise;
the omitted exterior remainder is included in the acceptance budget. These
are numerical diagnostics, not a proof of quadrature or rounding accuracy.
No automatic radius increase, tolerance relaxation, or term deletion occurs.

Only the supplied triple term is integrated; this is not full H6. No fields,
sources, contractions, integrals, checks or quadrature are constructed on import.
"""

from decimal import Decimal as D, InvalidOperation, Overflow, Underflow, localcontext
from fractions import Fraction as F
from functools import lru_cache
from math import isfinite
from numbers import Rational

from algebraic_heat_tail import TripleHeatTail
from run import gauss_nodes


MOMENT_GUARD_DIGITS = 12


def _integer(value, name, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError(name+" must be an integer >= "+str(minimum))
    return value


def _fraction_number(value, name):
    """Exact rational input, or the exact decimal represented by str(value)."""
    if isinstance(value, bool):
        raise TypeError(name+" must be a finite real number")
    if isinstance(value, Rational):
        return F(value)
    try:
        decimal = value if isinstance(value, D) else D(str(value))
    except (InvalidOperation, ValueError, TypeError) as exc:
        raise ValueError(name+" must be a finite real number") from exc
    if not decimal.is_finite():
        raise ValueError(name+" must be finite")
    return F(decimal)


def _decimal(value):
    return D(value.numerator)/D(value.denominator)


def _set_context(ctx, precision):
    ctx.prec = precision
    # A vanished positive seed/bound must fail rather than appear exact.
    ctx.traps[Underflow] = True
    ctx.traps[Overflow] = True
    ctx.traps[InvalidOperation] = True


def _display(value):
    converted = float(value)
    return None if not isfinite(converted) or (converted == 0 and value != 0) else converted


@lru_cache(maxsize=16)
def _decimal_pi(precision):
    """Machin's 16*atan(1/5)-4*atan(1/239); independent of binary64 pi."""
    with localcontext() as ctx:
        _set_context(ctx, precision+MOMENT_GUARD_DIGITS)
        threshold = D(10)**(-precision-6)

        def inverse_atan(denominator):
            x = D(1)/denominator
            x2, power, value = x*x, x, x
            for index in range(1, 4*precision+64):
                power *= -x2
                term = power/(2*index+1)
                value += term
                if abs(term) < threshold:
                    return value
            raise ArithmeticError("Decimal Machin series did not converge")

        value = 16*inverse_atan(5)-4*inverse_atan(239)
        ctx.prec = precision
        return +value


class _Work:
    def __init__(self, max_heat_steps, max_series_terms):
        self.heat_cap = max_heat_steps
        self.series_cap = max_series_terms
        self.heat_steps = self.series_terms = self.recurrence_steps = 0
        self.moment_evaluations = self.term_evaluations = self.nodes = 0
        self.max_seed_terms = 0

    def consume_heat(self, count, recurrence=False):
        if self.heat_steps+count > self.heat_cap:
            raise ValueError("Decimal HEAT work exceeds max_heat_steps")
        self.heat_steps += count
        if recurrence:
            self.recurrence_steps += count
        else:
            self.series_terms += count

    def statistics(self):
        return {"heat_steps": self.heat_steps, "kummer_series_terms": self.series_terms,
                "downward_recurrence_steps": self.recurrence_steps,
                "moment_evaluations": self.moment_evaluations,
                "term_evaluations": self.term_evaluations, "radial_nodes": self.nodes,
                "maximum_seed_terms": self.max_seed_terms}


def _heat_moments(z, maximum, precision, work):
    """Positive Kummer seed and downward recurrence, with a series-tail estimate.

    T_k=z**k/(m+3/2)_k, T_(k+1)/T_k=z/(m+3/2+k). Once that
    ratio q<1, the remaining positive sum is <=T_(k+1)/(1-q).
    A relative seed truncation bound also bounds every downward result,
    because the recurrence only multiplies/adds positive quantities.
    Decimal rounding and z evaluation are checked separately, not certified.
    """
    if not z.is_finite() or z < 0:
        raise ValueError("HEAT argument must be finite and nonnegative")
    with localcontext() as ctx:
        _set_context(ctx, precision+MOMENT_GUARD_DIGITS)
        work.moment_evaluations += 1
        if z == 0:
            work.consume_heat(maximum+1, recurrence=True)
            return tuple(D(1)/(2*m+1) for m in range(maximum+1)), D(0), 0
        exponential = (-z).exp()
        if not exponential.is_finite() or exponential <= 0:
            raise ArithmeticError("Positive Decimal Kummer exponential is unrepresentable")
        b = D(maximum)+D("1.5")
        series = term = D(1)
        count = 1
        work.consume_heat(1)
        target = D(10)**(-precision-6)
        index = 0
        while True:
            ratio = z/(b+index)
            following = term*ratio
            remainder = following/(1-ratio) if ratio < 1 else None
            if remainder is not None and remainder <= target*series:
                relative_tail = remainder/series
                break
            if count >= work.series_cap:
                raise ValueError("Positive Kummer series exceeds max_series_terms")
            work.consume_heat(1)
            series += following
            term = following
            index += 1
            count += 1
        work.max_seed_terms = max(work.max_seed_terms, count)
        values = [D(0)]*(maximum+1)
        values[maximum] = exponential*series/(2*maximum+1)
        for m in range(maximum, 0, -1):
            work.consume_heat(1, recurrence=True)
            values[m-1] = (2*z*values[m]+exponential)/(2*m-1)
        if any(not value.is_finite() or value <= 0 for value in values):
            raise ArithmeticError("Positive Decimal HEAT moment is unrepresentable")
        return tuple(values), relative_tail, count


def _prepare_exterior(tail, radius, max_products):
    """Reuse exact exterior coefficients and compile rational radius factors."""
    coefficients = dict(tail.exterior_coefficients(max_products=max_products))
    bad = {nu: str(value) for (k, nu), value in coefficients.items() if k == 3 and nu > -3}
    if bad:
        raise ValueError("Exact algebraic falloff gate failed: "+repr(bad))
    algebraic = tuple((nu, value*radius**(nu+1)/(-nu-1))
                      for (k, nu), value in sorted(coefficients.items()) if k == 3)
    bound_groups = {}
    for (k, nu), coefficient in coefficients.items():
        if k == 3:
            # Full/(4*pi) error from replacing erf**3 by one: 3*pi*|A|*J_(nu-1).
            family, rate, power, pi_half_power, prefactor = "erf_replacement", 1, nu-1, 2, 3*abs(coefficient)
        else:
            family, rate, power, pi_half_power, prefactor = "gaussian", 3-k, nu, k, abs(coefficient)
        denominator = 2*rate-F(max(power-1, 0), 1)/(radius*radius)
        if denominator <= 0:
            raise ValueError("Increase the radius: Gaussian remainder denominator is not positive")
        key = family, rate, pi_half_power
        weight = prefactor*radius**(power-1)/denominator
        bound_groups[key] = bound_groups.get(key, F(0))+weight
    return {"coefficients": coefficients, "algebraic": algebraic,
            "algebraic_sum": sum((term for _, term in algebraic), F(0)),
            "algebraic_absolute": sum((abs(term) for _, term in algebraic), F(0)),
            "bound_groups": tuple(sorted(bound_groups.items())), "radius_squared": radius*radius}


def _evaluate_exterior(prepared, precision):
    with localcontext() as ctx:
        _set_context(ctx, precision)
        decimal_pi = _decimal_pi(precision)
        sqrt_pi = decimal_pi.sqrt()
        pi_powers = (D(1), sqrt_pi, decimal_pi, decimal_pi*sqrt_pi)
        value = pi_powers[3]*_decimal(prepared["algebraic_sum"])
        absolute = pi_powers[3]*_decimal(prepared["algebraic_absolute"])
        terms = tuple(pi_powers[3]*_decimal(term) for _, term in prepared["algebraic"])
        r2 = _decimal(prepared["radius_squared"])
        exponentials = {}
        bounds = []
        for (_, rate, pi_power), weight in prepared["bound_groups"]:
            if rate not in exponentials:
                exponentials[rate] = (-rate*r2).exp()
                if exponentials[rate] <= 0:
                    raise ArithmeticError("Positive Decimal Gaussian remainder is unrepresentable")
            bound = _decimal(weight)*pi_powers[pi_power]*exponentials[rate]
            if not bound.is_finite() or bound <= 0:
                raise ArithmeticError("Positive Decimal remainder term is unrepresentable")
            bounds.append(bound)
        remainder = sum(bounds, D(0))
        if not value.is_finite() or not absolute.is_finite() or not remainder.is_finite():
            raise ArithmeticError("Nonfinite Decimal exterior evaluation")
        proxy = D(10)**(1-precision)*(absolute+remainder)*(1+len(terms)+len(bounds))
        return {"value": value, "absolute": absolute, "terms": terms,
                "remainder": remainder, "bound_terms": tuple(bounds), "pi": decimal_pi,
                "roundoff_proxy": proxy, "precision": precision}


def _evaluate_inner(entries, panels, order, precision, maximum_moment, work, progress):
    with localcontext() as ctx:
        _set_context(ctx, precision)
        coefficients = tuple(_decimal(coefficient) for _, coefficient in entries)
        totals = [D(0)]*len(entries)
        max_relative_tail, maximum_seed_terms, nodes = D(0), 0, 0
        degrees = sorted({key[0] for key, _ in entries})
        triples = sorted({key[1:] for key, _ in entries})
        # With normalized [-1,1] weights, panel mapping multiplies by its
        # FULL width, not half its width. No extra 4*pi or r**2 is applied.
        quadrature = gauss_nodes(order, precision) if entries else ()
        for panel_index, (left, right) in enumerate(panels if entries else ()):
            lower, width = _decimal(left), _decimal(right-left)
            for node, weight in quadrature:
                radius = lower+width*(node+1)/2
                if radius <= 0:
                    raise ArithmeticError("Positive mapped radial node is unrepresentable")
                moments, relative_tail, seed_terms = _heat_moments(radius*radius, maximum_moment, precision, work)
                max_relative_tail = max(max_relative_tail, relative_tail)
                maximum_seed_terms = max(maximum_seed_terms, seed_terms)
                radial = {degree: radius**(degree+2) for degree in degrees}
                products = {triple: moments[triple[0]]*moments[triple[1]]*moments[triple[2]] for triple in triples}
                mapped_weight = width*weight
                for index, ((degree, m, n, p), _) in enumerate(entries):
                    totals[index] += coefficients[index]*mapped_weight*radial[degree]*products[m, n, p]
                nodes += 1
                work.nodes += 1
                work.term_evaluations += len(entries)
            if progress is not None:
                partial = sum(totals, D(0))
                partial_absolute = sum((abs(value) for value in totals), D(0))
                progress({"stage": "triple_inner_panel", "order": order, "precision": precision,
                          "panel": panel_index+1, "panels": len(panels),
                          "left": str(left), "right": str(right),
                          "normalized_partial_decimal": str(partial),
                          "absolute_term_partial_decimal": str(partial_absolute), **work.statistics()})
        value = sum(totals, D(0))
        absolute = sum((abs(term) for term in totals), D(0))
        if not value.is_finite() or not absolute.is_finite():
            raise ArithmeticError("Nonfinite Decimal inner integral")
        # (1+rho)**3-1 evaluated without subtracting nearly equal numbers.
        product_relative_tail = max_relative_tail*(3+max_relative_tail*(3+max_relative_tail))
        series_indicator = absolute*product_relative_tail
        proxy = (D(10)**(1-precision))*absolute*(1+nodes+len(entries)+maximum_seed_terms+maximum_moment)
        return {"value": value, "absolute": absolute, "terms": tuple(totals), "order": order,
                "precision": precision, "nodes": nodes, "panels": len(panels) if entries else 0,
                "maximum_seed_terms": maximum_seed_terms, "maximum_seed_relative_tail": max_relative_tail,
                "series_indicator": series_indicator, "roundoff_proxy": proxy}


def _change(old, new, keys, precision):
    with localcontext() as ctx:
        _set_context(ctx, precision)
        changes = tuple(abs(a-b) for a, b in zip(old["terms"], new["terms"]))
        termwise = sum(changes, D(0))
        signed = abs(new["value"]-old["value"])
        worst_index = max(range(len(changes)), key=changes.__getitem__) if changes else None
        return {"termwise": termwise, "signed": signed, "indicator": max(termwise, signed),
                "maximum": changes[worst_index] if worst_index is not None else D(0),
                "worst_key": keys[worst_index] if worst_index is not None else None,
                "term_count": len(changes)}


def _json_change(change):
    key = change["worst_key"]
    return {"termwise_absolute_change_decimal": str(change["termwise"]),
            "signed_absolute_change_decimal": str(change["signed"]),
            "acceptance_change_indicator_decimal": str(change["indicator"]),
            "maximum_term_change_decimal": str(change["maximum"]),
            "worst_term": list(key) if isinstance(key, tuple) else key,
            "compared_terms": change["term_count"]}


def _json_inner(row):
    with localcontext() as ctx:
        _set_context(ctx, row["precision"])
        ratio = row["absolute"]/abs(row["value"]) if row["value"] else None
        return {"order": row["order"], "precision": row["precision"],
                "nodes": row["nodes"], "panels": row["panels"],
                "normalized_value_decimal": str(row["value"]),
                "absolute_term_integral_decimal": str(row["absolute"]),
                "cancellation_ratio_decimal": None if ratio is None else str(ratio),
                "zero_sum_with_nonzero_terms": row["value"] == 0 and row["absolute"] != 0,
                "maximum_seed_terms": row["maximum_seed_terms"],
                "maximum_seed_relative_tail_decimal": str(row["maximum_seed_relative_tail"]),
                "series_truncation_indicator_decimal": str(row["series_indicator"]),
                "roundoff_proxy_decimal": str(row["roundoff_proxy"])}


def integrate_triple_heat(tail, *, radius=12, panel_width=2, orders=(16, 24, 32),
                         precision=64, check_precision=96, absolute_tolerance=1e-10,
                         relative_tolerance=1e-9, max_term_evaluations=20_000_000,
                         max_heat_steps=5_000_000, max_series_terms=4096,
                         max_exterior_products=1_000_000, max_panels=128,
                         max_order=256, max_moment=256, max_degree=256, progress=None):
    """Integrate a supplied triple term; tolerances apply to physical/(4*pi).

    Returned total is inner integral plus the exact algebraic exterior. Its
    omitted Gaussian exterior is bounded by the displayed inequality.
    `converged` is a numerical gate, and `accuracy_certified` is always false.
    All planned term-node work is checked before requesting quadrature nodes.
    Exceptions/cap failures propagate; no hidden fallback changes the target.
    """
    if not isinstance(tail, TripleHeatTail):
        raise TypeError("Expected a supplied TripleHeatTail")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable or None")
    limits = {"max_term_evaluations": max_term_evaluations, "max_heat_steps": max_heat_steps,
              "max_series_terms": max_series_terms, "max_exterior_products": max_exterior_products,
              "max_panels": max_panels, "max_order": max_order,
              "max_moment": max_moment, "max_degree": max_degree}
    for key, value in limits.items():
        _integer(value, key)
    _integer(precision, "precision", 56)
    _integer(check_precision, "check_precision", precision+1)
    orders = tuple(orders)
    if (len(orders) < 2 or any(type(order) is not int or not 2 <= order <= max_order for order in orders)
            or any(high <= low for low, high in zip(orders, orders[1:]))):
        raise ValueError("orders must contain at least two strictly increasing integers in [2,max_order]")
    radius = _fraction_number(radius, "radius")
    panel_width = _fraction_number(panel_width, "panel_width")
    atol = _fraction_number(absolute_tolerance, "absolute_tolerance")
    rtol = _fraction_number(relative_tolerance, "relative_tolerance")
    if radius <= 0 or panel_width <= 0 or atol <= 0 or rtol < 0:
        raise ValueError("radius, panel_width and absolute_tolerance must be positive; relative_tolerance nonnegative")
    quotient = radius/panel_width
    panel_count = (quotient.numerator+quotient.denominator-1)//quotient.denominator
    if panel_count > max_panels:
        raise ValueError("Requested radial partition exceeds max_panels")
    panels = tuple((i*panel_width, min((i+1)*panel_width, radius)) for i in range(panel_count))
    entries = tuple(sorted(tail.coefficients.items()))
    if any(key[0] > max_degree or max(key[1:]) > max_moment for key, _ in entries):
        raise ValueError("Supplied triple bins exceed max_degree or max_moment")
    maximum_moment = max((max(key[1:]) for key, _ in entries), default=0)
    planned_nodes = panel_count*(sum(orders)+orders[-1]) if entries else 0
    planned_terms = planned_nodes*len(entries)
    if planned_terms > max_term_evaluations:
        raise ValueError("Planned order and precision grids exceed max_term_evaluations")
    if planned_nodes*(maximum_moment+1) > max_heat_steps:
        raise ValueError("Minimum planned HEAT work exceeds max_heat_steps")
    settings = {"radius": str(radius), "panel_width": str(panel_width), "orders": list(orders),
                "precision": precision, "check_precision": check_precision,
                "moment_guard_digits": MOMENT_GUARD_DIGITS,
                "absolute_tolerance_exact": str(atol), "relative_tolerance_exact": str(rtol),
                "ABSOLUTE_TOLERANCE_UNITS": "FULL physical integral divided by 4*pi",
                "normalization": "normalized_value_decimal = physical/(4*pi)",
                "gauss_weight_convention": "HALF [-1,1] weights multiplied by full panel width",
                "required_decay_power": -3, "panels": panel_count,
                "planned_radial_nodes": planned_nodes, "planned_term_evaluations": planned_terms,
                "source_bins": len(entries), "maximum_required_moment": maximum_moment,
                "moment_method": "Positive Decimal Kummer seed and downward recurrence",
                "pi_method": "Independent Decimal Machin formula", **limits}
    if progress is not None:
        progress({"stage": "triple_exterior_prepare", "settings": settings})
    prepared = _prepare_exterior(tail, radius, max_exterior_products)
    if progress is not None:
        progress({"stage": "triple_exterior_ready", "exterior_bins": len(prepared["coefficients"]),
                  "algebraic_terms": len(prepared["algebraic"]), "bound_groups": len(prepared["bound_groups"])})
    work = _Work(max_heat_steps, max_series_terms)
    levels = [_evaluate_inner(entries, panels, order, precision, maximum_moment, work, progress) for order in orders]
    keys = tuple(key for key, _ in entries)
    refinements = [_change(old, new, keys, check_precision) for old, new in zip(levels, levels[1:])]
    checked = _evaluate_inner(entries, panels, orders[-1], check_precision, maximum_moment, work, progress)
    inner_precision = _change(levels[-1], checked, keys, check_precision)
    exterior_working = _evaluate_exterior(prepared, precision)
    exterior_checked = _evaluate_exterior(prepared, check_precision)
    exterior_keys = tuple(nu for nu, _ in prepared["algebraic"])
    exterior_precision = _change(exterior_working, exterior_checked, exterior_keys, check_precision)
    with localcontext() as ctx:
        _set_context(ctx, check_precision)
        value = checked["value"]+exterior_checked["value"]
        absolute = checked["absolute"]+exterior_checked["absolute"]
        order_indicator = refinements[-1]["indicator"]
        bound_precision = sum((abs(a-b) for a, b in zip(exterior_working["bound_terms"], exterior_checked["bound_terms"])), D(0))
        precision_indicator = inner_precision["indicator"]+exterior_precision["indicator"]+bound_precision
        series_indicator = max(row["series_indicator"] for row in levels+[checked])
        roundoff_proxy = max(row["roundoff_proxy"] for row in levels+[checked])+max(
            exterior_working["roundoff_proxy"], exterior_checked["roundoff_proxy"])
        remainder_bound = max(exterior_working["remainder"], exterior_checked["remainder"])
        combined_indicator = order_indicator+precision_indicator+series_indicator+remainder_bound+roundoff_proxy
        acceptance = _decimal(atol)+_decimal(rtol)*abs(value)
        converged = combined_indicator <= acceptance
        ratio = absolute/abs(value) if value else None
        physical_factor = 4*exterior_checked["pi"]
        physical = physical_factor*value
        physical_bound = physical_factor*remainder_bound
        normalized_float, physical_float = _display(value), _display(physical)
        result = {
            "normalized_value_decimal": str(value), "normalized_value": normalized_float,
            "normalized_float_out_of_range": normalized_float is None,
            "physical_value_decimal": str(physical), "value": physical_float,
            "value_float_out_of_range": physical_float is None,
            "normalized_inner_decimal": str(checked["value"]),
            "normalized_algebraic_exterior_decimal": str(exterior_checked["value"]),
            "normalized_gaussian_remainder_bound_decimal": str(remainder_bound),
            "physical_gaussian_remainder_bound_decimal": str(physical_bound),
            "converged": converged, "accuracy_certified": False, "exact_zero_kernel": not entries,
            "acceptance_budget_normalized_decimal": str(acceptance),
            "combined_indicator_normalized_decimal": str(combined_indicator),
            "radial_termwise_change_decimal": str(order_indicator),
            "precision_change_indicator_decimal": str(precision_indicator),
            "kummer_series_truncation_indicator_decimal": str(series_indicator),
            "roundoff_proxy_decimal": str(roundoff_proxy),
            "absolute_inner_term_integral_decimal": str(checked["absolute"]),
            "absolute_algebraic_term_integral_decimal": str(exterior_checked["absolute"]),
            "absolute_term_sum_decimal": str(absolute),
            "cancellation_ratio_decimal": None if ratio is None else str(ratio),
            "zero_sum_with_nonzero_terms": value == 0 and absolute != 0,
            "inner_levels": [_json_inner(row) for row in levels],
            "radial_refinements": [{"coarse_order": low, "fine_order": high, **_json_change(change)}
                                   for low, high, change in zip(orders, orders[1:], refinements)],
            "precision_refinement": {"checked_inner": _json_inner(checked),
                                     "inner": _json_change(inner_precision),
                                     "algebraic_exterior": _json_change(exterior_precision),
                                     "remainder_bound_termwise_change_decimal": str(bound_precision)},
            "exterior": {
                "coefficient_bins": len(prepared["coefficients"]),
                "algebraic_coefficients": [{"power": nu, "coefficient": str(c)}
                                           for (k, nu), c in sorted(prepared["coefficients"].items()) if k == 3],
                "normalized_algebraic_rational_sum_before_pi": str(prepared["algebraic_sum"]),
                "normalized_absolute_rational_sum_before_pi": str(prepared["algebraic_absolute"]),
                "bound_groups": [{"family": key[0], "gaussian_rate": key[1], "pi_half_power": key[2],
                                  "exact_radius_prefactor": str(weight)} for key, weight in prepared["bound_groups"]],
                "working_remainder_bound_normalized_decimal": str(exterior_working["remainder"]),
                "checked_remainder_bound_normalized_decimal": str(exterior_checked["remainder"]),
                "pi_decimal": str(exterior_checked["pi"]),
                "falloff_gate": "All exact k=3 powers must be <=-3; no coefficient is dropped",
            },
            "settings": settings, "work": work.statistics(),
            "qualification": "Inner quadrature plus algebraic exterior with an explicit omitted Gaussian-remainder inequality; refinement, series and rounding diagnostics are not interval certification; only the supplied triple term",
        }
    if progress is not None:
        progress({"stage": "triple_integral_complete", "converged": result["converged"],
                  "normalized_value_decimal": result["normalized_value_decimal"],
                  "combined_indicator_normalized_decimal": result["combined_indicator_normalized_decimal"],
                  **work.statistics()})
    return result
