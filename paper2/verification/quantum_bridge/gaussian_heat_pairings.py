"""Exact angular collapse and refined radial quadrature for Gaussian HEAT pairings.

Integrates P(x)*exp(-beta*r^2) times zero, one, or two HeatFields. The
Gaussian rate beta is explicit and positive: DO NOT add the rates inside
the heat integrals to beta. No source construction or quadrature runs at
import. This is not a backend for beta=0 or products of >2 HEAT fields.

Public assembly:
    acc = PairingAccumulator()
    acc.add(poly, beta, fields=(field0, field1), factor=Fraction(-8))
    total.merge(acc)                 # no repeated angular contractions
    snapshots = total.collapsed()   # immutable exact bins
    result = integrate_pairings(snapshots, quadrature_orders=(16, 32, 64))

Each bin (degree,m,n) denotes 4*pi*c times the radial integral of
r**(degree+2)*exp(-beta*r*r)*I_m(a*r*r)*I_n(b*r*r).
A missing heat slot has index -1 and contributes one. Equal-rate slots
are interchangeable and their moment indices are sorted before merging.
Angular coefficients are sphere AVERAGES, not integrals; 4*pi is applied
only by the radial evaluator.

The angular average of xhat**powers is zero for any odd power, otherwise
prod((p_i-1)!!)/(sum(powers)+1)!!. Heat polynomials are grouped by moment,
degree and parity. For two fields P is convolved with the smaller field
once, then contracted with the other field without a full triple product.
Exact signed bins merge before floating evaluation.

The radial tail uses |I_m(z)| <= 1/(2*m+1) and an upper Gaussian moment
recurrence. Its analytic inequality is rigorous; the reported floating
bound is not an interval-arithmetic certificate. Gauss refinement and
cancellation diagnostics are required before `converged` can be true;
that flag is a numerical acceptance diagnostic, not a proved GL remainder.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import ceil, exp, fsum, isfinite, log, nextafter, pi, prod, sqrt
from numbers import Rational
from sys import float_info
from types import MappingProxyType
from typing import Mapping, Tuple

from gaussian_green_fields import HeatField, heat_integrals
from reduced_gaunt_exact import multiply, odd_double_factorial
from run import gauss_nodes

Bin = Tuple[int, int, int]
DEFAULT_MAX_PAIR_PRODUCTS = 20_000_000


def _fraction(value):
    if not isinstance(value, Rational):
        raise TypeError("Exact coefficients and rates must be rational, not float")
    return F(value)


def _positive_rate(value):
    value = _fraction(value)
    if value <= 0:
        raise ValueError("All Gaussian and HEAT rates must be positive")
    return value


def _powers(powers):
    powers = tuple(powers)
    if len(powers) != 3 or any(type(p) is not int or p < 0 for p in powers):
        raise ValueError("Expected three nonnegative Cartesian integer powers")
    return powers


def _polynomial(poly):
    result = {}
    for powers, coefficient in poly.items():
        powers, coefficient = _powers(powers), _fraction(coefficient)
        if coefficient:
            result[powers] = coefficient
    return result


def _accumulate(target, source, factor=F(1)):
    for key, coefficient in source.items():
        value = target.get(key, F(0)) + factor * coefficient
        if value:
            target[key] = value
        else:
            target.pop(key, None)


@lru_cache(maxsize=65536)
def angular_average(powers):
    """Exact normalized sphere average of the Cartesian unit-vector monomial."""
    powers = _powers(powers)
    if any(p % 2 for p in powers):
        return F(0)
    numerator = prod(odd_double_factorial(p - 1) for p in powers)
    return F(numerator, odd_double_factorial(sum(powers) + 1))


class _Budget:
    def __init__(self, maximum):
        if maximum is not None and (type(maximum) is not int or maximum < 0):
            raise ValueError("max_pair_products must be a nonnegative integer or None")
        self.maximum, self.used = maximum, 0

    def spend(self, count):
        if self.maximum is not None and self.used + count > self.maximum:
            raise ValueError("Exact angular-collapse pair-product budget exceeded: "
                             f"{self.used + count} > {self.maximum}")
        self.used += count


def _blocks(poly, include_degree=True):
    groups = {}
    for powers, coefficient in poly.items():
        parity = tuple(p % 2 for p in powers)
        key = (sum(powers), parity) if include_degree else parity
        groups.setdefault(key, {})[powers] = coefficient
    return groups


def _heat_polynomials(field):
    groups = {}
    for (powers, m), coefficient in field.terms.items():
        groups.setdefault(m, {})[powers] = coefficient
    return groups


def _angular_contract(left, right, budget):
    """Angular contraction of two polynomials, grouped by degree and parity."""
    left_blocks, right_blocks = _blocks(left), _blocks(right)
    result = {}
    for (left_degree, parity), a in left_blocks.items():
        for (right_degree, other_parity), b in right_blocks.items():
            if parity != other_parity:
                continue
            budget.spend(len(a) * len(b))
            coefficient = F(0)
            for pa, ca in a.items():
                for pb, cb in b.items():
                    powers = tuple(x + y for x, y in zip(pa, pb))
                    coefficient += ca * cb * angular_average(powers)
            _accumulate(result, {left_degree + right_degree: coefficient})
    return result


def _multiply_for_parities(left, right, wanted, budget):
    """Convolve only parity blocks that can survive the remaining contraction."""
    result = {}
    for pa, a in _blocks(left, False).items():
        for pb, b in _blocks(right, False).items():
            if tuple(x ^ y for x, y in zip(pa, pb)) not in wanted:
                continue
            budget.spend(len(a) * len(b))
            _accumulate(result, multiply(a, b))
    return result


@dataclass(frozen=True)
class CollapsedPairing:
    """Immutable exact angular bins; the factor 4*pi is not in coefficients."""

    beta: F
    rates: tuple
    coefficients: Mapping[Bin, F]
    statistics: Mapping = None

    def __post_init__(self):
        beta = _positive_rate(self.beta)
        rates = tuple(_positive_rate(rate) for rate in self.rates)
        if len(rates) > 2:
            raise ValueError("At most two HEAT rates are supported")
        coefficients = {}
        for key, value in self.coefficients.items():
            degree, m, n = key
            if type(degree) is not int or degree < 0 or degree % 2:
                raise ValueError("Collapsed degree must be even and nonnegative")
            indices = (m, n)
            if any(type(i) is not int for i in indices):
                raise ValueError("Moment indices must be integers")
            if any(indices[i] < 0 for i in range(len(rates))):
                raise ValueError("Present HEAT slots need nonnegative moment indices")
            if any(indices[i] != -1 for i in range(len(rates), 2)):
                raise ValueError("Missing HEAT slots must have moment index -1")
            value = _fraction(value)
            if len(rates) == 2 and rates[0] > rates[1]:
                m, n = n, m
            elif len(rates) == 2 and rates[0] == rates[1] and m > n:
                m, n = n, m
            _accumulate(coefficients, {(degree, m, n): value})
        object.__setattr__(self, "beta", beta)
        object.__setattr__(self, "rates", tuple(sorted(rates)))
        object.__setattr__(self, "coefficients", MappingProxyType(coefficients))
        object.__setattr__(self, "statistics", MappingProxyType(dict(self.statistics or {})))


def collapse_pairing(poly, beta, fields=(), factor=F(1),
                     max_pair_products=DEFAULT_MAX_PAIR_PRODUCTS):
    """Collapse one signed polynomial-Gaussian/HEAT product exactly.

    The pair-product budget counts rational convolution/contraction pairs,
    not time or integer bit complexity. None explicitly disables that guard.
    Inputs are copied/read only; cached source coefficients are not mutated.
    """
    poly, beta, factor = _polynomial(poly), _positive_rate(beta), _fraction(factor)
    fields = tuple(fields)
    if len(fields) > 2 or any(not isinstance(field, HeatField) for field in fields):
        raise ValueError("fields must contain zero, one, or two HeatField objects")
    fields = tuple(sorted(fields, key=lambda field: field.a))
    rates = tuple(field.a for field in fields)
    budget = _Budget(max_pair_products)
    naive_products = len(poly) * prod(len(field.terms) for field in fields)
    bins = {}
    if poly and factor and all(field.terms for field in fields):
        if not fields:
            for powers, coefficient in poly.items():
                _accumulate(bins, {(sum(powers), -1, -1): coefficient * angular_average(powers)})
        elif len(fields) == 1:
            for m, field_poly in _heat_polynomials(fields[0]).items():
                for degree, coefficient in _angular_contract(poly, field_poly, budget).items():
                    _accumulate(bins, {(degree, m, -1): coefficient})
        else:
            selected = 0 if len(fields[0].terms) <= len(fields[1].terms) else 1
            other = 1 - selected
            smaller = _heat_polynomials(fields[selected])
            remaining = _heat_polynomials(fields[other])
            wanted = {tuple(p % 2 for p in powers)
                      for field_poly in remaining.values() for powers in field_poly}
            for small_m, small_poly in smaller.items():
                combined = _multiply_for_parities(poly, small_poly, wanted, budget)
                if not combined:
                    continue
                for other_m, other_poly in remaining.items():
                    indices = [0, 0]
                    indices[selected], indices[other] = small_m, other_m
                    if rates[0] == rates[1]:
                        indices.sort()
                    for degree, coefficient in _angular_contract(combined, other_poly, budget).items():
                        _accumulate(bins, {(degree, *indices): coefficient})
    bins = {key: factor * value for key, value in bins.items() if factor * value}
    return CollapsedPairing(beta, rates, bins, {
        "source_polynomial_terms": len(poly),
        "heat_term_counts": [len(field.terms) for field in fields],
        "naive_cartesian_products": naive_products,
        "pair_products": budget.used,
        "collapsed_bins": len(bins),
    })


class PairingAccumulator:
    """Merge signed exact bins before quadrature, retaining explicit Gaussian rates."""

    def __init__(self, max_pair_products=DEFAULT_MAX_PAIR_PRODUCTS):
        _Budget(max_pair_products)
        self.max_pair_products = max_pair_products
        self._groups = {}
        self._add_calls = self._merge_calls = self._pair_products = 0

    def _merge_snapshot(self, snapshot):
        if not isinstance(snapshot, CollapsedPairing):
            raise TypeError("Expected a CollapsedPairing snapshot")
        key = (snapshot.beta, snapshot.rates)
        group = self._groups.setdefault(key, {})
        _accumulate(group, snapshot.coefficients)
        if not group:
            self._groups.pop(key, None)

    def add(self, poly, beta, fields=(), factor=F(1)):
        snapshot = collapse_pairing(poly, beta, fields, factor, self.max_pair_products)
        self._merge_snapshot(snapshot)
        self._add_calls += 1
        self._pair_products += snapshot.statistics["pair_products"]
        return self

    def merge(self, other_or_collapsed):
        """Merge another accumulator, one snapshot, or iterable of snapshots exactly."""
        if isinstance(other_or_collapsed, PairingAccumulator):
            snapshots = other_or_collapsed.collapsed()
        elif isinstance(other_or_collapsed, CollapsedPairing):
            snapshots = (other_or_collapsed,)
        else:
            snapshots = tuple(other_or_collapsed)
        if any(not isinstance(item, CollapsedPairing) for item in snapshots):
            raise TypeError("merge expects only CollapsedPairing snapshots")
        for snapshot in snapshots:
            self._merge_snapshot(snapshot)
        self._merge_calls += 1
        return self

    @property
    def statistics(self):
        return {"add_calls": self._add_calls, "merge_calls": self._merge_calls,
                "pair_products_executed_here": self._pair_products,
                "group_count": len(self._groups),
                "collapsed_bin_count": sum(len(group) for group in self._groups.values()),
                "per_add_pair_product_budget": self.max_pair_products}

    def collapsed(self):
        return tuple(CollapsedPairing(beta, rates, coefficients)
                     for (beta, rates), coefficients in sorted(self._groups.items()))

    def integrate(self, **options):
        result = integrate_pairings(self.collapsed(), **options)
        result["assembly_statistics"] = self.statistics
        return result


def _snapshots(pairings):
    if isinstance(pairings, PairingAccumulator):
        return pairings.collapsed()
    if isinstance(pairings, CollapsedPairing):
        return (pairings,)
    return tuple(pairings)


def _logadd(a, b):
    if a == float("-inf"):
        return b
    if b == float("-inf"):
        return a
    top = max(a, b)
    return top + log(exp(a - top) + exp(b - top))


def _gaussian_moment_log_bounds(beta, radius, max_power):
    """Upper log moments int_R^infty r^k exp(-beta*r^2) dr, for even k.

    At R>0 use M0 <= exp(-beta*R^2)/(2*beta*R), then
    Mk = R^(k-1)exp(-beta*R^2)/(2*beta) + (k-1) M_(k-2)/(2*beta).
    At R=0 the starting Gaussian moment is exact.
    """
    log_beta = log(beta)
    divisor = log(2) + log_beta
    if radius == 0:
        values = {0: 0.5 * log(pi) - log(2) - 0.5 * log_beta}
        for k in range(2, max_power + 1, 2):
            values[k] = log(k - 1) - divisor + values[k - 2]
        return values
    log_r = log(radius)
    damping = beta * radius * radius
    values = {0: -damping - divisor - log_r}
    for k in range(2, max_power + 1, 2):
        boundary = (k - 1) * log_r - damping - divisor
        inherited = log(k - 1) - divisor + values[k - 2]
        values[k] = _logadd(boundary, inherited)
    return values


def _reported_exp_upper(log_value):
    if log_value > log(float_info.max):
        return None
    value = exp(log_value)
    # Never silently report a zero tail merely because its estimate underflows.
    upper = nextafter(value, float("inf"))
    return upper if isfinite(upper) else None


def tail_bound(pairings, radius):
    """Analytic absolute radial tail inequality, evaluated in floating log arithmetic.

    Returns a JSON dictionary. None for absolute_bound means it overflows
    binary64; zero is reserved for an exactly empty assembled expression.
    This is not a directed-rounding interval certificate.
    """
    radius = float(radius)
    if not isfinite(radius) or radius < 0:
        raise ValueError("Tail radius must be finite and nonnegative")
    assembled = PairingAccumulator().merge(_snapshots(pairings)).collapsed()
    groups, log_total = [], float("-inf")
    for group in assembled:
        beta = float(group.beta)
        if not isfinite(beta) or beta <= 0:
            raise OverflowError("Gaussian rate is outside binary64 range")
        maximum = max(degree + 2 for degree, _, _ in group.coefficients)
        moments = _gaussian_moment_log_bounds(beta, radius, maximum)
        log_group = float("-inf")
        for (degree, m, n), coefficient in group.coefficients.items():
            term = (log(4 * pi) + log(abs(coefficient.numerator)) - log(coefficient.denominator)
                    + moments[degree + 2])
            for index in (m, n):
                if index >= 0:
                    term -= log(2 * index + 1)
            log_group = _logadd(log_group, term)
        log_total = _logadd(log_total, log_group)
        groups.append({"beta": str(group.beta), "rates": [str(a) for a in group.rates],
                       "absolute_bound": _reported_exp_upper(log_group),
                       "log_absolute_bound": log_group if isfinite(log_group) else None})
    return {"radius": radius,
            "absolute_bound": _reported_exp_upper(log_total) if assembled else 0.0,
            "log_absolute_bound": log_total if isfinite(log_total) else None,
            "groups": groups,
            "bound_type": "Im<=1/(2m+1), Mills upper bound for M0, positive Gaussian moment recurrence",
            "interval_certified": False}


@lru_cache(maxsize=32)
def _float_gauss_nodes(order, precision):
    # gauss_nodes returns half the standard [-1,1] weights (sum = 1).
    return tuple((float(node), float(weight)) for node, weight in gauss_nodes(order, precision))


def _finite(value, label):
    value = float(value)
    if not isfinite(value):
        raise OverflowError(label + " is outside binary64 range")
    return value


def _compile(groups):
    rates = sorted({rate for group in groups for rate in group.rates})
    rate_index = {rate: i for i, rate in enumerate(rates)}
    maxima = [0] * len(rates)
    rows = []
    for group_index, group in enumerate(groups):
        for (degree, m, n), coefficient in sorted(group.coefficients.items()):
            c = _finite(coefficient, "Collapsed coefficient")
            if c == 0:
                raise OverflowError("Nonzero exact collapsed coefficient underflows binary64")
            slots = tuple((rate_index[rate], moment)
                          for rate, moment in zip(group.rates, (m, n)))
            for index, moment in slots:
                maxima[index] = max(maxima[index], moment)
            rows.append((group_index, degree // 2 + 1, slots, _finite(4 * pi * c, "4*pi*coefficient")))
    return rows, tuple(float(rate) for rate in rates), tuple(maxima)


def _quadrature_level(groups, compiled, radius, panels, order, precision):
    rows, rates, maxima = compiled
    betas = tuple(float(group.beta) for group in groups)
    totals, compensation = [0.0] * len(rows), [0.0] * len(rows)
    max_power = max(row[1] for row in rows)
    nodes = _float_gauss_nodes(order, precision)
    width = radius / panels
    for panel in range(panels):
        midpoint = (panel + 0.5) * width
        for node, normalized_weight in nodes:
            r = midpoint + width * node / 2
            r2 = r * r
            powers = [1.0]
            for _ in range(max_power):
                powers.append(_finite(powers[-1] * r2, "Radial monomial"))
            gaussians = tuple(exp(-beta * r2) for beta in betas)
            moments = tuple(heat_integrals(rate * r2, mmax)
                            for rate, mmax in zip(rates, maxima))
            weight = width * normalized_weight
            for j, (group_index, power, slots, coefficient) in enumerate(rows):
                value = coefficient * powers[power] * gaussians[group_index] * weight
                for rate_index, moment in slots:
                    value *= moments[rate_index][moment]
                value = _finite(value, "Weighted radial bin")
                # Each bin has a fixed sign, so compensated accumulation precedes
                # the final cancellation between different exact bins.
                increment = value - compensation[j]
                updated = totals[j] + increment
                compensation[j] = (updated - totals[j]) - increment
                totals[j] = updated
    value = _finite(fsum(totals), "Pairing sum")
    term_l1 = _finite(fsum(abs(term) for term in totals), "Absolute bin integral")
    group_values = [_finite(fsum(term for term, row in zip(totals, rows) if row[0] == index),
                            "Gaussian group integral") for index in range(len(groups))]
    ratio = term_l1 / abs(value) if value else None
    ratio_overflow = ratio is not None and not isfinite(ratio)
    return {"order": order, "panels": panels, "nodes": panels * order,
            "value": value, "absolute_bin_integral": term_l1,
            "cancellation_ratio": None if ratio_overflow else ratio,
            "cancellation_ratio_overflow": ratio_overflow,
            "zero_sum_with_nonzero_bins": value == 0 and term_l1 != 0,
            "machine_epsilon_times_bin_l1": float_info.epsilon * term_l1,
            "group_values": group_values}, tuple(totals)


def integrate_pairings(pairings, quadrature_orders=(16, 32, 64), radius=None,
                       panel_width=None, absolute_tolerance=1e-10, relative_tolerance=1e-9,
                       decimal_precision=64, max_panels=256, max_nodes=65536,
                       max_bin_evaluations=10_000_000):
    """Integrate a merged expression with common nodes and HEAT tables.

    With radius=None, a deterministic Gaussian-tail search chooses R until
    the reported absolute tail bound is <= absolute_tolerance/4. Supplying
    R fixes it. Default panel width resolves the fastest Gaussian/HEAT rate.
    At least two increasing Gauss orders are mandatory. `converged` requires
    the last termwise absolute refinement change PLUS the tail bound to fit
    atol+rtol*abs(value), and epsilon*absolute_bin_integral to fit that budget.
    These diagnostics do not certify all quadrature or evaluation errors.
    """
    orders = tuple(quadrature_orders)
    if (len(orders) < 2 or any(type(n) is not int or n < 2 for n in orders)
            or any(a >= b for a, b in zip(orders, orders[1:]))):
        raise ValueError("Need at least two strictly increasing Gauss orders >=2")
    atol, rtol = float(absolute_tolerance), float(relative_tolerance)
    if not isfinite(atol) or atol <= 0 or not isfinite(rtol) or rtol < 0:
        raise ValueError("Require positive finite absolute and nonnegative finite relative tolerance")
    if type(decimal_precision) is not int or decimal_precision < 56:
        raise ValueError("decimal_precision must be >=56 for run.gauss_nodes")
    for label, limit in (("max_panels", max_panels), ("max_nodes", max_nodes),
                         ("max_bin_evaluations", max_bin_evaluations)):
        if type(limit) is not int or limit < 1:
            raise ValueError(label + " must be a positive integer")
    assembled = PairingAccumulator().merge(_snapshots(pairings))
    groups = assembled.collapsed()
    group_records = [{"beta": str(group.beta), "rates": [str(a) for a in group.rates],
                      "collapsed_bins": len(group.coefficients)} for group in groups]
    settings = {"quadrature_orders": list(orders), "absolute_tolerance": atol,
                "relative_tolerance": rtol, "decimal_precision_for_nodes": decimal_precision,
                "max_panels": max_panels, "max_nodes": max_nodes,
                "max_bin_evaluations_per_level": max_bin_evaluations}
    if not groups:
        return {"value": 0.0, "converged": True, "exact_angular_zero": True,
                "accuracy_certified": False, "tail": tail_bound((), 0.0),
                "refinement_absolute_change": 0.0, "termwise_refinement_change": 0.0,
                "absolute_bin_integral": 0.0, "cancellation_ratio": None,
                "cancellation_ratio_overflow": False,
                "machine_epsilon_times_bin_l1": 0.0,
                "levels": [], "groups": [], "settings": settings}
    betas = tuple(float(group.beta) for group in groups)
    all_rates = betas + tuple(float(rate) for group in groups for rate in group.rates)
    if any(not isfinite(rate) or rate <= 0 for rate in all_rates):
        raise OverflowError("A rate is outside positive binary64 range")
    degree = max(key[0] for group in groups for key in group.coefficients)
    radius_trials = []
    if radius is None:
        radius = (sqrt((degree + 2) / 2) + 6) / sqrt(min(betas))
        for _ in range(40):
            tail = tail_bound(groups, radius)
            radius_trials.append({"radius": radius, "absolute_tail_bound": tail["absolute_bound"]})
            if tail["absolute_bound"] is not None and tail["absolute_bound"] <= atol / 4:
                break
            radius *= 1.25
        else:
            raise ArithmeticError("Could not resolve the absolute Gaussian tail within 40 radius steps")
    else:
        radius = float(radius)
        if not isfinite(radius) or radius <= 0:
            raise ValueError("radius must be finite and positive")
        tail = tail_bound(groups, radius)
    width = 1 / sqrt(max(all_rates)) if panel_width is None else float(panel_width)
    if not isfinite(width) or width <= 0:
        raise ValueError("panel_width must be finite and positive")
    panels = ceil(radius / width)
    if panels > max_panels or panels * max(orders) > max_nodes:
        raise ValueError("Requested radial quadrature exceeds the panel/node budget")
    compiled = _compile(groups)
    bin_count = len(compiled[0])
    if bin_count * panels * max(orders) > max_bin_evaluations:
        raise ValueError("Requested radial quadrature exceeds max_bin_evaluations")
    levels, previous_bins = [], None
    for order in orders:
        level, bins = _quadrature_level(groups, compiled, radius, panels, order, decimal_precision)
        if previous_bins is not None:
            level["refinement_absolute_change"] = abs(level["value"] - levels[-1]["value"])
            level["termwise_refinement_change"] = _finite(
                fsum(abs(a - b) for a, b in zip(bins, previous_bins)), "Termwise refinement difference")
            level["group_refinement_changes"] = [abs(a - b) for a, b in
                                                  zip(level["group_values"], levels[-1]["group_values"])]
        levels.append(level)
        previous_bins = bins
    final = levels[-1]
    budget = _finite(atol + rtol * abs(final["value"]), "Acceptance budget")
    bound = tail["absolute_bound"]
    converged = (bound is not None and final["termwise_refinement_change"] + bound <= budget
                 and final["machine_epsilon_times_bin_l1"] <= budget)
    for record, value in zip(group_records, final["group_values"]):
        record["value"] = value
    settings.update({"radius": radius, "panel_width_requested": width,
                     "panel_width_actual": radius / panels, "panels": panels,
                     "collapsed_bin_count": bin_count})
    return {"value": final["value"], "converged": converged,
            "exact_angular_zero": False, "accuracy_certified": False,
            "acceptance_budget": budget, "tail": tail, "radius_search": radius_trials,
            "refinement_absolute_change": final["refinement_absolute_change"],
            "termwise_refinement_change": final["termwise_refinement_change"],
            "absolute_bin_integral": final["absolute_bin_integral"],
            "cancellation_ratio": final["cancellation_ratio"],
            "cancellation_ratio_overflow": final["cancellation_ratio_overflow"],
            "machine_epsilon_times_bin_l1": final["machine_epsilon_times_bin_l1"],
            "levels": levels, "groups": group_records, "settings": settings,
            "qualification": "Refinement and cancellation diagnostics plus an analytic tail inequality; not a rigorous quadrature error certificate"}
