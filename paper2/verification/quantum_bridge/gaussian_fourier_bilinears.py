"""Compiled scalar/York bilinears of localized Gaussian-HEAT sources.

Exact real coefficients only. The spatial functional is B(f,g), not a
sesquilinear extension to complex source coefficients. For real sources,
fhat(-k)=conj(fhat(k)), so Fourier real-real + imaginary-imaginary parts
implement the bilinear pairing. In the York k-dot numerator, retain these
ORIGINAL phase bins after multiplying by k; reclassifying from the new
polynomial degree would introduce the wrong sign.

A block is keyed by two endpoints ((beta,a),(beta,a)), a=None for a pure
Gaussian. Its exact term (m,n,i,j,k):c represents, WITHOUT pi**(3/2),
    c*t**(2*m)*u**(2*n)*A**i*B**j/(A+B)**k/sqrt(A+B),
A=beta_left+a_left*t*t and B=beta_right+a_right*u*u. Missing variables
use m/n=-1 and contribute one, with no integral. All widths stay positive.
The spatial and inverse integrations are analytic; only the heat variables
are integrated numerically. No heat_integrals, erf, or exponential routine
is used in the evaluator.

compile_scalar_bilinear returns <left,(-Delta)^-1 right/8>.
compile_york_bilinear returns <left,Delta_L^-1 right>, with multiplier
-delta_ij/k^2+k_i*k_j/(4*k^4). A constant k^-4 numerator is rejected.
Reciprocity swaps complete endpoints, including their heat weights and
Laurent powers. This preserves the integral; the canonical integrand need
not coincide pointwise with its pre-canonicalized expression.

Public API: compile_scalar_bilinear, compile_york_bilinear,
GaussianBilinearKernel.scale/add/merge/statistics/integrate,
integrate_heat_bilinear. All Decimal values and ABSOLUTE TOLERANCES are
normalized WITHOUT pi**(3/2). A full binary64 display value is also returned.
No compiler, quadrature, source construction, checks, or runs occur at import.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field as dataclass_field
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from functools import lru_cache
from math import isfinite, pi, prod
from numbers import Rational
from types import MappingProxyType

from gaussian_heat_sources import GaussianHeatSource
from reduced_gaunt_exact import odd_double_factorial
from run import gauss_nodes


DEFAULT_PAIR_CAP = 20_000_000
DEFAULT_FOURIER_CAP = 200_000
DEFAULT_TERM_CAP = 200_000


def _fraction(value, name="coefficient"):
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(name + " must be an exact real rational; float/Decimal/complex are unsupported")
    return F(value)


def _positive_int(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name + " must be a positive integer")
    return value


def _endpoint(value):
    beta, a = value
    beta = _fraction(beta, "beta")
    if beta <= 0:
        raise ValueError("Every Gaussian beta must be strictly positive")
    if a is not None:
        a = _fraction(a, "heat rate")
        if a <= 0:
            raise ValueError("A present heat rate must be strictly positive")
    return beta, a


def _endpoint_sort(endpoint):
    beta, a = endpoint
    return a is not None, beta, F(0) if a is None else a


def _block_sort(block):
    return tuple(_endpoint_sort(endpoint) for endpoint in block)


def _group_sort(group):
    beta, a, m = group
    return _endpoint_sort((beta, a)), -1 if m is None else m


def _canonical(block, term):
    left, right = block
    m, n, i, j, k = term
    if (_endpoint_sort(left) > _endpoint_sort(right)
            or (left == right and (m, i) > (n, j))):
        return (right, left), (n, m, j, i, k)
    return block, term


def _accumulate(out, key, coefficient):
    coefficient = out.get(key, F(0)) + coefficient
    if coefficient:
        out[key] = coefficient
    else:
        out.pop(key, None)


class _CompileBudget:
    def __init__(self, max_pair_products, max_fourier_terms, max_terms, max_degree, progress=None):
        self.pair_cap = _positive_int(max_pair_products, "max_pair_products")
        self.fourier_cap = _positive_int(max_fourier_terms, "max_fourier_terms")
        self.term_cap = _positive_int(max_terms, "max_terms")
        self.degree_cap = _positive_int(max_degree, "max_degree")
        self.pairs = self.emitted = self.transforms = self.group_pairs = 0
        self.peak_terms = self.live_terms = 0
        if progress is not None and not callable(progress):
            raise TypeError("progress must be callable or None")
        self.progress = progress
        self.active_block = None

    def enter_block(self, block, inverse_power):
        if self.progress is None:
            return
        block = tuple(sorted(block, key=_endpoint_sort))
        label = block, inverse_power
        if label != self.active_block:
            self.active_block = label
            self.progress({"stage": "fourier_compile_block", "endpoints": _json_block(block),
                           "inverse_power": inverse_power, **self.statistics()})

    def complete(self, kernel):
        if self.progress is not None:
            self.progress({"stage": "fourier_compile_complete", **kernel.statistics()})
        return kernel

    def emit(self, count):
        if self.emitted + count > self.fourier_cap:
            raise ValueError("Fourier expansion term budget exceeded")
        self.emitted += count

    def pair(self, count):
        if self.pairs + count > self.pair_cap:
            raise ValueError("Fourier contraction pair-product budget exceeded")
        self.pairs += count

    def put(self, blocks, block, term, coefficient):
        if not coefficient:
            return
        block, term = _canonical(block, term)
        terms = blocks.setdefault(block, {})
        existed = term in terms
        _accumulate(terms, term, coefficient)
        self.live_terms += int(term in terms) - int(existed)
        self.peak_terms = max(self.peak_terms, self.live_terms)
        if self.live_terms > self.term_cap:
            raise ValueError("Compiled Laurent kernel term budget exceeded")
        if not terms:
            blocks.pop(block, None)

    def statistics(self):
        return {"compile_pair_products": self.pairs,
                "fourier_terms_emitted": self.emitted,
                "polynomial_transforms": self.transforms,
                "source_group_pairs": self.group_pairs,
                "peak_kernel_terms": self.peak_terms,
                "max_pair_products": self.pair_cap,
                "max_source_group_pairs": self.pair_cap,
                "max_fourier_terms": self.fourier_cap,
                "max_degree": self.degree_cap}


@lru_cache(maxsize=8192)
def _width_hermite(powers):
    """Derivative polynomial of exp(-k^2/(4*A)), keyed by (k powers,A^-power).

    Differentiation lowers a Cartesian power, or appends -k_axis/(2*A).
    Returned immutable tuples contain no width-dependent numerical values.
    """
    if not any(powers):
        return (((0, 0, 0), 0, F(1)),)
    axis = next(i for i, degree in enumerate(powers) if degree)
    prior = list(powers)
    prior[axis] -= 1
    out = {}
    for gamma, inverse_power, coefficient in _width_hermite(tuple(prior)):
        if gamma[axis]:
            lower = list(gamma)
            lower[axis] -= 1
            _accumulate(out, (tuple(lower), inverse_power), gamma[axis] * coefficient)
        upper = list(gamma)
        upper[axis] += 1
        _accumulate(out, (tuple(upper), inverse_power + 1), -coefficient / 2)
    return tuple((gamma, power, c) for (gamma, power), c in sorted(out.items()))


def _fourier_polynomial(poly, budget, cache):
    key = tuple(sorted(poly.items()))
    if key in cache:
        return cache[key]
    real, imag = {}, {}
    for powers, coefficient in key:
        degree = sum(powers)
        if degree > budget.degree_cap:
            raise ValueError("Source polynomial degree exceeds max_degree")
        # The separable Hermite polynomial has exactly this many monomials.
        budget.emit(prod(power // 2 + 1 for power in powers))
        part = imag if degree % 2 else real
        phase_sign = (-1) ** (degree // 2)
        for gamma, inverse_power, h in _width_hermite(powers):
            _accumulate(part, (gamma, inverse_power), coefficient * phase_sign * h)
    result = (real, imag)
    cache[key] = result
    budget.transforms += 1
    return result


def _source_groups(source, budget, cache):
    if not isinstance(source, GaussianHeatSource):
        raise TypeError("Expected an immutable GaussianHeatSource")
    return {key: _fourier_polynomial(poly, budget, cache)
            for key, poly in source.grouped_terms().items()}


def _parity_groups(poly):
    groups = {}
    for (gamma, inverse_power), coefficient in poly.items():
        parity = tuple(power % 2 for power in gamma)
        groups.setdefault(parity, []).append((gamma, inverse_power, coefficient))
    return groups


@lru_cache(maxsize=65536)
def _moment_factor(powers, inverse_power):
    """Rational Fourier moment after pi^(3/2)/sqrt(A+B) is factored out."""
    if any(power % 2 for power in powers):
        return None
    n = sum(powers) // 2
    numerator = prod(odd_double_factorial(power - 1) for power in powers)
    if inverse_power == 1:
        return n, F(2) ** (n - 1) * numerator / (2 * n + 1)
    if inverse_power == 2:
        if n < 1:
            raise ValueError("Nonintegrable constant numerator multiplying k^-4")
        return n, F(2) ** (n - 2) * numerator / ((2 * n + 1) * (2 * n - 1))
    raise ValueError("Only inverse powers k^-2 and k^-4 are supported")


def _contract(left_key, left, right_key, right, inverse_power, factor, blocks, budget):
    beta_l, a_l, m_l = left_key
    beta_r, a_r, m_r = right_key
    block = ((beta_l, a_l), (beta_r, a_r))
    m, nheat = (-1 if m_l is None else m_l), (-1 if m_r is None else m_r)
    budget.enter_block(block, inverse_power)
    budget.group_pairs += 1
    if budget.group_pairs > budget.pair_cap:
        raise ValueError("Source-group traversal exceeds max_pair_products")
    for left_part, right_part in zip(left, right):
        lblocks, rblocks = _parity_groups(left_part), _parity_groups(right_part)
        for parity, lterms in lblocks.items():
            rterms = rblocks.get(parity, ())
            budget.pair(len(lterms) * len(rterms))
            for gamma, u, lc in lterms:
                for delta, v, rc in rterms:
                    powers = tuple(x + y for x, y in zip(gamma, delta))
                    degree, normalization = _moment_factor(powers, inverse_power)
                    shift = degree - inverse_power
                    term = (m, nheat, shift - u, shift - v, degree + 1 - inverse_power)
                    budget.put(blocks, block, term, factor * lc * rc * normalization)


def _k_dot(components, budget):
    """Multiply by k while retaining the ORIGINAL real/imaginary phase bins."""
    result = ({}, {})
    for axis, parts in enumerate(components):
        if parts is None:
            continue
        for output, source in zip(result, parts):
            budget.emit(len(source))
            for (gamma, inverse_power), coefficient in source.items():
                powers = list(gamma)
                powers[axis] += 1
                _accumulate(output, (tuple(powers), inverse_power), coefficient)
    return result


@dataclass(frozen=True)
class GaussianBilinearKernel:
    """Immutable exact integral representation, normalized without pi**(3/2)."""

    blocks: Mapping = dataclass_field(default_factory=dict)
    metadata: Mapping = dataclass_field(default_factory=dict)
    max_terms: int = DEFAULT_TERM_CAP

    def __post_init__(self):
        limit = _positive_int(self.max_terms, "max_terms")
        normalized, size = {}, 0
        if not isinstance(self.blocks, Mapping):
            raise TypeError("blocks must be a mapping")
        for raw_block, raw_terms in self.blocks.items():
            if len(raw_block) != 2:
                raise ValueError("Each block requires two endpoint descriptors")
            block = tuple(_endpoint(endpoint) for endpoint in raw_block)
            for term, coefficient in raw_terms.items():
                if len(term) != 5 or any(type(power) is not int for power in term):
                    raise ValueError("Kernel term must be five integer powers (m,n,i,j,k)")
                if term[4] < 0:
                    raise ValueError("The (A+B) denominator power must be nonnegative")
                for endpoint, moment in zip(block, term[:2]):
                    if (endpoint[1] is None and moment != -1) or (endpoint[1] is not None and moment < 0):
                        raise ValueError("Heat weights must stay attached to the corresponding endpoint")
                coefficient = _fraction(coefficient)
                key, powers = _canonical(block, tuple(term))
                group = normalized.setdefault(key, {})
                existed = powers in group
                _accumulate(group, powers, coefficient)
                size += int(powers in group) - int(existed)
                if size > limit:
                    raise ValueError("Kernel normalization exceeds max_terms")
                if not group:
                    normalized.pop(key, None)
        frozen = {key: MappingProxyType(dict(sorted(terms.items())))
                  for key, terms in sorted(normalized.items(), key=lambda item: _block_sort(item[0]))}
        object.__setattr__(self, "blocks", MappingProxyType(frozen))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    def scale(self, factor):
        factor = _fraction(factor, "scale factor")
        return type(self)({block: {term: factor * c for term, c in terms.items()}
                           for block, terms in self.blocks.items()}, self.metadata, self.max_terms)

    def add(self, other, factor=F(1)):
        if not isinstance(other, GaussianBilinearKernel):
            raise TypeError("Can only add another GaussianBilinearKernel")
        factor = _fraction(factor, "addition factor")
        blocks = {block: dict(terms) for block, terms in self.blocks.items()}
        limit = min(self.max_terms, other.max_terms)
        size = sum(len(terms) for terms in blocks.values())
        for block, terms in other.blocks.items():
            group = blocks.setdefault(block, {})
            for term, c in terms.items():
                existed = term in group
                _accumulate(group, term, factor * c)
                size += int(term in group) - int(existed)
                if size > limit:
                    raise ValueError("Merged kernel exceeds max_terms")
            if not group:
                blocks.pop(block, None)
        metadata = {"merged_kernel_count": self.metadata.get("merged_kernel_count", 1)
                    + other.metadata.get("merged_kernel_count", 1)}
        return type(self)(blocks, metadata, limit)

    def merge(self, other, factor=F(1)):
        """Immutable merge; equivalent to add, and returns a new kernel."""
        return self.add(other, factor)

    def statistics(self):
        return {**dict(self.metadata), "blocks": len(self.blocks),
                "terms": sum(len(terms) for terms in self.blocks.values()),
                "max_terms": self.max_terms,
                "coefficient_domain": "real Fraction only; bilinear functional",
                "normalization": "pi^(3/2) factored out",
                "block_details": [{"endpoints": _json_block(block),
                                   "dimension": sum(endpoint[1] is not None for endpoint in block),
                                   "terms": len(terms)} for block, terms in self.blocks.items()]}

    def integrate(self, **options):
        return integrate_heat_bilinear(self, **options)


def compile_scalar_bilinear(left, right, *, max_pair_products=DEFAULT_PAIR_CAP,
                            max_fourier_terms=DEFAULT_FOURIER_CAP,
                            max_terms=DEFAULT_TERM_CAP, max_degree=64, progress=None):
    """Compile <left, C right>, C=(-Delta)^-1/8, without evaluating sources."""
    budget = _CompileBudget(max_pair_products, max_fourier_terms, max_terms, max_degree, progress)
    cache, blocks = {}, {}
    left_groups = _source_groups(left, budget, cache)
    right_groups = _source_groups(right, budget, cache)
    for lk, lp in left_groups.items():
        for rk, rp in right_groups.items():
            _contract(lk, lp, rk, rp, 1, F(1, 8), blocks, budget)
    return budget.complete(GaussianBilinearKernel(
        blocks, {"operator": "C=(-Delta)^-1/8", **budget.statistics()}, max_terms))


def compile_york_bilinear(vector_left, vector_right, *, max_pair_products=DEFAULT_PAIR_CAP,
                          max_fourier_terms=DEFAULT_FOURIER_CAP,
                          max_terms=DEFAULT_TERM_CAP, max_degree=64, progress=None):
    """Compile general real vector cross bilinear <left, Delta_L^-1 right>."""
    left, right = tuple(vector_left), tuple(vector_right)
    if len(left) != 3 or len(right) != 3:
        raise ValueError("York bilinear needs three components on each side")
    budget = _CompileBudget(max_pair_products, max_fourier_terms, max_terms, max_degree, progress)
    cache, blocks = {}, {}
    lg = tuple(_source_groups(source, budget, cache) for source in left)
    rg = tuple(_source_groups(source, budget, cache) for source in right)
    for lcomponent, rcomponent in zip(lg, rg):
        for lk, lp in lcomponent.items():
            for rk, rp in rcomponent.items():
                _contract(lk, lp, rk, rp, 1, F(-1), blocks, budget)
    lkeys = set().union(*(component.keys() for component in lg))
    rkeys = set().union(*(component.keys() for component in rg))
    # Phase-bin labels are unchanged by multiplying the numerator by k.
    ldiv = {key: _k_dot(tuple(component.get(key) for component in lg), budget)
            for key in sorted(lkeys, key=_group_sort)}
    rdiv = {key: _k_dot(tuple(component.get(key) for component in rg), budget)
            for key in sorted(rkeys, key=_group_sort)}
    for lk, lp in ldiv.items():
        for rk, rp in rdiv.items():
            _contract(lk, lp, rk, rp, 2, F(1, 4), blocks, budget)
    return budget.complete(GaussianBilinearKernel(
        blocks, {"operator": "Delta_L^-1", **budget.statistics()}, max_terms))


def _json_block(block):
    return [{"beta": str(beta), "heat_rate": None if a is None else str(a)} for beta, a in block]


def _decimal_fraction(value):
    return D(value.numerator) / D(value.denominator)


def _tolerance(value, name, allow_zero=False):
    result = D(str(value))
    if not result.is_finite() or result < 0 or (not allow_zero and result == 0):
        raise ValueError(name + " must be finite and " + ("nonnegative" if allow_zero else "positive"))
    return result


def _shape(block, nleft, nright):
    return tuple(order if endpoint[1] is not None else 0
                 for endpoint, order in zip(block, (nleft, nright)))


def _planned_shapes(block, orders):
    shapes = {_shape(block, orders[0], orders[0])}
    for low, high in zip(orders, orders[1:]):
        shapes.update((_shape(block, high, low), _shape(block, low, high), _shape(block, high, high)))
    return tuple(sorted(shapes))


def _shape_nodes(shape):
    return prod(order or 1 for order in shape)


def _node_data(endpoint, order, precision, moments, powers):
    beta, a = endpoint
    beta = _decimal_fraction(beta)
    if a is None:
        nodes = ((D(0), D(1)),)
    else:
        # run.gauss_nodes weights are already half the standard [-1,1] weights.
        # The [0,1] map therefore uses them UNCHANGED.
        nodes = tuple(((node + 1) / 2, weight) for node, weight in gauss_nodes(order, precision))
        a = _decimal_fraction(a)
    result = []
    for t, weight in nodes:
        width = beta if a is None else beta + a * t * t
        if width <= 0:
            raise ArithmeticError("Nonpositive Gaussian width during evaluation")
        heat = {m: D(1) if m < 0 else t ** (2 * m) for m in moments}
        laurent = {power: width ** power for power in powers}
        result.append((weight, width, heat, laurent))
    return tuple(result)


def _evaluate_block(block, terms, shape, precision):
    """Evaluate one compiled grid. Every retained Laurent basis factor is positive."""
    with localcontext() as ctx:
        ctx.prec = precision
        entries = [(term, _decimal_fraction(c)) for term, c in terms.items()]
        left = _node_data(block[0], shape[0], precision,
                          {key[0] for key, _ in entries}, {key[2] for key, _ in entries})
        right = _node_data(block[1], shape[1], precision,
                           {key[1] for key, _ in entries}, {key[3] for key, _ in entries})
        denominator_powers = {key[4] for key, _ in entries}
        totals = [D(0)] * len(entries)
        for wl, a, heat_l, powers_l in left:
            for wr, b, heat_r, powers_r in right:
                width_sum = a + b
                common = wl * wr / width_sum.sqrt()
                denominators = {k: width_sum ** (-k) for k in denominator_powers}
                for index, ((m, n, i, j, k), coefficient) in enumerate(entries):
                    totals[index] += (coefficient * common * heat_l[m] * heat_r[n]
                                      * powers_l[i] * powers_r[j] * denominators[k])
        value = sum(totals, D(0))
        absolute = sum((abs(value) for value in totals), D(0))
        if not value.is_finite() or not absolute.is_finite():
            raise ArithmeticError("Nonfinite Decimal kernel integral")
        nodes = len(left) * len(right)
        proxy = (D(10) ** (1 - precision)) * absolute * (1 + nodes + len(entries))
        return {"value": value, "absolute": absolute, "terms": tuple(totals),
                "shape": shape, "nodes": nodes, "precision": precision, "roundoff_proxy": proxy}


def _difference(left, right, precision):
    with localcontext() as ctx:
        ctx.prec = precision
        return {"signed_absolute_change": abs(left["value"] - right["value"]),
                "termwise_absolute_change": sum((abs(a - b) for a, b in zip(left["terms"], right["terms"])), D(0))}


def _json_evaluation(row):
    with localcontext() as ctx:
        ctx.prec = row["precision"]
        ratio = row["absolute"] / abs(row["value"]) if row["value"] else None
        return {"orders_t_u": list(row["shape"]), "nodes": row["nodes"],
                "precision": row["precision"], "normalized_value_decimal": str(row["value"]),
                "absolute_term_integral_decimal": str(row["absolute"]),
                "cancellation_ratio_decimal": None if ratio is None else str(ratio),
                "zero_sum_with_nonzero_terms": row["value"] == 0 and row["absolute"] != 0,
                "roundoff_proxy_decimal": str(row["roundoff_proxy"])}


def _json_difference(row):
    return {key + "_decimal": str(value) for key, value in row.items()}


def _display(value):
    converted = float(value)
    if not isfinite(converted) or (converted == 0 and value != 0):
        return None
    return converted


def integrate_heat_bilinear(kernel, *, orders=(12, 24, 36), precision=64, check_precision=96,
                            absolute_tolerance=1e-10, relative_tolerance=1e-9,
                            max_term_evaluations=20_000_000, max_power=4096, max_order=256):
    """Independent t/u quadrature and final Decimal precision refinement.

    IMPORTANT: absolute_tolerance is in the normalized units WITHOUT
    pi**(3/2). relative_tolerance is unchanged. All acceptance arithmetic
    stays Decimal; the full binary64 value is for display/assembly only.
    Work is preflighted before any quadrature nodes or evaluations are made.
    `converged` is a diagnostic, never a rigorous quadrature certification.
    """
    if not isinstance(kernel, GaussianBilinearKernel):
        raise TypeError("Expected a compiled GaussianBilinearKernel")
    orders = tuple(orders)
    max_order = _positive_int(max_order, "max_order")
    if (len(orders) < 2 or any(type(n) is not int or not 2 <= n <= max_order for n in orders)
            or any(a >= b for a, b in zip(orders, orders[1:]))):
        raise ValueError("orders must contain >=2 increasing integers between 2 and max_order")
    if type(precision) is not int or precision < 56:
        raise ValueError("precision must be >=56 for run.gauss_nodes")
    if type(check_precision) is not int or check_precision <= precision:
        raise ValueError("check_precision must be an integer greater than precision")
    term_cap = _positive_int(max_term_evaluations, "max_term_evaluations")
    power_cap = _positive_int(max_power, "max_power")
    atol = _tolerance(absolute_tolerance, "absolute_tolerance")
    rtol = _tolerance(relative_tolerance, "relative_tolerance", True)
    planned = 0
    plans = []
    for block, terms in kernel.blocks.items():
        if any(abs(power) > power_cap for term in terms for power in term):
            raise ValueError("Compiled exponent exceeds max_power")
        shapes = _planned_shapes(block, orders)
        final_shape = _shape(block, orders[-1], orders[-1])
        work = len(terms) * (sum(_shape_nodes(shape) for shape in shapes) + _shape_nodes(final_shape))
        planned += work
        if planned > term_cap:
            raise ValueError("Planned grids and precision refinement exceed max_term_evaluations")
        plans.append((block, terms, shapes, final_shape))
    settings = {"orders": list(orders), "precision": precision, "check_precision": check_precision,
                "absolute_tolerance_normalized_decimal": str(atol),
                "relative_tolerance_decimal": str(rtol),
                "ABSOLUTE_TOLERANCE_UNITS": "value WITHOUT pi^(3/2)",
                "max_term_evaluations": term_cap, "planned_term_evaluations": planned,
                "max_power": power_cap, "max_order": max_order,
                "parameter_change_rule": "sum over blocks of max(final t+u termwise changes, final diagonal termwise change)",
                "gauss_weight_convention": "run.gauss_nodes HALF [-1,1] weights used unchanged on [0,1]"}
    block_records, working_values, checked_values = [], [], []
    working_absolute, checked_absolute = [], []
    parameter_errors, precision_errors, roundoff_proxies = [], [], []
    diagonal_changes = []
    for block, terms, shapes, final_shape in plans:
        evaluations = {shape: _evaluate_block(block, terms, shape, precision) for shape in shapes}
        refinement_rows = []
        finest_t = finest_u = {"signed_absolute_change": D(0), "termwise_absolute_change": D(0)}
        diagonal = {"signed_absolute_change": D(0), "termwise_absolute_change": D(0)}
        for low, high in zip(orders, orders[1:]):
            ll = evaluations[_shape(block, low, low)]
            hl = evaluations[_shape(block, high, low)]
            lh = evaluations[_shape(block, low, high)]
            hh = evaluations[_shape(block, high, high)]
            finest_t = _difference(hh, lh, check_precision)
            finest_u = _difference(hh, hl, check_precision)
            diagonal = _difference(hh, ll, check_precision)
            refinement_rows.append({"coarse_order": low, "fine_order": high,
                                    "t_present": block[0][1] is not None,
                                    "u_present": block[1][1] is not None,
                                    "t_refinement": _json_difference(finest_t),
                                    "u_refinement": _json_difference(finest_u),
                                    "diagonal_refinement": _json_difference(diagonal)})
        working = evaluations[final_shape]
        checked = _evaluate_block(block, terms, final_shape, check_precision)
        precision_change = _difference(checked, working, check_precision)
        with localcontext() as ctx:
            ctx.prec = check_precision
            parameter_errors.append(max(
                finest_t["termwise_absolute_change"] + finest_u["termwise_absolute_change"],
                diagonal["termwise_absolute_change"]))
        precision_errors.append(precision_change["termwise_absolute_change"])
        diagonal_changes.append(diagonal["termwise_absolute_change"])
        roundoff_proxies.append(working["roundoff_proxy"])
        working_values.append(working["value"])
        checked_values.append(checked["value"])
        working_absolute.append(working["absolute"])
        checked_absolute.append(checked["absolute"])
        block_records.append({"endpoints": _json_block(block), "terms": len(terms),
                              "dimension": sum(endpoint[1] is not None for endpoint in block),
                              "grids": [_json_evaluation(evaluations[shape]) for shape in shapes],
                              "refinements": refinement_rows,
                              "precision_refinement": {"checked": _json_evaluation(checked),
                                                       **_json_difference(precision_change)}})
    with localcontext() as ctx:
        ctx.prec = check_precision
        value = sum(checked_values, D(0))
        working_value = sum(working_values, D(0))
        absolute = sum(checked_absolute, D(0))
        parameter_change = sum(parameter_errors, D(0))
        precision_change = sum(precision_errors, D(0))
        roundoff_proxy = sum(roundoff_proxies, D(0))
        acceptance = atol + rtol * abs(value)
        converged = parameter_change + precision_change <= acceptance and roundoff_proxy <= acceptance
        ratio = absolute / abs(value) if value else None
        normalized_float = _display(value)
        full_float = None if normalized_float is None else normalized_float * pi ** 1.5
        if full_float is not None and not isfinite(full_float):
            full_float = None
        return {"normalized_value_decimal": str(value), "normalized_value": normalized_float,
                "value": full_float, "value_float_out_of_range": full_float is None,
                "converged": converged, "accuracy_certified": False,
                "exact_zero_kernel": not kernel.blocks,
                "acceptance_budget_normalized_decimal": str(acceptance),
                "parameter_termwise_change_decimal": str(parameter_change),
                "diagonal_termwise_change_decimal": str(sum(diagonal_changes, D(0))),
                "precision_termwise_change_decimal": str(precision_change),
                "precision_signed_change_decimal": str(abs(value - working_value)),
                "absolute_term_integral_decimal": str(absolute),
                "working_absolute_term_integral_decimal": str(sum(working_absolute, D(0))),
                "cancellation_ratio_decimal": None if ratio is None else str(ratio),
                "zero_sum_with_nonzero_terms": value == 0 and absolute != 0,
                "working_roundoff_proxy_decimal": str(roundoff_proxy),
                "blocks": block_records, "settings": settings, "kernel_statistics": kernel.statistics(),
                "qualification": "Independent parameter and arithmetic refinements, not rigorous quadrature or rounding certification; no spatial truncation"}
