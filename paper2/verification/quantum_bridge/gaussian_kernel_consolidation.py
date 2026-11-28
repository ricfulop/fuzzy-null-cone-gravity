"""Exact heat-weight elimination and common-denominator Gaussian kernels.

consolidate_kernel(kernel, *, max_expansion_products=50_000_000,
    max_terms=200_000, max_degree=128, max_coefficient_bits=4096,
    progress=None) -> (GaussianBilinearKernel, JSON-compatible statistics).

Only exact Fraction algebra is performed. No source/Fourier compiler,
quadrature, numerical evaluator or integration is called. Import is inert.
The returned kernel is consumed by the EXISTING integrator with unchanged
normalization, axis/precision grids, tolerances and evaluation-work cap.

For a present endpoint A=beta+a*t*t, eliminate t**(2*m) exactly as
(A-beta)**m/a**m. Absorb a pure endpoint's fixed beta**i before merging.
Merge the resulting Laurent coefficients by (i,j,k), choose I=min(i),
J=min(j), K=max(k), and form sum P_k(A,B)*(A+B)**(K-k) by sparse Horner.
Emit heat weights zero (pure sentinel -1), powers I+p,J+q, common k=K.
There is no change of integration variable and no Jacobian.

For equal endpoint descriptors, GaussianBilinearKernel canonicalization may
swap individual A/B monomials after heat elimination. Integral equality is
then audited by SYMMETRIZED rational integrands or square/paired grids, not
unsymmetrized point values. Coefficients and diagnostic basis change exactly;
cancellation and termwise convergence indicators can improve OR worsen.

stats['work']['total'] is the deterministic rewrite counter governed by
max_expansion_products. Its categories count generated binomial entries,
expanded heat products, dictionary scatters and fixed-width folds; they are
not processor instructions or a runtime estimate. Per-call caching makes
exact total/total-minus-one cap checks independent of previous calls.
"""

from collections import OrderedDict
from fractions import Fraction as F
from math import comb

from gaussian_fourier_bilinears import GaussianBilinearKernel


POWER_CACHE_COEFFICIENT_LIMIT = 4096


def _positive_int(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name+" must be a positive integer")
    return value


def _dimension(block):
    return sum(endpoint[1] is not None for endpoint in block)


def _endpoints(block):
    return [{"beta": str(beta), "heat_rate": None if a is None else str(a)} for beta, a in block]


class _Budget:
    def __init__(self, work_cap, term_cap, degree_cap, bit_cap):
        self.work_cap, self.term_cap = work_cap, term_cap
        self.degree_cap, self.bit_cap = degree_cap, bit_cap
        self.work = {"total": 0, "binomial_coefficients": 0, "heat_products": 0,
                     "laurent_scatters": 0, "numerator_scatters": 0,
                     "horner_scatters": 0, "output_scatters": 0,
                     "fixed_width_products": 0}
        self.live = self.peak_live = self.peak_degree = self.peak_bits = 0

    def charge(self, category, count=1):
        if self.work["total"]+count > self.work_cap:
            raise ValueError("Exact consolidation max_expansion_products exceeded; consumed="
                             +str(self.work["total"])+", attempted="+str(count)+", category="+category)
        self.work["total"] += count
        self.work[category] += count

    def coefficient(self, value):
        value = F(value)
        bits = max(abs(value.numerator).bit_length(), value.denominator.bit_length())
        self.peak_bits = max(self.peak_bits, bits)
        if bits > self.bit_cap:
            raise ValueError("Exact consolidation coefficient exceeds max_coefficient_bits")
        return value

    def degree(self, degree):
        self.peak_degree = max(self.peak_degree, degree)
        if degree > self.degree_cap:
            raise ValueError("Exact consolidation numerator exceeds max_degree")

    def put(self, mapping, key, value, category, polynomial=False):
        self.charge(category)
        value = self.coefficient(value)
        if polynomial:
            if min(key) < 0:
                raise ArithmeticError("Negative power in an exact numerator polynomial")
            self.degree(sum(key))
        existed = key in mapping
        value = self.coefficient(mapping.get(key, F(0))+value)
        if value:
            mapping[key] = value
        else:
            mapping.pop(key, None)
        self.live += int(key in mapping)-int(existed)
        self.peak_live = max(self.peak_live, self.live)
        if self.live > self.term_cap:
            raise ValueError("Exact consolidation live maps exceed max_terms")

    def take(self, mapping, key):
        value = mapping.pop(key)
        self.live -= 1
        return value

    def discard(self, mapping):
        self.live -= len(mapping)
        mapping.clear()


def _power(base, exponent, budget):
    """Bound coefficient size before constructing a possibly huge power."""
    base = budget.coefficient(base)
    if base not in (F(1), F(-1)):
        # This is a conservative upper bound, not a post-cancellation promise.
        prospective = abs(exponent)*max(abs(base.numerator).bit_length(), base.denominator.bit_length())
        if prospective > budget.bit_cap:
            raise ValueError("Exact consolidation power exceeds coefficient-bit preflight")
    return budget.coefficient(base**exponent)


class _HeatPowers:
    def __init__(self, budget):
        self.budget = budget
        self.cache = OrderedDict()
        self.coefficients = self.peak_coefficients = 0
        self.hits = self.misses = self.evictions = 0

    def get(self, beta, a, m):
        key = beta, a, m
        if key in self.cache:
            self.hits += 1
            self.cache.move_to_end(key)
            return self.cache[key]
        self.misses += 1
        self.budget.degree(m)
        rate_factor = _power(a, -m, self.budget)
        rows = []
        for power in range(m+1):
            self.budget.charge("binomial_coefficients")
            coefficient = self.budget.coefficient(comb(m, power)*_power(-beta, m-power, self.budget)*rate_factor)
            rows.append((power, coefficient))
        rows = tuple(rows)
        if len(rows) <= POWER_CACHE_COEFFICIENT_LIMIT:
            while self.cache and self.coefficients+len(rows) > POWER_CACHE_COEFFICIENT_LIMIT:
                _, old = self.cache.popitem(last=False)
                self.coefficients -= len(old)
                self.evictions += 1
            self.cache[key] = rows
            self.coefficients += len(rows)
            self.peak_coefficients = max(self.peak_coefficients, self.coefficients)
        return rows

    def statistics(self):
        return {"hits": self.hits, "misses": self.misses, "evictions": self.evictions,
                "peak_cached_coefficients": self.peak_coefficients,
                "coefficient_limit": POWER_CACHE_COEFFICIENT_LIMIT,
                "scope": "per call; no global cache or warming-dependent work counts"}


def _endpoint_expansion(endpoint, moment, exponent, powers, budget):
    beta, a = endpoint
    if a is None:
        budget.charge("fixed_width_products")
        return 0, ((0, _power(beta, exponent, budget)),)
    return exponent, powers.get(beta, a, moment)


def _emit(output, block, i, j, k, coefficient, budget):
    m = 0 if block[0][1] is not None else -1
    n = 0 if block[1][1] is not None else -1
    swapped = block[0] == block[1] and i > j
    if swapped:
        i, j = j, i
    budget.put(output, (m, n, i, j, k), coefficient, "output_scatters")
    return int(swapped)


def _zero_dimensional(block, terms, output, budget):
    a, b = block[0][0], block[1][0]
    total = F(0)
    for (_, _, i, j, k), coefficient in terms.items():
        budget.charge("fixed_width_products", 3)
        value = budget.coefficient(coefficient*_power(a, i, budget))
        value = budget.coefficient(value*_power(b, j, budget))
        value = budget.coefficient(value*_power(a+b, -k, budget))
        total = budget.coefficient(total+value)
    _emit(output, block, 0, 0, 0, total, budget)
    return {"after_heat_elimination_terms": int(bool(total)), "shift_i": 0, "shift_j": 0,
            "common_k": 0, "numerator_terms_before_canonicalization": int(bool(total)),
            "orientation_swaps": 0}


def _consolidate_block(block, terms, output, powers, budget):
    if _dimension(block) == 0:
        return _zero_dimensional(block, terms, output, budget)
    merged = {}
    for (m, n, i, j, k), coefficient in terms.items():
        coefficient = budget.coefficient(coefficient)
        base_i, left = _endpoint_expansion(block[0], m, i, powers, budget)
        base_j, right = _endpoint_expansion(block[1], n, j, powers, budget)
        for p, left_value in left:
            for q, right_value in right:
                budget.charge("heat_products")
                value = budget.coefficient(coefficient*left_value*right_value)
                budget.put(merged, (base_i+p, base_j+q, k), value, "laurent_scatters")
    merged_count = len(merged)
    if not merged:
        return {"after_heat_elimination_terms": 0, "shift_i": None, "shift_j": None,
                "common_k": None, "numerator_terms_before_canonicalization": 0, "orientation_swaps": 0}
    shift_i = min(key[0] for key in merged)
    shift_j = min(key[1] for key in merged)
    first_k, last_k = min(key[2] for key in merged), max(key[2] for key in merged)
    layers = {}
    for key in tuple(merged):
        i, j, k = key
        # Preflight the degree before expanding any power of S.
        budget.degree(i-shift_i+j-shift_j+last_k-k)
        value = budget.take(merged, key)
        budget.put(layers.setdefault(k, {}), (i-shift_i, j-shift_j), value,
                   "numerator_scatters", polynomial=True)
    s_terms = []
    fixed_sum = F(0)
    for axis, (beta, a) in enumerate(block):
        if a is None:
            fixed_sum += beta
        else:
            s_terms.append(((1, 0) if axis == 0 else (0, 1), F(1)))
    if fixed_sum:
        s_terms.append(((0, 0), budget.coefficient(fixed_sum)))
    polynomial = layers.pop(first_k)
    for k in range(first_k+1, last_k+1):
        product = {}
        for (p, q), coefficient in polynomial.items():
            for (u, v), factor in s_terms:
                budget.put(product, (p+u, q+v), coefficient*factor, "horner_scatters", polynomial=True)
        budget.discard(polynomial)
        polynomial = product
        layer = layers.pop(k, {})
        for key in tuple(layer):
            value = budget.take(layer, key)
            budget.put(polynomial, key, value, "numerator_scatters", polynomial=True)
    numerator_count, swaps = len(polynomial), 0
    for p, q in tuple(polynomial):
        coefficient = budget.take(polynomial, (p, q))
        swaps += _emit(output, block, shift_i+p, shift_j+q, last_k, coefficient, budget)
    return {"after_heat_elimination_terms": merged_count, "shift_i": shift_i, "shift_j": shift_j,
            "common_k": last_k, "numerator_terms_before_canonicalization": numerator_count,
            "orientation_swaps": swaps}


def consolidate_kernel(kernel, *, max_expansion_products=50_000_000, max_terms=200_000,
                       max_degree=128, max_coefficient_bits=4096, progress=None):
    """Return an equivalent kernel and exact rewrite-work statistics.

    No integration or scientific validation runs here. Input mappings remain
    untouched. All caps reject rather than prune. max_terms bounds concurrently
    live generated coefficient-map entries (including completed output blocks);
    a separate bounded cache holds at most 4096 binomial coefficients.
    """
    if not isinstance(kernel, GaussianBilinearKernel):
        raise TypeError("Expected an immutable GaussianBilinearKernel")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable or None")
    requested = {"max_expansion_products": max_expansion_products, "max_terms": max_terms,
                 "max_degree": max_degree, "max_coefficient_bits": max_coefficient_bits}
    for name, value in requested.items():
        _positive_int(value, name)
    effective_terms = min(max_terms, kernel.max_terms)
    budget = _Budget(max_expansion_products, effective_terms, max_degree, max_coefficient_bits)
    powers = _HeatPowers(budget)
    output_blocks, blocks = {}, []
    original_dimensions = {"0": 0, "1": 0, "2": 0}
    output_dimensions = {"0": 0, "1": 0, "2": 0}
    for index, (block, terms) in enumerate(kernel.blocks.items()):
        dimension = _dimension(block)
        original_dimensions[str(dimension)] += len(terms)
        before = budget.work["total"]
        if progress is not None:
            progress({"stage": "kernel_consolidation_block_start", "block": index+1,
                      "blocks": len(kernel.blocks), "endpoints": _endpoints(block),
                      "dimension": dimension, "original_terms": len(terms), "work": dict(budget.work)})
        output = output_blocks.setdefault(block, {})
        details = _consolidate_block(block, terms, output, powers, budget)
        output_dimensions[str(dimension)] += len(output)
        row = {"endpoints": _endpoints(block), "dimension": dimension, "original_terms": len(terms),
               "output_terms": len(output), "work_total": budget.work["total"]-before, **details}
        blocks.append(row)
        if not output:
            output_blocks.pop(block, None)
        if progress is not None:
            progress({"stage": "kernel_consolidation_block_complete", "block": index+1,
                      "block_statistics": row, "work": dict(budget.work)})
    output_count = sum(output_dimensions.values())
    if budget.live != output_count:
        raise ArithmeticError("Exact consolidation live-map accounting mismatch")
    statistics = {
        "original_terms": sum(original_dimensions.values()), "output_terms": output_count,
        "original_terms_by_dimension": original_dimensions,
        "output_terms_by_dimension": output_dimensions,
        "original_blocks": len(kernel.blocks), "output_blocks": len(output_blocks),
        "work": dict(budget.work), "limits": {**requested, "effective_max_terms": effective_terms},
        "peak_live_map_terms": budget.peak_live, "peak_numerator_degree": budget.peak_degree,
        "peak_coefficient_bits": budget.peak_bits, "binomial_cache": powers.statistics(),
        "blocks": blocks,
        "normalization": "pi^(3/2) factored out, unchanged",
        "diagnostic_basis_changed": True,
        "diagnostic_basis_note": "The existing integrator now refines compact Laurent terms; absolute-term and axis/precision changes need not match the raw representation",
        "identity_contract": "Equal endpoint blocks require symmetrized stripped rational integrands; canonical term swaps preserve the integral, not arbitrary unsymmetrized point values",
        "qualification": "Exact algebra only; no numerical convergence or performance claim; no source/evaluator caps or tolerances changed",
    }
    metadata = {**dict(kernel.metadata), "representation": "exact heat-eliminated common-denominator Laurent kernel",
                "consolidation": statistics}
    result = GaussianBilinearKernel(output_blocks, metadata, effective_terms)
    if progress is not None:
        progress({"stage": "kernel_consolidation_complete", "statistics": statistics})
    return result, statistics
