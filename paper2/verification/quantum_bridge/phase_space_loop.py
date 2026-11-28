"""Finite-mode bosonic one-loop 1PI derivatives, with hbar=1.

This is the ordered-partition prescription in section 7 of
reports/quantum-bridge/reduced-sixth-order-and-loop-kernel.md. It evaluates
supplied Hamiltonian Hessian derivatives; it does not construct gravitational
vertices or turn a 1PI coefficient into a scattering amplitude.

Only the Python standard library is used. Canonical coordinates are always
interleaved (Q0, P0, Q1, P1, ...). Missing blocks mean exactly zero vertices.
"""

from functools import lru_cache
from math import cos, exp, fsum, isfinite, log, pi, tan

try:  # Support both the existing script runner and package-style imports.
    from .run import gauss_nodes
except ImportError:
    from run import gauss_nodes


def _complex_sum(values):
    values = tuple(values)
    return complex(fsum(value.real for value in values),
                   fsum(value.imag for value in values))


@lru_cache(None)
def ordered_partitions(n):
    """All ordered partitions of labels 0,...,n-1, each exactly once.

    No cyclic quotient is taken: the 1/k in the logarithm is applied by
    integrate_vertex. Counts grouped by k are (1,2), (1,6,6), (1,14,36,24).
    """
    if isinstance(n, bool) or n not in (2, 3, 4):
        raise ValueError("the supported external leg counts are 2, 3 and 4")

    def visit(remaining):
        if not remaining:
            yield ()
            return
        subset = remaining
        while subset:
            block = frozenset(i for i in range(n) if subset & (1 << i))
            for tail in visit(remaining ^ subset):
                yield (block,) + tail
            subset = (subset - 1) & remaining

    return tuple(sorted(visit((1 << n) - 1), key=len))


@lru_cache(None)
def _line_nodes(order, scale):
    # run.gauss_nodes returns HALF the usual Legendre weights (sum=1).
    # With nu=scale*tan(pi*x/2), 2*w*dnu/dx/(2*pi)=w*scale/(2*cos^2).
    return tuple((scale * tan(pi * float(x) / 2),
                  float(weight) * scale / (2 * cos(pi * float(x) / 2) ** 2))
                 for x, weight in gauss_nodes(order, 60))


def _matmul(a, b):
    size = len(a)
    result = [[0j] * size for _ in range(size)]
    for row, out in zip(a, result):
        for k, value in enumerate(row):
            if value != 0:
                for j, other in enumerate(b[k]):
                    if other != 0:
                        out[j] += value * other
    return result


def _trace_product(factors):
    """Avoid constructing the last product when only its trace is needed."""
    product = factors[0]
    for factor in factors[1:-1]:
        product = _matmul(product, factor)
    last = factors[-1]
    return _complex_sum(product[i][j] * last[j][i]
                        for i in range(len(product))
                        for j in range(len(product))
                        if product[i][j] != 0 and last[j][i] != 0)


def _propagator_times_block(frequencies, nu, block):
    """Left-multiply by block-diagonal G, without building a dense G."""
    result = []
    for mode, omega in enumerate(frequencies):
        denominator = nu * nu + omega * omega
        qq, qp, pp = 1 / denominator, nu / denominator, omega * omega / denominator
        qrow, prow = block[2 * mode:2 * mode + 2]
        result.append([qq * q + qp * p for q, p in zip(qrow, prow)])
        result.append([-qp * q + pp * p for q, p in zip(qrow, prow)])
    return result


def _validated_inputs(frequencies, external_frequencies, blocks, quadrature_order):
    frequencies = tuple(float(value) for value in frequencies)
    if not frequencies or any(not isfinite(value) or value <= 0 for value in frequencies):
        raise ValueError("frequencies must be a nonempty list of finite positive oscillator frequencies")
    external_frequencies = tuple(float(value) for value in external_frequencies)
    n = len(external_frequencies)
    if n not in (2, 3, 4):
        raise ValueError("integrate_vertex supports n=2,3,4 external legs")
    if any(not isfinite(value) for value in external_frequencies):
        raise ValueError("external frequencies must be finite real numbers")
    if abs(fsum(external_frequencies)) > 1e-12 * max(1.0, fsum(abs(value) for value in external_frequencies)):
        raise ValueError("external frequencies must sum to zero (the conservation delta is stripped)")
    if isinstance(quadrature_order, bool) or not isinstance(quadrature_order, int) or quadrature_order < 2:
        raise ValueError("quadrature_order must be an integer at least 2")
    size = 2 * len(frequencies)
    labels = frozenset(range(n))
    nonzero = {}
    for subset, matrix in blocks.items():
        if (not isinstance(subset, frozenset) or not subset or not subset <= labels
                or any(type(label) is not int for label in subset)):
            raise ValueError("block keys must be nonempty frozensets of external labels 0,...,n-1")
        rows = [list(row) for row in matrix]
        if len(rows) != size or any(len(row) != size for row in rows):
            raise ValueError("each Hessian block must have shape (2*d, 2*d)")
        rows = [[complex(value) for value in row] for row in rows]
        if any(not isfinite(value.real) or not isfinite(value.imag) for row in rows for value in row):
            raise ValueError("Hessian block entries must be finite")
        for i in range(size):
            for j in range(i):
                a, b = rows[i][j], rows[j][i]
                if abs(a - b) > 1e-12 * max(1.0, abs(a), abs(b)):
                    raise ValueError("Hessian blocks must be symmetric under transpose, not conjugate transpose")
                # Remove only transpose asymmetry at floating-point tolerance.
                rows[i][j] = rows[j][i] = (a + b) / 2
        if any(value != 0 for row in rows for value in row):
            nonzero[subset] = tuple(tuple(row) for row in rows)
    return frequencies, external_frequencies, nonzero


def integrate_vertex(frequencies, external_frequencies, blocks, quadrature_order=64):
    """Return the complete supplied n=2,3,4 one-loop 1PI derivative.

    ``frequencies`` contains positive free oscillator frequencies. ``blocks``
    maps nonempty frozensets of external indices to complex symmetric 2d-by-2d
    Hessian *derivatives* B_S, already contracted with the external directions.
    There are no extra Taylor factorials to attach to these blocks. Only
    missing or exactly zero blocks are skipped; tiny nonzero entries are kept.
    Transpose asymmetries within 1e-12 scaled tolerance are averaged on input.

    External frequencies are real and energy conserving. For every ordered
    partition (S1,...,Sk), the factor is (-1)^(k+1)/(2*k), and propagator a
    carries nu + sum(Omega_r for r in preceding blocks). The free 2-by-2 block
    is [[1,nu],[-nu,omega**2]]/(nu**2+omega**2), in (Q,P) order.

    Gauss-Legendre quadrature integrates the whole real line through a tan
    map, scaled by the geometric mean of the free frequencies. This is a
    fixed-order double-precision integral, not an adaptive error estimate;
    compare orders for widely separated scales or large external frequencies.
    The result has hbar=1 and no added counterterms or vacuum normalization.
    """
    frequencies, external, nonzero = _validated_inputs(
        frequencies, external_frequencies, blocks, quadrature_order)
    if not nonzero:
        return 0j

    terms = []
    requests = set()
    for partition in ordered_partitions(len(external)):
        if not all(subset in nonzero for subset in partition):
            continue
        preceding = set()
        factors = []
        for subset in partition:
            shift = fsum(external[label] for label in sorted(preceding))
            factors.append((subset, shift))
            preceding.update(subset)
        k = len(partition)
        terms.append(((-1) ** (k + 1) / (2 * k), tuple(factors)))
        if k > 1:
            requests.update(factors)
    if not terms:
        return 0j

    scale = exp(fsum(log(omega) for omega in frequencies) / len(frequencies))
    samples = []
    for nu, weight in _line_nodes(quadrature_order, scale):
        products = {key: _propagator_times_block(frequencies, nu + key[1], nonzero[key[0]])
                    for key in requests}
        integrands = []
        for coefficient, factors in terms:
            if len(factors) == 1:
                block = nonzero[factors[0][0]]
                # Symmetric B cancels the off-diagonal, odd-in-nu trace exactly.
                trace = _complex_sum((block[2 * i][2 * i] + omega * omega * block[2 * i + 1][2 * i + 1])
                                     / (nu * nu + omega * omega)
                                     for i, omega in enumerate(frequencies))
            else:
                trace = _trace_product([products[key] for key in factors])
            integrands.append(coefficient * trace)
        samples.append(weight * _complex_sum(integrands))
    return _complex_sum(samples)
