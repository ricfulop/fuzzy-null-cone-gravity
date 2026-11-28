"""Independent rational Cartesian Laurent algebra and complete Kelvin traces.

CartesianLaurent.terms[(Cartesian_powers, radial_power, log_power)] = Fraction
means c*x**powers*r**radial_power*log(r/anchor)**log_power. This is exterior
algebra at positive radius, not a prescription for origin evaluation. The raw
dictionary is NOT a normal form modulo x*x+y*y+z*z=r*r; project_exterior
combines these identities exactly before checking the declared angular band.

project_exterior(field,ell_max,...)->(mode_mapping,anchor,statistics), with
mode_mapping[(l,m)][(p,q)] multiplying r**(-p)*log(r/anchor)**q*H_lm(n).
The algebra permits bounded/growing field powers; the radial inverse module
separately enforces its SOURCE-domain conditions. No mode is silently removed.
Different anchors cannot be combined if either operand has logarithms.

kelvin_exterior_basis(l,m,k,j) returns three rational CartesianLaurent fields
including the FULL normalization for the unnormalized-H source moment
  M_lmk = integral_0^infinity r**(l+2*k+2)*f_lm(r) dr,
where source component j is f_lm(r)*H_lm(n). The basis is
  (-1)**l/[8*2**k*k!*(2*l+2*k+1)!!] * H_lm(partial)*Delta**k*Kbar_ij,
  Kbar_ij = -7*delta_ij/r - x_i*x_j/r**3,  k=0 or 1 only.
This equals 4*pi*(-1)**l/[2**k*k!*(2*l+2*k+1)!!] times derivatives of
the physical Kelvin kernel Kbar/(32*pi). No extra h or pi belongs in M.
For normalized Y source coefficients first divide M by sqrt(4*pi*h_lm).

Only the algebraic exterior is represented; localized-source remainder and
moment uncertainties must be supplied/bounded separately by the caller.
No production inverse backend is imported. No scientific call runs on import.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial, fsum, isfinite, log, sqrt
from types import MappingProxyType

from solid_harmonic_reference import (
    MAX_DEGREE, MAX_ELL, _ExactWork, _accumulate, _add, _fischer,
    _integer, _mode, _multiply, _odd_factorial, _powers, _r2_power,
    _rational, solid_harmonic,
)


MAX_RADIAL_POWER = 128
MAX_LOG_POWER = 4


def _key(key):
    if not isinstance(key, tuple) or len(key) != 3:
        raise ValueError("Laurent keys must be (Cartesian powers, radial power, log power)")
    powers, radial, logarithm = key
    powers = _powers(powers)
    if sum(powers) > MAX_DEGREE:
        raise ValueError("Cartesian Laurent polynomial exceeds degree32")
    if type(radial) is not int or abs(radial) > MAX_RADIAL_POWER:
        raise ValueError("Cartesian Laurent radial exponent exceeds its declared range")
    _integer(logarithm, "log power", maximum=MAX_LOG_POWER)
    return powers, radial, logarithm


def _put(out, key, value, work):
    work.charge()
    key = _key(key)
    _accumulate(out, key, value)
    work.check(out)


def _derivative_terms(terms, axis, work):
    result = {}
    for (powers, radial, logarithm), coefficient in terms.items():
        if powers[axis]:
            lower = list(powers)
            lower[axis] -= 1
            _put(result, (tuple(lower), radial, logarithm), coefficient*powers[axis], work)
        raised = list(powers)
        raised[axis] += 1
        if radial:
            _put(result, (tuple(raised), radial-2, logarithm), coefficient*radial, work)
        if logarithm:
            _put(result, (tuple(raised), radial-2, logarithm-1), coefficient*logarithm, work)
    return result


def _differentiate(terms, orders, work):
    result = dict(terms)
    work.check(result)
    for axis, order in enumerate(orders):
        for _ in range(order):
            result = _derivative_terms(result, axis, work)
    return result


def _laplacian_terms(terms, work):
    result = {}
    for axis in range(3):
        differentiated = _derivative_terms(_derivative_terms(terms, axis, work), axis, work)
        for key, value in differentiated.items():
            _put(result, key, value, work)
    return result


@dataclass(frozen=True)
class CartesianLaurent:
    terms: Mapping
    anchor: F = 1

    def __post_init__(self):
        if not isinstance(self.terms, Mapping) or len(self.terms) > 100_000:
            raise ValueError("Laurent terms must be a mapping with at most100000 entries")
        anchor = _rational(self.anchor, "log anchor")
        if anchor <= 0:
            raise ValueError("Log anchor must be strictly positive")
        normalized = {}
        for key, coefficient in self.terms.items():
            _accumulate(normalized, _key(key), _rational(coefficient))
        object.__setattr__(self, "terms", MappingProxyType(normalized))
        object.__setattr__(self, "anchor", anchor)

    @property
    def has_logs(self):
        return any(key[2] for key in self.terms)

    def at_anchor(self, anchor):
        anchor = _rational(anchor, "log anchor")
        if anchor <= 0:
            raise ValueError("Log anchor must be strictly positive")
        if self.has_logs and anchor != self.anchor:
            raise ValueError("Changing a logarithmic anchor needs nonrational log-ratio coefficients")
        return CartesianLaurent(self.terms, anchor)

    def _common_anchor(self, other):
        if not isinstance(other, CartesianLaurent):
            raise TypeError("Expected another CartesianLaurent")
        if self.anchor != other.anchor and (self.has_logs or other.has_logs):
            raise ValueError("Cannot merge different log anchors with exact Fraction coefficients")
        return self.anchor

    def scale(self, factor):
        factor = _rational(factor)
        return CartesianLaurent({key: factor*c for key, c in self.terms.items()}, self.anchor)

    def add(self, other, factor=1, *, max_products=5_000_000, max_terms=100_000):
        anchor = self._common_anchor(other)
        factor = _rational(factor)
        work = _ExactWork(max_products, max_terms)
        result = dict(self.terms)
        work.check(result)
        for key, value in other.terms.items():
            _put(result, key, factor*value, work)
        return CartesianLaurent(result, anchor)

    def multiply(self, other, *, max_products=5_000_000, max_terms=100_000):
        anchor = self._common_anchor(other)
        work = _ExactWork(max_products, max_terms)
        result = {}
        for (p, nu, k), c in self.terms.items():
            for (q, mu, j), d in other.terms.items():
                _put(result, (tuple(a+b for a, b in zip(p, q)), nu+mu, k+j), c*d, work)
        return CartesianLaurent(result, anchor)

    def derivative(self, axis, order=1, *, max_products=5_000_000, max_terms=100_000):
        _integer(axis, "axis", maximum=2)
        _integer(order, "order", maximum=MAX_ELL+2)
        orders = [0, 0, 0]
        orders[axis] = order
        work = _ExactWork(max_products, max_terms)
        result = _differentiate(self.terms, orders, work)
        return CartesianLaurent(result, self.anchor)

    def laplacian(self, *, max_products=5_000_000, max_terms=100_000):
        work = _ExactWork(max_products, max_terms)
        return CartesianLaurent(_laplacian_terms(self.terms, work), self.anchor)

    def project(self, ell_max, *, max_products=5_000_000, max_terms=100_000):
        return project_exterior(self, ell_max, max_products=max_products, max_terms=max_terms)

    def evaluate(self, point):
        point = tuple(float(value) for value in point)
        if len(point) != 3 or any(not isfinite(value) for value in point):
            raise ValueError("Point must have three finite Cartesian coordinates")
        r2 = fsum(value*value for value in point)
        if not isfinite(r2) or r2 <= 0:
            raise ValueError("Laurent evaluation requires a positive finite radius")
        radius = sqrt(r2)
        logarithm = log(radius/float(self.anchor)) if self.has_logs else 0.0
        value = fsum(float(c)*point[0]**p[0]*point[1]**p[1]*point[2]**p[2]
                     *radius**nu*(logarithm**q if q else 1.0)
                     for (p, nu, q), c in self.terms.items())
        if not isfinite(value):
            raise ArithmeticError("Nonfinite Cartesian Laurent value")
        return value


class ExteriorSupportError(ValueError):
    def __init__(self, remainder, statistics):
        super().__init__("Nonzero exact exterior angular remainder outside ell_max")
        self.remainder = MappingProxyType({key: MappingProxyType(dict(poly)) for key, poly in remainder.items()})
        self.statistics = statistics


def project_exterior(field, ell_max, *, max_products=5_000_000, max_terms=100_000):
    """Return (H mode power/log mapping, exact anchor, work statistics).

    Terms with the same total radial power, log power and Cartesian parity
    are lifted to a common homogeneous polynomial before Fischer projection.
    This resolves r² identities and cross-representation cancellations before
    applying the exact angular-band gate. No high angular band is constructed.
    """
    if not isinstance(field, CartesianLaurent):
        raise TypeError("Expected a CartesianLaurent")
    _integer(ell_max, "ell_max", maximum=MAX_ELL)
    work = _ExactWork(max_products, max_terms)
    groups = {}
    for (powers, radial, logarithm), coefficient in field.terms.items():
        degree = sum(powers)
        group = groups.setdefault((degree+radial, logarithm, degree % 2), {})
        work.charge()
        _accumulate(group.setdefault(degree, {}), powers, coefficient)
    modes, remainders = {}, {}
    for group_key, homogeneous in sorted(groups.items()):
        radial_power, logarithm, _ = group_key
        degree = max(homogeneous)
        lifted = {}
        for current_degree, poly in homogeneous.items():
            product = _multiply(poly, _r2_power((degree-current_degree)//2), work)
            lifted = _add(lifted, product, work=work)
        coefficients, remainder = _fischer(lifted, ell_max, work)
        if remainder:
            remainders[group_key] = remainder
        for (ell, m, _), coefficient in coefficients.items():
            work.charge()
            _accumulate(modes.setdefault((ell, m), {}), (-radial_power, logarithm), coefficient)
            work.check(modes[ell, m])
    modes = {mode: MappingProxyType(dict(sorted(terms.items()))) for mode, terms in sorted(modes.items()) if terms}
    statistics = {"source_cartesian_terms": len(field.terms), "homogeneous_shell_groups": len(groups),
                  "ell_max": ell_max, "all_m": True, "mode_count": len(modes),
                  "mode_power_log_terms": sum(len(row) for row in modes.values()),
                  "remainder_groups": len(remainders), "anchor": str(field.anchor),
                  "coefficient_convention": "unnormalized H; c*r^(-p)*log(r/anchor)^q",
                  "work": work.statistics()}
    if remainders:
        raise ExteriorSupportError(remainders, statistics)
    return MappingProxyType(modes), field.anchor, statistics


def _kelvin_bar(i, j):
    terms = {}
    if i == j:
        terms[((0, 0, 0), -1, 0)] = F(-7)
    powers = [0, 0, 0]
    powers[i] += 1
    powers[j] += 1
    terms[tuple(powers), -3, 0] = F(-1)
    return terms


def kelvin_exterior_basis(ell, m, trace_level, source_component, *,
                          max_products=5_000_000, max_terms=100_000):
    """Three components for one UNIT unnormalized-H radial source moment.

    Multiply the returned fields ONLY by
      M_lmk=int r^(ell+2*k+2)*f_lm(r)dr.
    The full sign, 4*pi angular factor and physical Kelvin denominator have
    already been combined in the rational normalization. k=0 and1 are both
    required in general. No kernel for k>=2 is silently returned as zero.
    """
    _mode(ell, m)
    _integer(trace_level, "trace_level", maximum=1)
    _integer(source_component, "source_component", maximum=2)
    work = _ExactWork(max_products, max_terms)
    harmonic = solid_harmonic(ell, m).polynomial
    normalization = F((-1)**ell, 8*2**trace_level*factorial(trace_level)*_odd_factorial(2*ell+2*trace_level+1))
    output = []
    for component in range(3):
        base = _kelvin_bar(component, source_component)
        if trace_level:
            base = _laplacian_terms(base, work)
        terms = {}
        for orders, coefficient in harmonic.items():
            differentiated = _differentiate(base, orders, work)
            for key, value in differentiated.items():
                _put(terms, key, normalization*coefficient*value, work)
        output.append(CartesianLaurent(terms))
    return tuple(output)
