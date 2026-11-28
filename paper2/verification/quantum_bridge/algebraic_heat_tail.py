"""Exact exterior expansion of three rate-one HEAT factors.

The input is an already angular-averaged radial density
    4*pi sum c[D,m,n,p] r**(D+2) I_m(r*r) I_n(r*r) I_p(r*r).
No extra Gaussian multiplier is present. This is the new local H6 tail
class, not an extension that silently relaxes gaussian_heat_pairings.

The exact expansion has rational A[k,nu] multiplying
    4*pi*pi**(k/2) erf(r)**k exp(-(3-k)*r*r) r**nu.
Only the k=3 part has an algebraic tail. The scalar/York field falloffs
require its powers nu<=-3, after exact cancellation; the genuine nu=-3
coefficient must be retained. Exterior integrals include that whole
algebraic polynomial plus an explicit Gaussian-remainder inequality.

Point evaluation and bound evaluation use binary64, so the reported
inequality is not an interval-arithmetic certificate. No field, integral,
coefficient, or test is constructed/evaluated at import time.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import erf, exp, fsum, inf, isfinite, log, nextafter, pi
from numbers import Rational
from sys import float_info
from types import MappingProxyType

from gaussian_green_fields import heat_integrals


def _exact(value):
    if not isinstance(value, Rational):
        raise TypeError("Coefficients must be exact rational numbers")
    return F(value)


def _add(target, key, value):
    value = target.get(key, F(0)) + value
    if value:
        target[key] = value
    else:
        target.pop(key, None)


def _log_abs_exact(value):
    return log(abs(value.numerator))-log(value.denominator)


def _float_coefficient(value):
    converted = float(value)
    if not isfinite(converted) or (value and converted == 0):
        raise OverflowError("Nonzero exact coefficient is outside binary64 range")
    return converted


@lru_cache(maxsize=None)
def moment_exterior_coefficients(m):
    """Exact (erf-count, power)->coefficient for I_m(r*r).

    The factor sqrt(pi) accompanies erf-count=1; count=0 accompanies
    exp(-r*r). Recurrence construction avoids factorial overflow.
    """
    if type(m) is not int or m < 0:
        raise ValueError("The moment index must be a nonnegative integer")
    if m == 0:
        return MappingProxyType({(1, -1): F(1, 2)})
    previous = moment_exterior_coefficients(m-1)
    terms = {(k, nu-2): coefficient*F(2*m-1, 2)
             for (k, nu), coefficient in previous.items()}
    _add(terms, (0, -2), -F(1, 2))
    return MappingProxyType(terms)


def _positive_radius(radius):
    radius = float(radius)
    if not isfinite(radius) or radius <= 0:
        raise ValueError("Exterior radius must be finite and strictly positive")
    return radius


def _log_gaussian_tail_bound(power, rate, radius):
    """Log of the proved J_power(rate,R) upper-bound expression."""
    denominator = 2*rate-max(power-1, 0)/(radius*radius)
    if denominator <= 0:
        raise ValueError("Increase the radius: Gaussian tail denominator is not positive")
    return ((power-1)*log(radius)-rate*radius*radius-log(denominator))


def _log_sum(values):
    if not values:
        return None
    largest = max(values)
    return largest+log(fsum(exp(value-largest) for value in values))


def _from_log_bound(value):
    """Keep an underflowed positive bound visible, without claiming certification."""
    if value is None:
        return 0.0, False
    if value > log(float_info.max):
        raise OverflowError("Exterior remainder bound overflows binary64")
    smallest = nextafter(0.0, 1.0)
    if value < log(smallest):
        return smallest, True
    return nextafter(exp(value), inf), False


@dataclass(frozen=True)
class TripleHeatTail:
    """Immutable canonical angular bins for three rate-one HEAT fields."""

    coefficients: object

    def __post_init__(self):
        result = {}
        for key, value in self.coefficients.items():
            if (not isinstance(key, tuple) or len(key) != 4
                    or any(type(n) is not int or n < 0 for n in key)):
                raise ValueError("Bins must be (even degree, m, n, p) nonnegative integers")
            degree, *moments = key
            if degree % 2:
                raise ValueError("Odd Cartesian total degree has zero angular average")
            _add(result, (degree, *sorted(moments)), _exact(value))
        object.__setattr__(self, "coefficients", MappingProxyType(result))

    def add(self, other, factor=F(1)):
        if not isinstance(other, TripleHeatTail):
            raise TypeError("Expected TripleHeatTail")
        terms = dict(self.coefficients)
        factor = _exact(factor)
        for key, value in other.coefficients.items():
            _add(terms, key, factor*value)
        return TripleHeatTail(terms)

    def scale(self, factor):
        factor = _exact(factor)
        return TripleHeatTail({key: factor*value for key, value in self.coefficients.items()})

    def exterior_coefficients(self, max_products=1_000_000):
        """Return immutable exact A[k,nu], including all signed cancellations."""
        if type(max_products) is not int or max_products < 1:
            raise ValueError("max_products must be a positive integer")
        result = {}
        products = 0
        for (degree, *moments), coefficient in self.coefficients.items():
            partial = {(0, degree+2): coefficient}
            for m in moments:
                next_terms = {}
                expansion = moment_exterior_coefficients(m)
                products += len(partial)*len(expansion)
                if products > max_products:
                    raise RuntimeError("Exact exterior-expansion product budget exceeded")
                for (k, nu), value in partial.items():
                    for (other_k, other_nu), other_value in expansion.items():
                        _add(next_terms, (k+other_k, nu+other_nu), value*other_value)
                partial = next_terms
            for key, value in partial.items():
                _add(result, key, value)
        return MappingProxyType(result)

    def evaluate_radial(self, radius):
        """Regular moment representation; includes the radial measure and 4*pi."""
        radius = float(radius)
        if not isfinite(radius) or radius < 0:
            raise ValueError("Radius must be finite and nonnegative")
        if not self.coefficients:
            return 0.0
        moments = heat_integrals(radius*radius,
                                 max(max(key[1:]) for key in self.coefficients))
        return 4*pi*fsum(_float_coefficient(c)*radius**(degree+2)*moments[m]*moments[n]*moments[p]
                         for (degree, m, n, p), c in self.coefficients.items())

    def evaluate_exterior(self, radius, max_products=1_000_000):
        """Exact exterior identity evaluated in binary64, for radius>0 only."""
        radius = _positive_radius(radius)
        expansion = self.exterior_coefficients(max_products)
        er = erf(radius)
        return 4*pi*fsum(_float_coefficient(c)*pi**(k/2)*er**k
                         *exp(-(3-k)*radius*radius)*radius**nu
                         for (k, nu), c in expansion.items())

    def exterior_integral(self, radius, required_decay_power=-3,
                          max_products=1_000_000):
        """Algebraic exterior integral and Gaussian remainder inequality.

        required_decay_power<=-2 rejects nonintegrable or unexpectedly slow
        terms. The H6 application uses -3. No offending bin is discarded.
        Returned value is the algebraic integral, not the exact full tail.
        Its mathematical difference is bounded by the displayed inequality;
        floating evaluation of that inequality is a numerical diagnostic.
        """
        radius = _positive_radius(radius)
        if type(required_decay_power) is not int or required_decay_power > -2:
            raise ValueError("required_decay_power must be an integer <=-2")
        expansion = self.exterior_coefficients(max_products)
        algebraic = {nu: coefficient for (k, nu), coefficient in expansion.items() if k == 3}
        bad = {nu: str(c) for nu, c in algebraic.items() if nu > required_decay_power}
        if bad:
            raise ValueError("Exact algebraic falloff gate failed: "+repr(bad))
        terms = [4*pi**2.5*_float_coefficient(c)*radius**(nu+1)/(-nu-1)
                 for nu, c in algebraic.items()]
        if not all(isfinite(value) for value in terms):
            raise OverflowError("Algebraic tail contribution overflows binary64")
        value, absolute = fsum(terms), fsum(abs(term) for term in terms)
        bounds = []
        for (k, nu), coefficient in expansion.items():
            if k == 3:
                bounds.append(log(12)+2*log(pi)+_log_abs_exact(coefficient)
                              +_log_gaussian_tail_bound(nu-1, 1, radius))
            else:
                bounds.append(log(4)+(1+k/2)*log(pi)+_log_abs_exact(coefficient)
                              +_log_gaussian_tail_bound(nu, 3-k, radius))
        log_bound = _log_sum(bounds)
        bound, underflow = _from_log_bound(log_bound)
        return {
            "radius": radius,
            "algebraic_integral": value,
            "absolute_algebraic_term_sum": absolute,
            "algebraic_cancellation_ratio": absolute/abs(value) if value else None,
            "machine_epsilon_times_absolute_sum": float_info.epsilon*absolute,
            "absolute_remainder_bound": bound,
            "log_absolute_remainder_bound": log_bound,
            "positive_bound_below_binary64_range": underflow,
            "required_decay_power": required_decay_power,
            "largest_algebraic_power": max(algebraic, default=None),
            "algebraic_coefficients": [{"power": nu, "coefficient": str(c)}
                                       for nu, c in sorted(algebraic.items())],
            "input_bins": len(self.coefficients),
            "exterior_bins": len(expansion),
            "interpretation": "Algebraic tail plus a Gaussian remainder inequality; binary64 evaluation is not interval certified",
        }
