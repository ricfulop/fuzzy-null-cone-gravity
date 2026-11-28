"""Exact localized sources: Gaussian polynomials times zero or one heat moment.

GaussianHeatSource.terms is an immutable mapping
    (beta, a, m, powers) -> Fraction
whose entry represents c*x**powers*exp(-beta*r^2)*I_m(a*r^2).
Here I_m(z)=integral_0^1 t**(2*m)*exp(-z*t*t) dt, beta>0 and a>0.
The PURE GAUSSIAN sentinel is exactly (a,m)=(None,None); that entry instead
represents c*x**powers*exp(-beta*r^2), with no heat factor or integration.
Zero a, negative m, or a lone None are not alternate sentinels.

All coefficients and rates are exact rational numbers (ints or Rational,
normalized to Fraction); float, Decimal, bool, infinite and NaN input are
not accepted as exact scalars. Cartesian powers and derivative orders are
nonnegative built-in integers. The external beta exponent is NOT the
Gaussian derivative rate 2*beta used by some polynomial helper routines.

Construction is available as class methods and same-named module functions:
    from_polynomial(poly, beta)
    from_heat(field, poly, beta)  # existing gaussian_green_fields.HeatField
No Green inverse, source construction, evaluation or checks run on import.
This module supplies algebra only: no numerical evaluation or integration.

grouped_terms() returns a read-only mapping (beta,a,m) -> read-only polynomial.
A heat group explicitly denotes
    integral_0^1 t**(2*m) P(x) exp[-(beta+a*t*t)*r^2] dt.
A pure group denotes P(x) exp(-beta*r^2), without an auxiliary variable.
Groups merge equal keys only; no recurrence identities are used to identify
different heat groups. Products of two heat sources are outside this class.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field as dataclass_field
from fractions import Fraction as F
from numbers import Rational
from types import MappingProxyType
from typing import Optional, Tuple

from gaussian_green_fields import HeatField


Powers = Tuple[int, int, int]
SourceKey = Tuple[F, Optional[F], Optional[int], Powers]
GroupKey = Tuple[F, Optional[F], Optional[int]]
ZERO_POWERS = (0, 0, 0)


def _rational(value, name):
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(name + " must be an exact integer or rational, not float/Decimal")
    # Fraction has finite numerator/denominator; zero denominator and any
    # nonfinite or unsupported representation cannot enter the stored form.
    return F(value)


def _positive_rate(value, name):
    value = _rational(value, name)
    if value <= 0:
        raise ValueError(name + " must be strictly positive")
    return value


def _powers(value, name="Cartesian powers"):
    try:
        value = tuple(value)
    except TypeError as exc:
        raise TypeError(name + " must be a sequence of three nonnegative integers") from exc
    if len(value) != 3 or any(type(n) is not int or n < 0 for n in value):
        raise ValueError(name + " must contain three nonnegative integers")
    return value


def _accumulate(out, key, coefficient):
    coefficient = out.get(key, F(0)) + coefficient
    if coefficient:
        out[key] = coefficient
    else:
        out.pop(key, None)


def _polynomial(poly):
    if not isinstance(poly, Mapping):
        raise TypeError("poly must map Cartesian powers to exact coefficients")
    result = {}
    for powers, coefficient in poly.items():
        powers = _powers(powers)
        coefficient = _rational(coefficient, "polynomial coefficient")
        _accumulate(result, powers, coefficient)
    return result


def _key(key):
    if not isinstance(key, tuple) or len(key) != 4:
        raise ValueError("source keys must be (beta,a,m,powers) tuples")
    beta, a, m, powers = key
    beta = _positive_rate(beta, "beta")
    powers = _powers(powers)
    if a is None:
        if m is not None:
            raise ValueError("pure Gaussian terms require both a=None and m=None")
    else:
        a = _positive_rate(a, "heat rate a")
        if type(m) is not int or m < 0:
            raise ValueError("heat terms require a nonnegative integer m")
    return beta, a, m, powers


def _group_sort_key(key):
    beta, a, m = key
    # Pure groups precede heat groups at equal beta; never compare None to F.
    return beta, a is not None, F(0) if a is None else a, -1 if m is None else m


@dataclass(frozen=True)
class GaussianHeatSource:
    """Immutable finite sum with a positive external Gaussian in every term.

    Construction defensively copies terms and normalizes their exact keys.
    Returned operations create new sources without modifying either operand.
    The empty mapping represents zero and is available as GaussianHeatSource().
    """

    terms: Mapping[SourceKey, F] = dataclass_field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.terms, Mapping):
            raise TypeError("terms must be a mapping from source keys to exact coefficients")
        normalized = {}
        for key, coefficient in self.terms.items():
            key = _key(key)
            coefficient = _rational(coefficient, "source coefficient")
            _accumulate(normalized, key, coefficient)
        object.__setattr__(self, "terms", MappingProxyType(normalized))

    @classmethod
    def from_polynomial(cls, poly, beta):
        """P(x)*exp(-beta*r^2), using the explicit (None,None) pure sentinel."""
        beta = _positive_rate(beta, "beta")
        return cls({(beta, None, None, powers): c for powers, c in _polynomial(poly).items()})

    @classmethod
    def from_heat(cls, field, poly, beta):
        """P(x)*exp(-beta*r^2)*field(x), preserving the HeatField's exact a,m.

        beta is the added, positive external Gaussian exponent; field.a is
        inside I_m and must not be folded into beta or counted twice.
        """
        if not isinstance(field, HeatField):
            raise TypeError("field must be an existing gaussian_green_fields.HeatField")
        beta = _positive_rate(beta, "beta")
        a = _positive_rate(field.a, "heat rate a")
        poly = _polynomial(poly)
        terms = {}
        for p, c in poly.items():
            for (q, m), d in field.terms.items():
                powers = tuple(left + right for left, right in zip(p, q))
                _accumulate(terms, (beta, a, m, powers), c * d)
        return cls(terms)

    def add(self, other, factor=F(1)):
        """Return self + factor*other; no common rate is required."""
        if not isinstance(other, GaussianHeatSource):
            raise TypeError("other must be a GaussianHeatSource")
        factor = _rational(factor, "addition factor")
        terms = dict(self.terms)
        for key, coefficient in other.terms.items():
            _accumulate(terms, key, factor * coefficient)
        return type(self)(terms)

    def scale(self, coefficient):
        coefficient = _rational(coefficient, "scale factor")
        return type(self)({key: coefficient * value for key, value in self.terms.items()})

    def derivative(self, axis, order=1):
        """Exact partial derivative: dP-2*beta*x_i*P, plus the heat derivative.

        For heat terms the additional contribution is -2*a*x_i*P*I_(m+1).
        Pure Gaussian terms never acquire a heat factor. Zero-order requests
        return this immutable source unchanged.
        """
        if type(axis) is not int or not 0 <= axis < 3:
            raise ValueError("axis must be 0, 1, or 2")
        if type(order) is not int or order < 0:
            raise ValueError("order must be a nonnegative integer")
        if order == 0:
            return self
        terms = dict(self.terms)
        for _ in range(order):
            differentiated = {}
            for (beta, a, m, powers), coefficient in terms.items():
                if powers[axis]:
                    lower = list(powers)
                    lower[axis] -= 1
                    _accumulate(differentiated, (beta, a, m, tuple(lower)), coefficient * powers[axis])
                upper = list(powers)
                upper[axis] += 1
                upper = tuple(upper)
                _accumulate(differentiated, (beta, a, m, upper), -2 * beta * coefficient)
                if a is not None:
                    _accumulate(differentiated, (beta, a, m + 1, upper), -2 * a * coefficient)
            terms = differentiated
            if not terms:
                break
        return type(self)(terms)

    def diff(self, multi=ZERO_POWERS):
        """Apply the Cartesian derivative multi-index, without numerical differences."""
        multi = _powers(multi, "derivative multi-index")
        result = self
        for axis, order in enumerate(multi):
            if order:
                result = result.derivative(axis, order)
        return result

    def multiply_polynomial(self, poly, beta=0):
        """Multiply by P(x)*exp(-beta*r^2), preserving at most one heat factor.

        beta is an exact additive exponent shift. Zero means an ordinary
        polynomial multiplier. A negative shift is allowed only if every
        resulting source exponent stays strictly positive; this condition
        is checked even when the supplied polynomial is zero.
        """
        shift = _rational(beta, "Gaussian exponent shift")
        poly = _polynomial(poly)
        for original, _, _, _ in self.terms:
            if original + shift <= 0:
                raise ValueError("multiplication must leave every total beta strictly positive")
        terms = {}
        for (original, a, m, powers), coefficient in self.terms.items():
            for q, c in poly.items():
                combined = tuple(left + right for left, right in zip(powers, q))
                _accumulate(terms, (original + shift, a, m, combined), coefficient * c)
        return type(self)(terms)

    def grouped_terms(self):
        """Read-only (beta,a,m) -> read-only P for Gaussian-mixture integration.

        Heat groups have weight t**(2*m) and exponent beta+a*t*t on [0,1].
        Pure (beta,None,None) groups have no auxiliary integration variable.
        Both mapping layers are defensive snapshots with deterministic order.
        """
        groups = {}
        for (beta, a, m, powers), coefficient in self.terms.items():
            groups.setdefault((beta, a, m), {})[powers] = coefficient
        return MappingProxyType({key: MappingProxyType(dict(sorted(groups[key].items())))
                                 for key in sorted(groups, key=_group_sort_key)})

    @property
    def max_degree(self):
        """Maximum Cartesian polynomial degree; zero for an empty source."""
        return max((sum(powers) for _, _, _, powers in self.terms), default=0)

    @property
    def max_powers(self):
        """Maximum exponent in each Cartesian coordinate; (0,0,0) for zero."""
        return tuple(max((powers[i] for _, _, _, powers in self.terms), default=0) for i in range(3))

    @property
    def max_m(self):
        """Largest heat-moment index, or None if there are no heat terms."""
        return max((m for _, a, m, _ in self.terms if a is not None), default=None)

    def statistics(self):
        """JSON-serializable structural information, not numerical diagnostics."""
        pure = sum(a is None for _, a, _, _ in self.terms)
        groups = {(beta, a, m) for beta, a, m, _ in self.terms}
        return {
            "terms": len(self.terms), "pure_gaussian_terms": pure,
            "heat_terms": len(self.terms) - pure, "groups": len(groups),
            "is_zero": not self.terms,
            "min_degree": min((sum(powers) for _, _, _, powers in self.terms), default=None),
            "max_degree": self.max_degree, "max_powers": list(self.max_powers),
            "max_m": self.max_m,
            "gaussian_rates": [str(beta) for beta in sorted({key[0] for key in self.terms})],
            "heat_rates": [str(a) for a in sorted({key[1] for key in self.terms if key[1] is not None})],
            "pure_gaussian_sentinel": {"a": None, "m": None},
        }

    def __bool__(self):
        return bool(self.terms)

    def __len__(self):
        return len(self.terms)


def from_polynomial(poly, beta):
    """Module-level constructor for a pure polynomial times exp(-beta*r^2)."""
    return GaussianHeatSource.from_polynomial(poly, beta)


def from_heat(field, poly, beta):
    """Module-level constructor preserving the existing HeatField's rate/moments."""
    return GaussianHeatSource.from_heat(field, poly, beta)
