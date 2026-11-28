"""Exact first Green inverses of polynomial-Gaussian sources on R^3.

A source is P(x) exp(-a |x|^2), with positive rational a. Polynomial keys
are Cartesian triples, as in reduced_gaunt_exact. A HEAT polynomial maps
(powers, m) to a Fraction and denotes sum c*x**powers*I_m(a |x|^2), where
I_m(z) = integral_0^1 t**(2*m) exp(-z*t*t) dt. Algebraic Green tails are
retained; there is no outer-radius cutoff or quadrature over this inverse.

The Hermite basis gaussian_multi(alpha, 2*a) represents derivatives of
exp(-a |x|^2). Constant-coefficient Green inverses commute with these
source derivatives. The bases are

    C[g] = I_0/(16*a),                         C = (-Delta)^(-1)/8,
    W_ij[g] = -((7*I_0+I_1)*delta_ij
                 + 2*a*x_i*x_j*(I_1-I_2))/(16*a),
    Delta_L W_i = Delta W_i + (1/3)*partial_i partial_j W_j.

To check the latter, write W_ij = u*delta_ij + v*x_i*x_j. The identity
(2*m+1)*I_m - 2*z*I_(m+1) = exp(-z) gives

    partial_i W_ij = (3/4)*x_j*I_1,
    Delta W_ij = delta_ij*(exp(-z)-I_1/4) + (a/2)*x_i*x_j*I_2.

The remaining York term cancels both I terms. Similarly Delta I_0 =
-2*a*exp(-z). The bases are regular at zero and decay at infinity.

Public construction: scalar_field(P, a), york_field((P0, P1, P2), a).
Both accept an additional derivative multi-index. HeatField.diff,
gradient, hessian and laplacian expose analytic derivative fields.
evaluate_fields shares integral tables; evaluate_source evaluates the
original Gaussian source independently. scalar_residual/york_residual
return exact recurrence-normalized residuals (two empty dictionaries
mean zero), for optional validation by the tracked caller.

No calculations, residual checks, files, or runs are initiated at import.
Coefficients remain exact; point evaluation uses binary64 and may lose
relative accuracy for strongly cancelling high derivatives or far tails.
This module is a first-inverse building block, not an inverse for products
of HEAT fields or a complete H5/H6 spatial integration backend.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import erf, exp, fsum, isfinite, log, log1p, pi, sqrt
from numbers import Rational
from types import MappingProxyType
from typing import Mapping, Tuple

from reduced_gaunt_exact import add, gaussian_derivative, gaussian_multi, multiply

Powers = Tuple[int, int, int]
HeatKey = Tuple[Powers, int]
ZERO_POWERS = (0, 0, 0)
R2 = {(2, 0, 0): F(1), (0, 2, 0): F(1), (0, 0, 2): F(1)}


def _rational(value):
    if not isinstance(value, Rational):
        raise TypeError("Exact coefficients and rates must be integers or rationals")
    return F(value)


def _rate(a):
    a = _rational(a)
    if a <= 0:
        raise ValueError("Gaussian exponent a must be positive")
    return a


def _orders(powers):
    powers = tuple(powers)
    if len(powers) != 3 or any(type(n) is not int or n < 0 for n in powers):
        raise ValueError("Expected three nonnegative integer Cartesian powers")
    return powers


def _poly(poly):
    out = {}
    for powers, coefficient in poly.items():
        powers = _orders(powers)
        coefficient = _rational(coefficient)
        if coefficient:
            out[powers] = coefficient
    return out


def _accumulate(target, terms, factor=F(1)):
    """Accumulate either polynomial dictionary format, removing exact zeros."""
    for key, value in terms.items():
        value = target.get(key, F(0)) + factor * value
        if value:
            target[key] = value
        else:
            target.pop(key, None)


def hermite_decomposition(poly, a=F(1)):
    """P = sum c_alpha gaussian_multi(alpha, 2*a), with exact reconstruction.

    Descending degree is triangular because the leading coefficient of
    gaussian_multi(alpha, 2*a) is (-2*a)**sum(alpha). Cached polynomials
    belonging to reduced_gaunt_exact are never mutated.
    """
    a = _rate(a)
    original = _poly(poly)
    remainder = dict(original)
    coefficients = {}
    while remainder:
        alpha = max(remainder, key=lambda powers: (sum(powers), powers))
        coefficient = remainder[alpha] / ((-2 * a) ** sum(alpha))
        coefficients[alpha] = coefficient
        remainder = add(remainder, gaussian_multi(alpha, 2 * a), -coefficient)
    reconstructed = {}
    for alpha, coefficient in coefficients.items():
        reconstructed = add(reconstructed, gaussian_multi(alpha, 2 * a), coefficient)
    if reconstructed != original:
        raise AssertionError("Exact Gaussian-Hermite reconstruction failed")
    return coefficients


def heat_integrals(z, mmax):
    """Return (I_0(z), ..., I_mmax(z)) without subtractive small-z recurrence.

    Below mmax+5, seed the highest moment by the positive Kummer series
    exp(-z)/(2*m+1) * sum_k z**k/(m+3/2)_k, then recurse downward.
    Above that threshold, erf supplies I_0 and upward recursion is used.
    A logarithmic downward variant avoids an underflowed seed when both
    z and the requested degree are large. Underflow of an individual
    correctly tiny final moment to zero is ordinary binary64 behavior.
    """
    z = float(z)
    if not isfinite(z) or z < 0:
        raise ValueError("z must be finite and nonnegative")
    if type(mmax) is not int or mmax < 0:
        raise ValueError("mmax must be a nonnegative integer")
    if z == 0:
        return tuple(1.0 / (2 * m + 1) for m in range(mmax + 1))
    values = [0.0] * (mmax + 1)
    if z > mmax + 5:
        values[0] = sqrt(pi) * erf(sqrt(z)) / (2 * sqrt(z))
        boundary = exp(-z)
        for m in range(1, mmax + 1):
            values[m] = ((2 * m - 1) * values[m - 1] - boundary) / z / 2
        return tuple(values)

    terms = [1.0]
    term = total = 1.0
    # Ratios decrease, so the next ratio bounds the whole positive tail.
    for k in range(1, 100001):
        term *= z / (mmax + 1.5 + k - 1)
        terms.append(term)
        total += term
        ratio = z / (mmax + 1.5 + k)
        if ratio < 1 and term * ratio / (1 - ratio) <= (2.0 ** -54) * total:
            break
    else:
        raise ArithmeticError("Positive Kummer series did not converge")
    series = fsum(terms) / (2 * mmax + 1)
    if z <= 600:
        boundary = exp(-z)
        values[mmax] = boundary * series
        for m in range(mmax, 0, -1):
            values[m - 1] = (2 * z * values[m] + boundary) / (2 * m - 1)
    else:
        log_value = -z + log(series)
        values[mmax] = exp(log_value)
        for m in range(mmax, 0, -1):
            left, right = log(2 * z) + log_value, -z
            high, low = max(left, right), min(left, right)
            log_value = high + log1p(exp(low - high)) - log(2 * m - 1)
            values[m - 1] = exp(log_value)
    return tuple(values)


@dataclass(frozen=True)
class HeatField:
    """Immutable single-rate HEAT field with exact rational coefficients.

    ``terms[(powers, m)]`` multiplies x**powers I_m(a*|x|^2).
    Scalar multiples and sums at the same rate remain in this class;
    multiplication of two fields generally does not.
    """

    a: F
    terms: Mapping[HeatKey, F]

    def __post_init__(self):
        object.__setattr__(self, "a", _rate(self.a))
        terms = {}
        for (powers, m), coefficient in self.terms.items():
            powers = _orders(powers)
            if type(m) is not int or m < 0:
                raise ValueError("HEAT index m must be a nonnegative integer")
            coefficient = _rational(coefficient)
            if coefficient:
                terms[powers, m] = coefficient
        object.__setattr__(self, "terms", MappingProxyType(terms))

    @property
    def max_m(self):
        return max((m for _, m in self.terms), default=0)

    def scale(self, coefficient):
        coefficient = _rational(coefficient)
        return HeatField(self.a, {key: coefficient * c for key, c in self.terms.items()})

    def add(self, other, factor=F(1)):
        if not isinstance(other, HeatField) or other.a != self.a:
            raise ValueError("Addition requires HeatFields with the same exact rate")
        terms = dict(self.terms)
        _accumulate(terms, other.terms, _rational(factor))
        return HeatField(self.a, terms)

    def derivative(self, axis, order=1):
        if type(axis) is not int or not 0 <= axis < 3:
            raise ValueError("axis must be 0, 1, or 2")
        if type(order) is not int or order < 0:
            raise ValueError("order must be a nonnegative integer")
        terms = dict(self.terms)
        for _ in range(order):
            result = {}
            for (powers, m), coefficient in terms.items():
                if powers[axis]:
                    lower = list(powers)
                    lower[axis] -= 1
                    _accumulate(result, {(tuple(lower), m): powers[axis] * coefficient})
                upper = list(powers)
                upper[axis] += 1
                _accumulate(result, {(tuple(upper), m + 1): -2 * self.a * coefficient})
            terms = result
        return HeatField(self.a, terms)

    def diff(self, orders=ZERO_POWERS):
        result = self
        for axis, order in enumerate(_orders(orders)):
            if order:
                result = result.derivative(axis, order)
        return result

    def gradient(self):
        return tuple(self.derivative(axis) for axis in range(3))

    def hessian(self):
        gradient = self.gradient()
        return tuple(tuple(gradient[i].derivative(j) for j in range(3)) for i in range(3))

    def laplacian(self):
        terms = {}
        for axis in range(3):
            _accumulate(terms, self.derivative(axis, 2).terms)
        return HeatField(self.a, terms)

    def evaluate(self, x):
        return evaluate_fields((self,), x)[0]

    __call__ = evaluate


def scalar_base(a=F(1)):
    a = _rate(a)
    return HeatField(a, {(ZERO_POWERS, 0): 1 / (16 * a)})


def york_base(i, j, a=F(1)):
    """Component i of Delta_L^-1 applied to exp(-a*r^2) in source slot j."""
    a = _rate(a)
    if any(type(axis) is not int or not 0 <= axis < 3 for axis in (i, j)):
        raise ValueError("York component indices must be 0, 1, or 2")
    terms = {}
    if i == j:
        terms[ZERO_POWERS, 0] = -F(7) / (16 * a)
        terms[ZERO_POWERS, 1] = -F(1) / (16 * a)
    powers = [0, 0, 0]
    powers[i] += 1
    powers[j] += 1
    terms[tuple(powers), 1] = -F(1, 8)
    terms[tuple(powers), 2] = F(1, 8)
    return HeatField(a, terms)


@lru_cache(maxsize=4096)
def _basis_derivative(kind, i, j, a, orders):
    if orders == ZERO_POWERS:
        return scalar_base(a) if kind == "scalar" else york_base(i, j, a)
    axis = next(axis for axis, order in enumerate(orders) if order)
    previous = list(orders)
    previous[axis] -= 1
    return _basis_derivative(kind, i, j, a, tuple(previous)).derivative(axis)


def scalar_field(source, a=F(1), derivatives=ZERO_POWERS):
    """Return partial**derivatives C[P exp(-a*r^2)] as a HeatField."""
    a, derivatives = _rate(a), _orders(derivatives)
    terms = {}
    for alpha, coefficient in hermite_decomposition(source, a).items():
        orders = tuple(n + k for n, k in zip(alpha, derivatives))
        basis = _basis_derivative("scalar", 0, 0, a, orders)
        _accumulate(terms, basis.terms, coefficient)
    return HeatField(a, terms)


def york_field(sources, a=F(1), derivatives=ZERO_POWERS):
    """Return three analytic fields solving Delta_L W = vector Gaussian source.

    ``sources`` contains three polynomial mappings at the same rate a.
    Requested spatial derivatives act on every output component.
    """
    a, derivatives = _rate(a), _orders(derivatives)
    sources = tuple(sources)
    if len(sources) != 3:
        raise ValueError("York inversion needs three source polynomials")
    hermite = tuple(hermite_decomposition(poly, a) for poly in sources)
    fields = []
    for i in range(3):
        terms = {}
        for j, expansion in enumerate(hermite):
            for alpha, coefficient in expansion.items():
                orders = tuple(n + k for n, k in zip(alpha, derivatives))
                basis = _basis_derivative("york", i, j, a, orders)
                _accumulate(terms, basis.terms, coefficient)
        fields.append(HeatField(a, terms))
    return tuple(fields)


def _point(x):
    point = tuple(float(value) for value in x)
    if len(point) != 3 or not all(isfinite(value) for value in point):
        raise ValueError("Point must have three finite Cartesian coordinates")
    return point


def _powers_at(point, maxima):
    tables = []
    for value, degree in zip(point, maxima):
        row = [1.0]
        for _ in range(degree):
            row.append(row[-1] * value)
            if not isfinite(row[-1]):
                raise OverflowError("Cartesian monomial exceeds binary64 range")
        tables.append(row)
    return tables


def evaluate_fields(fields, x):
    """Evaluate a flat sequence of fields at one point, sharing tables per a.

    Flatten gradient/Hessian component arrays before passing them here.
    Rates may differ. This is a scalar-point backend, not a NumPy batch API.
    """
    fields, point = tuple(fields), _point(x)
    if not fields:
        return ()
    maxima = [0, 0, 0]
    by_rate = {}
    for field in fields:
        by_rate[field.a] = max(by_rate.get(field.a, 0), field.max_m)
        for powers, _ in field.terms:
            maxima = [max(old, new) for old, new in zip(maxima, powers)]
    powers = _powers_at(point, maxima)
    r2 = fsum(value * value for value in point)
    integrals = {}
    for a, mmax in by_rate.items():
        z = float(a) * r2
        if not isfinite(z):
            raise OverflowError("a*|x|^2 exceeds binary64 range")
        integrals[a] = heat_integrals(z, mmax)
    result = []
    for field in fields:
        moments = integrals[field.a]
        value = fsum(float(c) * powers[0][p[0]] * powers[1][p[1]]
                     * powers[2][p[2]] * moments[m]
                     for (p, m), c in field.terms.items())
        if not isfinite(value):
            raise OverflowError("HEAT evaluation exceeds binary64 range")
        result.append(value)
    return tuple(result)


def source_derivative(source, a=F(1), derivatives=ZERO_POWERS):
    """Polynomial of partial**derivatives[P exp(-a*r^2)], with the same a."""
    a, derivatives = _rate(a), _orders(derivatives)
    poly = _poly(source)
    for axis, order in enumerate(derivatives):
        for _ in range(order):
            poly = gaussian_derivative(poly, axis, 2 * a)
    return poly


def evaluate_source(source, x, a=F(1), derivatives=ZERO_POWERS):
    """Evaluate the Gaussian source directly, without invoking a Green field."""
    a, point = _rate(a), _point(x)
    poly = source_derivative(source, a, derivatives)
    maxima = tuple(max((p[i] for p in poly), default=0) for i in range(3))
    powers = _powers_at(point, maxima)
    z = float(a) * fsum(value * value for value in point)
    if not isfinite(z):
        raise OverflowError("a*|x|^2 exceeds binary64 range")
    value = fsum(float(c) * powers[0][p[0]] * powers[1][p[1]] * powers[2][p[2]]
                 for p, c in poly.items()) * exp(-z)
    if not isfinite(value):
        raise OverflowError("Gaussian source evaluation exceeds binary64 range")
    return value


def _divide_r2(poly):
    """Exact polynomial division by x0^2+x1^2+x2^2, remainder degree_x0 < 2."""
    quotient, remainder = {}, dict(poly)
    while True:
        candidates = [powers for powers in remainder if powers[0] >= 2]
        if not candidates:
            return quotient, remainder
        powers = max(candidates)
        reduced = (powers[0] - 2, powers[1], powers[2])
        coefficient = remainder[powers]
        _accumulate(quotient, {reduced: coefficient})
        remainder = add(remainder, multiply(R2, {reduced: coefficient}), -1)


def heat_normal_form(field):
    """Return (HEAT remainder, Gaussian polynomial) using exact I recurrences.

    The pair represents the original field as HEAT + P*exp(-a*r^2).
    At each m>=1 divide its polynomial by r^2 and replace
    r^2 I_m = ((2*m-1) I_(m-1) - exp(-a*r^2))/(2*a).
    Nonzero remainders are retained rather than declared identities.
    """
    groups = {}
    for (powers, m), coefficient in field.terms.items():
        groups.setdefault(m, {})[powers] = coefficient
    exponential = {}
    for m in range(field.max_m, 0, -1):
        quotient, groups[m] = _divide_r2(groups.get(m, {}))
        if quotient:
            groups[m - 1] = add(groups.get(m - 1, {}), quotient, F(2 * m - 1) / (2 * field.a))
            exponential = add(exponential, quotient, -1 / (2 * field.a))
    heat = {(powers, m): c for m, poly in groups.items() for powers, c in poly.items()}
    return heat, exponential


def scalar_residual(field, source):
    """Exact Delta(field)+source/8 residual; return ({}, {}) iff reduced to zero.

    For a differentiated Green field supply the corresponding polynomial
    from source_derivative, since its Gaussian exponent is unchanged.
    """
    heat, exponential = heat_normal_form(field.laplacian())
    return heat, add(exponential, _poly(source), F(1, 8))


def york_operator(fields):
    """Return three analytic fields Delta W_i + partial_i div(W)/3."""
    fields = tuple(fields)
    if len(fields) != 3 or len({field.a for field in fields}) != 1:
        raise ValueError("York operator needs three fields at one exact rate")
    divergence = HeatField(fields[0].a, {})
    for j in range(3):
        divergence = divergence.add(fields[j].derivative(j))
    return tuple(fields[i].laplacian().add(divergence.derivative(i), F(1, 3))
                 for i in range(3))


def york_residual(fields, sources):
    """Exact Delta_L W - source residuals; each zero component is ({}, {})."""
    sources = tuple(sources)
    if len(sources) != 3:
        raise ValueError("York residual needs three source polynomials")
    residuals = []
    for field, source in zip(york_operator(fields), sources):
        heat, exponential = heat_normal_form(field)
        residuals.append((heat, add(exponential, _poly(source), -1)))
    return tuple(residuals)
