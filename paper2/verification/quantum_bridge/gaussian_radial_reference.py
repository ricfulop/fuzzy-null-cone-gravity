"""Independent polynomial-Gaussian source and exterior-response companion.

GaussianRadialSource(polynomial, beta=Fraction(1), *, ell_max, ...) constructs
an exact Fischer representation and GaussianProfile, never a Green solve.
Physical source coefficients multiply the UNNORMALIZED H_lm(n). Beta in this
module means exp(-beta*r*r); gravitational R2,s2,b1 have beta=1, not the
half-Gaussian rate 2 stored by the local-geometry metadata.

.sample(radii) -> physical shell samples; .origin_data(l,m,degree,patch_radius,
derivative_order=4) -> exact source Taylor coefficients and uniform bounds.
.asymptotic_moments() -> exact Fraction coefficients of sqrt(pi/beta), for
radial moments int_0^infinity r^(l+2*trace+2) f_l(r) dr, trace=0,1.
.finite_moments(radius) -> enclosed finite-core/tail moments, not FD targets.
.exterior(radius) -> GaussianExteriorTail with the protocol agreed with the
separate spherical_york_reference wrapper:
  mode_keys, radius, radius_exact (Fraction anchor; wrapper must verify it)
  outside_response_regular_jet(l,m,r,alpha,order=3,side=None)
  origin_shift(l,m,alpha)
  algebraic_response_terms(l,m,alpha)
  response_error_bounds(l,m,r,alpha,order=3,side=None)
  remainder_response_regular_jet(l,m,r,alpha,order=3)
  source_regular_jet(l,m,r,order=1)
  weighted_outer_bounds(), provenance()
  response_remainder_derivative(axis,alpha=1) -> GaussianDivergenceTail

Regular jets are derivatives in s=r*r. Alpha is 1 for P=(-Delta)^-1 and
1/8 for C. The derivative provider is the exterior SOURCE derivative of
[outside response minus its algebraic response], NOT of the full correction.
It supplies both angular channels of div(V)'s Gaussian remainder. Its source
is defined only at r>=R>0; its outside-only inverse is regular for r<R.
Deeper residual differentiation is explicitly unsupported. No inverse
quadrature, production HEAT/Fourier/pairing code, or radial solver is used.

Decimal interval arithmetic encloses moment evaluation and conversion errors.
Gaussian moment intervals use positive recurrences; U0 uses a bounded positive
Kummer series or integration-by-parts (Mills) brackets, not libm erfc. These
analytic/arithmetic bounds do NOT bound FD, interpolation, angular projection,
matched a0, or final energy errors. A caller must preserve those separate
uncertainties and matching diagnostics. Exact all-space source moments must
not replace numerically matched homogeneous continuation coefficients.

Import is inert. All source construction/evaluation belongs to tracked calls.
"""

from collections import OrderedDict
from collections.abc import Mapping
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction as F
from math import factorial, isfinite, nextafter, inf
from numbers import Rational
from types import MappingProxyType

from solid_harmonic_reference import GaussianProfile, fischer_project, solid_harmonic


def _integer(value, name, lower=0, upper=None):
    if type(value) is not int or value < lower or (upper is not None and value > upper):
        raise ValueError(name+" is outside its declared integer range")
    return value


def _exact(value, name):
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(name+" must be an exact real rational")
    return F(value)


def _radius(value):
    if isinstance(value, bool):
        raise TypeError("radius must be a finite nonnegative real number")
    if isinstance(value, Rational):
        result = F(value)
    elif isinstance(value, float) and isfinite(value):
        result = F.from_float(value)
    else:
        raise TypeError("radius must be rational or an explicitly represented binary64 value")
    if result < 0:
        raise ValueError("radius must be nonnegative")
    return result


def _alpha(value):
    value = _radius(value)
    if value <= 0:
        raise ValueError("inverse alpha must be positive")
    return value


def _mode(ell, m):
    _integer(ell, "ell", upper=16)
    if type(m) is not int or abs(m) > ell:
        raise ValueError("m must satisfy -ell<=m<=ell")
    return ell, m


def _put(result, key, value):
    value = result.get(key, F(0))+value
    if value:
        result[key] = value
    else:
        result.pop(key, None)


class _Interval:
    __slots__ = ("lo", "hi", "engine")

    def __init__(self, lo, hi, engine):
        if not lo.is_finite() or not hi.is_finite() or lo > hi:
            raise ArithmeticError("Invalid finite Decimal enclosure")
        self.lo, self.hi, self.engine = lo, hi, engine

    def _other(self, other):
        return self.engine.number(other)

    def __neg__(self):
        return _Interval(self.hi.copy_negate(), self.lo.copy_negate(), self.engine)

    def __add__(self, other):
        other = self._other(other)
        return self.engine.bounds(lambda: self.lo+other.lo, lambda: self.hi+other.hi)

    __radd__ = __add__

    def __sub__(self, other):
        return self+-self._other(other)

    def __rsub__(self, other):
        return self._other(other)+-self

    def __mul__(self, other):
        other = self._other(other)
        return self.engine.bounds(lambda: min(a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)),
                                  lambda: max(a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self._other(other)
        if other.lo <= 0 <= other.hi:
            raise ArithmeticError("Interval division crosses zero")
        return self.engine.bounds(lambda: min(a/b for a in (self.lo, self.hi) for b in (other.lo, other.hi)),
                                  lambda: max(a/b for a in (self.lo, self.hi) for b in (other.lo, other.hi)))

    def __rtruediv__(self, other):
        return self._other(other)/self

    def __pow__(self, power):
        _integer(abs(power), "absolute integer power", upper=256)
        if type(power) is not int:
            raise TypeError("Only integer interval powers are supported")
        if power < 0:
            return self.engine.number(1)/(self**(-power))
        if power == 0:
            return self.engine.number(1)
        value, base, remaining = self.engine.number(1), self, power
        while remaining:
            if remaining % 2:
                value = value*base
            base, remaining = base*base, remaining//2
        return value

    def exp(self):
        with localcontext() as context:
            context.prec = self.engine.precision
            context.rounding = ROUND_HALF_EVEN
            lo, hi = self.lo.exp(), self.hi.exp()
            return _Interval(max(D(0), lo.next_minus()), hi.next_plus(), self.engine)

    def sqrt(self):
        if self.lo < 0:
            raise ArithmeticError("Negative square-root enclosure")
        with localcontext() as context:
            context.prec = self.engine.precision
            lo, hi = self.lo.sqrt(), self.hi.sqrt()
            return _Interval(max(D(0), lo.next_minus()), hi.next_plus(), self.engine)

    def nonnegative(self):
        if self.hi < 0:
            raise ArithmeticError("A known nonnegative moment has negative upper enclosure")
        return _Interval(max(D(0), self.lo), self.hi, self.engine)


class _Engine:
    def __init__(self, *, precision, max_moment_steps, max_series_terms,
                 max_point_evaluations, max_cached_radii, max_z):
        self.precision = _integer(precision, "precision", 32, 256)
        self.max_steps = _integer(max_moment_steps, "max_moment_steps", 1)
        self.max_series = _integer(max_series_terms, "max_series_terms", 1)
        self.max_points = _integer(max_point_evaluations, "max_point_evaluations", 1)
        self.max_cache = _integer(max_cached_radii, "max_cached_radii", 1, 4096)
        self.max_z = _exact(max_z, "max_z")
        if self.max_z <= 0:
            raise ValueError("max_z must be positive")
        self.steps = self.points = self.underflow_displays = 0
        self.cache, self._pi = OrderedDict(), None

    def step(self, count=1):
        if self.steps+count > self.max_steps:
            raise ValueError("Gaussian moments exceed max_moment_steps")
        self.steps += count

    def point(self):
        if self.points >= self.max_points:
            raise ValueError("Gaussian response exceeds max_point_evaluations")
        self.points += 1

    def bounds(self, lower, upper):
        with localcontext() as context:
            context.prec, context.rounding = self.precision, ROUND_FLOOR
            lo = lower()
            context.rounding = ROUND_CEILING
            hi = upper()
        return _Interval(lo, hi, self)

    def number(self, value):
        if isinstance(value, _Interval):
            if value.engine is not self:
                raise TypeError("Cannot combine unrelated interval engines")
            return value
        value = _exact(value, "interval scalar")
        return self.bounds(lambda: D(value.numerator)/D(value.denominator),
                           lambda: D(value.numerator)/D(value.denominator))

    def pi(self):
        if self._pi is None:
            def arctangent(q):
                total, power = F(0), q
                target = F(1, 10**(self.precision+8))
                for k in range(self.max_series):
                    self.step()
                    total += F((-1)**k, (2*k+1)*power)
                    power *= q*q
                    remainder = F(1, (2*k+3)*power)
                    if remainder <= target:
                        lower, upper = (total-remainder, total) if k % 2 == 0 else (total, total+remainder)
                        return _Interval(self.number(lower).lo, self.number(upper).hi, self)
                raise ValueError("Machin pi exceeds max_series_terms")
            self._pi = 16*arctangent(5)-4*arctangent(239)
        return self._pi

    def _u0(self, beta, radius, z, envelope):
        complete = (self.pi()/beta).sqrt()/2
        if radius == 0:
            return complete
        # Positive Kummer expansion of int_0^R exp(-beta*r^2)dr. All
        # subsequent ratios decrease; next/(1-q) bounds the positive tail.
        if z <= self.precision:
            term = total = self.number(1)
            zi = self.number(z)
            tolerance = D(10)**(-(self.precision-8))
            for k in range(self.max_series):
                self.step()
                q = zi/F(2*k+3, 2)
                following = term*q
                if q.hi < 1:
                    geometric = following/(1-q)
                    tail = _Interval(D(0), geometric.hi, self)
                    if tail.hi <= tolerance*total.lo:
                        lower_integral = radius*envelope*(total+tail)
                        answer = complete-lower_integral
                        coarse_upper = envelope/(2*beta*radius)
                        return _Interval(max(D(0), answer.lo), min(answer.hi, coarse_upper.hi), self)
                total, term = total+following, following
            raise ValueError("Positive Gaussian moment series exceeds max_series_terms")
        # Alternating integration-by-parts brackets are valid even though
        # the Mills series is asymptotic. Stop at sufficient width or before
        # the next term grows; never assume convergence of this series.
        prefactor = envelope/(2*beta*radius)
        best = _Interval(D(0), prefactor.hi, self)
        term = total = self.number(1)
        zi = self.number(z)
        tolerance = D(10)**(-(self.precision-8))
        for k in range(self.max_series):
            self.step()
            following = term*F(2*k+1, 2)/zi
            adjacent = total-following if k % 2 == 0 else total+following
            lower, upper = (adjacent, total) if k % 2 == 0 else (total, adjacent)
            bracket = prefactor*_Interval(lower.lo, upper.hi, self)
            best = _Interval(max(best.lo, bracket.lo, D(0)), min(best.hi, bracket.hi), self)
            if best.hi-best.lo <= tolerance*best.hi or following.lo >= term.hi:
                return best
            total, term = adjacent, following
        raise ValueError("Gaussian Mills brackets exceed max_series_terms")

    def moments(self, beta, radius, maximum):
        if beta*radius*radius > self.max_z:
            raise ValueError("Gaussian moment argument exceeds max_z")
        key = beta, radius, maximum
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        _integer(maximum, "maximum moment", upper=68)
        envelope = self.number(-beta*radius*radius).exp()
        rows = [self._u0(beta, radius, beta*radius*radius, envelope)]
        if maximum:
            rows.append(envelope/(2*beta))
        for n in range(2, maximum+1):
            self.step()
            rows.append((radius**(n-1)*envelope+(n-1)*rows[n-2])/(2*beta))
        result = tuple(rows)
        self.cache[key] = result
        if len(self.cache) > self.max_cache:
            self.cache.popitem(last=False)
        return result

    def value(self, interval):
        with localcontext() as context:
            context.prec = self.precision
            midpoint = (interval.lo+interval.hi)/2
        result = float(midpoint)
        if not isfinite(result):
            raise ArithmeticError("Gaussian response is outside finite binary64 display range")
        if result == 0 and midpoint != 0:
            self.underflow_displays += 1
        return result

    def upper(self, value):
        if value < 0 or not value.is_finite():
            raise ArithmeticError("Error bound must be finite and nonnegative")
        result = float(value)
        if not isfinite(result):
            raise ArithmeticError("Gaussian bound exceeds binary64 range")
        if D.from_float(result) < value:
            result = nextafter(result, inf)
        if not isfinite(result):
            raise ArithmeticError("Upward Gaussian bound conversion exceeds binary64 range")
        return result

    def error(self, interval):
        value = D.from_float(self.value(interval))
        with localcontext() as context:
            context.prec, context.rounding = self.precision, ROUND_CEILING
            error = max(D(0), value-interval.lo, interval.hi-value)
        return self.upper(error)

    def record(self, interval):
        return {"value": self.value(interval), "absolute_error_bound": self.error(interval),
                "lower_decimal": str(interval.lo), "upper_decimal": str(interval.hi)}

    def statistics(self):
        return {"precision": self.precision, "moment_steps": self.steps, "max_moment_steps": self.max_steps,
                "point_evaluations": self.points, "max_point_evaluations": self.max_points,
                "max_series_terms": self.max_series, "cached_radii": len(self.cache),
                "max_cached_radii": self.max_cache, "max_z": str(self.max_z),
                "central_float_underflows_recorded": self.underflow_displays}


class GaussianRadialSource:
    """Exact PG source data plus a separately evaluated exterior companion.

    Constructor/projection work is charged cumulatively to max_products.
    Every nonzero Fischer mode is retained; an insufficient ell_max raises.
    Numerical moment/point caps are cumulative across derived providers.
    """
    def __init__(self, polynomial, beta=F(1), *, ell_max, precision=64,
                 max_products=5_000_000, max_terms=100_000,
                 max_moment_steps=2_000_000, max_series_terms=2048,
                 max_point_evaluations=100_000, max_cached_radii=256,
                 max_z=10_000):
        self.beta = _exact(beta, "beta")
        if self.beta <= 0:
            raise ValueError("beta must be positive")
        self.ell_max = _integer(ell_max, "ell_max", upper=16)
        self.max_products = _integer(max_products, "max_products", 1)
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.exact_work = 0
        coefficients, stats = fischer_project(polynomial, ell_max, max_products=max_products, max_terms=max_terms)
        self.exact_work += stats["work"]["products_and_scatters"]
        self.fischer_statistics = stats
        self.coefficients = coefficients
        widths = {}
        for ell, m, k in coefficients:
            widths[ell, m] = max(widths.get((ell, m), 0), k+1)
        self._charge_exact(len(coefficients)+6*sum(widths.values()))
        self.profile = GaussianProfile(coefficients, self.beta)
        self.mode_keys = self.profile.modes
        self._engine = _Engine(precision=precision, max_moment_steps=max_moment_steps,
                               max_series_terms=max_series_terms, max_point_evaluations=max_point_evaluations,
                               max_cached_radii=max_cached_radii, max_z=max_z)
        self._derivatives = {}
        for mode, row in self.profile.profiles.items():
            derivatives = [dict(row)]
            for _ in range(8):
                derivatives.append(self._differentiate_regular(derivatives[-1]))
            self._derivatives[mode] = tuple(MappingProxyType(row) for row in derivatives)
        self._derivatives = MappingProxyType(self._derivatives)
        self.maximum_moment = max((2*l+2*k+4 for l, m, k in coefficients), default=4)

    def _charge_exact(self, count):
        if self.exact_work+count > self.max_products:
            raise ValueError("PG exact construction exceeds max_products")
        self.exact_work += count

    def _differentiate_regular(self, polynomial):
        self._charge_exact(2*len(polynomial))
        result = {}
        for k, coefficient in polynomial.items():
            _put(result, k, -self.beta*coefficient)
            if k:
                _put(result, k-1, k*coefficient)
        return result

    def sample(self, radii, *, max_values=2_000_000):
        """Audited binary64 GaussianProfile samples; no FD error is implied.

        For enclosed point values use source_regular_jet/error methods on the
        exterior companion. This sampling method preserves the tested profile
        implementation and its explicit regular r^l factor.
        """
        return self.profile.sample(radii, max_values=max_values)

    def origin_data(self, ell, m, degree, patch_radius, derivative_order=4):
        _mode(ell, m)
        _integer(degree, "origin Taylor degree", upper=32)
        _integer(derivative_order, "source derivative_order", upper=8)
        patch_radius = _exact(patch_radius, "origin patch_radius")
        if patch_radius <= 0:
            raise ValueError("origin patch_radius must be positive")
        patch_float = float(patch_radius)
        if not isfinite(patch_float) or patch_float <= 0:
            raise ValueError("origin patch_radius must remain positive and finite in binary64")
        bound_radius = max(patch_radius, F.from_float(patch_float))
        row = self.profile.profiles.get((ell, m), {})
        # Reserve an explicit finite upper bound before the existing exact
        # profile constructor; higher-order derivative polynomials are cached.
        width = max(row, default=0)+1
        self._charge_exact((degree+1)*len(row)+12*width+(derivative_order+1)*width)
        base = self.profile.origin_data(ell, m, degree, patch_radius, order=min(3, derivative_order))
        bounds, s0 = [], bound_radius*bound_radius
        derivatives = self._derivatives.get((ell, m), ({},)*9)
        for j in range(derivative_order+1):
            n = degree-j
            bound = F(0)
            for k, coefficient in derivatives[j].items():
                if k <= n:
                    bound += abs(coefficient)*self.beta**(n-k+1)*s0**(n+1)/factorial(n-k+1)
                else:
                    bound += abs(coefficient)*s0**k
            bounds.append(bound)
        self._charge_exact((degree+1)*(derivative_order+1))
        conversion = []
        converted = []
        for coefficient in base["coefficients"]:
            value = float(coefficient)
            if not isfinite(value):
                raise ArithmeticError("Origin coefficient exceeds binary64 range")
            converted.append(abs(F.from_float(value)-coefficient))
        for j in range(derivative_order+1):
            conversion.append(sum((error*F(factorial(k), factorial(k-j))*s0**(k-j)
                                   for k, error in enumerate(converted) if k >= j), F(0)))
        combined = tuple(a+b for a, b in zip(bounds, conversion))
        return {"coefficients": base["coefficients"], "patch_radius": patch_radius,
                "patch_radius_float": patch_float, "bound_patch_radius_exact": bound_radius,
                "error_bounds_s": combined,
                "error_bounds_float_s": tuple(self._engine.upper(self._engine.number(value).hi) for value in combined),
                "source_truncation_error_bounds_s": tuple(bounds),
                "coefficient_conversion_error_bounds_s": tuple(conversion),
                "derivative_order": derivative_order,
                "coefficient_status": "exact source Fractions; no fitted a0",
                "float_conversion_error_included": True,
                "conversion_rule": "error_bounds_float_s rounds the combined exact bounds UP for OriginSeries",
                "bound_method": "Exact polynomial times exponential Taylor remainder, derivatives0..8",
                "global_fd_a0_uncertainty": "separate; not bounded here"}

    def asymptotic_moments(self):
        """Immutable (l,m,trace)->Fraction multiplying sqrt(pi/beta).

        These are exact ALL-SPACE SOURCE moments. Trace0/1 are both retained;
        they are held-out checks, never replacements for FD matching data.
        """
        result = {}
        for (ell, m, k), coefficient in self.coefficients.items():
            for trace in (0, 1):
                self._charge_exact(1)
                n = ell+k+trace
                odd_factorial = factorial(2*n+1)//(2**n*factorial(n))
                _put(result, (ell, m, trace), coefficient*F(odd_factorial, 2**(n+2))/self.beta**(n+1))
        return MappingProxyType(result)

    def finite_moments(self, radius):
        """Enclosed finite-core/tail radial moments, with exact source symbols."""
        radius = _radius(radius)
        self._engine.point()
        outer = self._engine.moments(self.beta, radius, self.maximum_moment)
        all_space = self._engine.moments(self.beta, F(0), self.maximum_moment)
        records = []
        for mode, polynomial in self.profile.profiles.items():
            ell, m = mode
            for trace in (0, 1):
                total = sum((c*all_space[2*ell+2*k+2*trace+2] for k, c in polynomial.items()), self._engine.number(0))
                tail = sum((c*outer[2*ell+2*k+2*trace+2] for k, c in polynomial.items()), self._engine.number(0))
                records.append({"mode": [ell, m], "trace": trace,
                                "core": self._engine.record(total-tail), "tail": self._engine.record(tail),
                                "all_space": self._engine.record(total)})
        return {"radius": str(radius), "moments": records,
                "meaning": "source-exact integrals with arithmetic enclosures; not matched numerical coefficients"}

    def moment_enclosures(self, radius, powers):
        """Public JSON records for U_n=int_radius^infinity r^n exp(-beta*r²)dr."""
        radius = _radius(radius)
        powers = tuple(powers)
        if not powers or len(powers) > 69 or len(set(powers)) != len(powers):
            raise ValueError("powers must contain1..69 distinct moment indices")
        for power in powers:
            _integer(power, "moment power", upper=68)
        self._engine.point()
        rows = self._engine.moments(self.beta, radius, max(powers))
        return {"beta": str(self.beta), "radius": str(radius),
                "moments": [{"power": n, **self._engine.record(rows[n])} for n in powers],
                "statistics": self.statistics(), "scope": "analytic/arithmetic moment enclosure only"}

    def statistics(self):
        return {"exact_products_counted": self.exact_work, "max_products": self.max_products,
                "max_terms": self.max_terms, "numerical_work": self._engine.statistics()}

    def exterior(self, radius):
        return GaussianExteriorTail(self, radius)

    def _F(self, mode, radius, order):
        _integer(order, "source regular derivative order", upper=8)
        if self.beta*radius*radius > self._engine.max_z:
            raise ValueError("Gaussian source argument exceeds max_z")
        si = self._engine.number(radius*radius)
        envelope = self._engine.number(-self.beta*radius*radius).exp()
        rows = self._derivatives.get(mode, ({},)*9)
        return tuple(envelope*sum((c*si**k for k, c in row.items()), self._engine.number(0)) for row in rows[:order+1])

    def _raw(self, mode, radius, absolute=False):
        ell, m = mode
        moments = self._engine.moments(self.beta, radius, self.maximum_moment)
        row = self.profile.profiles.get(mode, {})
        coefficients = {k: abs(c) if absolute else c for k, c in row.items()}
        zero = self._engine.number(0)
        return {
            "B": sum((c*moments[2*k+1] for k, c in coefficients.items()), zero),
            "B1": sum((c*moments[2*k+3] for k, c in coefficients.items()), zero),
            "T": sum((c*moments[2*ell+2*k+2] for k, c in coefficients.items()), zero),
            "T1": sum((c*moments[2*ell+2*k+4] for k, c in coefficients.items()), zero),
        }

    def provenance(self):
        return {"kind": "exact_polynomial_Gaussian_source", "beta": str(self.beta),
                "envelope": "exp(-beta*r^2)", "harmonic_convention": "unnormalized H_lm",
                "mode_keys": [list(mode) for mode in self.mode_keys], "ell_max": self.ell_max,
                "fischer": self.fischer_statistics, "exact_products_counted": self.exact_work,
                "max_products": self.max_products, "max_terms": self.max_terms,
                "all_space_moment_prefactor": "sqrt(pi/"+str(self.beta)+")",
                "exact_source_moments_are_not_fd_boundary_data": True,
                "numerical_work": self._engine.statistics(), "fd_error_bound": None,
                "total_error_certified": False}


def _regular_green_jet(engine, ell, radius, alpha, A, B, F0, F1, order):
    """A_s=s^(l+1/2)F/2; cancellation is performed symbolically first."""
    s = engine.number(radius*radius)
    lam = F(2*ell+1, 2)
    # s^(-lambda) = r^(-2l-1); every remaining exponent is integral.
    power = engine.number(radius)**(-2*ell-1)
    rows = [alpha*F(1, 2*ell+1)*(power*A+B)]
    if order >= 1:
        rows.append(-alpha*F(1, 2)*power*A/s)
    if order >= 2:
        rows.append(alpha*F(1, 2)*(lam+1)*power*A/s**2-alpha*F(1, 4)*F0/s)
    if order >= 3:
        rows.append(-alpha*F(1, 2)*(lam+1)*(lam+2)*power*A/s**3
                    +alpha*F(1, 4)*(lam+2)*F0/s**2-alpha*F(1, 4)*F1/s)
    return tuple(rows)


class _ResidualProvider:
    def _arguments(self, ell, m, radius, alpha, order):
        mode = _mode(ell, m)
        _integer(order, "response derivative order", upper=3)
        radius, alpha = _radius(radius), _alpha(alpha)
        self._engine.point()
        return mode, radius, alpha

    def _at_R(self, mode):
        if mode not in self._boundary_data:
            self._boundary_data[mode] = self._tail_data(mode, self._radius)
        return self._boundary_data[mode]

    def _response(self, mode, radius, alpha, order, side):
        if side not in (None, "left", "right"):
            raise ValueError("side must be None, left, or right")
        if radius == self._radius and order and side is None:
            raise ValueError("Response jets at R require an explicit side")
        if radius < self._radius or (radius == self._radius and side != "right"):
            value = alpha*self._at_R(mode)["B"]/(2*mode[0]+1)
            return (value,)+(self._engine.number(0),)*order
        data = self._tail_data(mode, radius)
        A = self._engine.number(0) if radius == self._radius else self._at_R(mode)["T"]-data["T"]
        return _regular_green_jet(self._engine, mode[0], radius, alpha, A, data["B"], data["F"][0], data["F"][1], order)

    def _remainder(self, mode, radius, alpha, order):
        if radius < self._radius:
            raise ValueError("The subtracted Gaussian remainder is exterior-only and is singular at the origin")
        data = self._tail_data(mode, radius)
        return _regular_green_jet(self._engine, mode[0], radius, alpha, -data["T"], data["B"], data["F"][0], data["F"][1], order)

    def outside_response_regular_jet(self, ell, m, r, alpha, order=3, side=None):
        mode, radius, alpha = self._arguments(ell, m, r, alpha, order)
        return tuple(self._engine.value(value) for value in self._response(mode, radius, alpha, order, side))

    def response_error_bounds(self, ell, m, r, alpha, order=3, side=None):
        mode, radius, alpha = self._arguments(ell, m, r, alpha, order)
        return tuple(self._engine.error(value) for value in self._response(mode, radius, alpha, order, side))

    def remainder_response_regular_jet(self, ell, m, r, alpha, order=3):
        mode, radius, alpha = self._arguments(ell, m, r, alpha, order)
        return tuple(self._engine.value(value) for value in self._remainder(mode, radius, alpha, order))

    def remainder_response_error_bounds(self, ell, m, r, alpha, order=3):
        mode, radius, alpha = self._arguments(ell, m, r, alpha, order)
        return tuple(self._engine.error(value) for value in self._remainder(mode, radius, alpha, order))

    def remainder_response_absolute_bounds(self, ell, m, r, alpha, order=3):
        """Bounds on the remainder itself, distinct from arithmetic errors."""
        mode, radius, alpha = self._arguments(ell, m, r, alpha, order)
        return tuple(self._engine.upper(max(value.lo.copy_abs(), value.hi.copy_abs()))
                     for value in self._remainder(mode, radius, alpha, order))

    def response_enclosures(self, ell, m, r, alpha, order=3, side=None):
        mode, radius, alpha = self._arguments(ell, m, r, alpha, order)
        rows = self._response(mode, radius, alpha, order, side)
        return {"mode": [ell, m], "radius": str(radius), "alpha": str(alpha), "side": side,
                "regular_s_jets": [self._engine.record(row) for row in rows],
                "meaning": "outside-source correction; FD and matching errors excluded"}

    def origin_shift(self, ell, m, alpha):
        return self.outside_response_regular_jet(ell, m, 0, alpha, order=0)[0]

    def origin_shift_error_bound(self, ell, m, alpha):
        return self.response_error_bounds(ell, m, 0, alpha, order=0)[0]

    def algebraic_response_terms(self, ell, m, alpha):
        mode, _, alpha = self._arguments(ell, m, self._radius, alpha, 0)
        value = alpha*self._at_R(mode)["T"]/(2*ell+1)
        if value.lo == value.hi == 0:
            return {}
        return {(ell+1, 0): self._engine.value(value)}

    def algebraic_response_error_bounds(self, ell, m, alpha):
        mode, _, alpha = self._arguments(ell, m, self._radius, alpha, 0)
        value = alpha*self._at_R(mode)["T"]/(2*ell+1)
        return {(ell+1, 0): self._engine.error(value)}

    def source_regular_jet(self, ell, m, r, order=1):
        mode = _mode(ell, m)
        _integer(order, "exterior source derivative order", upper=2)
        radius = _radius(r)
        self._engine.point()
        if radius < self._radius:
            raise ValueError("Residual source is defined only for r>=R")
        return tuple(self._engine.value(value) for value in self._tail_data(mode, radius)["F"][:order+1])

    def source_error_bounds(self, ell, m, r, order=1):
        mode = _mode(ell, m)
        _integer(order, "exterior source derivative order", upper=2)
        radius = _radius(r)
        self._engine.point()
        if radius < self._radius:
            raise ValueError("Residual source is defined only for r>=R")
        return tuple(self._engine.error(value) for value in self._tail_data(mode, radius)["F"][:order+1])

    def _provenance(self, kind):
        return {"kind": kind, "radius": self.radius, "radius_exact_binary_or_rational": str(self._radius),
                "mode_keys": [list(mode) for mode in self.mode_keys],
                "policy": "explicit full outside-source correction, plus separate algebraic-response contribution",
                "central_exterior_source_omitted": False,
                "bounds_cover": "Gaussian moment enclosures and conversion of returned numbers only",
                "fd_a0_interpolation_projection_errors_included": False,
                "matched_homogeneous_coefficients_supplied": False,
                "arithmetic_uncertainty": "response/source/algebraic error-bound methods",
                "source_origin_remainder": "separate source OriginSeries declaration",
                "fd_uncertainty": "unknown here; preserve independent refinement and matching diagnostics",
                "numerical_work": self._engine.statistics(), "total_error_certified": False}


class GaussianExteriorTail(_ResidualProvider):
    """Full noncompact PG exterior source, with an explicit scalar response."""
    def __init__(self, source, radius):
        if not isinstance(source, GaussianRadialSource):
            raise TypeError("source must be a GaussianRadialSource")
        self.source, self._engine = source, source._engine
        self._radius = _radius(radius)
        self.radius_exact = self._radius
        self.radius = float(self._radius)
        if not isfinite(self.radius) or self.radius <= 0:
            raise ValueError("Exterior radius must be positive and representable")
        if source.beta*self._radius*self._radius > self._engine.max_z:
            raise ValueError("Exterior radius exceeds max_z")
        self.mode_keys = source.mode_keys
        self._boundary_data = {}

    def _tail_data(self, mode, radius):
        result = self.source._raw(mode, radius)
        result["F"] = self.source._F(mode, radius, 2)
        return result

    def weighted_outer_bounds(self):
        """Absolute source-tail moments int_R∞r^(1-l)|f_l|dr, not zero."""
        self._engine.point()
        return {mode: self._engine.upper(self.source._raw(mode, self._radius, True)["B"].hi)
                for mode in self.mode_keys}

    def response_remainder_derivative(self, axis, alpha=1):
        """Derivative of the Gaussian RESPONSE remainder as a new SOURCE.

        The lower branch is absent for ell0. Every angular coefficient is
        projected exactly, so neither magnetic modes nor trace terms are fit.
        """
        _integer(axis, "axis", upper=2)
        alpha = _alpha(alpha)
        channels = {}
        for ell, m in self.mode_keys:
            harmonic = solid_harmonic(ell, m)
            orders = tuple(int(i == axis) for i in range(3))
            lower = dict(harmonic.derivative(orders))
            self.source._charge_exact(len(harmonic.poly)+3*len(lower))
            upper = {}
            for powers, coefficient in harmonic.poly.items():
                shifted = list(powers)
                shifted[axis] += 1
                _put(upper, tuple(shifted), coefficient)
            for powers, coefficient in lower.items():
                for direction in range(3):
                    shifted = list(powers)
                    shifted[direction] += 2
                    _put(upper, tuple(shifted), -coefficient/F(2*ell+1))
            for branch, polynomial, output_ell in (("lower", lower, ell-1), ("upper", upper, ell+1)):
                if not polynomial:
                    continue
                if output_ell > 16:
                    raise ValueError("Derivative needs ell>16; no angular mode may be dropped")
                cap = self.source.max_products-self.source.exact_work
                if cap < 1:
                    raise ValueError("Residual angular projection exhausted max_products")
                try:
                    projected, stats = fischer_project(polynomial, output_ell, max_products=cap,
                                                       max_terms=self.source.max_terms)
                except Exception:
                    self.source._charge_exact(cap)
                    raise
                self.source._charge_exact(stats["work"]["products_and_scatters"])
                for (l_out, m_out, trace), coefficient in projected.items():
                    if trace or l_out != output_ell:
                        raise ArithmeticError("Solid-harmonic derivative channel is not harmonic")
                    _put(channels, ((l_out, m_out), (ell, m), branch), alpha*coefficient)
                    if len(channels) > self.source.max_terms:
                        raise ValueError("Derivative residual channels exceed max_terms")
        return GaussianDivergenceTail(self, channels, axis, alpha)

    def provenance(self):
        return {**self._provenance("PG_Gaussian_exterior_source"), "source": self.source.provenance(),
                "outside_response": "A=T(R)-T(r), B=B(r); inside shift alpha*B(R)/(2l+1)",
                "algebraic_response": "additional alpha*T(R)/(2l+1)*r^(-l-1), included exactly once"}


class GaussianDivergenceTail(_ResidualProvider):
    """One Cartesian derivative of the scalar PG response's Gaussian part.

    Exact Fubini moment identities avoid inverse quadrature. With d=2l+1,
    B=int_r∞rho^(1-l)f, T=int_r∞rho^(l+2)f, B1=int rho^(3-l)f,
    T1=int rho^(l+4)f, the lower/upper SOURCE moments are:
      lower: Bh=(B1-r²B)/(2d), Mh=(T-r^d B)/d²;
      upper: Bh=(r^-d T-B)/d, Mh=(T1-r²T)/2.
    Angular derivative coefficients, including the parent's alpha, are in
    channels. These formulas apply at positive exterior radius ONLY.
    """
    def __init__(self, parent, channels, axis, alpha):
        self.parent, self.source, self._engine = parent, parent.source, parent._engine
        self.radius, self._radius = parent.radius, parent._radius
        self.radius_exact = self._radius
        self.axis, self.parent_alpha = axis, alpha
        self.channels = MappingProxyType(dict(channels))
        self.mode_keys = tuple(sorted({key[0] for key in channels}))
        self._boundary_data = {}
        grouped = {}
        for (output, original, branch), coefficient in channels.items():
            grouped.setdefault(output, []).append((original, branch, coefficient))
        self._grouped = MappingProxyType({mode: tuple(rows) for mode, rows in grouped.items()})

    def _channel_moments(self, original, branch, radius, absolute=False):
        ell, m = original
        raw = self.source._raw(original, radius, absolute)
        d = 2*ell+1
        if branch == "lower":
            B = (raw["B1"]-radius*radius*raw["B"])/(2*d)
            T = (raw["T"]-radius**d*raw["B"])/(d*d)
        else:
            B = (radius**(-d)*raw["T"]-raw["B"])/d
            T = (raw["T1"]-radius*radius*raw["T"])/2
        if absolute:
            # For nonnegative source coefficients the Fubini kernels are
            # positive. Clamp only an enclosure's lower endpoint, not data.
            B, T = B.nonnegative(), T.nonnegative()
        return B, T, raw

    def _tail_data(self, mode, radius):
        zero = self._engine.number(0)
        B, T, derivatives = zero, zero, [zero, zero, zero]
        for original, branch, coefficient in self._grouped.get(mode, ()):
            ell, m = original
            b, t, raw = self._channel_moments(original, branch, radius)
            source_F = self.source._F(original, radius, 1)
            B, T = B+coefficient*b, T+coefficient*t
            d, lam = 2*ell+1, F(2*ell+1, 2)
            si = self._engine.number(radius*radius)
            if branch == "lower":
                # F_lower=B/d. Cancel the two T terms algebraically before
                # interval arithmetic; the following derivatives are exact.
                source_jets = (raw["B"]/d, -source_F[0]/(2*d), -source_F[1]/(2*d))
            else:
                power = self._engine.number(radius)**(-2*ell-3)
                source_jets = (power*raw["T"],
                               -(lam+1)*power*raw["T"]/si-source_F[0]/(2*si),
                               (lam+1)*(lam+2)*power*raw["T"]/si**2
                               +(lam+2)*source_F[0]/(2*si**2)-source_F[1]/(2*si))
            for j in range(3):
                derivatives[j] = derivatives[j]+coefficient*source_jets[j]
        return {"B": B, "T": T, "F": tuple(derivatives)}

    def weighted_outer_bounds(self):
        self._engine.point()
        result = {}
        for mode, channels in self._grouped.items():
            bound = self._engine.number(0)
            for original, branch, coefficient in channels:
                moment, _, _ = self._channel_moments(original, branch, self._radius, True)
                bound = bound+abs(coefficient)*moment
            result[mode] = self._engine.upper(bound.hi)
        return result

    def response_remainder_derivative(self, axis, alpha=1):
        raise ValueError("Deeper Gaussian residual differentiation is unsupported; cannot claim closure")

    def provenance(self):
        return {**self._provenance("derivative_of_PG_response_Gaussian_remainder"),
                "cartesian_axis": self.axis, "parent_inverse_alpha": str(self.parent_alpha),
                "exact_angular_channels": len(self.channels),
                "source_domain": "r>=R>0 only; not a regular-origin source polynomial",
                "moments": "Closed lower/upper Fubini identities including trace1; no inverse quadrature",
                "independent_of_numerical_fd_homogeneous_coefficients": True,
                "localized_h00_may_be_nonzero": True,
                "deeper_derivative_supported": False, "source": self.source.provenance()}
