"""Independent radial scalar inverse and Cartesian York-jet combination.

This module uses only the standard library. Import constructs no fields and
performs no solve. It does not import the production HEAT, Fourier, pairing,
tail, or geometry implementations. Angular projection belongs to the caller.

CONVENTIONS / PUBLIC API
  RadialGrid(radius, intervals, breakpoints=())
      .radii is the immutable tuple (0,h,...,R). All joins are explicit knots.
  OriginSeries(radius, coefficients, error_bounds)
      f_l(r)=r**l F(r*r); coefficients[k] multiplies s**k in F(s).
      error_bounds[j] declares a uniform bound on the j-th s derivative of
      F minus that polynomial, on the origin patch. At least j=0,1,2 are
      required. These are caller declarations, NOT certified by this module.
  ExteriorSource(radius, modes, remainder_outer_bounds)
      modes[(l,m)][(p,q)] multiplies r**(-p)*log(r/R)**q in f_l(r).
      Each remainder bound is integral_R^infinity r**(1-l)*|epsilon_l| dr.
      Explicit zero bounds describe exact tails. Allowed source powers obey
      p>=2, p+l>2; p=2 permits only l>=1 and q=0 (bounded inverse domain).
  solve_scalar_modes(grid, source_modes, *, inverse='C', exterior,
                     origin_series, source_limits=None, progress=None)
      source_modes[(l,m)] contains len(grid.radii) real shell samples.
      inverse='P' means (-Delta)^-1; 'C' means (-Delta)^-1/8.
      source_limits[(l,m)][join_radius]=(left_f,right_f) optionally supplies
      distinct traces. A jump row uses their arithmetic average; interpolation
      and core quadrature retain the separate traces. Interpolation stencils
      stay within a smooth segment; the FD transmission row couples both sides.
  ScalarModes.radial_jet(l,m,r,order=3,side=None)
      -> (u, d_r u, d_r^2 u, d_r^3 u), truncated to order+1.
  ScalarModes.regular_jet(l,m,r,order=3,side=None)
      -> (G, G_s, G_ss, G_sss), G(s)=u_l(r)/r**l, s=r*r.
  combine_york_jets(vector_jets, chi_jet)
      V=P f, chi=P div(V); returns W=-V-grad(chi)/4 and dW,ddW.

All supplied magnetic modes are kept. Shell coefficients multiply the real
UNNORMALIZED solid harmonics H_lm(n); any H/Y conversion is explicit upstream.
The radial equation itself does not depend on harmonic normalization. The
origin bridge consumes source Taylor data, never a fit to f_l/r**l. Its single
homogeneous coefficient a0 is matched to the global FD value at the declared
patch edge. Its discrepancies from other FD values/derivatives are reported.

The scalar solve is second-order radial FD with second-order inhomogeneous
Robin data. Away from the origin patch, u,u' use local cubic interpolation;
u'',u''' use the ODE and a differentiated source interpolant. Thus a residual
formed with those same ODE derivatives is NOT an independent PDE test.
Interior jets at a source join or R require side='left'/'right'. At r=0 the
regular series is used directly, including for Cartesian reconstruction.
Formal second-order accuracy is not a uniform error theorem for large l or
shrinking origin patches. Refine with a fixed physical origin patch and check
its Taylor error and a0 matching separately; a patch of only a fixed number
of cells can retain an l-dependent bias in a0 as h decreases. Conversion to
regular derivatives outside the patch can also amplify cancellation.

Analytic power/log continuation matches numerical u(R), not an archived mass
target. An independently integrated finite-core moment discrepancy is exposed.
Unknown exterior remainders are omitted from the central value and their input
bounds are reported separately. No total solution/jet/energy error certificate
is inferred from those inputs or from the diagnostic rows.
"""

from bisect import bisect_right
from dataclasses import dataclass, field
from math import factorial, fsum, isfinite, log
from sys import float_info
from types import MappingProxyType


def _number(value, name):
    if isinstance(value, (bool, str, bytes)):
        raise TypeError(name + " must be a real number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(name + " must be a finite real number") from exc
    if not isfinite(result):
        raise ValueError(name + " must be finite")
    return result


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(name + " must be an integer >= " + str(minimum))
    return value


def _mode(key):
    if not isinstance(key, tuple) or len(key) != 2:
        raise ValueError("mode keys must be (l,m) tuples")
    ell = _integer(key[0], "l")
    m = key[1]
    if isinstance(m, bool) or not isinstance(m, int) or abs(m) > ell:
        raise ValueError("m must be an integer with |m|<=l")
    return ell, m


def _order(order):
    _integer(order, "order")
    if order > 3:
        raise ValueError("jets are implemented only through order 3")
    return order


def _finite_tuple(values, name):
    return tuple(_number(value, name) for value in values)


def _finite_tree(value, name):
    if isinstance(value, dict):
        for item in value.values():
            _finite_tree(item, name)
    elif isinstance(value, (tuple, list)):
        for item in value:
            _finite_tree(item, name)
    elif isinstance(value, float) and not isfinite(value):
        raise ArithmeticError(name + " contains nonfinite arithmetic")


@dataclass(frozen=True)
class RadialGrid:
    radius: float
    intervals: int
    breakpoints: tuple = ()
    h: float = field(init=False)
    radii: tuple = field(init=False)
    breakpoint_indices: tuple = field(init=False)

    def __post_init__(self):
        radius = _number(self.radius, "radius")
        intervals = _integer(self.intervals, "intervals", 3)
        if radius <= 0:
            raise ValueError("radius must be positive")
        h = radius / intervals
        if h <= 0 or not isfinite(h):
            raise ValueError("grid spacing must be finite and positive")
        indices = []
        for point in self.breakpoints:
            point = _number(point, "breakpoint")
            index = round(point / h)
            if not 0 < index < intervals or abs(point - index*h) > 32*float_info.epsilon*radius:
                raise ValueError("each interior breakpoint must be a grid knot")
            indices.append(index)
        if len(set(indices)) != len(indices):
            raise ValueError("duplicate breakpoint")
        indices.sort()
        endpoints = (0, *indices, intervals)
        if any(b-a < 3 for a, b in zip(endpoints, endpoints[1:])):
            raise ValueError("each smooth grid segment needs at least three intervals")
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "radii", tuple(i*h for i in range(intervals)) + (radius,))
        object.__setattr__(self, "breakpoint_indices", tuple(indices))
        object.__setattr__(self, "breakpoints", tuple(i*h for i in indices))


@dataclass(frozen=True)
class OriginSeries:
    radius: float
    coefficients: tuple
    error_bounds: tuple

    def __post_init__(self):
        radius = _number(self.radius, "origin radius")
        coefficients = _finite_tuple(self.coefficients, "origin coefficient")
        bounds = _finite_tuple(self.error_bounds, "origin error bound")
        if radius <= 0 or not coefficients:
            raise ValueError("origin radius must be positive and coefficients nonempty")
        if len(bounds) < 3 or any(value < 0 for value in bounds):
            raise ValueError("nonnegative origin error bounds for F,F_s,F_ss are required")
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "error_bounds", bounds)


@dataclass(frozen=True)
class ExteriorSource:
    radius: float
    modes: object
    remainder_outer_bounds: object

    def __post_init__(self):
        radius = _number(self.radius, "exterior radius")
        if radius <= 0:
            raise ValueError("exterior radius must be positive")
        modes = {}
        for key, terms in self.modes.items():
            ell, m = _mode(key)
            cleaned = {}
            for pair, coefficient in terms.items():
                if not isinstance(pair, tuple) or len(pair) != 2:
                    raise ValueError("exterior keys must be (p,log_power)")
                p = _integer(pair[0], "exterior p", 2)
                q = _integer(pair[1], "exterior log_power")
                coefficient = _number(coefficient, "exterior coefficient")
                if not coefficient:
                    continue
                if p+ell <= 2 or (p == 2 and (ell < 1 or q != 0)):
                    raise ValueError("source tail is outside the bounded/decaying inverse domain")
                cleaned[p, q] = coefficient
            modes[ell, m] = MappingProxyType(cleaned)
        bounds = {}
        for key, value in self.remainder_outer_bounds.items():
            key = _mode(key)
            value = _number(value, "exterior remainder bound")
            if value < 0:
                raise ValueError("exterior remainder bounds must be nonnegative")
            bounds[key] = value
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "modes", MappingProxyType(modes))
        object.__setattr__(self, "remainder_outer_bounds", MappingProxyType(bounds))

    def outer_moment(self, ell, m):
        """Algebraic integral_R^infinity r^(1-l) f_l(r) dr, without 1/8."""
        key = _mode((ell, m))
        terms = self.modes.get(key, {})
        return _number(fsum(c*factorial(q)*self.radius**(2-ell-p) /
                            (p+ell-2)**(q+1) for (p, q), c in terms.items()),
                       "exterior outer moment")


def _polynomial_jet(coefficients, s, order):
    result = []
    current = tuple(coefficients)
    for _ in range(order+1):
        value = 0.0
        for coefficient in reversed(current):
            value = value*s + coefficient
        result.append(value)
        current = tuple(i*current[i] for i in range(1, len(current)))
    return tuple(result)


def _regular_radial_jet(coefficients, ell, r, order):
    result = []
    for derivative in range(order+1):
        terms = []
        for k, coefficient in enumerate(coefficients):
            power = ell+2*k
            if power < derivative or not coefficient:
                continue
            factor = 1
            for j in range(derivative):
                factor *= power-j
            terms.append(coefficient*factor*r**(power-derivative))
        result.append(fsum(terms))
    return tuple(result)


def _power_log_jet(terms, r, radius, order):
    """Here keys (p,q) multiply r^-p log(r/R)^q; p may be negative."""
    current = dict(terms)
    logarithm = log(r/radius)
    result = []
    for _ in range(order+1):
        result.append(fsum(c*r**(-p)*logarithm**q for (p, q), c in current.items()))
        following = {}
        for (p, q), coefficient in current.items():
            if p:
                pair = p+1, q
                following[pair] = following.get(pair, 0.0)-p*coefficient
            if q:
                pair = p+1, q-1
                following[pair] = following.get(pair, 0.0)+q*coefficient
        current = following
    return tuple(result)


def _particular_terms(ell, terms, alpha):
    """Solve the power/log radial ODE; resonances raise the log degree."""
    result = {}
    for (p, q), coefficient in terms.items():
        exponent = 2-p
        constant = exponent*(exponent+1)-ell*(ell+1)
        first = 2*exponent+1
        polynomial = [0.0]*(q+3)
        if constant:
            for j in range(q, -1, -1):
                rhs = -alpha*coefficient if j == q else 0.0
                polynomial[j] = (rhs-first*(j+1)*polynomial[j+1] -
                                 (j+2)*(j+1)*polynomial[j+2])/constant
        else:
            # P(0)=0 fixes only the particular solution's homogeneous freedom.
            for j in range(q, -1, -1):
                rhs = -alpha*coefficient if j == q else 0.0
                polynomial[j+1] = (rhs-(j+2)*(j+1)*polynomial[j+2])/(first*(j+1))
        for j, value in enumerate(polynomial):
            if value:
                pair = p-2, j
                result[pair] = result.get(pair, 0.0)+value
    return {pair: value for pair, value in result.items() if value}


def _factorize(grid, ell):
    n, h, radius = grid.intervals, grid.h, grid.radius
    lower, diagonal, upper = ([0.0]*(n+1) for _ in range(3))
    if ell == 0:
        diagonal[0], upper[0] = -6/h**2, 6/h**2
    else:
        diagonal[0] = 1.0
    for i in range(1, n):
        r = grid.radii[i]
        lower[i] = 1/h**2-1/(r*h)
        diagonal[i] = -2/h**2-ell*(ell+1)/r**2
        upper[i] = 1/h**2+1/(r*h)
    elimination = (1/(2*h))/lower[-2]
    lower[-1] = -2/h-elimination*diagonal[-2]
    diagonal[-1] = 3/(2*h)+(ell+1)/radius-elimination*upper[-2]
    pivots, modified = [], []
    for i in range(n+1):
        pivot = diagonal[i]-(lower[i]*modified[i-1] if i else 0.0)
        if not isfinite(pivot) or pivot == 0:
            raise ArithmeticError("nonfinite or singular radial factorization for l=" + str(ell))
        pivots.append(1/pivot)
        modified.append(upper[i]/pivot)
    return tuple(lower), tuple(pivots), tuple(modified), elimination


def _solve_fd(grid, ell, samples, limits, alpha, outer, factors):
    forcing = list(samples)
    for index, (left, right) in limits.items():
        forcing[index] = (left+right)/2
    rhs = [-alpha*value for value in forcing]
    if ell:
        rhs[0] = 0.0
    lower, pivots, modified, elimination = factors
    rhs[-1] = alpha*grid.radius**(ell-1)*outer-elimination*rhs[-2]
    answer = []
    for i, value in enumerate(rhs):
        answer.append((value-(lower[i]*answer[i-1] if i else 0.0))*pivots[i])
    for i in range(grid.intervals-1, -1, -1):
        answer[i] -= modified[i]*answer[i+1]
    return _finite_tuple(answer, "radial FD solution")


def _segment(grid, r, side, require_side=False):
    if side not in (None, "left", "right"):
        raise ValueError("side must be None, 'left', or 'right'")
    boundaries = (0, *grid.breakpoint_indices, grid.intervals)
    radii = tuple(grid.radii[i] for i in boundaries)
    if r in radii[1:] and side is None and require_side:
        raise ValueError("an explicit side is required at a join or R")
    segment = bisect_right(radii, r)-1
    if r in radii[1:] and side != "right":
        segment -= 1
    segment = min(max(segment, 0), len(boundaries)-2)
    return boundaries[segment], boundaries[segment+1]


def _cubic_jet(values, t, h):
    value_terms, derivative_terms = [], []
    for j in range(4):
        other = [k for k in range(4) if k != j]
        denominator = 1
        numerator = 1.0
        for k in other:
            denominator *= j-k
            numerator *= t-k
        derivative = 0.0
        for omitted in other:
            product = 1.0
            for k in other:
                if k != omitted:
                    product *= t-k
            derivative += product
        value_terms.append(values[j]*numerator/denominator)
        derivative_terms.append(values[j]*derivative/(denominator*h))
    return fsum(value_terms), fsum(derivative_terms)


def _sample_jet(grid, samples, limits, r, side, is_source=False):
    lo, hi = _segment(grid, r, side)
    start = min(max(int(r/grid.h)-1, lo), hi-3)
    values = []
    for index in range(start, start+4):
        value = samples[index]
        if is_source and index == lo and index in limits:
            value = limits[index][1]
        if is_source and index == hi and index in limits:
            value = limits[index][0]
        values.append(value)
    return _cubic_jet(values, r/grid.h-start, grid.h)


def _interior_jet(grid, ell, values, samples, limits, alpha, r, side):
    value, first = _sample_jet(grid, values, limits, r, side)
    source, dsource = _sample_jet(grid, samples, limits, r, side, True)
    second = ell*(ell+1)*value/r**2-2*first/r-alpha*source
    third = ((ell*(ell+1)+2)*first/r**2-2*ell*(ell+1)*value/r**3 -
             2*second/r-alpha*dsource)
    return value, first, second, third


def _core_moment(grid, ell, samples, limits):
    """Separate composite Simpson 1/3 + 3/8 integration on each smooth piece."""
    boundaries = (0, *grid.breakpoint_indices, grid.intervals)
    pieces = []
    for lo, hi in zip(boundaries, boundaries[1:]):
        values = []
        for i in range(lo, hi+1):
            value = samples[i]
            if i == lo and i in limits:
                value = limits[i][1]
            if i == hi and i in limits:
                value = limits[i][0]
            values.append(value*grid.radii[i]**(ell+2))
        count = hi-lo
        end = count if count % 2 == 0 else count-3
        for i in range(0, end, 2):
            pieces.append(grid.h*(values[i]+4*values[i+1]+values[i+2])/3)
        if count % 2:
            pieces.append(3*grid.h*(values[-4]+3*values[-3]+3*values[-2]+values[-1])/8)
    return fsum(pieces)


@dataclass(frozen=True)
class _ModeSolution:
    samples: tuple
    values: tuple
    limits: object
    origin: OriginSeries
    origin_solution: tuple
    exterior_solution: object
    diagnostics: object


@dataclass(frozen=True)
class ScalarModes:
    grid: RadialGrid
    inverse: str
    exterior: ExteriorSource
    _solutions: object

    @property
    def mode_keys(self):
        return tuple(self._solutions)

    @property
    def values(self):
        """Raw FD values; the declared origin patch uses the matched series."""
        return MappingProxyType({key: value.values for key, value in self._solutions.items()})

    def _arguments(self, ell, m, r, order, side):
        key, order = _mode((ell, m)), _order(order)
        r = _number(r, "evaluation radius")
        if r < 0:
            raise ValueError("evaluation radius must be nonnegative")
        if side not in (None, "left", "right"):
            raise ValueError("side must be None, 'left', or 'right'")
        if key not in self._solutions:
            raise KeyError("mode was not supplied: " + str(key))
        if r <= self.grid.radius:
            _segment(self.grid, r, side, require_side=order > 0)
        return key, r, order, self._solutions[key]

    def radial_jet(self, ell, m, r, order=3, side=None):
        key, r, order, mode = self._arguments(ell, m, r, order, side)
        if r > self.grid.radius or (r == self.grid.radius and side == "right"):
            result = _power_log_jet(mode.exterior_solution, r, self.grid.radius, order)
        elif r < mode.origin.radius or (r == mode.origin.radius and side != "right"):
            result = _regular_radial_jet(mode.origin_solution, ell, r, order)
        else:
            alpha = 0.125 if self.inverse == "C" else 1.0
            result = _interior_jet(self.grid, ell, mode.values, mode.samples, mode.limits,
                                   alpha, r, side)[:order+1]
        return _finite_tuple(result, "radial jet " + str(key))

    def regular_jet(self, ell, m, r, order=3, side=None):
        key, r, order, mode = self._arguments(ell, m, r, order, side)
        if r < mode.origin.radius or (r == mode.origin.radius and side != "right"):
            return _finite_tuple(_polynomial_jet(mode.origin_solution, r*r, order),
                                 "regular origin jet " + str(key))
        jet = self.radial_jet(ell, m, r, order, side)
        result = [jet[0]/r**ell]
        if order >= 1:
            result.append((jet[1]-ell*jet[0]/r)/(2*r**(ell+1)))
        if order >= 2:
            result.append((jet[2]-(2*ell+1)*jet[1]/r +
                           ell*(ell+2)*jet[0]/r**2)/(4*r**(ell+2)))
        if order >= 3:
            result.append((jet[3]-3*(ell+1)*jet[2]/r +
                           3*(ell*ell+3*ell+1)*jet[1]/r**2 -
                           ell*(ell+2)*(ell+4)*jet[0]/r**3)/(8*r**(ell+3)))
        return _finite_tuple(result, "regular jet " + str(key))

    def source_jet(self, ell, m, r, order=1, side=None):
        key, r, order, mode = self._arguments(ell, m, r, order, side)
        if order > 1:
            raise ValueError("source jets are exposed only through first radial derivative")
        if r in self.grid.breakpoints and side is None:
            raise ValueError("source traces at a join require an explicit side")
        if r == self.grid.radius and side is None:
            raise ValueError("source trace at R requires an explicit side")
        if r > self.grid.radius or (r == self.grid.radius and side == "right"):
            result = _power_log_jet(self.exterior.modes.get(key, {}), r, self.grid.radius, order)
        elif r < mode.origin.radius or (r == mode.origin.radius and side != "right"):
            result = _regular_radial_jet(mode.origin.coefficients, ell, r, order)
        else:
            result = _sample_jet(self.grid, mode.samples, mode.limits, r, side, True)[:order+1]
        return _finite_tuple(result, "source jet " + str(key))

    def exterior_terms(self, ell, m):
        """Read-only response coefficients (p,q): c*r^-p*log(r/R)^q."""
        return self._solutions[_mode((ell, m))].exterior_solution

    def diagnostics(self):
        """JSON-serializable provenance and errors; not an acceptance decision."""
        rows = [dict(mode=key, **dict(value.diagnostics))
                for key, value in self._solutions.items()]
        return {
            "inverse": self.inverse, "harmonic_convention": "unnormalized_real_H_lm",
            "radius": self.grid.radius, "intervals": self.grid.intervals,
            "breakpoints": self.grid.breakpoints, "mode_count": len(rows), "rows": rows,
            "work": {"source_samples": len(rows)*(self.grid.intervals+1),
                     "fd_unknowns": len(rows)*(self.grid.intervals+1),
                     "factorizations": len({key[0] for key in self._solutions})},
            "origin_error_bounds_are_caller_declarations": True,
            "ode_derived_derivatives_are_not_independent_pde_checks": True,
            "total_error_certified": False,
        }


def solve_scalar_modes(grid, source_modes, *, inverse="C", exterior,
                       origin_series, source_limits=None, progress=None):
    """Solve all supplied (l,m) modes, with no implicit angular truncation."""
    if not isinstance(grid, RadialGrid) or not isinstance(exterior, ExteriorSource):
        raise TypeError("grid and exterior must be RadialGrid and ExteriorSource")
    if inverse not in ("C", "P"):
        raise ValueError("inverse must be 'C' or 'P'")
    if exterior.radius != grid.radius:
        raise ValueError("exterior logarithm anchor/radius must equal grid.radius")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable or None")
    keys = tuple(sorted(_mode(key) for key in source_modes))
    keyset = set(keys)
    if set(origin_series) != keyset:
        raise ValueError("origin_series must declare exactly every source mode")
    if set(exterior.remainder_outer_bounds) != keyset or not set(exterior.modes) <= keyset:
        raise ValueError("exterior bounds must declare every source mode; no extra tail modes")
    source_limits = {} if source_limits is None else source_limits
    if not set(source_limits) <= keyset:
        raise ValueError("source_limits contains an unsupplied mode")
    alpha = 0.125 if inverse == "C" else 1.0
    factors, solutions = {}, {}
    for completed, key in enumerate(keys, 1):
        ell, m = key
        samples = _finite_tuple(source_modes[key], "source sample " + str(key))
        if len(samples) != grid.intervals+1:
            raise ValueError("source sample count must equal len(grid.radii)")
        origin = origin_series[key]
        if not isinstance(origin, OriginSeries):
            raise TypeError("each origin declaration must be OriginSeries")
        patch = round(origin.radius/grid.h)
        first_join = grid.breakpoint_indices[0] if grid.breakpoint_indices else grid.intervals
        if not 2 <= patch <= first_join or abs(origin.radius-grid.radii[patch]) > 32*float_info.epsilon*grid.radius:
            raise ValueError("origin patch edge must be a grid knot >=2h and no later than first join")
        limits = {}
        for radius, traces in source_limits.get(key, {}).items():
            radius = _number(radius, "source join")
            if radius not in grid.breakpoints:
                raise ValueError("source_limits keys must be declared interior breakpoints")
            pair = _finite_tuple(traces, "source trace")
            if len(pair) != 2:
                raise ValueError("each source join needs (left,right) traces")
            limits[grid.breakpoint_indices[grid.breakpoints.index(radius)]] = pair
        if ell not in factors:
            factors[ell] = _factorize(grid, ell)
        outer = exterior.outer_moment(ell, m)
        values = _solve_fd(grid, ell, samples, limits, alpha, outer, factors[ell])
        origin_coefficients = [0.0] + [
            -alpha*d/(2*(k+1)*(2*ell+2*k+3)) for k, d in enumerate(origin.coefficients)]
        patch_radius = grid.radii[patch]
        particular_patch = _regular_radial_jet(origin_coefficients, ell, patch_radius, 0)[0]
        origin_coefficients[0] = (values[patch]-particular_patch)/patch_radius**ell
        origin_coefficients = _finite_tuple(origin_coefficients, "origin solution coefficient")
        exterior_terms = _particular_terms(ell, exterior.modes.get(key, {}), alpha)
        particular_boundary = _power_log_jet(exterior_terms, grid.radius, grid.radius, 0)[0]
        homogeneous = (values[-1]-particular_boundary)*grid.radius**(ell+1)
        pair = ell+1, 0
        exterior_terms[pair] = exterior_terms.get(pair, 0.0)+homogeneous
        exterior_terms = {key: _number(value, "exterior solution coefficient")
                          for key, value in exterior_terms.items() if value}
        core = _core_moment(grid, ell, samples, limits)
        boundary_core = ((2*ell+1)*grid.radius**(ell+1)*values[-1]/alpha -
                         grid.radius**(2*ell+1)*outer)
        source_mismatches, source_excesses, patch_differences, a0_spread = [], [], [], []
        for i in range(patch+1):
            radius = grid.radii[i]
            source = limits[i][0] if i in limits else samples[i]
            source_series = _regular_radial_jet(origin.coefficients, ell, radius, 0)[0]
            mismatch = abs(source-source_series)
            source_mismatches.append(mismatch)
            rounding = 100*float_info.epsilon*max(1.0, abs(source), abs(source_series))
            source_excesses.append(max(0.0, mismatch-radius**ell*origin.error_bounds[0]-rounding))
            predicted = _regular_radial_jet(origin_coefficients, ell, radius, 0)[0]
            patch_differences.append(abs(values[i]-predicted))
            if i >= max(1, (patch+1)//2):
                a0_spread.append(abs(values[i]-predicted)/radius**ell)
        left_patch = _interior_jet(grid, ell, values, samples, limits, alpha, patch_radius, "left")
        series_patch = _regular_radial_jet(origin_coefficients, ell, patch_radius, 3)
        derivative_matching = tuple(series_patch[i]-left_patch[i] for i in range(1, 4))
        fd_boundary_derivative = (3*values[-1]-4*values[-2]+values[-3])/(2*grid.h)
        exterior_boundary = _power_log_jet(exterior_terms, grid.radius, grid.radius, 1)
        cubic_boundary = _sample_jet(grid, values, limits, grid.radius, "left")[1]
        robin_rhs = alpha*grid.radius**ell*outer
        joins = tuple((grid.radii[i],
                       _sample_jet(grid, values, limits, grid.radii[i], "right")[1] -
                       _sample_jet(grid, values, limits, grid.radii[i], "left")[1])
                      for i in grid.breakpoint_indices)
        bound = exterior.remainder_outer_bounds[key]
        diagnostics = {
            "origin_radius": patch_radius, "origin_source_error_bounds_s": origin.error_bounds,
            "origin_source_sample_max_error": max(source_mismatches),
            "origin_source_sample_max_bound_excess": max(source_excesses),
            "origin_source_samples_within_declared_bound": max(source_excesses) == 0.0,
            "origin_derivative_error_bounds_independently_checked": False,
            "origin_a0": origin_coefficients[0], "origin_a0_sample_spread": max(a0_spread),
            "origin_patch_max_fd_difference": max(patch_differences),
            "origin_patch_derivative_matching_differences": derivative_matching,
            "origin_patch_value_matching_is_constructed": True,
            "core_moment_quadrature": core, "core_moment_from_boundary": boundary_core,
            "core_moment_difference": boundary_core-core,
            "exterior_outer_moment": outer, "exterior_remainder_outer_bound": bound,
            "exterior_omitted_value_bound_at_R": alpha*grid.radius**ell*bound/(2*ell+1),
            "robin_rhs": robin_rhs,
            "robin_fd_residual": grid.radius*fd_boundary_derivative+(ell+1)*values[-1]-robin_rhs,
            "robin_cubic_residual": grid.radius*cubic_boundary+(ell+1)*values[-1]-robin_rhs,
            "exterior_derivative_matching_difference": exterior_boundary[1]-cubic_boundary,
            "exterior_value_matching_is_constructed": True,
            "join_derivative_jumps": joins,
            "source_jump_rows_averaged": len(limits), "exterior_response_terms": len(exterior_terms),
        }
        _finite_tree(diagnostics, "scalar diagnostics " + str(key))
        solutions[key] = _ModeSolution(samples, values, MappingProxyType(limits), origin,
                                        origin_coefficients, MappingProxyType(exterior_terms),
                                        MappingProxyType(diagnostics))
        if progress is not None:
            progress({"stage": "scalar_mode", "inverse": inverse, "mode": key,
                      "completed": completed, "total": len(keys), "intervals": grid.intervals})
    return ScalarModes(grid, inverse, exterior, MappingProxyType(solutions))


def _jet_array(value, shape, name):
    if not shape:
        return _number(value, name)
    if not isinstance(value, (tuple, list)) or len(value) != shape[0]:
        raise ValueError(name + " must have shape " + str(shape))
    return tuple(_jet_array(item, shape[1:], name) for item in value)


def combine_york_jets(vector_jets, chi_jet):
    """Combine independently supplied Cartesian potential jets, without solves.

    vector_jets contains three dictionaries, each with value, gradient[a], and
    hessian[a][b]. chi_jet additionally supplies third[a][b][i]. Output indices
    are W[i], dW[a][i], ddW[a][b][i], matching york_operator_fields. No assumption
    that chi decays is made; a localized div(V) monopole is never deleted.
    """
    if not isinstance(vector_jets, (tuple, list)) or len(vector_jets) != 3:
        raise ValueError("vector_jets must contain three scalar Cartesian jets")
    values = tuple(_number(jet["value"], "V value") for jet in vector_jets)
    gradients = tuple(_jet_array(jet["gradient"], (3,), "V gradient") for jet in vector_jets)
    hessians = tuple(_jet_array(jet["hessian"], (3, 3), "V Hessian") for jet in vector_jets)
    gradient = _jet_array(chi_jet["gradient"], (3,), "chi gradient")
    hessian = _jet_array(chi_jet["hessian"], (3, 3), "chi Hessian")
    third = _jet_array(chi_jet["third"], (3, 3, 3), "chi third derivative")
    result = {
        "W": tuple(-values[i]-gradient[i]/4 for i in range(3)),
        "dW": tuple(tuple(-gradients[i][a]-hessian[a][i]/4 for i in range(3)) for a in range(3)),
        "ddW": tuple(tuple(tuple(-hessians[i][a][b]-third[a][b][i]/4 for i in range(3))
                           for b in range(3)) for a in range(3)),
        "div_V": fsum(gradients[i][i] for i in range(3)),
    }
    _finite_tree(result, "York combination")
    return result
