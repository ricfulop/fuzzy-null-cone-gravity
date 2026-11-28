"""General independent scalar/York field adapters for the unreduced H6 chain.

No scientific objects are constructed on import. No fixture targets or
production HEAT/Fourier/pairing implementation is imported. Public operators
are C=(-Delta)^-1/8, P=(-Delta)^-1 and Y=Delta_L^-1. Shell coefficients multiply
unnormalized real H_lm. Spatial derivatives have derivative axes first.

PUBLIC CONTRACT
  ScalarSourceSpec(samples, origin_series, exterior, *, residual_policy,
                   source_limits=None, residual_provider=None, provenance=None)
      samples[(l,m)] is an immutable copy of the radial sample sequence.
      origin_series[(l,m)] is an explicit OriginSeries, never a fitted source.
      residual_policy is 'zero', 'exact_correction', or 'bounded_omission'.
      Bounds in ExteriorSource never modify a central Robin datum.
  WorkBudget(max_fd_unknowns=1000000, max_point_evaluations=5000000,
             max_exact_products=20000000, max_stored_values=5000000,
             max_cartesian_calls=20000000, max_projection_products=50000000)
      One shared instance must cover the complete caller's computation.
  solve_scalar_source(grid, spec, *, inverse, budget, progress=None)
  solve_york_source(grid, vector_specs, angular, *, h_lmax, guard_lmax,
                    support_tolerance, budget, progress=None)
      vector_specs contains three ScalarSourceSpecs. The angular grid covers
      the complete guard band; no fixture-specific parity/m restriction exists.
      Origin patches must be declared grid joins, common to the York inputs.

ScalarField and YorkField provide jet(point,order=2,side=None), shell_jets,
origin_series_data(), exterior_data(), exterior_remainder_jet(), diagnostics().
Jets use value/gradient/hessian (and scalar third if requested), matching the
pointwise unreduced algebra. Taylor data and Laurent objects are mathematical
objects; diagnostics/checkpoint payloads alone are JSON-serializable.

EXACT-CORRECTION PROVIDER (duck typed; no circular import)
  mode_keys, radius, radius_exact (Fraction)
  outside_response_regular_jet(l,m,r,alpha,order=3,side=None)
  origin_shift(l,m,alpha)
  algebraic_response_terms(l,m,alpha) -> {(p,q): coefficient}
  response_error_bounds(l,m,r,alpha,order=3,side=None)
  remainder_response_regular_jet(l,m,r,alpha,order=3), r>=R
  response_remainder_derivative(axis,alpha=1) -> a SOURCE residual provider
  source_regular_jet(l,m,r,order=1), weighted_outer_bounds(), provenance()
The derivative provider differentiates only the response's nonalgebraic
remainder. Differentiating its full outside-source response would count the
additional algebraic moment twice. SumProviders composes this protocol.

All source Taylor coefficients plus computed a0 are retained. The computed h
origin polynomial is an identity for the nominal, polynomial origin bridge;
its inherited source-Taylor, correction-arithmetic and FD/a0 uncertainties
remain separate. Neither numerical a0 nor finite-difference error is certified
by a single solve. Bounded omissions retain their supplied errors and are not
silently upgraded to exact corrections or total field/energy bounds.
Inside that origin patch, h samples use the computed polynomial after a
separate comparison with projected div(V); projection noise is never divided
by r^l or fitted into a source series. Right-hand patch traces remain explicit.

Provider internals have their own declared work caps; this budget counts
provider calls and adapter operations, not unreported provider internal cost.
Laurent derivative scatter reservations are conservative upper work counts.
Every restriction is preceded by a finite/support/domain gate and checkpoint.
"""

from collections.abc import Mapping
from fractions import Fraction
import json
from math import fsum, inf, isfinite, nextafter, sqrt
from types import MappingProxyType

from spherical_green_reference import ExteriorSource, OriginSeries, RadialGrid, solve_scalar_modes
from solid_harmonic_reference import AngularGrid, cartesian_jet, fischer_project, solid_harmonic
from kelvin_exterior_reference import CartesianLaurent, project_exterior


def _real(value, name, nonnegative=False):
    if isinstance(value, (bool, str, bytes)):
        raise TypeError(name + " must be real")
    result = float(value)
    if not isfinite(result) or (nonnegative and result < 0):
        raise ValueError(name + " must be finite" + (" and nonnegative" if nonnegative else ""))
    return result


def _integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name + " must be an integer >= " + str(minimum))
    return value


def _upper_fraction(value):
    """Convert a nonnegative exact rational upward to finite binary64."""
    if value < 0:
        raise ValueError("a bound must be nonnegative")
    result = _real(float(value), "converted bound", True)
    if Fraction.from_float(result) < value:
        result = nextafter(result, inf)
    return _real(result, "upward-rounded bound", True)


def _upper_sum(values):
    return _upper_fraction(sum((Fraction.from_float(_real(value, "bound", True))
                                for value in values), Fraction(0)))


def _mode(value):
    if not isinstance(value, tuple) or len(value) != 2:
        raise ValueError("mode must be an (l,m) tuple")
    ell = _integer(value[0], "ell")
    if type(value[1]) is not int or abs(value[1]) > ell:
        raise ValueError("invalid magnetic index")
    return value


def _jet_order(order, maximum=3):
    _integer(order, "jet order")
    if order > maximum:
        raise ValueError("requested jet order exceeds " + str(maximum))
    return order


def _point(point):
    point = tuple(_real(value, "coordinate") for value in point)
    if len(point) != 3:
        raise ValueError("point needs three coordinates")
    return point


def _json_copy(value):
    return json.loads(json.dumps(value, allow_nan=False))


def _freeze_json(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze_json(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze_json(item) for item in value)
    return value


def _plain(value):
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    return value


def _emit(progress, stage, **payload):
    if progress is not None:
        progress({"stage": stage, **_json_copy(payload)})


def _exact_anchor(provider):
    anchor = getattr(provider, "radius_exact", None)
    if isinstance(anchor, bool) or not isinstance(anchor, (int, Fraction)):
        raise TypeError("provider must expose its exact rational radius_exact")
    anchor = Fraction(anchor)
    if anchor <= 0:
        raise ValueError("provider exact radius must be positive")
    return anchor


def _provider_validate(provider, radius, radius_exact=None):
    names = ("outside_response_regular_jet", "origin_shift", "algebraic_response_terms",
             "response_error_bounds", "remainder_response_regular_jet",
             "response_remainder_derivative", "source_regular_jet",
             "weighted_outer_bounds", "provenance")
    if _real(provider.radius, "provider radius") != radius:
        raise ValueError("provider radius does not match the source exterior anchor")
    expected_anchor = Fraction.from_float(radius) if radius_exact is None else radius_exact
    if _exact_anchor(provider) != expected_anchor:
        raise ValueError("provider exact radius differs from the requested FD/exterior anchor")
    keys = tuple(_mode(key) for key in provider.mode_keys)
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate provider mode")
    if any(not callable(getattr(provider, name, None)) for name in names):
        raise TypeError("residual provider does not implement the complete protocol")
    return keys


class WorkBudget:
    """Cumulative bounded counters; no reset is performed by a solve or query."""
    def __init__(self, max_fd_unknowns=1_000_000, max_point_evaluations=5_000_000,
                 max_exact_products=20_000_000, max_stored_values=5_000_000,
                 max_cartesian_calls=20_000_000, max_projection_products=50_000_000):
        self.limits = {"fd_unknowns": max_fd_unknowns, "point_evaluations": max_point_evaluations,
                       "exact_products": max_exact_products, "stored_values": max_stored_values,
                       "cartesian_calls": max_cartesian_calls,
                       "projection_products": max_projection_products}
        self.limits = MappingProxyType({key: _integer(value, key, 1) for key, value in self.limits.items()})
        self.counts = {key: 0 for key in self.limits}
        self.provider_calls = 0

    def remaining(self, key):
        return self.limits[key]-self.counts[key]

    def charge(self, key, count=1):
        _integer(count, "work count")
        if count > self.remaining(key):
            raise RuntimeError("cumulative spherical work cap exceeded: " + key)
        self.counts[key] += count

    def snapshot(self):
        return {"counts": dict(self.counts), "limits": dict(self.limits),
                "provider_calls": self.provider_calls,
                "exact_derivative_work_includes_upper_reservations": True,
                "provider_internal_work_not_inferred": True}


class ScalarSourceSpec:
    """Immutable copies of source data; provider objects must themselves be immutable."""
    __slots__ = ("samples", "origin_series", "exterior", "source_limits",
                 "residual_policy", "residual_provider", "provider_outer_bounds", "provenance", "_locked")

    def __setattr__(self, name, value):
        if getattr(self, "_locked", False):
            raise AttributeError("ScalarSourceSpec is immutable")
        object.__setattr__(self, name, value)

    def __init__(self, samples, origin_series, exterior, *, residual_policy,
                 source_limits=None, residual_provider=None, provenance=None):
        if not isinstance(exterior, ExteriorSource):
            raise TypeError("exterior must be ExteriorSource")
        if residual_policy not in ("zero", "exact_correction", "bounded_omission"):
            raise ValueError("an explicit residual policy is required")
        copied = {_mode(mode): tuple(_real(x, "source sample") for x in values)
                  for mode, values in samples.items()}
        keys = set(copied)
        if set(origin_series) != keys or set(exterior.remainder_outer_bounds) != keys:
            raise ValueError("source, origin and exterior-bound mode sets must agree")
        if not set(exterior.modes) <= keys:
            raise ValueError("algebraic source exterior has an undeclared mode")
        if any(not isinstance(value, OriginSeries) for value in origin_series.values()):
            raise TypeError("every origin declaration must be OriginSeries")
        provider_bounds = {}
        if residual_policy == "exact_correction":
            if residual_provider is None:
                raise ValueError("exact_correction requires its provider")
            provider_keys = set(_provider_validate(residual_provider, exterior.radius))
            if not provider_keys <= keys:
                raise ValueError("residual provider has an undeclared source mode")
            raw_bounds = residual_provider.weighted_outer_bounds()
            if set(raw_bounds) != provider_keys:
                raise ValueError("provider weighted bounds must declare exactly its mode_keys")
            provider_bounds = {mode: _real(value, "provider weighted outer bound", True)
                               for mode, value in raw_bounds.items()}
            uncovered = [mode for mode, bound in exterior.remainder_outer_bounds.items()
                         if bound and mode not in provider_keys]
            if uncovered:
                raise ValueError("exact_correction leaves nonzero source remainders uncovered: " + str(uncovered))
            if any(exterior.remainder_outer_bounds[mode] < bound for mode, bound in provider_bounds.items()):
                raise ValueError("declared source omission bounds must cover the provider's stated outer bounds")
        elif residual_provider is not None:
            raise ValueError("only exact_correction may attach a central correction provider")
        if residual_policy == "zero" and any(exterior.remainder_outer_bounds.values()):
            raise ValueError("zero residual policy cannot discard nonzero residual bounds")
        limits = {}
        for mode, entries in (source_limits or {}).items():
            if mode not in keys:
                raise ValueError("source_limits has an undeclared mode")
            row = {}
            for radius, traces in entries.items():
                pair = tuple(_real(x, "source trace") for x in traces)
                if len(pair) != 2:
                    raise ValueError("source_limits entries require left/right traces")
                row[_real(radius, "join radius")] = pair
            limits[mode] = MappingProxyType(row)
        self.samples = MappingProxyType(copied)
        self.origin_series = MappingProxyType(dict(origin_series))
        self.exterior = exterior
        self.source_limits = MappingProxyType(limits)
        self.residual_policy = residual_policy
        self.residual_provider = residual_provider
        self.provider_outer_bounds = MappingProxyType(provider_bounds)
        self.provenance = _freeze_json(_json_copy(provenance or {}))
        self._locked = True


class SumProviders:
    """Linear sum of independently supplied exterior residual providers."""
    def __init__(self, providers, radius=None):
        self.providers = tuple(providers)
        if radius is None:
            if not self.providers:
                raise ValueError("empty provider sum requires an explicit radius")
            radius = _exact_anchor(self.providers[0])
        self.radius = _real(radius, "provider-sum radius")
        self.radius_exact = Fraction(radius) if isinstance(radius, (int, Fraction)) else Fraction.from_float(self.radius)
        self.mode_keys = tuple(sorted({key for provider in self.providers
                                      for key in _provider_validate(provider, self.radius, self.radius_exact)}))

    def _sum_jet(self, method, ell, m, radius, order, **kwargs):
        _integer(order, "provider jet order")
        terms = [[] for _ in range(order+1)]
        for provider in self.providers:
            if (ell, m) not in provider.mode_keys:
                continue
            jet = tuple(getattr(provider, method)(ell, m, radius, order=order, **kwargs))
            if len(jet) != order+1:
                raise ValueError("provider returned the wrong jet length")
            for i, value in enumerate(jet):
                terms[i].append(_real(value, "provider jet", method == "response_error_bounds"))
        if method == "response_error_bounds":
            return tuple(_upper_sum(values) for values in terms)
        return tuple(_real(fsum(values), "summed provider jet") for values in terms)

    def outside_response_regular_jet(self, ell, m, radius, alpha, order=3, side=None):
        _jet_order(order)
        return self._sum_jet("outside_response_regular_jet", ell, m, radius, order, alpha=alpha, side=side)

    def response_error_bounds(self, ell, m, radius, alpha, order=3, side=None):
        _jet_order(order)
        errors, centers = ([[] for _ in range(order+1)] for _ in range(2))
        for provider in self.providers:
            if (ell, m) not in provider.mode_keys:
                continue
            bounds = tuple(provider.response_error_bounds(ell, m, radius, alpha, order=order, side=side))
            central = tuple(provider.outside_response_regular_jet(ell, m, radius, alpha, order=order, side=side))
            if len(bounds) != order+1 or len(central) != order+1:
                raise ValueError("provider returned the wrong error/central jet length")
            for i, (bound, value) in enumerate(zip(bounds, central)):
                errors[i].append(Fraction.from_float(_real(bound, "component error bound", True)))
                centers[i].append(_real(value, "component central jet"))
        result = []
        for bounds, central in zip(errors, centers):
            exact_center_sum = sum((Fraction.from_float(value) for value in central), Fraction(0))
            returned_center = _real(fsum(central), "summed central jet")
            rounding = abs(Fraction.from_float(returned_center)-exact_center_sum)
            result.append(_upper_fraction(sum(bounds, Fraction(0))+rounding))
        return tuple(result)

    def remainder_response_regular_jet(self, ell, m, radius, alpha, order=3):
        _jet_order(order)
        return self._sum_jet("remainder_response_regular_jet", ell, m, radius, order, alpha=alpha)

    def source_regular_jet(self, ell, m, radius, order=1):
        return self._sum_jet("source_regular_jet", ell, m, radius, order)

    def origin_shift(self, ell, m, alpha):
        return _real(fsum(_real(provider.origin_shift(ell, m, alpha), "origin shift")
                    for provider in self.providers if (ell, m) in provider.mode_keys), "summed origin shift")

    def algebraic_response_terms(self, ell, m, alpha):
        entries = {}
        for provider in self.providers:
            if (ell, m) in provider.mode_keys:
                for pair, value in provider.algebraic_response_terms(ell, m, alpha).items():
                    entries.setdefault(pair, []).append(_real(value, "algebraic correction"))
        result = {pair: _real(fsum(values), "summed algebraic correction") for pair, values in entries.items()}
        return {pair: value for pair, value in result.items() if value}

    def response_remainder_derivative(self, axis, alpha=1):
        return SumProviders(tuple(provider.response_remainder_derivative(axis, alpha=alpha)
                                  for provider in self.providers), self.radius_exact)

    def weighted_outer_bounds(self):
        entries = {mode: [] for mode in self.mode_keys}
        for provider in self.providers:
            bounds = provider.weighted_outer_bounds()
            if set(bounds) != set(provider.mode_keys):
                raise ValueError("provider weighted bounds must declare every mode")
            for mode, value in bounds.items():
                entries[mode].append(_real(value, "provider outer bound", True))
        return {mode: _upper_sum(values) for mode, values in entries.items()}

    def provenance(self):
        return {"kind": "sum_of_exterior_residuals", "radius": self.radius,
                "radius_exact": str(self.radius_exact),
                "response_bounds_include_rounding_of_the_returned_sum": True,
                "providers": [_json_copy(provider.provenance()) for provider in self.providers]}


def _accumulate(out, key, value, budget):
    budget.charge("exact_products")
    value = out.get(key, Fraction(0))+value
    if value:
        out[key] = value
    else:
        out.pop(key, None)


def _fischer(poly, ell_max, budget):
    remaining = budget.remaining("exact_products")
    if remaining <= 0:
        raise RuntimeError("no exact projection work remains")
    try:
        result, stats = fischer_project(poly, ell_max, max_products=remaining)
    except BaseException as exc:
        used = getattr(exc, "statistics", {}).get("work", {}).get("products_and_scatters", remaining)
        budget.charge("exact_products", min(used, remaining))
        raise
    budget.charge("exact_products", stats["work"]["products_and_scatters"])
    return result


def _laurent_derivative(field, axis, budget):
    reservation = 3*len(field.terms)
    budget.charge("exact_products", reservation)
    return field.derivative(axis, max_products=max(1, reservation))


def _laurent_project(field, ell_max, budget):
    remaining = budget.remaining("exact_products")
    if remaining <= 0:
        raise RuntimeError("no exterior projection work remains")
    try:
        result, anchor, stats = project_exterior(field, ell_max, max_products=remaining)
    except BaseException as exc:
        used = getattr(exc, "statistics", {}).get("work", {}).get("products_and_scatters", remaining)
        budget.charge("exact_products", min(used, remaining))
        raise
    budget.charge("exact_products", stats["work"]["products_and_scatters"])
    return result, anchor, stats


def _laurent_from_modes(modes, radius, budget):
    terms = {}
    for (ell, m), row in modes.items():
        harmonic = solid_harmonic(ell, m).polynomial
        for pair, coefficient in row.items():
            if not isinstance(pair, tuple) or len(pair) != 2 or type(pair[0]) is not int:
                raise ValueError("response powers must be integer (p,q) pairs")
            p, q = pair
            _integer(q, "response logarithm power")
            coefficient = Fraction.from_float(_real(coefficient, "response coefficient"))
            for powers, value in harmonic.items():
                _accumulate(terms, (powers, -p-ell, q), coefficient*value, budget)
    return CartesianLaurent(terms, Fraction.from_float(radius))


def _sum_laurent(fields, factors, radius, budget):
    terms = {}
    anchor = Fraction.from_float(radius)
    for field, factor in zip(fields, factors):
        if field.anchor != anchor:
            raise ValueError("all Laurent logarithm anchors must agree")
        for key, value in field.terms.items():
            _accumulate(terms, key, factor*value, budget)
    return CartesianLaurent(terms, anchor)


def _origin_derivative(coefficients, axis, ell_max, budget):
    """Project small angular factors, retaining the entire radial series."""
    out = {}
    for (ell, m), series in coefficients.items():
        harmonic = solid_harmonic(ell, m)
        derivative_order = tuple(1 if i == axis else 0 for i in range(3))
        derivative = dict(harmonic.derivative(derivative_order))
        plus = {}
        for powers, value in harmonic.polynomial.items():
            shifted = list(powers)
            shifted[axis] += 1
            _accumulate(plus, tuple(shifted), value, budget)
        for powers, value in derivative.items():
            for i in range(3):
                shifted = list(powers)
                shifted[i] += 2
                _accumulate(plus, tuple(shifted), -value/Fraction(2*ell+1), budget)
        lower_map = _fischer(derivative, ell_max, budget)
        upper_map = _fischer(plus, ell_max, budget)
        for mapping, upper in ((lower_map, False), (upper_map, True)):
            for (degree, magnetic, trace), angular in mapping.items():
                if trace:
                    raise ArithmeticError("solid-harmonic derivative decomposition is not harmonic")
                bucket = out.setdefault((degree, magnetic), {})
                for k, coefficient in enumerate(series):
                    if upper and k == 0:
                        continue
                    factor = 2*k if upper else Fraction(2*ell+1+2*k, 2*ell+1)
                    radial = k-1 if upper else k
                    _accumulate(bucket, radial, angular*factor*Fraction.from_float(coefficient), budget)
    return {mode: tuple(float(row.get(k, 0)) for k in range(max(row, default=0)+1))
            for mode, row in out.items() if row}


def _sum_origins(origins, factors, budget):
    rows = {}
    for origin, factor in zip(origins, factors):
        for mode, coefficients in origin.items():
            bucket = rows.setdefault(mode, {})
            for k, value in enumerate(coefficients):
                _accumulate(bucket, k, factor*Fraction.from_float(value), budget)
    return {mode: tuple(float(row.get(k, 0)) for k in range(max(row, default=0)+1))
            for mode, row in rows.items() if row}


class ScalarField:
    def __init__(self, grid, spec, inverse, base, budget):
        self.grid, self.spec, self.inverse, self.base, self.budget = grid, spec, inverse, base, budget
        self.alpha = .125 if inverse == "C" else 1.0
        self.provider = spec.residual_provider
        self.mode_keys = tuple(sorted(spec.samples))
        diagnostics = {tuple(row["mode"]): row for row in base.diagnostics()["rows"]}
        origins, exteriors = {}, {}
        for mode in self.mode_keys:
            ell, m = mode
            coefficients = [diagnostics[mode]["origin_a0"]] + [
                -self.alpha*d/(2*(k+1)*(2*ell+2*k+3))
                for k, d in enumerate(spec.origin_series[mode].coefficients)]
            exterior = dict(base.exterior_terms(ell, m))
            if self.provider is not None and mode in self.provider.mode_keys:
                budget.provider_calls += 2
                coefficients[0] += _real(self.provider.origin_shift(ell, m, self.alpha), "origin correction")
                for pair, value in self.provider.algebraic_response_terms(ell, m, self.alpha).items():
                    exterior[pair] = exterior.get(pair, 0.0)+_real(value, "exterior correction")
            origins[mode] = tuple(_real(x, "origin response coefficient") for x in coefficients)
            exteriors[mode] = MappingProxyType(exterior)
        budget.charge("stored_values", sum(map(len, origins.values()))+sum(map(len, exteriors.values())))
        self.origin_coefficients = MappingProxyType(origins)
        self.exterior_modes = MappingProxyType(exteriors)
        self.algebraic = _laurent_from_modes(exteriors, grid.radius, budget)

    def regular_jet(self, ell, m, radius, order=3, side=None):
        result = self.base.regular_jet(ell, m, radius, order, side)
        if self.provider is not None and (ell, m) in self.provider.mode_keys:
            self.budget.provider_calls += 1
            correction = tuple(self.provider.outside_response_regular_jet(
                ell, m, radius, self.alpha, order=order, side=side))
            if len(correction) != order+1:
                raise ValueError("provider returned the wrong correction jet length")
            result = tuple(fsum((a, _real(b, "correction jet"))) for a, b in zip(result, correction))
        return tuple(_real(x, "corrected regular jet") for x in result)

    def correction_error_bounds(self, ell, m, radius, order=3, side=None):
        _jet_order(order)
        if self.provider is None or (ell, m) not in self.provider.mode_keys:
            return (0.0,)*(order+1)
        self.budget.provider_calls += 1
        values = tuple(self.provider.response_error_bounds(ell, m, radius, self.alpha, order=order, side=side))
        if len(values) != order+1:
            raise ValueError("provider returned the wrong error-bound length")
        return tuple(_real(value, "correction error bound", True) for value in values)

    def _regular(self, radius, order, side, remainder=False):
        if remainder:
            if radius < self.grid.radius:
                raise ValueError("exterior remainder is defined only at r>=R")
            result = {}
            if self.provider is not None:
                for mode in self.provider.mode_keys:
                    self.budget.provider_calls += 1
                    jet = tuple(self.provider.remainder_response_regular_jet(*mode, radius, self.alpha, order=order))
                    if len(jet) != order+1:
                        raise ValueError("provider returned the wrong remainder jet length")
                    result[mode] = tuple(_real(value, "remainder jet") for value in jet)
            return result
        return {mode: self.regular_jet(*mode, radius, order, side) for mode in self.mode_keys}

    def _cartesian(self, regular, point, order):
        self.budget.charge("cartesian_calls")
        return cartesian_jet(regular, point, order=order)

    def jet(self, point, order=2, side=None):
        _jet_order(order)
        point = _point(point)
        radius = sqrt(fsum(x*x for x in point))
        self.budget.charge("point_evaluations")
        return self._cartesian(self._regular(radius, order, side), point, order)

    def shell_jets(self, radius, angular, order=2, side=None):
        _jet_order(order)
        radius = _real(radius, "shell radius", True)
        regular = self._regular(radius, order, side)
        values = []
        for direction in angular.directions:
            self.budget.charge("point_evaluations")
            values.append(self._cartesian(regular, tuple(radius*x for x in direction), order))
        return tuple(values)

    def exterior_remainder_jet(self, point, order=2):
        _jet_order(order)
        point = _point(point)
        radius = sqrt(fsum(x*x for x in point))
        self.budget.charge("point_evaluations")
        return self._cartesian(self._regular(radius, order, "right", True), point, order)

    def origin_series_data(self):
        return {"coefficients": self.origin_coefficients,
                "patch_radii": MappingProxyType({mode: origin.radius for mode, origin in self.spec.origin_series.items()}),
                "nominal_polynomial_bridge": True, "a0_error_certified": False,
                "input_source_error_bounds": MappingProxyType({mode: origin.error_bounds
                    for mode, origin in self.spec.origin_series.items()}),
                "source_induced_regular_derivative_bounds": MappingProxyType({mode: tuple(
                    _upper_fraction(Fraction.from_float(self.alpha)*Fraction.from_float(bound)/(4*mode[0]+4*j+2))
                    for j, bound in enumerate(origin.error_bounds, 1))
                    for mode, origin in self.spec.origin_series.items()})}

    def exterior_data(self):
        return {"algebraic": self.algebraic, "modes": self.exterior_modes,
                "residual_policy": self.spec.residual_policy, "provider": self.provider,
                "omitted_outer_bounds": self.spec.exterior.remainder_outer_bounds}

    def diagnostics(self):
        provider = None
        if self.provider is not None:
            self.budget.provider_calls += 1
            provider = _json_copy(self.provider.provenance())
        return {"inverse": self.inverse, "base": self.base.diagnostics(),
                "residual_policy": self.spec.residual_policy, "provider": provider,
                "source_provenance": _plain(self.spec.provenance),
                "origin_series": [{"mode": mode, "coefficients": coefficients,
                                   "source_error_bounds": self.spec.origin_series[mode].error_bounds}
                                  for mode, coefficients in self.origin_coefficients.items()],
                "exterior": [{"mode": mode, "terms": [{"p": p, "q": q, "coefficient": c}
                               for (p, q), c in row.items()]} for mode, row in self.exterior_modes.items()],
                "source_residual_outer_bounds": [{"mode": mode, "bound": value}
                    for mode, value in self.spec.exterior.remainder_outer_bounds.items()],
                "bounds_do_not_modify_base_robin": True,
                "base_omission_bounds_superseded_by_explicit_correction": self.provider is not None,
                "exact_correction_modes": list(self.spec.provider_outer_bounds),
                "fd_and_a0_uncertainty": "requires separate measured refinement",
                "correction_arithmetic_bounds": "query correction_error_bounds separately",
                "total_error_certified": False, "work": self.budget.snapshot()}


def solve_scalar_source(grid, spec, *, inverse, budget, progress=None):
    if not isinstance(grid, RadialGrid) or not isinstance(spec, ScalarSourceSpec) or not isinstance(budget, WorkBudget):
        raise TypeError("expected RadialGrid, ScalarSourceSpec, and shared WorkBudget")
    if inverse not in ("C", "P"):
        raise ValueError("inverse must be C or P")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable")
    for origin in spec.origin_series.values():
        if origin.radius not in grid.breakpoints:
            raise ValueError("the origin patch edge must be a declared interior grid join")
    count = len(spec.samples)*(grid.intervals+1)
    budget.charge("fd_unknowns", count)
    budget.charge("stored_values", 2*count)
    base = solve_scalar_modes(grid, spec.samples, inverse=inverse, exterior=spec.exterior,
                              origin_series=spec.origin_series, source_limits=spec.source_limits,
                              progress=progress)
    failures = [row for row in base.diagnostics()["rows"]
                if not row["origin_source_samples_within_declared_bound"]]
    _emit(progress, "scalar_source_origin_gate", passed=not failures,
          failures=failures, work=budget.snapshot())
    _emit(progress, "scalar_source_base_result", inverse=inverse,
          residual_policy=spec.residual_policy, base=base.diagnostics(), work=budget.snapshot())
    if failures:
        raise ValueError("source samples violate their declared origin-series bounds")
    result = ScalarField(grid, spec, inverse, base, budget)
    _emit(progress, "scalar_source_result", result=result.diagnostics())
    return result


class YorkField:
    def __init__(self, vfields, chi, origin, exterior, source_record, budget):
        self.vfields, self.chi, self.budget = tuple(vfields), chi, budget
        self.V = self.vfields
        self.h_source = chi.spec
        self.grid = chi.grid
        self.origin_coefficients = tuple(MappingProxyType(row) for row in origin)
        self.algebraic = tuple(exterior)
        self.source_record = source_record

    def _regular(self, radius, order, side, remainder=False):
        return (tuple(field._regular(radius, order, side, remainder) for field in self.vfields),
                self.chi._regular(radius, order+1, side, remainder))

    def _combine(self, regular, point, order):
        vregular, chiregular = regular
        self.budget.charge("cartesian_calls", 4)
        v = tuple(cartesian_jet(item, point, order=order) for item in vregular)
        chi = cartesian_jet(chiregular, point, order=order+1)
        return {"value": tuple(_real(-v[i]["value"]-chi["gradient"][i]/4, "York value") for i in range(3)),
                "gradient": tuple(tuple(_real(-v[i]["gradient"][a]-chi["hessian"][a][i]/4, "York gradient")
                    for i in range(3)) for a in range(3)) if order >= 1 else None,
                "hessian": tuple(tuple(tuple(_real(-v[i]["hessian"][a][b]-chi["third"][a][b][i]/4, "York Hessian")
                    for i in range(3)) for b in range(3)) for a in range(3)) if order >= 2 else None}

    def jet(self, point, order=2, side=None):
        _jet_order(order, 2)
        point = _point(point)
        radius = sqrt(fsum(x*x for x in point))
        self.budget.charge("point_evaluations")
        return self._combine(self._regular(radius, order, side), point, order)

    def shell_jets(self, radius, angular, order=2, side=None):
        _jet_order(order, 2)
        radius = _real(radius, "shell radius", True)
        regular = self._regular(radius, order, side)
        result = []
        for direction in angular.directions:
            self.budget.charge("point_evaluations")
            result.append(self._combine(regular, tuple(radius*x for x in direction), order))
        return tuple(result)

    def divergence_jet(self, point, order=0, side=None):
        """Computed h=div(V), without substituting the second solve's source."""
        _jet_order(order, 2)
        point = _point(point)
        radius = sqrt(fsum(x*x for x in point))
        self.budget.charge("point_evaluations")
        self.budget.charge("cartesian_calls", 3)
        jets = tuple(cartesian_jet(field._regular(radius, order+1, side), point, order=order+1)
                     for field in self.vfields)
        return {"value": _real(fsum(jet["gradient"][i] for i, jet in enumerate(jets)), "computed divergence"),
                "gradient": tuple(_real(fsum(jet["hessian"][a][i] for i, jet in enumerate(jets)), "divergence gradient")
                                  for a in range(3)) if order >= 1 else None,
                "hessian": tuple(tuple(_real(fsum(jet["third"][a][b][i] for i, jet in enumerate(jets)), "divergence Hessian")
                                       for b in range(3)) for a in range(3)) if order >= 2 else None}

    def exterior_remainder_jet(self, point, order=2):
        _jet_order(order, 2)
        point = _point(point)
        radius = sqrt(fsum(x*x for x in point))
        self.budget.charge("point_evaluations")
        return self._combine(self._regular(radius, order, "right", True), point, order)

    def origin_series_data(self):
        return {"components": self.origin_coefficients, "nominal_polynomial_bridge": True,
                "inherited_errors": self.source_record["inherited_errors"],
                "total_error_certified": False}

    def exterior_data(self):
        return {"components": self.algebraic,
                "remainder_available": any(field.provider is not None for field in (*self.vfields, self.chi)),
                "inherited_errors": self.source_record["inherited_errors"],
                "total_error_certified": False}

    def diagnostics(self):
        return {"operator": "Delta_L_inverse", "source": self.source_record,
                "first_inverse": [field.diagnostics() for field in self.vfields],
                "second_inverse": self.chi.diagnostics(), "work": self.budget.snapshot(),
                "derivative_axes_first": True, "fixture_targets_used": False,
                "full_H6_validation": False, "total_error_certified": False}


def solve_york_source(grid, vector_specs, angular, *, h_lmax, guard_lmax,
                      support_tolerance, budget, progress=None):
    if not isinstance(grid, RadialGrid) or not isinstance(angular, AngularGrid) or not isinstance(budget, WorkBudget):
        raise TypeError("expected grid, angular grid and shared WorkBudget")
    specs = tuple(vector_specs)
    if len(specs) != 3 or any(not isinstance(spec, ScalarSourceSpec) for spec in specs):
        raise ValueError("York source must contain three ScalarSourceSpecs")
    _integer(h_lmax, "h_lmax")
    _integer(guard_lmax, "guard_lmax")
    tolerance = _real(support_tolerance, "support tolerance", True)
    if guard_lmax <= h_lmax or angular.ell_max != guard_lmax:
        raise ValueError("angular grid must cover a strictly larger complete guard band")
    patches = {origin.radius for spec in specs for origin in spec.origin_series.values()}
    if len(patches) != 1:
        raise ValueError("York source modes require one common explicit origin patch")
    patch = next(iter(patches))
    if patch not in grid.breakpoints:
        raise ValueError("origin patch must be an explicit grid join")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable")
    vfields = tuple(solve_scalar_source(grid, spec, inverse="P", budget=budget, progress=progress) for spec in specs)
    # All central origins and algebraic exteriors come from the computed first
    # solves. No analytical source mass or a0 is imposed on those responses.
    derivatives = [_origin_derivative(field.origin_coefficients, axis, guard_lmax, budget)
                   for axis, field in enumerate(vfields)]
    origin = _sum_origins(derivatives, (Fraction(1),)*3, budget)
    exterior_field = _sum_laurent([_laurent_derivative(field.algebraic, axis, budget)
                                  for axis, field in enumerate(vfields)],
                                 (Fraction(1),)*3, grid.radius, budget)
    exterior_modes, _, exterior_stats = _laurent_project(exterior_field, guard_lmax, budget)
    providers = []
    for axis, field in enumerate(vfields):
        if field.provider is not None:
            budget.provider_calls += 1
            providers.append(field.provider.response_remainder_derivative(axis, alpha=1))
    residual = SumProviders(providers, grid.radius) if providers else None
    if residual is not None:
        budget.provider_calls += 1
    provider_bounds = residual.weighted_outer_bounds() if residual is not None else {}
    allowed = tuple((ell, m) for ell in range(h_lmax+1) for m in range(-ell, ell+1))
    allowed_set = set(allowed)
    unsupported_origin = [{"mode": mode, "coefficients": row} for mode, row in origin.items()
                          if mode not in allowed_set and any(row)]
    unsupported_exterior = [{"mode": mode, "terms": [{"p": p, "q": q, "coefficient": str(c)}
                             for (p, q), c in row.items()]} for mode, row in exterior_modes.items()
                            if mode not in allowed_set and row]
    unsupported_provider = [mode for mode in provider_bounds if mode not in allowed_set]
    domain_failures = [{"mode": mode, "p": p, "q": q, "coefficient": str(c)}
        for mode, row in exterior_modes.items() for (p, q), c in row.items()
        if c and (p < 2 or p+mode[0] <= 2 or (p == 2 and (mode[0] < 1 or q != 0)))]
    structural = {"origin_support_failures": unsupported_origin,
                  "exterior_support_failures": unsupported_exterior,
                  "provider_support_failures": unsupported_provider, "domain_failures": domain_failures,
                  "exterior_projection": exterior_stats}
    _emit(progress, "york_exact_support_gate", passed=not any((unsupported_origin, unsupported_exterior,
          unsupported_provider, domain_failures)), diagnostics=structural, work=budget.snapshot())
    if any((unsupported_origin, unsupported_exterior, unsupported_provider, domain_failures)):
        raise ValueError("computed divergence exceeds declared support or scalar inverse domain")
    budget.charge("stored_values", len(angular.modes)*(grid.intervals+1))
    samples = {mode: [] for mode in angular.modes}
    limits = {mode: {} for mode in angular.modes}
    guard_maxima = {mode: 0.0 for mode in angular.modes if mode not in allowed_set}
    reconstruction_max = 0.0
    origin_projection_max = 0.0

    def project_shell(radius, side):
        regular = [field._regular(radius, 1, side) for field in vfields]
        values = []
        for direction in angular.directions:
            budget.charge("point_evaluations")
            budget.charge("cartesian_calls", 3)
            point = tuple(radius*x for x in direction)
            values.append(fsum(cartesian_jet(regular[i], point, order=1)["gradient"][i] for i in range(3)))
        budget.charge("projection_products", 2*len(values)*len(angular.modes))
        projected = angular.project(values)
        reconstructed = angular.reconstruct(projected)
        error = max((abs(a-b)/max(1.0, abs(a), abs(b)) for a, b in zip(values, reconstructed)), default=0.0)
        for mode in guard_maxima:
            guard_maxima[mode] = max(guard_maxima[mode], abs(projected[mode]))
        return projected, error

    for i, radius in enumerate(grid.radii):
        side = "left" if radius in grid.breakpoints or radius == grid.radius else None
        left, error = project_shell(radius, side)
        reconstruction_max = max(reconstruction_max, error)
        # This is algebraically computed source data, not a fitted angular
        # coefficient or an expected fixture. Check its independent projection
        # before using it to avoid spurious nonregular coefficients at r=0.
        if radius <= patch:
            for mode in allowed:
                polynomial_value = fsum(c*radius**(mode[0]+2*k)
                                        for k, c in enumerate(origin.get(mode, (0.0,))))
                origin_projection_max = max(origin_projection_max,
                    abs(left[mode]-polynomial_value)/max(1.0, abs(left[mode]), abs(polynomial_value)))
                left[mode] = polynomial_value
        for mode in angular.modes:
            samples[mode].append(left[mode])
        if radius in grid.breakpoints:
            right, error = project_shell(radius, "right")
            reconstruction_max = max(reconstruction_max, error)
            for mode in angular.modes:
                limits[mode][radius] = (left[mode], right[mode])
        if i % 16 == 0 or i == grid.intervals:
            _emit(progress, "york_divergence_shells", completed=i+1, total=grid.intervals+1,
                  radius=radius, work=budget.snapshot())
    passed = (reconstruction_max <= tolerance and origin_projection_max <= tolerance
              and all(value <= tolerance for value in guard_maxima.values()))
    angular_gate = {"passed": passed, "tolerance": tolerance,
                    "reconstruction_max": reconstruction_max,
                    "origin_polynomial_projection_max": origin_projection_max,
                    "guard_maxima": [{"mode": mode, "value": value} for mode, value in guard_maxima.items()]}
    inherited = [{"component": axis, "residual_policy": field.spec.residual_policy,
                  "source_outer_bounds": [{"mode": mode, "bound": bound}
                      for mode, bound in field.spec.exterior.remainder_outer_bounds.items()],
                  "source_origin_bounds": [{"mode": mode, "bounds": row.error_bounds}
                      for mode, row in field.spec.origin_series.items()],
                  "fd_and_a0_error_certified": False,
                  "provider_arithmetic_bounds_require_separate_query": field.provider is not None}
                 for axis, field in enumerate(vfields)]
    if residual is not None:
        budget.provider_calls += 1
    source_record = {"h_lmax": h_lmax, "guard_lmax": guard_lmax, "exact_support": structural,
                     "angular_gate": angular_gate, "inherited_errors": inherited,
                     "nominal_origin_is_computed_polynomial_not_exact_physical_source": True,
                     "origin_samples_from_checked_computed_polynomial": True,
                     "origin": [{"mode": mode, "coefficients": row} for mode, row in origin.items()],
                     "exterior": [{"mode": mode, "terms": [{"p": p, "q": q, "coefficient": str(c)}
                                    for (p, q), c in row.items()]} for mode, row in exterior_modes.items()],
                     "residual": residual.provenance() if residual is not None else None}
    _emit(progress, "york_source_result", source=source_record,
          first_inverse=[field.diagnostics() for field in vfields], work=budget.snapshot())
    if not passed:
        raise ValueError("computed divergence failed the complete angular guard/reconstruction gate")
    # Only now restrict from the explicitly tested guard band. Preserve every m
    # in the accepted band, including zeros and permitted numerical errors.
    origin_specs = {mode: OriginSeries(patch, origin.get(mode, (0.0,)), (0.0, 0.0, 0.0)) for mode in allowed}
    exterior = ExteriorSource(grid.radius,
        {mode: dict(exterior_modes.get(mode, {})) for mode in allowed},
        {mode: provider_bounds.get(mode, 0.0) for mode in allowed})
    spec = ScalarSourceSpec({mode: tuple(samples[mode]) for mode in allowed}, origin_specs, exterior,
        residual_policy="exact_correction" if residual is not None else "zero",
        residual_provider=residual, source_limits={mode: limits[mode] for mode in allowed},
        provenance={"kind": "divergence_of_computed_first_inverse", "inherited_errors": inherited,
                    "zero_nominal_Taylor_remainder_is_not_zero_physical_error": True})
    chi = solve_scalar_source(grid, spec, inverse="P", budget=budget, progress=progress)
    w_origin, w_exterior = [], []
    final_lmax = max(h_lmax+1, max((mode[0] for field in vfields for mode in field.mode_keys), default=0))
    for axis, field in enumerate(vfields):
        gradient = _origin_derivative(chi.origin_coefficients, axis, final_lmax, budget)
        w_origin.append(_sum_origins((field.origin_coefficients, gradient), (Fraction(-1), Fraction(-1, 4)), budget))
        w_exterior.append(_sum_laurent((field.algebraic, _laurent_derivative(chi.algebraic, axis, budget)),
                                      (Fraction(-1), Fraction(-1, 4)), grid.radius, budget))
    result = YorkField(vfields, chi, w_origin, w_exterior, source_record, budget)
    _emit(progress, "york_field_result", result=result.diagnostics())
    return result
