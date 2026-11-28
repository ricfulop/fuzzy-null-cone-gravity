"""Two-stage independent York inverse controls, evaluated only on request.

run_york_controls(config, progress=None) solves V=P f, constructs h=div(V)
from those computed fields (origin AND exterior included), solves chi=P h,
and combines W=-V-grad(chi)/4. No production inverse/HEAT/Fourier/tail module
is imported. The only numerical backend dependencies are the new spherical
scalar and solid-harmonic modules. Import runs no construction or experiment.

Fixtures: a nonaxial uniform vector ball; and f=grad(1-r^2)^2 inside the unit
ball, zero outside. Exact targets are in spherical-green-control-design.md.
For the latter, independent polynomial/power integration additionally gives
<f,W>=-128*pi/1155, ||LW||^2=256*pi/1155, and exterior norm
128*pi/(1225*R^3). These targets never determine a computed source or boundary.

Configuration defaults: york_intervals=[128,256,512], york_radius=4,
york_origin_radius=.5, york_component_tolerance=.003,
york_value_tolerance=.0003, york_integral_tolerance=.003,
york_guard_tolerance=1e-10, york_n_mu=8, york_n_phi=16,
york_ell_max=3, york_tail_order=12, york_max_point_evaluations=5000000.
The fixed physical origin patch is .5. Both .5 and 1 are declared joins:
the first preserves the computed bridge/FD derivative discrepancy in h.
Separate radial/angle/radius configurations are the driver's responsibility.
Exterior radial GL uses its own order, held fixed when angular order changes.

Returned keys: passed, rows, levels, refinements, settings, maximum,
limitations, work, scope, full_H6_validation. Every completed comparison row
is sent immediately through progress, including failed guard rows. Values
must remain finite. A work-cap exception preserves preceding callback rows.
source_result callbacks preserve computed h coefficient maps and both inverse
diagnostics before point/integral work; level_result preserves each complete
level immediately, including a stopped support-gate level.
Raw row passed flags are never relabeled: coarse physical failures remain
visible. Acceptance requires structural gates at EVERY level, physical targets
at the finest supplied level, and the latest radial refinement. A single-level
call requires all its physical comparisons. maximum includes diagnostic rows;
maximum_required and required/diagnostic failure counts distinguish acceptance.
Fraction.from_float is used only to perform exact polynomial algebra on
computed binary64 coefficient surrogates; it does not certify their accuracy.
"""

from fractions import Fraction
from math import factorial, fsum, isfinite, pi, sqrt

from spherical_green_reference import (
    ExteriorSource, OriginSeries, RadialGrid, combine_york_jets,
    solve_scalar_modes,
)
from solid_harmonic_reference import AngularGrid, cartesian_jet, fischer_project, solid_harmonic


def _real(value, name, positive=False):
    if isinstance(value, bool):
        raise ValueError(name + " must be real")
    value = float(value)
    if not isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value


def _integer(value, name, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError(name + " must be an integer >= " + str(minimum))
    return value


def _settings(config):
    levels = tuple(_integer(n, "york_intervals", 12)
                   for n in config.get("york_intervals", (128, 256, 512)))
    if not levels or any(a >= b for a, b in zip(levels, levels[1:])):
        raise ValueError("york_intervals must be nonempty and strictly increasing")
    settings = {
        "intervals": levels, "radius": _real(config.get("york_radius", 4), "radius", True),
        "origin_radius": _real(config.get("york_origin_radius", .5), "origin radius", True),
        "component_tolerance": _real(config.get("york_component_tolerance", .003), "component tolerance", True),
        "value_tolerance": _real(config.get("york_value_tolerance", .0003), "value tolerance", True),
        "integral_tolerance": _real(config.get("york_integral_tolerance", .003), "integral tolerance", True),
        "guard_tolerance": _real(config.get("york_guard_tolerance", 1e-10), "guard tolerance", True),
        "n_mu": _integer(config.get("york_n_mu", 8), "n_mu", 5),
        "n_phi": _integer(config.get("york_n_phi", 16), "n_phi", 9),
        "ell_max": _integer(config.get("york_ell_max", 3), "ell_max", 3),
        "tail_order": _integer(config.get("york_tail_order", 12), "tail order", 4),
        "max_point_evaluations": _integer(config.get("york_max_point_evaluations", 5_000_000), "point cap"),
        "uniform_ball_vector": (1.0, 2.0, -1.0),
    }
    if settings["radius"] < 2 or settings["origin_radius"] != .5:
        raise ValueError("these controls require R>=2 and a fixed physical origin patch .5")
    return settings


def _flatten(value):
    if isinstance(value, (tuple, list)):
        return [entry for item in value for entry in _flatten(item)]
    return [_real(value, "comparison value")]


def _row(rows, progress, case, n, name, actual, target, tolerance,
         acceptance_required=True, role="structural", **extra):
    left, right = _flatten(actual), _flatten(target)
    if len(left) != len(right):
        raise ValueError("comparison shapes differ: " + name)
    residual = max((abs(a-b)/max(1.0, abs(a), abs(b)) for a, b in zip(left, right)), default=0.0)
    row = {"case": case, "intervals": n, "check": name, "actual": actual,
           "target": target, "scaled_residual": residual, "tolerance": tolerance,
           "passed": residual <= tolerance, "acceptance_required": acceptance_required,
           "role": role, **extra}
    rows.append(row)
    if progress is not None:
        progress({"stage": "york_control_row", "row": row})
    return row


def _charge(work, settings, stage, cartesian_calls):
    if work["point_evaluations"] >= settings["max_point_evaluations"]:
        raise RuntimeError("York controls exceeded declared point-evaluation cap")
    work["point_evaluations"] += 1
    work[stage] += 1
    work["cartesian_jet_calls"] += cartesian_calls


def _add(out, powers, value):
    value = out.get(powers, Fraction(0))+value
    if value:
        out[powers] = value
    else:
        out.pop(powers, None)


def _multiply(left, right):
    out = {}
    for p, a in left.items():
        for q, b in right.items():
            _add(out, tuple(x+y for x, y in zip(p, q)), a*b)
    return out


def _derivative(poly, axis):
    out = {}
    for powers, value in poly.items():
        if powers[axis]:
            reduced = list(powers)
            reduced[axis] -= 1
            _add(out, tuple(reduced), value*powers[axis])
    return out


def _times_coordinate(poly, axis):
    out = {}
    for powers, value in poly.items():
        shifted = list(powers)
        shifted[axis] += 1
        _add(out, tuple(shifted), value)
    return out


def _computed_origin(vfields, ell_max):
    """Differentiate computed V Taylor polynomials, including their a0."""
    r2 = {(2, 0, 0): Fraction(1), (0, 2, 0): Fraction(1), (0, 0, 2): Fraction(1)}
    total = {}
    for axis, field in enumerate(vfields):
        for mode in field.mode_keys:
            harmonic = solid_harmonic(*mode).polynomial
            power = {(0, 0, 0): Fraction(1)}
            poly = {}
            for k, derivative in enumerate(field.regular_jet(*mode, 0, order=3)):
                coefficient = Fraction.from_float(derivative)/factorial(k)
                for powers, value in _multiply(harmonic, power).items():
                    _add(poly, powers, coefficient*value)
                power = _multiply(power, r2)
            for powers, value in _derivative(poly, axis).items():
                _add(total, powers, value)
    coefficients, stats = fischer_project(total, ell_max)
    modes = {}
    for (ell, m, k), value in coefficients.items():
        modes.setdefault((ell, m), {})[k] = value
    return modes, stats


def _computed_exterior(vfields, ell_max):
    """Exact derivative of computed homogeneous scalar exteriors; no target moments."""
    modes, projection_stats = {}, []
    for axis, field in enumerate(vfields):
        for mode in field.mode_keys:
            ell, _ = mode
            harmonic = solid_harmonic(*mode).polynomial
            for (p, q), number in field.exterior_terms(*mode).items():
                if q or p != ell+1:
                    raise ValueError("compact fixture V must have homogeneous, log-free scalar exterior")
                coefficient = Fraction.from_float(number)
                exponent = -p-ell
                pieces = ((_derivative(harmonic, axis), exponent, coefficient),
                          (_times_coordinate(harmonic, axis), exponent-2, exponent*coefficient))
                for poly, radial_power, factor in pieces:
                    projected, stats = fischer_project(poly, ell_max)
                    projection_stats.append(stats)
                    for (degree, m, k), value in projected.items():
                        tail_power = -(radial_power+2*k+degree)
                        bucket = modes.setdefault((degree, m), {})
                        _add(bucket, (tail_power, 0), factor*value)
    return {mode: terms for mode, terms in modes.items() if terms}, projection_stats


def _source(case, point, side, vector, radius=None):
    # Use the declared shell radius at a join, not a rounded unit-vector norm.
    r2 = fsum(x*x for x in point) if radius is None else radius*radius
    inside = r2 < 1 or (r2 == 1 and side == "left")
    if not inside:
        return (0.0, 0.0, 0.0)
    if case == "uniform_ball":
        return vector
    return tuple(-4*x*(1-r2) for x in point)


def _first_solve(case, grid, settings, progress, work):
    fields = []
    for axis in range(3):
        if case == "uniform_ball":
            mode, factor, coefficients = (0, 0), settings["uniform_ball_vector"][axis], (1.0,)
        else:
            mode = ((1, 1), (1, -1), (1, 0))[axis]
            factor, coefficients = (-1.0 if axis < 2 else 1.0), (-4.0, 4.0)
        samples = tuple(factor*r**mode[0]*fsum(c*r**(2*k) for k, c in enumerate(coefficients))
                        if r <= 1 else 0.0 for r in grid.radii)
        at_patch = factor*.5**mode[0]*fsum(c*.5**(2*k) for k, c in enumerate(coefficients))
        at_join = factor*fsum(coefficients)
        origin = OriginSeries(.5, tuple(factor*c for c in coefficients), (0.0, 0.0, 0.0))
        field = solve_scalar_modes(grid, {mode: samples}, inverse="P",
                    exterior=ExteriorSource(grid.radius, {}, {mode: 0.0}), origin_series={mode: origin},
                    source_limits={mode: {.5: (at_patch, at_patch), 1.0: (at_join, 0.0)}}, progress=progress)
        work["fd_unknowns"] += grid.intervals+1
        fields.append(field)
    return tuple(fields)


def _regular_fields(fields, radius, side, order):
    return tuple({mode: field.regular_jet(*mode, radius, order=order, side=side)
                  for mode in field.mode_keys} for field in fields)


def _h_shell(vfields, radius, side, angular, settings, work):
    regular = _regular_fields(vfields, radius, side, 1)
    values = []
    for direction in angular.directions:
        _charge(work, settings, "source_sampling_points", 3)
        point = tuple(radius*x for x in direction)
        values.append(fsum(cartesian_jet(regular[i], point, order=1)["gradient"][i] for i in range(3)))
    work["angular_projection_products"] += len(values)*len(angular.modes)
    projected = angular.project(values)
    reconstructed = angular.reconstruct(projected)
    residual = max(abs(a-b)/max(1.0, abs(a), abs(b)) for a, b in zip(values, reconstructed))
    return projected, residual


def _second_solve(case, grid, vfields, angular, settings, rows, progress, work):
    allowed = tuple(mode for mode in angular.modes if mode[0] in ((1,) if case == "uniform_ball" else (0, 2)))
    samples = {mode: [] for mode in angular.modes}
    limits = {mode: {} for mode in angular.modes}
    maxima = {mode: 0.0 for mode in angular.modes}
    reconstruction_max = 0.0
    for index, radius in enumerate(grid.radii):
        side = "left" if radius in grid.breakpoints or radius == grid.radius else None
        projected, residual = _h_shell(vfields, radius, side, angular, settings, work)
        reconstruction_max = max(reconstruction_max, residual)
        for mode, value in projected.items():
            samples[mode].append(value)
            maxima[mode] = max(maxima[mode], abs(value))
        if radius in grid.breakpoints:
            right, residual = _h_shell(vfields, radius, "right", angular, settings, work)
            reconstruction_max = max(reconstruction_max, residual)
            for mode in angular.modes:
                limits[mode][radius] = projected[mode], right[mode]
                maxima[mode] = max(maxima[mode], abs(right[mode]))
        if progress is not None and (index % 32 == 0 or index == grid.intervals):
            progress({"stage": "york_h_shells", "case": case, "intervals": grid.intervals,
                      "completed": index+1, "total": grid.intervals+1,
                      "point_evaluations": work["point_evaluations"]})
    origin, origin_stats = _computed_origin(vfields, settings["ell_max"])
    exterior, exterior_stats = _computed_exterior(vfields, settings["ell_max"])
    guards = []
    for mode in angular.modes:
        if mode not in allowed:
            guards.append(_row(rows, progress, case, grid.intervals, "h_guard_"+str(mode),
                               maxima[mode], 0.0, settings["guard_tolerance"], mode=mode))
    guards.append(_row(rows, progress, case, grid.intervals, "h_angular_reconstruction",
                       reconstruction_max, 0.0, settings["guard_tolerance"]))
    unsupported_origin = tuple((mode, k, str(value)) for mode, terms in origin.items()
                               for k, value in terms.items() if mode not in allowed and value)
    unsupported_tail = tuple((mode, pair, str(value)) for mode, terms in exterior.items()
                             for pair, value in terms.items()
                             if value and (mode not in allowed or pair[0] < 2 or pair[0]+mode[0] <= 2
                                           or (pair[0] == 2 and (mode[0] < 1 or pair[1]))))
    guards.append(_row(rows, progress, case, grid.intervals, "computed_origin_support",
                       len(unsupported_origin), 0, 0, unsupported=unsupported_origin))
    guards.append(_row(rows, progress, case, grid.intervals, "computed_exterior_support_and_monopole",
                       len(unsupported_tail), 0, 0, unsupported=unsupported_tail))
    metadata = {"h_guard_channel_maxima": [{"mode": mode, "value": value} for mode, value in maxima.items()],
                "origin_projection": origin_stats, "exterior_projections": exterior_stats,
                "computed_origin_h": [{"mode": mode, "coefficients": {str(k): str(c) for k, c in terms.items()}}
                                      for mode, terms in origin.items()],
                "computed_exterior_h": [{"mode": mode, "terms": [{"p": p, "q": q, "coefficient": str(c)}
                                         for (p, q), c in terms.items()]} for mode, terms in exterior.items()],
                "exact_surrogates_are_computed_binary64_coefficients": True}
    if not all(row["passed"] for row in guards):
        return None, metadata
    # Restriction happens only after every unsupported channel above is gated.
    origins = {}
    for mode in allowed:
        terms = origin.get(mode, {})
        coefficients = tuple(float(terms.get(k, 0)) for k in range(max(terms, default=0)+1))
        origins[mode] = OriginSeries(.5, coefficients, (0.0, 0.0, 0.0))
    tail = {mode: dict(exterior.get(mode, {})) for mode in allowed}
    chi = solve_scalar_modes(grid, {mode: tuple(samples[mode]) for mode in allowed}, inverse="P",
              exterior=ExteriorSource(grid.radius, tail, {mode: 0.0 for mode in allowed}),
              origin_series=origins, source_limits={mode: limits[mode] for mode in allowed}, progress=progress)
    work["fd_unknowns"] += len(allowed)*(grid.intervals+1)
    return chi, metadata


def _potentials(vfields, chi, radius, side):
    return (_regular_fields(vfields, radius, side, 2),
            {mode: chi.regular_jet(*mode, radius, order=3, side=side) for mode in chi.mode_keys})


def _evaluate(potentials, point, settings, work):
    _charge(work, settings, "field_evaluation_points", 4)
    vregular, chiregular = potentials
    vjets = tuple(cartesian_jet(regular, point, order=2) for regular in vregular)
    chijet = cartesian_jet(chiregular, point, order=3)
    result = combine_york_jets(vjets, chijet)
    result["V"] = tuple(jet["value"] for jet in vjets)
    result["chi"] = chijet["value"]
    return result


def _power(s, exponent):
    return s**exponent, exponent*s**(exponent-1), exponent*(exponent-1)*s**(exponent-2)


def _target(case, point, side, vector):
    """Direct differentiated analytic fixture, independent of harmonic/FD jets."""
    s = fsum(x*x for x in point)
    inside = s < 1 or (s == 1 and side == "left")
    if case == "uniform_ball":
        dot = fsum(x*c for x, c in zip(point, vector))
        if inside:
            a, da, dda = -11/24+19*s/120, 19/120, 0.0
            b, db, ddb = -1/60, 0.0, 0.0
            chi, h = dot*(-1/6+s/30), -dot/3
            v = tuple(c*(1/2-s/6) for c in vector)
        else:
            p1, p3, p5 = _power(s, -.5), _power(s, -1.5), _power(s, -2.5)
            a, da, dda = tuple(-7*x/24-y/120 for x, y in zip(p1, p3))
            b, db, ddb = tuple(-x/24+y/40 for x, y in zip(p3, p5))
            chi, h = dot*(-p1[0]/6+p3[0]/30), -dot*p3[0]/3
            v = tuple(c*p1[0]/3 for c in vector)
        w = tuple(a*vector[i]+b*point[i]*dot for i in range(3))
        dw = tuple(tuple(2*point[j]*da*vector[i]+2*point[j]*db*point[i]*dot +
                         b*((i == j)*dot+point[i]*vector[j]) for i in range(3)) for j in range(3))
        ddw = tuple(tuple(tuple((2*(j == k)*da+4*point[j]*point[k]*dda)*vector[i] +
                     (2*(j == k)*db+4*point[j]*point[k]*ddb)*point[i]*dot +
                     2*point[j]*db*((i == k)*dot+point[i]*vector[k]) +
                     2*point[k]*db*((i == j)*dot+point[i]*vector[j]) +
                     b*((i == j)*vector[k]+(i == k)*vector[j]) for i in range(3))
                     for k in range(3)) for j in range(3))
    else:
        if inside:
            a, da, dda = 1/4-3*s/10+3*s*s/28, -3/10+3*s/14, 3/14
            chi, h = -(1/6-s/6+s*s/10-s**3/42), -(1-s)**2
            v = tuple(x*(-1/3+2*s/5-s*s/7) for x in point)
        else:
            a, da, dda = tuple(2*x/35 for x in _power(s, -1.5))
            chi, h = -8/(105*sqrt(s)), 0.0
            v = tuple(-8*x/(105*s**1.5) for x in point)
        w = tuple(a*x for x in point)
        dw = tuple(tuple(a*(i == j)+2*da*point[i]*point[j] for i in range(3)) for j in range(3))
        ddw = tuple(tuple(tuple(2*da*((j == k)*point[i]+(i == j)*point[k]+(i == k)*point[j]) +
                                 4*dda*point[i]*point[j]*point[k] for i in range(3))
                          for k in range(3)) for j in range(3))
    return {"W": w, "dW": dw, "ddW": ddw, "chi": chi, "div_V": h, "V": v}


def _longitudinal(dw):
    divergence = fsum(dw[i][i] for i in range(3))
    return tuple(tuple(dw[i][j]+dw[j][i]-(2/3)*divergence*(i == j)
                       for j in range(3)) for i in range(3))


def _radial_weights(grid):
    boundaries = (0, *grid.breakpoint_indices, grid.intervals)
    for lo, hi in zip(boundaries, boundaries[1:]):
        weights = {i: 0.0 for i in range(lo, hi+1)}
        count = hi-lo
        end = hi if count % 2 == 0 else hi-3
        for i in range(lo, end, 2):
            for j, factor in enumerate((1, 4, 1)):
                weights[i+j] += grid.h*factor/3
        if count % 2:
            for j, factor in enumerate((1, 3, 3, 1)):
                weights[hi-3+j] += 3*grid.h*factor/8
        for i, weight in weights.items():
            side = "right" if i == lo and lo else "left" if i == hi else None
            yield grid.radii[i], side, weight


def _integrals(case, grid, vfields, chi, angular, tail_rule, settings, work, progress=None):
    pairing_terms, norm_terms = [], []
    for shell, (radius, side, radial_weight) in enumerate(_radial_weights(grid)):
        if radius == 0:
            continue
        potentials = _potentials(vfields, chi, radius, side)
        for direction, weight in zip(angular.directions, angular.weights):
            point = tuple(radius*x for x in direction)
            result = _evaluate(potentials, point, settings, work)
            tensor = _longitudinal(result["dW"])
            measure = radial_weight*radius*radius*weight
            source = _source(case, point, side, settings["uniform_ball_vector"], radius)
            pairing_terms.append(measure*fsum(a*b for a, b in zip(source, result["W"])))
            norm_terms.append(measure*fsum(x*x for row in tensor for x in row))
        if progress is not None and (shell % 16 == 0 or radius == grid.radius):
            progress({"stage": "york_volume_integrals", "case": case, "intervals": grid.intervals,
                      "completed_shells": shell+1, "radius": radius,
                      "point_evaluations": work["point_evaluations"]})
    outside_terms = []
    # This independent one-dimensional order is fixed across angular slices.
    for t, radial_weight in tail_rule:
        radius = grid.radius/t
        potentials = _potentials(vfields, chi, radius, None)
        for direction, weight in zip(angular.directions, angular.weights):
            result = _evaluate(potentials, tuple(radius*x for x in direction), settings, work)
            tensor = _longitudinal(result["dW"])
            outside_terms.append(radial_weight*grid.radius/t**2*radius**2*weight*
                                 fsum(x*x for row in tensor for x in row))
    potentials = _potentials(vfields, chi, grid.radius, "right")
    flux_terms = []
    for direction, weight in zip(angular.directions, angular.weights):
        result = _evaluate(potentials, tuple(grid.radius*x for x in direction), settings, work)
        tensor = _longitudinal(result["dW"])
        flux_terms.append(-2*grid.radius**2*weight*fsum(result["W"][i]*tensor[i][j]*direction[j]
                                                      for i in range(3) for j in range(3)))
    inner, outside = fsum(norm_terms), fsum(outside_terms)
    return {"source_pairing": fsum(pairing_terms), "norm_inside": inner,
            "norm_outside": outside, "norm_total": inner+outside,
            "exterior_flux_completion": fsum(flux_terms)}


def _point_checks(case, grid, vfields, chi, settings, rows, progress, work):
    points = ((0., 0., 0.), (.2, -.15, .1), (.6, .2, -.1), (.5, 0., 0.),
              (1., 0., 0.), (1.25, -.5, .25), (.3*grid.radius, -.4*grid.radius, .5*grid.radius),
              (1.2*grid.radius, .2*grid.radius, -.1*grid.radius))
    snapshots = {}
    for point in points:
        radius = sqrt(fsum(x*x for x in point))
        sides = ("left", "right") if radius in grid.breakpoints or radius == grid.radius else (None,)
        for side in sides:
            actual = _evaluate(_potentials(vfields, chi, radius, side), point, settings, work)
            target = _target(case, point, side, settings["uniform_ball_vector"])
            label = repr((point, side))
            snapshots[label] = actual
            for key in ("W", "dW", "ddW", "chi", "div_V", "V"):
                tolerance = settings["component_tolerance"] if key in ("dW", "ddW", "div_V") else settings["value_tolerance"]
                _row(rows, progress, case, grid.intervals, "point_"+key, actual[key], target[key], tolerance,
                     acceptance_required=grid.intervals == settings["intervals"][-1],
                     role="physical_target", point=point, side=side)
    return snapshots


def run_york_controls(config, progress=None):
    """Run only when explicitly invoked by the tracked driver; strict finite rows."""
    settings = _settings(config)
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable or None")
    angular = AngularGrid(settings["ell_max"], settings["n_mu"], settings["n_phi"])
    # Reuse only the independent GL implementation; this table has its own
    # order and one azimuthal point, unrelated to the integration sphere grid.
    tail_grid = AngularGrid(0, settings["tail_order"], 1)
    tail_rule = tuple(((1+direction[2])/2, weight/(4*pi))
                      for direction, weight in zip(tail_grid.directions, tail_grid.weights))
    rows, levels, refinements, previous = [], [], [], {}
    work = {"point_evaluations": 0, "source_sampling_points": 0,
            "field_evaluation_points": 0, "cartesian_jet_calls": 0,
            "angular_projection_products": 0, "fd_unknowns": 0}
    for n in settings["intervals"]:
        grid = RadialGrid(settings["radius"], n, (.5, 1.0))
        for case in ("uniform_ball", "compact_gradient"):
            start_row = len(rows)
            def stage(update):
                if progress is not None:
                    progress({"case": case, "intervals": n, **update})
            vfields = _first_solve(case, grid, settings, stage, work)
            chi, metadata = _second_solve(case, grid, vfields, angular, settings, rows, progress, work)
            level = {"case": case, "intervals": n, "first_inverse": [field.diagnostics() for field in vfields],
                     "computed_h": metadata, "row_start": start_row}
            if chi is not None:
                level["second_inverse"] = chi.diagnostics()
            if progress is not None:
                progress({"stage": "source_result", "case": case, "intervals": n,
                          "first_inverse": level["first_inverse"],
                          "second_inverse": level.get("second_inverse"),
                          "computed_h": metadata, "work": dict(work)})
            if chi is None:
                level.update(passed=False, raw_passed=False,
                             stopped_at="computed_h_support_gate", row_stop=len(rows))
                levels.append(level)
                if progress is not None:
                    progress({"stage": "level_result", "level": level, "work": dict(work)})
                continue
            for label, field in [("V"+str(i), field) for i, field in enumerate(vfields)] + [("chi", chi)]:
                for diagnostics in field.diagnostics()["rows"]:
                    _row(rows, progress, case, n, "origin_source_bound_"+label,
                         diagnostics["origin_source_sample_max_bound_excess"], 0.0, 0.0,
                         mode=diagnostics["mode"], source_bounds=diagnostics["origin_source_error_bounds_s"])
            points = _point_checks(case, grid, vfields, chi, settings, rows, progress, work)
            integrals = _integrals(case, grid, vfields, chi, angular, tail_rule, settings, work, progress)
            radius = grid.radius
            if case == "uniform_ball":
                c2 = fsum(c*c for c in settings["uniform_ball_vector"])
                pair = -22*pi*c2/45
                norm = 44*pi*c2/45
                outside = pi*c2*(5/(6*radius)-1/(45*radius**3)+1/(150*radius**5))
            else:
                pair, norm, outside = -128*pi/1155, 256*pi/1155, 128*pi/(1225*radius**3)
            targets = {"source_pairing": pair, "norm_total": norm, "norm_outside": outside,
                       "exterior_flux_completion": outside}
            finest = n == settings["intervals"][-1]
            for key, target in targets.items():
                _row(rows, progress, case, n, key, integrals[key], target, settings["integral_tolerance"],
                     acceptance_required=finest, role="physical_target")
            _row(rows, progress, case, n, "independent_exterior_norm_vs_flux", integrals["norm_outside"],
                 integrals["exterior_flux_completion"], settings["integral_tolerance"],
                 acceptance_required=finest, role="physical_target")
            _row(rows, progress, case, n, "norm_vs_negative_twice_pairing", integrals["norm_total"],
                 -2*integrals["source_pairing"], settings["integral_tolerance"],
                 acceptance_required=finest, role="physical_target")
            if case in previous:
                old_n, old_points, old_integrals = previous[case]
                for key in ("W", "dW", "ddW", "chi", "div_V", "V"):
                    change = max(abs(a-b)/max(1.0, abs(a), abs(b)) for label in points
                                 for a, b in zip(_flatten(points[label][key]), _flatten(old_points[label][key])))
                    tolerance = settings["component_tolerance"] if key in ("dW", "ddW", "div_V") else settings["value_tolerance"]
                    refinements.append(_row(rows, progress, case, n, "radial_refinement_"+key,
                                             change, 0.0, tolerance, acceptance_required=finest,
                                             role="radial_refinement", previous_intervals=old_n))
                for key in integrals:
                    refinements.append(_row(rows, progress, case, n, "radial_refinement_"+key,
                                             integrals[key], old_integrals[key], settings["integral_tolerance"],
                                             acceptance_required=finest, role="radial_refinement",
                                             previous_intervals=old_n))
            previous[case] = n, points, integrals
            level.update(points=points, integrals=integrals, row_stop=len(rows),
                         raw_passed=all(row["passed"] for row in rows[start_row:]),
                         passed=all(row["passed"] for row in rows[start_row:] if row["acceptance_required"]))
            levels.append(level)
            if progress is not None:
                progress({"stage": "level_result", "level": level, "work": dict(work)})
    required = [row for row in rows if row["acceptance_required"]]
    return {"passed": bool(required) and all(row["passed"] for row in required), "rows": rows,
            "raw_passed": bool(rows) and all(row["passed"] for row in rows),
            "required_failures": sum(not row["passed"] for row in required),
            "diagnostic_failures": sum(not row["passed"] for row in rows if not row["acceptance_required"]),
            "levels": levels, "refinements": refinements, "settings": settings,
            "maximum": max((row["scaled_residual"] for row in rows), default=0.0),
            "maximum_required": max((row["scaled_residual"] for row in required), default=0.0),
            "work": {**work, "angular_grid": angular.statistics(),
                     "tail_quadrature": {"order": settings["tail_order"], "points": len(tail_rule),
                                         "weight_sum_on_unit_interval": fsum(weight for _, weight in tail_rule),
                                         "independent_of_angular_order": True}},
            "scope": "independent_two_stage_York_inverse_primitives", "full_H6_validation": False,
            "limitations": [
                "Formal second-order radial FD; jump and origin errors require measured refinement.",
                "Source Taylor data for h are derived from computed V, including its numerical a0.",
                "Zero origin/exterior remainders describe computed finite polynomial/power surrogates, not true solution errors.",
                "Unsupported angular channels are checked before restriction; permitted channels retain numerical errors.",
                "ODE-derived jets do not provide an independent finite-difference PDE residual.",
                "The total norm and exterior flux use distinct quadratures but share the computed fields.",
                "This invocation refines radius resolution only; separate angular and exterior-radius slices require driver configurations.",
                "No gravitational source, W2/u4 mixed functional, or four-point amplitude is validated.",
            ]}
