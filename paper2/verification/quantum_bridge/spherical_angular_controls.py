"""Tracked, independent polynomial controls for angular/Kelvin primitives.

run_angular_controls(config, progress=None) returns JSON-compatible evidence.
Import is inert. No production inverse, HEAT, Fourier, or pairing is imported.
All exact reference differentiation, multiplication, sphere moments and radial
trace normalization below are written independently of the tested modules.

The normalized Gram and roundtrips cover EVERY m through at least ell=10.
Numerical tolerances use abs(actual-target)/max(1,abs(target)); normalized
coefficient comparisons use Y coefficients, although module inputs use H.
These are bounded numerical controls, not accuracy certification or validation
of the full gravitational functional. Actual two-stage York inversions and
their localized h00 response belong to the separate radial/York controls.

Progress events angular_control_completed and angular_control_exception carry
a complete control_result, including each completed grid/profile row. An
unexpected exception is emitted before being re-raised; it is not a pass.
Work counts exact products/scatters and numerical summands. Small analytic
upper bounds replace unavailable derivative/Kelvin counters. Basis-cache
construction, Python bookkeeping, and scalar elementary-function costs are
excluded; degree, grid, routine and total caps bound these separately.
"""

from collections.abc import Mapping
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations
from math import comb, exp, factorial, fsum, isfinite, pi, sqrt

from solid_harmonic_reference import (
    AngularGrid, AngularSupportError, GaussianProfile, cartesian_jet,
    fischer_project, solid_harmonic,
)
from kelvin_exterior_reference import CartesianLaurent, kelvin_exterior_basis


def _settings(config):
    if not isinstance(config, Mapping):
        raise TypeError("config must be a mapping")
    defaults = {
        "angular_ell_max": 10,
        "angular_grids": ((12, 24), (16, 32)),
        "angular_tolerance": 2e-10,
        "angular_cartesian_tolerance": 2e-10,
        "angular_profile_tolerance": 2e-12,
        "angular_max_exact_products": 5_000_000,
        "angular_max_terms": 100_000,
        "angular_max_grid_points": 10_000,
        "angular_max_routine_work": 5_000_000,
        "angular_max_total_work": 50_000_000,
    }
    result = {key: config.get(key, value) for key, value in defaults.items()}
    ell = result["angular_ell_max"]
    if type(ell) is not int or not 10 <= ell <= 16:
        raise ValueError("angular_ell_max must be an integer from 10 through 16")
    for key in ("angular_tolerance", "angular_cartesian_tolerance", "angular_profile_tolerance"):
        if isinstance(result[key], bool):
            raise ValueError(key+" must be finite and positive")
        result[key] = float(result[key])
        if not isfinite(result[key]) or result[key] <= 0:
            raise ValueError(key+" must be finite and positive")
    for key in ("angular_max_exact_products", "angular_max_terms", "angular_max_grid_points",
                "angular_max_routine_work", "angular_max_total_work"):
        if type(result[key]) is not int or result[key] < 1:
            raise ValueError(key+" must be a positive integer")
    grids = tuple(tuple(row) for row in result["angular_grids"])
    if len(grids) != 2:
        raise ValueError("angular_grids must contain exactly two refinement grids")
    for row in grids:
        if (len(row) != 2 or any(type(n) is not int for n in row)
                or not ell+1 <= row[0] <= 256 or not 2*ell+1 <= row[1] <= 512
                or row[0]*row[1] > result["angular_max_grid_points"]):
            raise ValueError("angular_grids violates the complete-band or bounded-grid floors")
    if not (grids[1][0] > grids[0][0] and grids[1][1] > grids[0][1]):
        raise ValueError("Both angular grid axes must increase independently of the target tolerance")
    result["angular_grids"] = [list(row) for row in grids]
    result.update({"coefficient_convention": "unnormalized H; numerical coefficient gates in normalized Y",
                   "error_scale": "abs(actual-target)/max(1,abs(target))",
                   "profile_reference_precision": 64,
                   "profile_bound_roundoff_allowance": "1e-50*max(1,abs(reference derivative))",
                   "accuracy_certified": False})
    return result


class _Work:
    def __init__(self, settings):
        self.cap = settings["angular_max_total_work"]
        self.routine_cap = settings["angular_max_routine_work"]
        self.total = 0
        self.by_kind = {}

    @property
    def remaining(self):
        return self.cap-self.total

    def ensure(self, count):
        if count < 0 or count > self.remaining:
            raise ValueError("Angular controls exceed angular_max_total_work")

    def charge(self, count, kind="independent_exact"):
        self.ensure(count)
        self.total += count
        self.by_kind[kind] = self.by_kind.get(kind, 0)+count

    def routine(self, count, kind="numerical_summands"):
        if count > self.routine_cap:
            raise ValueError("Angular control routine exceeds angular_max_routine_work")
        self.charge(count, kind)

    def statistics(self):
        return {"total": self.total, "max_total": self.cap, "max_routine": self.routine_cap,
                "by_kind": dict(self.by_kind),
                "unit": "Exact products/scatters or numerical summands; derivative/Kelvin bounds explicit"}


class _Controls:
    def __init__(self, settings, progress):
        self.settings, self.progress = settings, progress
        self.work = _Work(settings)
        self.rows = []
        self.exception_emitted = False

    def run(self, name, function):
        before = self.work.total
        try:
            row = function()
            row = {"name": name, **row, "counted_work": self.work.total-before}
            if type(row.get("passed")) is not bool:
                raise TypeError("Control must return an explicit boolean passed")
        except Exception as exc:
            row = {"name": name, "passed": False, "exception": type(exc).__name__,
                   "message": str(exc), "counted_work": self.work.total-before,
                   "work": self.work.statistics(),
                   "completed_control_names": [item["name"] for item in self.rows]}
            self.rows.append(row)
            self.exception_emitted = True
            if self.progress is not None:
                self.progress({"stage": "angular_control_exception", "control_result": row})
            raise
        self.rows.append(row)
        if self.progress is not None:
            self.progress({"stage": "angular_control_completed", "control_result": row})
        return row

    def exact(self, function):
        cap = min(self.settings["angular_max_exact_products"], self.work.remaining)
        if cap < 1:
            raise ValueError("Angular controls exhausted angular_max_total_work")
        try:
            result = function({"max_products": cap, "max_terms": self.settings["angular_max_terms"]})
        except AngularSupportError as exc:
            self.work.charge(exc.statistics["work"]["products_and_scatters"], "primitive_reported_exact")
            raise
        except Exception:
            # Only an exceptional call has unknown partial work; reserve its
            # entire already-bounded allowance before propagating the failure.
            self.work.charge(cap, "exceptional_exact_upper_bound")
            raise
        self.work.charge(result[-1]["work"]["products_and_scatters"], "primitive_reported_exact")
        return result

    def bounded(self, count, function):
        if count > self.settings["angular_max_exact_products"]:
            raise ValueError("Analytic exact-call bound exceeds angular_max_exact_products")
        self.work.charge(count, "primitive_exact_upper_bounds")
        return function({"max_products": max(1, count), "max_terms": self.settings["angular_max_terms"]})


def _put(poly, powers, value):
    total = poly.get(powers, F(0))+value
    if total:
        poly[powers] = total
    else:
        poly.pop(powers, None)


def _add(left, right, work, factor=1):
    work.charge(len(right))
    result = dict(left)
    for powers, coefficient in right.items():
        _put(result, powers, factor*coefficient)
    return result


def _multiply(left, right, work):
    work.charge(len(left)*len(right))
    result = {}
    for powers, coefficient in left.items():
        for other, value in right.items():
            _put(result, tuple(a+b for a, b in zip(powers, other)), coefficient*value)
    return result


def _derivative(poly, axis, work):
    work.charge(len(poly))
    result = {}
    for powers, coefficient in poly.items():
        if powers[axis]:
            lower = list(powers)
            lower[axis] -= 1
            _put(result, tuple(lower), coefficient*powers[axis])
    return result


def _differentiate(poly, indices, work):
    result = dict(poly)
    for axis in indices:
        result = _derivative(result, axis, work)
    return result


def _r2_power(power, work):
    result = {(0, 0, 0): F(1)}
    r2 = {(2, 0, 0): F(1), (0, 2, 0): F(1), (0, 0, 2): F(1)}
    for _ in range(power):
        result = _multiply(result, r2, work)
    return result


def _evaluate(poly, point, work, exact=False):
    work.charge(len(poly), "independent_exact" if exact else "independent_numeric_summands")
    if exact:
        return sum((c*point[0]**p[0]*point[1]**p[1]*point[2]**p[2] for p, c in poly.items()), F(0))
    return fsum(float(c)*point[0]**p[0]*point[1]**p[1]*point[2]**p[2] for p, c in poly.items())


def _sphere_monomial(powers):
    """Independent factorial sphere AVERAGE, not a production moment call."""
    if any(power % 2 for power in powers):
        return F(0)
    halves = tuple(power//2 for power in powers)
    n = sum(halves)
    numerator, denominator = factorial(n), factorial(2*n+1)
    for power, half in zip(powers, halves):
        numerator *= factorial(power)
        denominator *= factorial(half)
    return F(numerator, denominator)


def _norm(poly, work):
    work.charge(len(poly)**2)
    return sum((c*d*_sphere_monomial(tuple(a+b for a, b in zip(p, q)))
                for p, c in poly.items() for q, d in poly.items()), F(0))


def _scaled(actual, target):
    actual, target = float(actual), float(target)
    if not isfinite(actual) or not isfinite(target):
        raise ArithmeticError("Nonfinite angular-control comparison")
    return abs(actual-target)/max(1.0, abs(target))


def _harmonics(ell, controls):
    rows = []
    for m in range(-ell, ell+1):
        harmonic = solid_harmonic(ell, m)
        laplacian = {}
        for axis in range(3):
            laplacian = _add(laplacian, _differentiate(harmonic.poly, (axis, axis), controls.work), controls.work)
        norm = _norm(harmonic.poly, controls.work)
        homogeneous = all(sum(p) == ell for p in harmonic.poly)
        rows.append({"mode": [ell, m], "terms": len(harmonic.poly), "homogeneous": homogeneous,
                     "laplacian_terms": len(laplacian), "independent_norm_over_4pi": str(norm),
                     "declared_norm_over_4pi": str(harmonic.norm_over_4pi),
                     "passed": homogeneous and not laplacian and norm == harmonic.norm_over_4pi and norm > 0})
    return {"passed": all(row["passed"] for row in rows), "all_m": True, "modes": rows}


def _azimuthal(power, imaginary=False):
    result = {}
    for y_power in range(power+1):
        if y_power % 2 == int(imaginary):
            result[(power-y_power, y_power, 0)] = F(comb(power, y_power)*(-1)**(y_power//2))
    return result


def _nonaxial_fixture(work):
    poly = {(0, 1, 0): F(3), (1, 0, 1): F(-2)}
    poly = _add(poly, _azimuthal(8), work, F(1, 17))
    poly = _add(poly, _azimuthal(8, True), work, F(-2, 19))
    double_factorial = 1
    for value in range(1, 16, 2):
        double_factorial *= value
    coefficients = {(1, -1): F(-3), (2, 1): F(2, 3),
                    (8, 8): F(1, 17*double_factorial), (8, -8): F(-2, 19*double_factorial)}
    return poly, coefficients


def _grid_control(shape, controls, fixture, projections):
    settings, work = controls.settings, controls.work
    ell = settings["angular_ell_max"]
    modes = tuple((l, m) for l in range(ell+1) for m in range(-l, l+1))
    harmonics = tuple(solid_harmonic(*mode) for mode in modes)
    points, count = shape[0]*shape[1], len(modes)
    table_terms = points*sum(len(harmonic.poly) for harmonic in harmonics)
    gram_terms = points*count*(count+1)//2
    # Actual per-basis calls each use a complete dense normalized table.
    plan = 2*table_terms+gram_terms+2*points*count*count+8*points*count+40*shape[0]**2
    work.ensure(plan)
    work.routine(table_terms, "angular_table_monomials")
    work.charge(40*shape[0]**2, "angular_gauss_root_upper_bound")
    grid = AngularGrid(ell, *shape, max_points=settings["angular_max_grid_points"],
                       max_work=settings["angular_max_routine_work"])
    work.routine(gram_terms, "normalized_gram_summands")
    gram_error, worst = 0.0, None
    for i in range(count):
        for j in range(i, count):
            value = fsum(weight*row[i]*row[j] for weight, row in zip(grid.weights, grid.normalized_table))
            error = _scaled(value, int(i == j))
            if error > gram_error:
                gram_error, worst = error, {"left": list(modes[i]), "right": list(modes[j]), "value": value}
    tolerance = settings["angular_tolerance"]
    basis_rows = []
    for index, (mode, harmonic) in enumerate(zip(modes, harmonics)):
        factor = sqrt(4*pi*float(harmonic.norm_over_4pi))
        work.routine(points*count)
        values = grid.reconstruct({mode: 1/factor}, max_work=settings["angular_max_routine_work"])
        direct_error = max((_scaled(value, _evaluate(harmonic.poly, point, work)/factor)
                            for value, point in zip(values, grid.directions)), default=0.0)
        work.routine(points*count)
        recovered = grid.project(values, max_work=settings["angular_max_routine_work"])
        errors = [_scaled(recovered[other]*grid.factors[j], int(j == index)) for j, other in enumerate(modes)]
        basis_rows.append({"mode": list(mode), "direct_polynomial_error": direct_error,
                           "normalized_roundtrip_error": max(errors),
                           "passed": direct_error <= tolerance and max(errors) <= tolerance})
    mixture_rows, observed = [], {}
    all_mode_y = {mode: F((-1)**(l+m)*(abs(m) % 3+1), (l+1)**2) for mode in modes for l, m in (mode,)}
    all_mode_h = {mode: float(all_mode_y[mode])/grid.factors[index] for index, mode in enumerate(modes)}
    for name, coefficients, target_poly in (("all_modes_signed", all_mode_h, None),
                                             ("y_xz_real_imag_m8", fixture[1], fixture[0])):
        work.routine(points*count)
        reconstructed = grid.reconstruct(coefficients, max_work=settings["angular_max_routine_work"])
        values = reconstructed if target_poly is None else tuple(_evaluate(target_poly, point, work) for point in grid.directions)
        direct_error = max((_scaled(a, b) for a, b in zip(reconstructed, values)), default=0.0)
        work.routine(points*count)
        recovered = grid.project(values, max_work=settings["angular_max_routine_work"])
        normalized = [recovered[mode]*grid.factors[index] for index, mode in enumerate(modes)]
        target = [float(coefficients.get(mode, 0))*grid.factors[index] for index, mode in enumerate(modes)]
        errors = [_scaled(a, b) for a, b in zip(normalized, target)]
        observed[name] = normalized
        mixture_rows.append({"name": name, "direct_polynomial_error": direct_error,
                             "normalized_projection_error": max(errors), "all_modes_retained": len(recovered) == count,
                             "normalized_coefficients": normalized,
                             "passed": direct_error <= tolerance and max(errors) <= tolerance and len(recovered) == count})
    projections.append(observed)
    weight_error = _scaled(fsum(grid.weights), 4*pi)
    return {"passed": gram_error <= tolerance and weight_error <= tolerance
            and all(row["passed"] for row in basis_rows+mixture_rows),
            "grid": list(shape), "grid_statistics": grid.statistics(), "planned_core_work_bound": plan,
            "normalized_gram_max_error": gram_error, "normalized_gram_worst": worst,
            "gram_pairs_checked": count*(count+1)//2, "sphere_weight_error": weight_error,
            "basis_roundtrips": basis_rows, "mixtures": mixture_rows, "tolerance": tolerance}


def _reconstruct(coefficients, work):
    result = {}
    for (ell, m, k), coefficient in coefficients.items():
        term = _multiply(_r2_power(k, work), solid_harmonic(ell, m).poly, work)
        result = _add(result, term, work, coefficient)
    return result


def _fischer_control(poly, ell, controls):
    coefficients, stats = controls.exact(lambda limits: fischer_project(poly, ell, **limits))
    reconstructed = _reconstruct(coefficients, controls.work)
    return {"passed": reconstructed == poly, "source_terms": len(poly), "ell_max": ell,
            "coefficients": [{"lmk": list(key), "coefficient": str(value)} for key, value in sorted(coefficients.items())],
            "exact_reconstruction": reconstructed == poly, "statistics": stats}


def _expected_rejection(function, words):
    try:
        function()
    except (ValueError, TypeError) as exc:
        return {"passed": any(word in str(exc) for word in words), "exception": type(exc).__name__, "message": str(exc)}
    return {"passed": False, "message": "Expected bounded-domain rejection did not occur"}


def _fischer_rejection(controls):
    poly = _azimuthal(4)
    try:
        controls.exact(lambda limits: fischer_project(poly, 3, **limits))
    except AngularSupportError as exc:
        return {"passed": dict(exc.remainder) == poly and bool(exc.remainder),
                "exception": type(exc).__name__, "remainder_is_entire_ell4_source": dict(exc.remainder) == poly,
                "statistics": exc.statistics}
    return {"passed": False, "message": "An exact ell4 source was silently accepted in ell<=3"}


def _regular_polynomial_jets(coefficients, s, work):
    result = {}
    for (ell, m, k), coefficient in coefficients.items():
        row = result.setdefault((ell, m), [F(0)]*4)
        for order in range(4):
            work.charge(1)
            if k >= order:
                row[order] += coefficient*F(factorial(k), factorial(k-order))*s**(k-order)
    return result


def _jet_item(jet, indices):
    value = jet[("value", "gradient", "hessian", "third")[len(indices)]]
    for axis in indices:
        value = value[axis]
    return value


def _compare_jets(poly, regular, point, controls, beta=None):
    work = controls.work
    bound = 20*sum(len(solid_harmonic(*mode).poly) for mode in regular)
    work.routine(bound, "cartesian_jet_monomial_upper_bound")
    jet = cartesian_jet(regular, point, order=3, max_work=controls.settings["angular_max_routine_work"])
    errors, worst = [0.0]*4, [None]*4
    for order in range(4):
        for indices in combinations_with_replacement(range(3), order):
            differentiated = dict(poly)
            for axis in indices:
                derivative = _derivative(differentiated, axis, work)
                if beta is not None:
                    x = tuple(int(i == axis) for i in range(3))
                    derivative = _add(derivative, _multiply(differentiated, {x: -2*beta}, work), work)
                differentiated = derivative
            target = float(_evaluate(differentiated, point, work, exact=True))
            if beta is not None:
                target *= exp(-float(beta*sum(x*x for x in point)))
            for actual_indices in sorted(set(permutations(indices))):
                value = _jet_item(jet, actual_indices)
                error = _scaled(value, target)
                if error > errors[order]:
                    errors[order], worst[order] = error, {"indices": list(actual_indices), "value": value, "target": target}
    tolerance = controls.settings["angular_cartesian_tolerance"] if beta is None else controls.settings["angular_profile_tolerance"]
    return {"point": [str(x) for x in point], "maximum_errors_by_order": errors, "worst_by_order": worst,
            "tensor_entries_checked_by_order": [1, 3, 9, 27],
            "tolerance": tolerance, "passed": max(errors) <= tolerance}


def _polynomial_jet_control(controls):
    poly = {(0, 0, 0): F(2), (1, 0, 0): F(-3), (0, 2, 0): F(5, 7), (1, 0, 1): F(-2),
            (1, 1, 1): F(3), (2, 0, 2): F(-2, 5), (3, 1, 0): F(4, 9), (0, 0, 6): F(1, 11)}
    coefficients, stats = controls.exact(lambda limits: fischer_project(poly, 6, **limits))
    rows = []
    for point in ((F(0),)*3, (F(1, 3), F(-2, 5), F(1, 4)), (F(-1, 2), F(2, 3), F(1, 5))):
        jets = _regular_polynomial_jets(coefficients, sum(x*x for x in point), controls.work)
        rows.append(_compare_jets(poly, jets, point, controls))
    return {"passed": all(row["passed"] for row in rows), "points": rows, "fischer_statistics": stats,
            "reference": "Independent direct Cartesian differentiation of a signed degree-six polynomial"}


def _decimal_fraction(value):
    value = F(value)
    return D(value.numerator)/D(value.denominator)


def _one_variable_derivative(poly, beta, work):
    result = {}
    for power, coefficient in poly.items():
        work.charge(2)
        _put(result, power, -beta*coefficient)
        if power:
            _put(result, power-1, power*coefficient)
    return result


def _profile_fixture(work):
    poly = {(0, 0, 0): F(2), (0, 1, 0): F(-2), (1, 0, 0): F(3, 5), (1, 0, 1): F(7, 6)}
    poly = _add(poly, _multiply(_r2_power(1, work),
                {(0, 0, 0): F(1), (0, 1, 0): F(-2, 3), (1, 0, 1): F(5, 8)}, work), work)
    poly = _add(poly, _r2_power(2, work), work, F(3, 7))
    profiles = {(0, 0): {0: F(2), 1: F(1), 2: F(3, 7)}, (1, -1): {0: F(2), 1: F(2, 3)},
                (1, 1): {0: F(-3, 5)}, (2, 1): {0: F(-7, 18), 1: F(-5, 24)}}
    return poly, profiles


def _profile_control(beta, fixture, controls):
    poly, expected = fixture
    coefficients, stats = controls.exact(lambda limits: fischer_project(poly, 2, **limits))
    expected_coefficients = {(l, m, k): c for (l, m), row in expected.items() for k, c in row.items()}
    controls.work.routine(len(coefficients)+6*sum(max(row)+1 for row in expected.values()),
                          "profile_construction_upper_bound")
    profile = GaussianProfile(coefficients, beta)
    radii = (F(0), F(1, 8), F(1, 3), F(3, 4), F(5, 4))
    controls.work.routine(len(radii)*sum(1+max(row) for row in expected.values()), "profile_sampling_upper_bound")
    sampled = profile.sample(radii, max_values=controls.settings["angular_max_routine_work"])
    tolerance = controls.settings["angular_profile_tolerance"]
    rows = []
    with localcontext() as ctx:
        ctx.prec = 64
        for mode, source_poly in expected.items():
            samples = []
            for index, radius in enumerate(radii):
                s = radius*radius
                controls.work.charge(len(source_poly), "decimal_reference_summands")
                radial_factor = D(1) if mode[0] == 0 else _decimal_fraction(radius)**mode[0]
                target = (radial_factor*(-_decimal_fraction(beta*s)).exp()
                          *sum((_decimal_fraction(c*s**k) for k, c in source_poly.items()), D(0)))
                error = _scaled(sampled[mode][index], target)
                samples.append({"radius": str(radius), "value": sampled[mode][index], "reference_decimal": str(target),
                                "scaled_error": error, "passed": error <= tolerance})
            controls.work.charge(7*len(source_poly)+12*(max(source_poly)+1), "profile_origin_upper_bound")
            data = profile.origin_data(*mode, degree=6, patch_radius=F(1, 2), order=3)
            independent_taylor = []
            for n in range(7):
                value = F(0)
                for k, c in source_poly.items():
                    controls.work.charge(1)
                    if n >= k:
                        value += c*(-beta)**(n-k)/factorial(n-k)
                independent_taylor.append(value)
            bound_rows, derivative_poly = [], dict(source_poly)
            for order in range(4):
                for radius in (F(0), F(1, 6), F(1, 3), F(1, 2)):
                    s = radius*radius
                    controls.work.charge(len(derivative_poly)+7, "decimal_reference_summands")
                    target = (-_decimal_fraction(beta*s)).exp()*sum(
                        (_decimal_fraction(c*s**k) for k, c in derivative_poly.items()), D(0))
                    taylor = sum((_decimal_fraction(c*F(factorial(k), factorial(k-order))*s**(k-order))
                                  for k, c in enumerate(data["coefficients"]) if k >= order), D(0))
                    error = abs(target-taylor)
                    bound = _decimal_fraction(data["error_bounds_s"][order])
                    allowance = D("1e-50")*max(D(1), abs(target))
                    bound_rows.append({"derivative_order": order, "radius": str(radius),
                                       "observed_absolute_remainder_decimal": str(error),
                                       "declared_uniform_bound_decimal": str(bound),
                                       "decimal_roundoff_allowance": str(allowance), "passed": error <= bound+allowance})
                derivative_poly = _one_variable_derivative(derivative_poly, beta, controls.work)
            coefficient_match = tuple(independent_taylor) == data["coefficients"]
            rows.append({"mode": list(mode), "samples": samples, "origin_coefficients_exact": coefficient_match,
                         "origin_coefficients": [str(c) for c in data["coefficients"]],
                         "origin_patch_radius": str(data["patch_radius"]), "origin_bound_samples": bound_rows,
                         "passed": coefficient_match and all(row["passed"] for row in samples+bound_rows)})
    jet_rows = []
    for point in ((F(0),)*3, (F(1, 4), F(-1, 3), F(1, 5))):
        radius = sqrt(float(sum(x*x for x in point)))
        controls.work.charge(4*sum(len(row)+4 for row in expected.values()), "profile_jet_upper_bound")
        jets = {mode: profile.regular_jet(*mode, radius, order=3) for mode in profile.modes}
        jet_rows.append(_compare_jets(poly, jets, point, controls, beta))
    exact_match = dict(coefficients) == expected_coefficients
    return {"passed": exact_match and all(row["passed"] for row in rows+jet_rows), "beta": str(beta),
            "independent_h_profile_coefficients_exact": exact_match, "modes": rows, "cartesian_gaussian_jets": jet_rows,
            "fischer_statistics": stats, "tolerance": tolerance, "uniform_bound_test": "Sampled diagnostics, not a certification"}


def _laurent_normal(terms, work):
    """Independent exact common-degree Cartesian trace normal form."""
    groups = {}
    for (powers, radial, log_power), coefficient in terms.items():
        degree = sum(powers)
        groups.setdefault((degree+radial, log_power, degree % 2), []).append((powers, coefficient))
    normalized = {}
    for key, rows in groups.items():
        degree = max(sum(powers) for powers, _ in rows)
        poly = {}
        for powers, coefficient in rows:
            poly = _add(poly, _multiply({powers: coefficient}, _r2_power((degree-sum(powers))//2, work), work), work)
        if poly:
            normalized[key] = poly
    return normalized


def _laurent_difference(left, right, work):
    work.charge(len(right))
    result = dict(left)
    for key, value in right.items():
        _put(result, key, -value)
    return result


def _laurent_derivative(terms, axis, work):
    result = {}
    for (powers, radial, log_power), c in terms.items():
        work.charge(3)
        if powers[axis]:
            lower = list(powers)
            lower[axis] -= 1
            _put(result, (tuple(lower), radial, log_power), c*powers[axis])
        higher = list(powers)
        higher[axis] += 1
        _put(result, (tuple(higher), radial-2, log_power), c*radial)
        if log_power:
            _put(result, (tuple(higher), radial-2, log_power-1), c*log_power)
    return result


def _field_add(left, right, factor, controls):
    return controls.bounded(len(right.terms), lambda limits: left.add(right, factor, **limits))


def _field_derivative(field, axis, controls):
    branches = 3 if field.has_logs else 2
    return controls.bounded(branches*len(field.terms), lambda limits: field.derivative(axis, **limits))


def _field_laplacian(field, controls):
    bound = (63 if field.has_logs else 30)*len(field.terms)
    return controls.bounded(bound, lambda limits: field.laplacian(**limits))


def _field_scale(field, factor, controls):
    controls.work.charge(len(field.terms), "primitive_scale_products")
    return field.scale(factor)


def _project(field, ell, controls):
    return controls.exact(lambda limits: field.project(ell, **limits))


def _ball_kelvin(component, controls):
    leading = controls.bounded(4, lambda limits: kelvin_exterior_basis(0, 0, 0, component, **limits))
    traced = controls.bounded(168, lambda limits: kelvin_exterior_basis(0, 0, 1, component, **limits))
    fields, targets, rows = [], [], []
    for i in range(3):
        field = _field_add(_field_scale(leading[i], F(1, 3), controls), traced[i], F(1, 5), controls)
        powers = [0, 0, 0]
        powers[i] += 1
        powers[component] += 1
        target = {(tuple(powers), -3, 0): F(-1, 24), (tuple(powers), -5, 0): F(1, 40)}
        if i == component:
            target[((0, 0, 0), -1, 0)] = F(-7, 24)
            target[((0, 0, 0), -3, 0)] = F(-1, 120)
        difference = _laurent_normal(_laurent_difference(field.terms, target, controls.work), controls.work)
        projected, _, stats = _project(field, 2, controls)
        target_projected, _, _ = _project(CartesianLaurent(target), 2, controls)
        derivative_rows = []
        for axis in range(3):
            derivative = _field_derivative(field, axis, controls)
            reference = _laurent_derivative(target, axis, controls.work)
            exact = not _laurent_normal(_laurent_difference(derivative.terms, reference, controls.work), controls.work)
            derivative_rows.append({"axis": axis, "exact_direct_reference_match": exact, "passed": exact})
        rows.append({"component": i, "exact_cartesian_trace_identity": not difference,
                     "exact_projected_identity": projected == target_projected,
                     "exterior_modes": [{"mode": list(mode), "powers": [{"pq": list(key), "coefficient": str(c)}
                                         for key, c in terms.items()]} for mode, terms in projected.items()],
                     "derivatives": derivative_rows, "statistics": stats,
                     "passed": not difference and projected == target_projected and all(row["passed"] for row in derivative_rows)})
        fields.append(field)
        targets.append(target)
    divergence = CartesianLaurent({})
    for i in range(3):
        divergence = _field_add(divergence, _field_derivative(fields[i], i, controls), F(1), controls)
    residual_rows = []
    for i in range(3):
        residual = _field_add(_field_laplacian(fields[i], controls),
                              _field_derivative(divergence, i, controls), F(1, 3), controls)
        modes, _, stats = _project(residual, 4, controls)
        independent = {}
        for axis in range(3):
            second = _laurent_derivative(_laurent_derivative(targets[i], axis, controls.work), axis, controls.work)
            independent = _laurent_difference(independent, {key: -c for key, c in second.items()}, controls.work)
        for j in range(3):
            term = _laurent_derivative(_laurent_derivative(targets[j], j, controls.work), i, controls.work)
            independent = _laurent_difference(independent, {key: -c/3 for key, c in term.items()}, controls.work)
        direct_zero = not _laurent_normal(independent, controls.work)
        residual_rows.append({"component": i, "projected_delta_L_modes": len(modes),
                              "independent_direct_delta_L_zero": direct_zero,
                              "statistics": stats, "passed": not modes and direct_zero})
    return {"passed": all(row["passed"] for row in rows+residual_rows), "source_component": component,
            "radial_moments": ["1/3", "1/5"], "components": rows, "exterior_delta_L": residual_rows,
            "target": "-[7c+n(n.c)]/(24r)+[3n(n.c)-c]/(120r^3); complete trace level retained"}


def _log_algebra(controls):
    terms = {((1, 0, 0), -3, 2): F(1), ((0, 1, 0), -2, 1): F(-3, 2), ((0, 0, 0), -1, 0): F(2)}
    other = {((0, 0, 1), -2, 1): F(2), ((0, 0, 0), -1, 0): F(-1)}
    field, right = CartesianLaurent(terms, 2), CartesianLaurent(other, 2)
    dx_target = {((0, 0, 0), -3, 2): F(1), ((2, 0, 0), -5, 2): F(-3), ((2, 0, 0), -5, 1): F(2),
                 ((1, 1, 0), -4, 1): F(3), ((1, 1, 0), -4, 0): F(-3, 2), ((1, 0, 0), -3, 0): F(-2)}
    product_target = {((1, 0, 1), -5, 3): F(2), ((1, 0, 0), -4, 2): F(-1),
                      ((0, 1, 1), -4, 2): F(-3), ((0, 1, 0), -3, 1): F(3, 2),
                      ((0, 0, 1), -3, 1): F(4), ((0, 0, 0), -2, 0): F(-2)}
    derivative = _field_derivative(field, 0, controls)
    product = controls.bounded(len(terms)*len(other), lambda limits: field.multiply(right, **limits))
    derivative_product = _field_derivative(product, 0, controls)
    first = controls.bounded(len(derivative.terms)*len(other), lambda limits: derivative.multiply(right, **limits))
    d_right = _field_derivative(right, 0, controls)
    second = controls.bounded(len(terms)*len(d_right.terms), lambda limits: field.multiply(d_right, **limits))
    leibniz = _field_add(first, second, F(1), controls)
    cancellation = _field_add(field, field, F(-1), controls)
    log_free = CartesianLaurent({((0, 0, 0), -1, 0): F(1)}, 2)
    reanchored = log_free.at_anchor(3)
    free_cancellation = _field_add(log_free, reanchored, F(-1), controls)
    checks = {"handwritten_log_derivative": dict(derivative.terms) == dx_target,
              "handwritten_same_anchor_product": dict(product.terms) == product_target,
              "exact_product_rule": dict(derivative_product.terms) == dict(leibniz.terms),
              "signed_cancellation": not cancellation.terms, "log_free_anchor_cancellation": not free_cancellation.terms}
    failures = []
    for name, function, words in (
        ("log_anchor_merge", lambda: field.add(CartesianLaurent(other, 3)), ("anchor",)),
        ("log_anchor_change", lambda: field.at_anchor(3), ("anchor",)),
        ("zero_order_copy_term_cap", lambda: field.derivative(0, order=0, max_terms=1), ("max_terms",)),
        ("multiplication_work_cap", lambda: field.multiply(right, max_products=1), ("max_products",)),
        ("unsupported_kelvin_trace", lambda: kelvin_exterior_basis(0, 0, 2, 0), ("trace_level",)),
        ("cartesian_degree_cap", lambda: CartesianLaurent({((33, 0, 0), -34, 0): F(1)}), ("degree",)),
    ):
        controls.work.charge(1, "expected_rejection_upper_bound")
        failures.append({"name": name, **_expected_rejection(function, words)})
    return {"passed": all(checks.values()) and all(row["passed"] for row in failures),
            "exact_identities": checks, "expected_rejections": failures, "anchor": "2"}


def _origin_exterior_distinction(controls):
    bounded = CartesianLaurent({((0, 1, 0), -1, 0): F(1)})
    modes, _, stats = _project(bounded, 1, controls)
    bounded_ok = modes == {(1, -1): {(0, 0): F(-1)}}
    monopole = GaussianProfile({(0, 0, 0): F(1)}, F(1))
    controls.work.charge(4, "profile_sampling_upper_bound")
    samples = monopole.sample((0, F(1, 2)))
    localized_ok = monopole.modes == ((0, 0),) and samples[(0, 0)][0] == 1.0 and samples[(0, 0)][1] > 0
    return {"passed": bounded_ok and localized_ok, "bounded_dipole_p0_field_retained": bounded_ok,
            "localized_source_monopole_retained": localized_ok, "statistics": stats,
            "scope": "Representation only. No deletion of localized h00; inverse source admissibility and two-stage York checks are separate."}


def _run_controls(controls):
    settings = controls.settings
    for ell in range(settings["angular_ell_max"]+1):
        controls.run("exact_harmonics_ell="+str(ell), lambda ell=ell: _harmonics(ell, controls))
    fixture = _nonaxial_fixture(controls.work)
    projections = []
    for shape in settings["angular_grids"]:
        controls.run("angular_grid="+str(shape), lambda shape=shape: _grid_control(shape, controls, fixture, projections))
    def refinement():
        rows = []
        for name in projections[0]:
            error = max(_scaled(a, b) for a, b in zip(projections[0][name], projections[1][name]))
            rows.append({"mixture": name, "normalized_coefficient_grid_change": error,
                         "passed": error <= settings["angular_tolerance"]})
        return {"passed": all(row["passed"] for row in rows), "mixtures": rows,
                "tolerance": settings["angular_tolerance"]}
    controls.run("angular_grid_refinement", refinement)
    small = {(0, 0, 0): F(2), (1, 0, 0): F(-3, 2), (0, 1, 0): F(1, 3),
             (1, 0, 1): F(5, 7), (2, 0, 0): F(2), (0, 2, 0): F(-1), (0, 0, 2): F(4, 3)}
    controls.run("fischer_small_signed", lambda: _fischer_control(small, 2, controls))
    radial = _multiply(_r2_power(2, controls.work), {(0, 1, 0): F(2), (1, 0, 1): F(-3)}, controls.work)
    radial = _add(radial, _r2_power(3, controls.work), controls.work, F(-5, 7))
    controls.run("fischer_high_cartesian_low_angular", lambda: _fischer_control(radial, 2, controls))
    controls.run("fischer_nonaxial_m8", lambda: _fischer_control(fixture[0], 8, controls))
    controls.run("fischer_out_of_band_rejection", lambda: _fischer_rejection(controls))
    controls.run("cartesian_jets_through_three", lambda: _polynomial_jet_control(controls))
    profile_fixture = _profile_fixture(controls.work)
    for beta in (F(1, 2), F(3, 2)):
        controls.run("gaussian_profile_beta="+str(beta), lambda beta=beta: _profile_control(beta, profile_fixture, controls))
    for component in range(3):
        controls.run("kelvin_uniform_ball_component="+str(component), lambda component=component: _ball_kelvin(component, controls))
    controls.run("kelvin_log_product_and_caps", lambda: _log_algebra(controls))
    controls.run("bounded_field_localized_source_distinction", lambda: _origin_exterior_distinction(controls))
    return {"passed": all(row["passed"] for row in controls.rows), "settings": settings,
            "controls": controls.rows, "work": controls.work.statistics(), "accuracy_certified": False,
            "counted_preparation_work": controls.work.total-sum(row["counted_work"] for row in controls.rows),
            "scope": "Independent angular, polynomial jet, profile and exterior Kelvin primitives; no full gravitational or two-stage inverse claim"}


def run_angular_controls(config, progress=None):
    """Run only when invoked by a tracked driver; return complete JSON rows."""
    controls = None
    try:
        settings = _settings(config)
        controls = _Controls(settings, progress)
        return _run_controls(controls)
    except Exception as exc:
        # Also preserve failures between controls (fixture construction or
        # settings), while avoiding a second copy of a row already emitted.
        if controls is None or not controls.exception_emitted:
            row = {"name": "settings" if controls is None else "control_preparation",
                   "passed": False, "exception": type(exc).__name__, "message": str(exc),
                   "work": None if controls is None else controls.work.statistics()}
            if progress is not None:
                progress({"stage": "angular_control_exception", "control_result": row})
        raise
