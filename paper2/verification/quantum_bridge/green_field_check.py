"""Bounded tracked checks for the polynomial-Gaussian first Green inverses.

run_green_checks receives the actual scalar and vector source polynomials
from its caller; it does not construct gravitational sources or launch runs.
All returned objects are JSON-compatible. Import is inert. Actual source
identities use exact Fractions; numeric PDE residuals use separately
computed analytic Hessian components, never finite differences.

Independent heat-integral references use Decimal Gauss-Legendre integration,
not the production recurrence. run.gauss_nodes weights sum to ONE, so an
interval of width L multiplies each normalized weight by L. For z>1 use
u=sqrt(z)*t, composite panels, and a bounded exponentially small tail.
Coarse/fine agreement is required, but is a convergence diagnostic rather
than a rigorous Gauss quadrature remainder bound. The omitted tail has a
separate analytic bound.

Controls check scalar erf normalization, the York value at the origin,
its monopole AND r^-3 tail, and exact Laplacian/longitudinal source
identities. Large-degree stress cases compare only low moments to avoid
expensive high-degree Decimal quadrature. No calculation runs at import.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from functools import lru_cache
from math import ceil, erf, fsum, isfinite, nextafter, pi, sqrt
from time import perf_counter

from run import gauss_nodes
from reduced_gaunt_exact import add
from gaussian_green_fields import (
    scalar_field, york_field, scalar_residual, york_residual,
    heat_normal_form, heat_integrals, evaluate_fields, evaluate_source,
)


DEFAULT_POINTS = ((0.0, 0.0, 0.0), (0.37, -0.41, 0.29), (8.0, -6.0, 4.0))
CONTROL_RATES = (F(1, 2), F(1), F(3, 2))
ONE = {(0, 0, 0): F(1)}


def _settings(config):
    defaults = {
        "green_pde_scaled_tolerance": 1e-8,
        "green_moment_relative_tolerance": 2e-12,
        "green_high_moment_relative_tolerance": 2e-10,
        "green_quadrature_relative_tolerance": 1e-22,
        "green_control_relative_tolerance": 2e-12,
        "green_tail_ratio_tolerance": 1e-6,
    }
    settings = {}
    for key, default in defaults.items():
        value = float(config.get(key, default))
        if not isfinite(value) or value <= 0:
            raise ValueError(key + " must be finite and positive")
        settings[key] = value
    orders = tuple(config.get("green_quadrature_orders", (24, 48)))
    if (len(orders) < 2 or any(type(n) is not int or n < 8 for n in orders)
            or any(a >= b for a, b in zip(orders, orders[1:]))):
        raise ValueError("green_quadrature_orders must contain >=2 increasing integers >=8")
    precision = config.get("green_decimal_precision", 64)
    if type(precision) is not int or precision < 56:
        raise ValueError("green_decimal_precision must be >=56 for run.gauss_nodes")
    mmax = config.get("green_moment_mmax", 16)
    if type(mmax) is not int or mmax < 0:
        raise ValueError("green_moment_mmax must be a nonnegative integer")
    points = tuple(tuple(float(x) for x in point)
                   for point in config.get("green_points", DEFAULT_POINTS))
    if (len(points) < 3 or any(len(p) != 3 or not all(isfinite(x) for x in p) for p in points)
            or (0.0, 0.0, 0.0) not in points):
        raise ValueError("green_points must include the origin and at least two finite 3D points")
    extra_z = tuple(float(z) for z in config.get("green_moment_z_values", ()))
    if any(not isfinite(z) or z < 0 for z in extra_z):
        raise ValueError("green_moment_z_values must be finite and nonnegative")
    settings.update({
        "green_quadrature_orders": list(orders),
        "green_decimal_precision": precision,
        "green_moment_mmax": mmax,
        "green_points": [list(p) for p in points],
        "green_moment_z_values": list(extra_z),
        "quadrature_max_panel_width": 1.5,
        "control_rates": [str(a) for a in CONTROL_RATES],
        "dimensionless_tail_radii": [12.0, 24.0],
    })
    return settings


def _finite(value):
    value = float(value)
    if not isfinite(value):
        raise ArithmeticError("Nonfinite check result")
    return value


def _residual_summary(pair):
    heat, exponential = pair
    samples = []
    for (powers, m), coefficient in list(sorted(heat.items()))[:3]:
        samples.append({"kind": "HEAT", "powers": list(powers), "m": m,
                        "coefficient": str(coefficient)})
    for powers, coefficient in list(sorted(exponential.items()))[:3]:
        samples.append({"kind": "Gaussian", "powers": list(powers),
                        "coefficient": str(coefficient)})
    return {"passed": not heat and not exponential,
            "heat_terms": len(heat), "gaussian_terms": len(exponential),
            "nonzero_samples": samples}


def _hessian_components(field):
    return {(i, j): field.derivative(i).derivative(j)
            for i in range(3) for j in range(i, 3)}


def _scalar_check(name, source, a, settings):
    field = scalar_field(source, a)
    exact = _residual_summary(scalar_residual(field, source))
    hessian = _hessian_components(field)
    keys = tuple(hessian)
    rows = []
    for point in settings["green_points"]:
        values = dict(zip(keys, evaluate_fields(tuple(hessian.values()), point)))
        terms = [values[i, i] for i in range(3)]
        lhs = fsum(terms)
        target = -evaluate_source(source, point, a) / 8
        error = _finite(abs(lhs - target) / (1 + abs(target)))
        rows.append({"point": point, "laplacian": lhs, "target": target,
                     "absolute_error": _finite(abs(lhs - target)), "scaled_error": error,
                     "scaled_term_l1": _finite(fsum(abs(x) for x in terms) / (1 + abs(target)))})
    largest = max(row["scaled_error"] for row in rows)
    return {"name": str(name), "rate": str(a), "source_terms": len(source),
            "field_terms": len(field.terms), "exact": exact, "points": rows,
            "max_m": max(f.max_m for f in hessian.values()),
            "max_scaled_error": largest,
            "passed": exact["passed"] and largest <= settings["green_pde_scaled_tolerance"]}


def _vector_check(name, sources, a, settings):
    sources = tuple(sources)
    fields = york_field(sources, a)
    exact = [_residual_summary(pair) for pair in york_residual(fields, sources)]
    hessians = tuple(_hessian_components(field) for field in fields)
    indexed = [(component, i, j, field) for component, hessian in enumerate(hessians)
               for (i, j), field in hessian.items()]
    rows = []
    for point in settings["green_points"]:
        values = {(component, i, j): value for (component, i, j, _), value in
                  zip(indexed, evaluate_fields(tuple(item[3] for item in indexed), point))}
        residuals, targets, lhs_values, l1_values = [], [], [], []
        for i in range(3):
            terms = [values[i, k, k] for k in range(3)]
            terms += [values[j, min(i, j), max(i, j)] / 3 for j in range(3)]
            lhs = fsum(terms)
            target = evaluate_source(sources[i], point, a)
            residuals.append(_finite(abs(lhs - target) / (1 + abs(target))))
            targets.append(target)
            lhs_values.append(lhs)
            l1_values.append(_finite(fsum(abs(x) for x in terms) / (1 + abs(target))))
        rows.append({"point": point, "york_operator": lhs_values, "target": targets,
                     "component_scaled_errors": residuals, "max_scaled_error": max(residuals),
                     "component_scaled_term_l1": l1_values})
    largest = max(row["max_scaled_error"] for row in rows)
    return {"name": str(name), "rate": str(a), "source_terms": [len(p) for p in sources],
            "field_terms": [len(f.terms) for f in fields], "exact": exact, "points": rows,
            "max_m": max(item[3].max_m for item in indexed),
            "max_scaled_error": largest,
            "passed": all(row["passed"] for row in exact)
                      and largest <= settings["green_pde_scaled_tolerance"]}


@lru_cache(maxsize=64)
def _decimal_u_quadrature(upper_string, mmax, order, precision, gaussian_rate_string):
    """Integrals int_0^upper u^(2m) exp(-rate*u^2) du, all m in one pass."""
    with localcontext() as ctx:
        ctx.prec = precision
        upper, rate = D(upper_string), D(gaussian_rate_string)
        panels = max(1, ceil(float(upper) / 1.5))
        width = upper / panels
        nodes = gauss_nodes(order, precision)
        totals = [D(0)] * (mmax + 1)
        for panel in range(panels):
            midpoint = width * (D(panel) + D("0.5"))
            for node, weight in nodes:
                u = midpoint + width * node / 2
                u2 = u * u
                contribution = width * weight * (-rate * u2).exp()
                for m in range(mmax + 1):
                    totals[m] += contribution
                    contribution *= u2
        return tuple(totals), panels


def _decimal_reference(z_float, mmax, order, precision):
    """Return positive moments, absolute truncation bounds, and panel count."""
    with localcontext() as ctx:
        ctx.prec = precision
        # Reference the exact binary64 input used by the production function.
        z = D.from_float(z_float)
        if z <= 1:
            moments, panels = _decimal_u_quadrature("1", mmax, order, precision, str(z))
            return moments, (D(0),) * (mmax + 1), panels
        root = z.sqrt()
        cap = max(D(12), D(mmax).sqrt() + 8)
        upper = min(root, cap)
        integrals, panels = _decimal_u_quadrature(str(upper), mmax, order, precision, "1")
        factor = 1 / root
        moments, tails = [], []
        u2, boundary_power = upper * upper, upper * (-(upper * upper)).exp()
        for m, integral in enumerate(integrals):
            moments.append(integral * factor)
            # log derivative of u^(2m)e^-u^2 is <= -2(U-m/U) for u>=U.
            bound = boundary_power / (2 * (u2 - m)) if upper < root else D(0)
            tails.append(bound * factor)
            boundary_power *= u2
            factor /= z
        return tuple(moments), tuple(tails), panels


def _moment_check(z, requested_mmax, selected, settings, stress=False, expect_seed_zero=False):
    z = float(z)
    estimates = heat_integrals(z, requested_mmax)
    precision = settings["green_decimal_precision"]
    references = [_decimal_reference(z, max(selected), n, precision)
                  for n in settings["green_quadrature_orders"]]
    relative_errors, convergence, tail_errors, comparisons = [], [], [], []
    with localcontext() as ctx:
        ctx.prec = precision
        for m in selected:
            fine = references[-1][0][m]
            coarse = references[-2][0][m]
            if fine <= 0:
                raise ArithmeticError("Independent integral is not positive")
            estimate = _finite(estimates[m])
            error = _finite(abs(D.from_float(estimate) - fine) / fine)
            delta = _finite(abs(fine - coarse) / fine)
            tail = _finite(references[-1][1][m] / fine)
            relative_errors.append(error)
            convergence.append(delta)
            tail_errors.append(tail)
            comparisons.append({"m": m, "estimate": estimate, "reference_decimal": str(fine),
                                "relative_error": error, "quadrature_relative_difference": delta,
                                "relative_tail_bound": tail})
    tolerance_key = "green_high_moment_relative_tolerance" if stress else "green_moment_relative_tolerance"
    quad_tol = settings["green_quadrature_relative_tolerance"]
    seed_ok = not expect_seed_zero or estimates[-1] == 0.0
    passed = (max(relative_errors) <= settings[tolerance_key]
              and max(convergence) <= quad_tol and max(tail_errors) <= quad_tol
              and all(estimates[m] > 0 for m in selected) and seed_ok)
    return {"name": "heat_moments", "z": z, "requested_mmax": requested_mmax,
            "selected_moments": list(selected), "stress_case": stress,
            "relative_tolerance": settings[tolerance_key],
            "quadrature_orders": settings["green_quadrature_orders"],
            "quadrature_panels": [reference[2] for reference in references],
            "max_relative_error": max(relative_errors),
            "max_quadrature_relative_difference": max(convergence),
            "max_relative_tail_bound": max(tail_errors),
            "highest_requested_moment": _finite(estimates[-1]),
            "expect_highest_underflow_to_zero": expect_seed_zero,
            "highest_underflow_check_passed": seed_ok,
            "comparisons": comparisons, "passed": passed}


def _known_controls(a, settings):
    scalar = scalar_field(ONE, a)
    v = (F(1), F(-2), F(3))
    vector = york_field(tuple({(0, 0, 0): c} for c in v), a)
    af = float(a)
    direction = (2 / 3, -1 / 3, 2 / 3)
    vector_origin = evaluate_fields(vector, (0, 0, 0))
    origin_target = tuple(-11 * float(c) / (24 * af) for c in v)
    origin_error = max(abs(x - y) for x, y in zip(vector_origin, origin_target)) / max(map(abs, origin_target))
    scalar_rows = []
    for dimensionless_radius in (0.0, 0.6, 12.0, 24.0):
        radius = dimensionless_radius / sqrt(af)
        point = tuple(radius * n for n in direction)
        target = (1 / (16 * af) if radius == 0 else
                  sqrt(pi) * erf(sqrt(af) * radius) / (32 * af ** 1.5 * radius))
        value = scalar.evaluate(point)
        scalar_rows.append({"radius": radius, "value": value, "erf_target": target,
                            "relative_error": _finite(abs(value - target) / abs(target))})
    dot = fsum(n * float(c) for n, c in zip(direction, v))
    leading = tuple(-sqrt(pi) * (7 * float(c) + n * dot) / (32 * af ** 1.5)
                    for n, c in zip(direction, v))
    subleading = tuple(-sqrt(pi) * (float(c) - 3 * n * dot) / (64 * af ** 2.5)
                       for n, c in zip(direction, v))
    tail_rows = []
    for dimensionless_radius in (12.0, 24.0):
        radius = dimensionless_radius / sqrt(af)
        point = tuple(radius * n for n in direction)
        r_values = tuple(radius * value for value in evaluate_fields(vector, point))
        corrected = tuple(l + s / radius ** 2 for l, s in zip(leading, subleading))
        relative = max(abs(x - y) for x, y in zip(r_values, corrected)) / max(map(abs, corrected))
        leading_defect = max(abs(x - y) for x, y in zip(r_values, leading)) / max(map(abs, leading))
        scalar_coefficient = radius * scalar.evaluate(point)
        scalar_target = sqrt(pi) / (32 * af ** 1.5)
        tail_rows.append({"radius": radius, "r_times_W": list(r_values),
                          "leading_coefficient": list(leading), "corrected_coefficient": list(corrected),
                          "corrected_relative_error": _finite(relative),
                          "leading_relative_defect": _finite(leading_defect),
                          "scalar_monopole_relative_error": _finite(abs(scalar_coefficient - scalar_target) / scalar_target)})
    if tail_rows[1]["leading_relative_defect"] == 0:
        raise ArithmeticError("Resolved York subleading tail unexpectedly vanished")
    ratio = tail_rows[0]["leading_relative_defect"] / tail_rows[1]["leading_relative_defect"]
    ratio_error = _finite(abs(ratio - 4) / 4)
    errors = [origin_error] + [row["relative_error"] for row in scalar_rows]
    errors += [row[key] for row in tail_rows
               for key in ("corrected_relative_error", "scalar_monopole_relative_error")]
    maximum = _finite(max(errors))
    return {"name": "known_gaussian_controls", "rate": str(a),
            "vector_origin": list(vector_origin), "vector_origin_target": list(origin_target),
            "vector_origin_relative_error": _finite(origin_error),
            "scalar_erf": scalar_rows, "tails": tail_rows,
            "tail_defect_ratio": ratio, "expected_tail_defect_ratio": 4,
            "tail_ratio_relative_error": ratio_error, "max_relative_error": maximum,
            "passed": maximum <= settings["green_control_relative_tolerance"]
                      and ratio_error <= settings["green_tail_ratio_tolerance"]}


def _exact_controls(a):
    laplacian_gaussian = {(0, 0, 0): -6 * a}
    gradient_gaussian = []
    for axis in range(3):
        powers = [0, 0, 0]
        powers[axis] = 2
        laplacian_gaussian[tuple(powers)] = 4 * a * a
        powers[axis] = 1
        gradient_gaussian.append({tuple(powers): -2 * a})
    heat, exponential = heat_normal_form(scalar_field(laplacian_gaussian, a))
    scalar_check = _residual_summary((heat, add(exponential, ONE, F(1, 8))))
    scalar_gradient = scalar_field(ONE, a).gradient()
    longitudinal = york_field(tuple(gradient_gaussian), a)
    vector_checks = [_residual_summary(heat_normal_form(w.add(g, F(6))))
                     for w, g in zip(longitudinal, scalar_gradient)]
    return {"name": "exact_laplacian_and_longitudinal_controls", "rate": str(a),
            "C_Delta_gaussian_plus_gaussian_over_8": scalar_check,
            "York_gradient_plus_6_gradient_C": vector_checks,
            "passed": scalar_check["passed"] and all(row["passed"] for row in vector_checks)}


def run_green_checks(config, scalar_sources, vector_sources, progress=None):
    """Run the bounded Green suite inside the caller's tracked experiment.

    Inputs: scalar_sources=[(name, polynomial, Fraction_rate)], and
    vector_sources=[(name, tuple_of_three_polynomials, Fraction_rate)].
    Both collections must be nonempty for the actual-source coverage gate.
    A progress callback, if present, receives a small JSON-compatible dict
    after each case. No source construction, output files or run management
    is performed here. Invalid configuration raises ValueError; case errors
    are retained as failed JSON rows so the caller can report the evidence.
    """
    start = perf_counter()
    settings = _settings(config)
    scalar_sources, vector_sources = list(scalar_sources), list(vector_sources)
    failures = []

    def perform(stage, label, function):
        try:
            row = function()
        except (ArithmeticError, ValueError, TypeError, AssertionError) as error:
            row = {"name": str(label), "passed": False,
                   "exception": type(error).__name__, "message": str(error)}
        if not row["passed"]:
            failures.append({"stage": stage, "name": str(label),
                             "exception": row.get("exception"), "message": row.get("message")})
        if progress:
            progress({"stage": stage, "name": str(label), "passed": row["passed"]})
        return row

    scalars = [perform("scalar_sources", name,
                       lambda source=source, a=a, name=name: _scalar_check(name, source, a, settings))
               for name, source, a in scalar_sources]
    vectors = [perform("vector_sources", name,
                       lambda sources=sources, a=a, name=name: _vector_check(name, sources, a, settings))
               for name, sources, a in vector_sources]
    exact_controls = [perform("exact_controls", str(a), lambda a=a: _exact_controls(a))
                      for a in CONTROL_RATES]
    controls = [perform("known_controls", str(a), lambda a=a: _known_controls(a, settings))
                for a in CONTROL_RATES]
    mmax = max([settings["green_moment_mmax"]]
               + [row.get("max_m", 0) for row in scalars + vectors])
    threshold = mmax + 5.0
    z_values = list(dict.fromkeys([0.0, 1e-12, 0.125, 2.0, nextafter(threshold, 0.0),
                                  threshold, nextafter(threshold, float("inf")), max(100.0, 2 * threshold),
                                  10000.0] + settings["green_moment_z_values"]))
    moments = [perform("heat_moments", f"z={z}:mmax={mmax}",
                       lambda z=z: _moment_check(z, mmax, tuple(range(mmax + 1)), settings))
               for z in z_values]
    stress_specs = ((599.75, 640, False), (nextafter(600.0, 0.0), 640, False),
                    (600.0, 640, False), (nextafter(600.0, float("inf")), 640, False),
                    (600.25, 640, False), (800.0, 820, True))
    stress = [perform("high_degree_heat_moments", f"z={z}:mmax={n}",
                      lambda z=z, n=n, seed=seed: _moment_check(z, n, (0, 1, 2, 3, 4), settings,
                                                              stress=True, expect_seed_zero=seed))
              for z, n, seed in stress_specs]
    coverage = bool(scalars) and bool(vectors)
    if not coverage:
        failures.append({"stage": "source_coverage", "name": "nonempty_scalar_and_vector_inputs",
                         "exception": None, "message": "Both actual-source collections are required"})

    def largest(rows, key):
        values = [row[key] for row in rows if key in row]
        return max(values) if values else None

    exact_failure_count = sum(not row["exact"]["passed"] for row in scalars if "exact" in row)
    exact_failure_count += sum(not component["passed"] for row in vectors
                               for component in row.get("exact", ()))
    exact_failure_count += sum(not row["passed"] for row in exact_controls)
    all_rows = scalars + vectors + exact_controls + controls + moments + stress
    return {
        "passed": coverage and all(row["passed"] for row in all_rows),
        "scope": "First scalar/York inverse checks on caller-supplied polynomial-Gaussian sources; not H5/H6 integrals",
        "settings": settings, "actual_moment_mmax": mmax,
        "error_definitions": {
            "pde_scaled": "abs(operator-target)/(1+abs(target)), scalar target=-source/8, York target=source",
            "moments": "abs(binary64-Decimal_reference)/Decimal_reference, including small nonzero tail values",
            "quadrature_convergence": "last two orders' relative difference; resolution diagnostic, not rigorous GL remainder",
            "quadrature_tail": "analytic omitted rescaled-Gaussian tail bound divided by positive reference",
            "stress_tolerance": "Separate a-priori binary64 allowance for hundreds of logarithmic recurrence steps",
            "tail_control": "r*York field includes both known constant and r^-2 coefficients; leading defect ratio must approach four",
        },
        "source_counts": {"scalar": len(scalars), "vector": len(vectors)},
        "summary": {
            "exact_failure_count": exact_failure_count,
            "exception_count": sum("exception" in row for row in all_rows),
            "max_pde_scaled_error": largest(scalars + vectors, "max_scaled_error"),
            "max_moment_relative_error": largest(moments, "max_relative_error"),
            "max_high_moment_relative_error": largest(stress, "max_relative_error"),
            "max_quadrature_relative_difference": largest(moments + stress, "max_quadrature_relative_difference"),
            "max_quadrature_relative_tail_bound": largest(moments + stress, "max_relative_tail_bound"),
            "max_control_relative_error": largest(controls, "max_relative_error"),
            "max_tail_ratio_relative_error": largest(controls, "tail_ratio_relative_error"),
        },
        "scalar_sources": scalars, "vector_sources": vectors,
        "exact_controls": exact_controls, "known_controls": controls,
        "heat_moments": moments, "high_degree_heat_moments": stress,
        "failures": failures, "elapsed_seconds": perf_counter() - start,
    }
