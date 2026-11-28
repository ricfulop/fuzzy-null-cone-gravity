"""Closed-form controls for the independent radial inverse, inert on import.

All numerical evaluations belong to a tracked run. These controls test an
inverse implementation, not a gravitational Hamiltonian or scattering limit.
The reference formulas integrate the piecewise source analytically; they do
not use the numerical solver's exterior recurrence or ODE-derived jets.
"""

from fractions import Fraction as F
from math import fsum, isfinite, log

from spherical_green_reference import (
    ExteriorSource, OriginSeries, RadialGrid, solve_scalar_modes,
)


def _finite_positive(config, key, default):
    value = float(config.get(key, default))
    if not isfinite(value) or value <= 0:
        raise ValueError(key + " must be positive and finite")
    return value


def _power_derivative(power, order):
    value = 1
    for j in range(order):
        value *= power-j
    return value


def _power_log_derivative(power, radius, order):
    # d^j(r^a log r) = r^(a-j)(A_j log r+B_j).
    a, b = 1, 0
    for j in range(order):
        a, b = (power-j)*a, a+(power-j)*b
    return radius**(power-order)*(a*log(radius)+b)


def exact_scalar_jet(ell, p, radius, order=3, side="left", inverse="C"):
    """Analytic jet for f_l=r^l inside 1, r^-p outside, in H convention."""
    scale = 1 if inverse == "C" else 8
    if radius < 1 or (radius == 1 and side == "left"):
        a = (F(1, 2)+F(1, p+ell-2))/ (8*(2*ell+1))
        b = -F(1, 16*(2*ell+3))
        result = []
        for j in range(order+1):
            terms = []
            for degree, coefficient in ((ell, a), (ell+2, b)):
                if degree >= j:
                    terms.append(float(coefficient)*_power_derivative(degree, j)*radius**(degree-j))
            result.append(scale*fsum(terms))
        return tuple(result)
    divisor = 8*(2*ell+1)
    k = ell+3-p
    if k:
        terms = ((-ell-1, F(1, 2*ell+3)-F(1, k)),
                 (2-p, F(1, k)+F(1, p+ell-2)))
        return tuple(scale*fsum(float(c)*_power_derivative(a, j)*radius**(a-j)
                                for a, c in terms)/divisor for j in range(order+1))
    return tuple(scale*(float(F(1, 2*ell+3)+F(1, p+ell-2))*
                        _power_derivative(-ell-1, j)*radius**(-ell-1-j)
                        +_power_log_derivative(-ell-1, radius, j))/divisor
                 for j in range(order+1))


def _reference_core(ell, p, radius):
    k = ell+3-p
    return 1/(2*ell+3)+(log(radius) if k == 0 else (radius**k-1)/k)


def _solve(case, radius, intervals, patch, *, omit_exterior=False):
    ell, m, p, inverse = case
    grid = RadialGrid(radius, intervals, breakpoints=(1,))
    source = tuple(r**ell if r <= 1 else r**(-p) for r in grid.radii)
    key = ell, m
    exterior = ExteriorSource(radius, {key: {} if omit_exterior else {(p, 0): 1}}, {key: 0})
    solution = solve_scalar_modes(grid, {key: source}, inverse=inverse,
        exterior=exterior, origin_series={key: OriginSeries(patch, (1,), (0, 0, 0))})
    return grid, solution


def _case_row(case, radius, intervals, patch, settings):
    ell, m, p, inverse = case
    grid, solution = _solve(case, radius, intervals, patch)
    diagnostic = solution.diagnostics()["rows"][0]
    probes = [(0.0, None), (.25, None), (patch, "left"), (patch, "right"),
              (.75, None), (1.0, "left"), (1.0, "right"), (1.5, None),
              (radius, "left"), (radius, "right"), (1.5*radius, None)]
    points, value_error, jet_error, regular_error = [], 0.0, 0.0, 0.0
    reference_a0 = float((F(1, 2)+F(1, p+ell-2))/(8*(2*ell+1)))*(1 if inverse == "C" else 8)
    reference_b = -(1 if inverse == "C" else 8)/(16*(2*ell+3))
    for r, side in probes:
        actual = solution.radial_jet(ell, m, r, order=3, side=side)
        expected = exact_scalar_jet(ell, p, r, order=3, side=side or "left", inverse=inverse)
        errors = [abs(a-b)/max(1.0, abs(b)) for a, b in zip(actual, expected)]
        value_error = max(value_error, errors[0])
        jet_error = max(jet_error, *errors)
        point = {"radius": r, "side": side, "actual": actual, "reference": expected,
                 "scaled_errors": errors}
        if r < 1 or (r == 1 and side == "left"):
            regular = solution.regular_jet(ell, m, r, order=3, side=side)
            regular_target = (reference_a0+reference_b*r*r, reference_b, 0.0, 0.0)
            regular_errors = [abs(a-b)/max(1.0, abs(b)) for a, b in zip(regular, regular_target)]
            regular_error = max(regular_error, *regular_errors)
            point.update(regular_jet=regular, regular_reference=regular_target,
                         regular_scaled_errors=regular_errors)
        points.append(point)
    alpha = .125 if inverse == "C" else 1.0
    a0_target = reference_a0
    a0_relative_error = abs(diagnostic["origin_a0"]-a0_target)/abs(a0_target)
    # Use raw unknowns and a distinct fourth-order stencil. Do not substitute
    # radial_jet's ODE-derived second/third derivatives in this residual.
    values = solution.values[ell, m]
    residual_rows = []
    for fraction in (.16, .22, .34, .43, .57, .71, .83):
        index = round(fraction*intervals)
        r = grid.radii[index]
        if index < 2 or index > intervals-2 or r < patch or abs(r-1) <= 2*grid.h:
            continue
        first = (values[index-2]-8*values[index-1]+8*values[index+1]-values[index+2])/(12*grid.h)
        second = (-values[index-2]+16*values[index-1]-30*values[index]
                  +16*values[index+1]-values[index+2])/(12*grid.h**2)
        source = r**ell if r < 1 else r**(-p)
        residual = second+2*first/r-ell*(ell+1)*values[index]/r**2+alpha*source
        residual_rows.append({"radius": r, "residual": residual,
                              "scaled_error": abs(residual)/max(1, abs(alpha*source))})
    if not residual_rows:
        raise ValueError("No admissible independent PDE comparison points")
    pde_error = max(row["scaled_error"] for row in residual_rows)
    core = _reference_core(ell, p, radius)
    core_error = abs(diagnostic["core_moment_from_boundary"]-core)/max(1, abs(core))
    core_quadrature_error = abs(diagnostic["core_moment_quadrature"]-core)/max(1, abs(core))
    expected_outer = radius**(2-ell-p)/(p+ell-2)
    outer_error = abs(diagnostic["exterior_outer_moment"]-expected_outer)/max(1, abs(expected_outer))
    patch_error = max(abs(v) for v in diagnostic["origin_patch_derivative_matching_differences"])
    final_gates = {
        "value": value_error <= settings["value_tolerance"],
        "jet": jet_error <= settings["jet_tolerance"],
        "regular_jet": regular_error <= settings["jet_tolerance"],
        "origin_a0": a0_relative_error <= settings["origin_relative_tolerance"],
        "pde": pde_error <= settings["pde_tolerance"],
        "core_moment": core_error <= settings["core_tolerance"],
        "core_quadrature": core_quadrature_error <= settings["core_quadrature_tolerance"],
        "outer_moment": outer_error <= 1e-12,
        "origin_source": diagnostic["origin_source_samples_within_declared_bound"],
        "patch_derivatives": patch_error <= settings["jet_tolerance"],
    }
    return {"case": {"ell": ell, "m": m, "p": p, "inverse": inverse},
            "radius": radius, "intervals": intervals, "h": grid.h, "origin_patch": patch,
            "points": points, "independent_pde": residual_rows, "diagnostics": diagnostic,
            "value_error": value_error, "jet_error": jet_error, "regular_jet_error": regular_error,
            "origin_a0_target": a0_target, "origin_a0_relative_error": a0_relative_error,
            "pde_error": pde_error, "core_error": core_error,
            "core_quadrature_error": core_quadrature_error, "outer_error": outer_error,
            "patch_derivative_error": patch_error, "final_resolution_gates": final_gates,
            "within_final_tolerances": all(final_gates.values())}


def run_scalar_controls(config, progress=None):
    """Run fixed-patch radial and independent fixed-h outer-radius slices."""
    settings = {name: _finite_positive(config, "scalar_"+name, default) for name, default in (
        ("value_tolerance", 1e-4), ("jet_tolerance", .002),
        ("origin_relative_tolerance", .003), ("pde_tolerance", .001),
        ("core_tolerance", .003), ("core_quadrature_tolerance", 1e-7),
        ("minimum_error_reduction", 2.5))}
    radius = _finite_positive(config, "scalar_radius", 4)
    patch = _finite_positive(config, "scalar_origin_radius", .5)
    intervals = tuple(config.get("scalar_intervals", (256, 512, 1024)))
    radii = tuple(float(r) for r in config.get("scalar_exterior_radii", (2, 4, 6)))
    if (len(intervals) < 3 or any(type(n) is not int or n <= 0 for n in intervals)
            or any(b != 2*a for a, b in zip(intervals, intervals[1:]))
            or len(radii) < 3 or len(set(radii)) != len(radii)
            or any(not isfinite(r) or r <= 1 for r in radii)
            or radius <= 1 or patch != .5):
        raise ValueError("Declare three doubled radial grids, three distinct R>1, and fixed origin patch .5")
    cases = ((0, 0, 4, "C"), (1, -1, 4, "C"), (2, 1, 4, "C"),
             (1, 1, 2, "C"), (8, -3, 4, "C"), (10, 7, 4, "C"), (2, -1, 4, "P"))
    max_unknowns = config.get("scalar_max_unknowns", 200_000)
    h = radius/intervals[-1]
    radius_counts = tuple(round(r/h) for r in radii)
    planned = len(cases)*(sum(n+1 for n in intervals)+sum(n+1 for n in radius_counts))+intervals[-1]+1
    if type(max_unknowns) is not int or max_unknowns < 1 or planned > max_unknowns:
        raise ValueError("Planned scalar unknowns exceed scalar_max_unknowns")
    if any(abs(r-n*h) > 1e-12 for r, n in zip(radii, radius_counts)):
        raise ValueError("Exterior-radius slice must preserve the exact requested h")
    rows = []

    def record(row):
        rows.append(row)
        if progress is not None:
            progress({"stage": "scalar_control_complete", "control": row["name"],
                      "passed": row["passed"], "control_result": row})

    current = None
    try:
        for case in cases:
            current = "scalar_%s_%s_p%s_%s" % case
            levels = []
            for n in intervals:
                row = _case_row(case, radius, n, patch, settings)
                levels.append(row)
                if progress is not None:
                    progress({"stage": "scalar_level_complete", "control": current,
                              "slice": "radial", "level_result": row})
            reductions = {}
            for field in ("value_error", "origin_a0_relative_error"):
                reductions[field] = [None if b[field] < 1e-10
                    else a[field]/b[field] for a, b in zip(levels, levels[1:])]
            reduction_passed = all(r is None or r >= settings["minimum_error_reduction"]
                                   for values in reductions.values() for r in values)
            exterior_levels = []
            for r, n in zip(radii, radius_counts):
                row = _case_row(case, r, n, patch, settings)
                exterior_levels.append(row)
                if progress is not None:
                    progress({"stage": "scalar_level_complete", "control": current,
                              "slice": "exterior_radius", "level_result": row})
            structural_passed = all(row["final_resolution_gates"][key]
                for row in levels+exterior_levels for key in ("origin_source", "outer_moment"))
            record({"name": current, "radial_levels": levels, "error_reductions": reductions,
                    "exterior_radius_levels": exterior_levels,
                    "radial_reduction_passed": reduction_passed,
                    "all_level_structural_passed": structural_passed,
                    "passed": levels[-1]["within_final_tolerances"] and reduction_passed and structural_passed
                              and all(row["within_final_tolerances"] for row in exterior_levels)})
        current = "exterior_source_omission_ablation"
        case = (1, 1, 2, "C")
        _, wrong = _solve(case, radius, intervals[-1], patch, omit_exterior=True)
        actual = wrong.radial_jet(1, 1, .75, order=0)[0]
        expected = exact_scalar_jet(1, 2, .75, order=0)[0]
        discrepancy = abs(actual-expected)
        record({"name": current, "actual": actual, "reference": expected,
                "absolute_discrepancy": discrepancy, "passed": discrepancy > 10*settings["value_tolerance"],
                "scope": "Omitting a legitimate exterior dipole must visibly fail"})
        current = "domain_and_join_guards"
        guards = []
        for ell, p, q in ((0, 2, 0), (2, 1, 0), (1, 2, 1)):
            try:
                ExteriorSource(4, {(ell, 0): {(p, q): 1}}, {(ell, 0): 0})
            except (ValueError, ArithmeticError):
                guards.append(True)
            else:
                guards.append(False)
        try:
            wrong.radial_jet(1, 1, 1, order=1)
        except ValueError:
            guards.append(True)
        else:
            guards.append(False)
        record({"name": current, "rejections": guards, "passed": all(guards)})
    except BaseException as exc:
        record({"name": current or "preparation", "passed": False,
                "exception": type(exc).__name__+": "+str(exc)})
        raise
    return {"passed": bool(rows) and all(row["passed"] for row in rows), "rows": rows,
            "settings": dict(settings, radius=radius, intervals=intervals, origin_radius=patch,
                             exterior_radii=radii, radius_slice_h=h, max_unknowns=max_unknowns),
            "planned_unknowns": planned, "accuracy_certified": False,
            "limitations": ["Fixed finite harmonic fixtures only; no gravitational functional validation",
                "ODE-derived jets are compared to analytic solutions, not reused as independent PDE residuals",
                "Fourth-order raw-grid residual and refinement remain numerical diagnostics",
                "Artificial origin patch is fixed in physical radius and its a0 and derivatives are tested"]}
