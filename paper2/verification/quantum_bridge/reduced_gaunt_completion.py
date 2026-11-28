"""Constraint reconstruction for the explicitly revised finite-free-data model.

The free angular rank is N=3 throughout. Auxiliary J, angular quadrature and
radial grid refine dependent geometry; they are NOT a physical large-N limit.
This calculation supplies initial-data and Hamiltonian ingredients, not an
interacting quantum amplitude or a dynamical corner-state count.
"""

from functools import lru_cache
from math import comb, erf, erfc, exp, factorial, fsum, pi, sqrt
from fractions import Fraction as F
import json
import os
import platform
import sys
import time

from run import ROOT, gauss_nodes, report
from reduced_gaunt_exact import exact_mode_basis, axisymmetric_source


@lru_cache(None)
def quadrature(order):
    return tuple((float(x), float(w)) for x, w in gauss_nodes(order, 60))


def legendre_values(maximum, x):
    values = [1.0]
    if maximum:
        values.append(x)
    for ell in range(2, maximum + 1):
        values.append(((2 * ell - 1) * x * values[-1] - (ell - 1) * values[-2]) / ell)
    return values


def lower_gaussian_moment(n, radius):
    if radius == 0:
        return 0.0
    if radius <= 2:
        term = 1.0
        summands = []
        for k in range(120):
            value = term / (n + 2 * k + 1)
            summands.append(value)
            if k > 10 and abs(value) < 1e-18 * abs(fsum(summands)):
                break
            term *= -radius * radius / (k + 1)
        else:
            raise ArithmeticError("Lower Gaussian moment did not converge")
        return radius ** (n + 1) * fsum(summands)
    exponential = exp(-radius * radius)
    if n % 2 == 0:
        value, start = sqrt(pi) * erf(radius) / 2, 2
    else:
        value, start = (1 - exponential) / 2, 3
    for power in range(start, n + 1, 2):
        value = ((power - 1) * value - radius ** (power - 1) * exponential) / 2
    return value


def upper_gaussian_moment(n, radius):
    exponential = exp(-radius * radius)
    if n % 2 == 0:
        value, start = sqrt(pi) * erfc(radius) / 2, 2
    else:
        value, start = exponential / 2, 3
    for power in range(start, n + 1, 2):
        value = ((power - 1) * value + radius ** (power - 1) * exponential) / 2
    return value


def source_radial(ell, radius, source):
    return exp(-radius * radius) * fsum(
        float(value) * radius ** power for power, value in source.get(ell, {}).items())


def source_point(radius, x):
    """Independent positive square expression in an orthonormal spherical frame."""
    r2 = radius * radius
    p2 = (3 * x * x - 1) / 2
    weight = 1 - x * x
    return exp(-r2) * (54 * p2 * p2 + 18 * (3 - r2) ** 2 * x * x * weight +
                      9 / 8 * (r2 * r2 - 8 * r2 + 6) ** 2 * weight * weight)


def linear_green(ell, radius, source):
    radial = source[ell]
    if radius == 0:
        if ell:
            return 0.0
        return fsum(float(value) * factorial(power // 2) / 16
                    for power, value in radial.items())
    interior = fsum(float(value) * lower_gaussian_moment(power + ell + 2, radius)
                    for power, value in radial.items())
    exterior = fsum(float(value) * upper_gaussian_moment(power + 1 - ell, radius)
                    for power, value in radial.items())
    return (interior / radius ** (ell + 1) + radius ** ell * exterior) / (8 * (2 * ell + 1))


def radial_majorant_tails(radius, source):
    potential = fsum(abs(float(value)) * upper_gaussian_moment(power + 1, radius) / 8
                     for radial in source.values() for power, value in radial.items())
    mass = 4 * pi * fsum(abs(float(value)) * upper_gaussian_moment(power + 2, radius)
                        for radial in source.values() for power, value in radial.items())
    return {"uniform_linear_potential_tail_bound": potential,
            "source_integral_tail_bound": mass}


class RadialPoisson:
    """Second-order radial Laplacian, regular origin and exterior Robin boundary.

The last Robin row uses a second-order one-sided derivative. Its second
subdiagonal is eliminated with the preceding interior row before a Thomas
factorization. This is separate from the analytic Gaussian Green function.
"""

    def __init__(self, ell, intervals, radius):
        self.ell, self.intervals, self.radius = ell, intervals, radius
        self.h = radius / intervals
        h = self.h
        self.lower = [0.0] * (intervals + 1)
        diagonal = [0.0] * (intervals + 1)
        upper = [0.0] * (intervals + 1)
        if ell == 0:
            diagonal[0], upper[0] = -6 / h ** 2, 6 / h ** 2
        else:
            diagonal[0] = 1.0
        for i in range(1, intervals):
            r = i * h
            self.lower[i] = 1 / h ** 2 - 1 / (r * h)
            diagonal[i] = -2 / h ** 2 - ell * (ell + 1) / r ** 2
            upper[i] = 1 / h ** 2 + 1 / (r * h)
        self.boundary_elimination = (1 / (2 * h)) / self.lower[-2]
        self.lower[-1] = -2 / h - self.boundary_elimination * diagonal[-2]
        diagonal[-1] = 3 / (2 * h) + (ell + 1) / radius - self.boundary_elimination * upper[-2]
        self.inverse_pivots, self.modified_upper = [], []
        for i in range(intervals + 1):
            pivot = diagonal[i] - (self.lower[i] * self.modified_upper[i - 1] if i else 0)
            if abs(pivot) < 1e-20:
                raise ArithmeticError("Singular radial Poisson solve")
            self.inverse_pivots.append(1 / pivot)
            self.modified_upper.append(upper[i] / pivot)

    def solve(self, source, epsilon_squared=1.0):
        rhs = [-epsilon_squared * value / 8 for value in source]
        if self.ell:
            rhs[0] = 0.0
        rhs[-1] = -self.boundary_elimination * rhs[-2]
        forward = []
        for i, value in enumerate(rhs):
            forward.append((value - (self.lower[i] * forward[i - 1] if i else 0)) *
                           self.inverse_pivots[i])
        result = list(forward)
        for i in range(self.intervals - 1, -1, -1):
            result[i] -= self.modified_upper[i] * result[i + 1]
        return result


def linear_grid_check(intervals, radius, source):
    h = radius / intervals
    errors = {}
    for ell in source:
        forcing = [source_radial(ell, i * h, source) for i in range(intervals + 1)]
        numerical = RadialPoisson(ell, intervals, radius).solve(forcing)
        target = [linear_green(ell, i * h, source) for i in range(intervals + 1)]
        errors[str(ell)] = max(abs(a - b) for a, b in zip(numerical, target))
    return {"radial_intervals": intervals, "h": h, "errors_by_ell": errors,
            "maximum_linear_green_error": max(errors.values())}


def energy_integrals(radial_order, angular_order, radius, source):
    second, third = [], []
    angles = [(x, w, legendre_values(4, x)) for x, w in quadrature(angular_order)]
    for node, weight in quadrature(radial_order):
        r = radius * (node + 1) / 2
        potentials = {ell: linear_green(ell, r, source) for ell in source}
        q2 = fsum(source_radial(ell, r, source) * potentials[ell] / (2 * ell + 1)
                  for ell in source)
        q3 = fsum(w * source_point(r, x) *
                  fsum(potentials[ell] * values[ell] for ell in source) ** 2
                  for x, w, values in angles)
        second.append(radius * weight * 4 * pi * r * r * q2)
        third.append(radius * weight * 4 * pi * r * r * q3)
    return {"radial_order": radial_order, "angular_order": angular_order,
            "source_times_u2": fsum(second), "source_times_u2_squared": fsum(third)}


def simpson(values, h):
    intervals = len(values) - 1
    assert intervals % 2 == 0
    return h / 3 * (values[0] + values[-1] +
                    4 * fsum(values[1:-1:2]) + 2 * fsum(values[2:-1:2]))


def within_free_area_budget(electric, magnetic, budget):
    """Exact comparison sqrt(11/2)*ne + sqrt(13/2)*nm <= rational B."""
    budget = F(budget)
    remainder = budget * budget - F(11, 2) * electric ** 2 - F(13, 2) * magnetic ** 2
    return budget >= 0 and remainder >= 0 and remainder ** 2 >= 143 * electric ** 2 * magnetic ** 2


def free_corner_count(budget):
    """Leading area spectrum only, with vacuum subtraction explicitly fixed."""
    budget = F(str(budget))
    levels = []
    dimension = 0
    for electric in range(int(budget) + 1):
        for magnetic in range(int(budget) + 1):
            if within_free_area_budget(electric, magnetic, budget):
                degeneracy = comb(electric + 4, 4) * comb(magnetic + 4, 4)
                dimension += degeneracy
                levels.append({"electric_quanta": electric, "magnetic_quanta": magnetic,
                               "degeneracy": degeneracy})

    @lru_cache(None)
    def count_individual_modes(index, electric, magnetic):
        if index == 10:
            return 1
        total = 0
        for occupation in range(int(budget) + 1):
            ne = electric + (occupation if index < 5 else 0)
            nm = magnetic + (occupation if index >= 5 else 0)
            if not within_free_area_budget(ne, nm, budget):
                break
            total += count_individual_modes(index + 1, ne, nm)
        return total

    independent = count_individual_modes(0, 0, 0)
    assert dimension == independent
    return {
        "dimensionless_area_budget_B": str(budget), "dimension": dimension,
        "independent_mode_recursion_count": independent, "levels": levels,
        "budget_definition": "B=(A_max-A_0)*a/(8*pi*G*hbar*R), at quadratic order and a far sphere",
        "scope": "bulk-induced uniform-area matching sector of ten free modes; not a full corner Hilbert space or an entropy law",
        "renormalization": "normal-ordered quadratic ADM energy; Minkowski corner area set to A_0",
    }


def fourth_order_laplacian(values, ell, h, i):
    r = i * h
    first = (values[i - 2] - 8 * values[i - 1] + 8 * values[i + 1] - values[i + 2]) / (12 * h)
    second = (-values[i + 2] + 16 * values[i + 1] - 30 * values[i] +
              16 * values[i - 1] - values[i - 2]) / (12 * h * h)
    return second + 2 * first / r - ell * (ell + 1) * values[i] / (r * r)


def nonlinear_case(case, source, potential_bound, energy, config):
    intervals, auxiliary = case["radial_intervals"], case["auxiliary_J"]
    angular_order = case["angular_nodes"]
    radius = config["outer_radius"]
    h = radius / intervals
    epsilon_squared = case["rho"] / float(potential_bound)
    contraction_bound = 7 * epsilon_squared * float(potential_bound)
    assert contraction_bound < 1
    radii = [i * h for i in range(intervals + 1)]
    ells = list(range(0, auxiliary + 1, 2))
    nodes = quadrature(angular_order)
    legendre = [[legendre_values(auxiliary, x)[ell] for ell in ells] for x, _ in nodes]
    weights = [[(2 * ell + 1) * w * values[j] for (x, w), values in zip(nodes, legendre)]
               for j, ell in enumerate(ells)]
    raw_source = [[source_point(r, x) for x, _ in nodes] for r in radii]
    solvers = [RadialPoisson(ell, intervals, radius) for ell in ells]

    def project_source(fields):
        coefficients = [[0.0] * (intervals + 1) for _ in ells]
        minimum, maximum = float("inf"), 0.0
        for i in range(intervals + 1):
            psi = [1 + fsum(fields[j][i] * values[j] for j in range(len(ells)))
                   for values in legendre]
            minimum, maximum = min(minimum, min(psi)), max(maximum, max(psi))
            if min(psi) <= 0:
                raise ArithmeticError("Conformal factor left its positive domain")
            values = [density * conformal ** -7 for density, conformal in zip(raw_source[i], psi)]
            for j in range(len(ells)):
                coefficients[j][i] = fsum(w * value for w, value in zip(weights[j], values))
        return coefficients, minimum, maximum

    fields = [[0.0] * (intervals + 1) for _ in ells]
    history = []
    for iteration in range(1, config["maximum_iterations"] + 1):
        coefficients, _, _ = project_source(fields)
        updated = [solver.solve(coeff, epsilon_squared) for solver, coeff in zip(solvers, coefficients)]
        # |P_l|<=1 converts this coefficient norm to a pointwise upper bound.
        change = max(fsum(abs(updated[j][i] - fields[j][i]) for j in range(len(ells)))
                     for i in range(intervals + 1))
        fields = updated
        history.append(change)
        report(row="iteration_history", m={"case": case["name"], "iteration": iteration,
                                          "fixed_point_change": change})
        if change < config["iteration_tolerance"]:
            break
    else:
        raise ArithmeticError("Conformal fixed-point iteration did not converge")
    coefficients, minimum, maximum = project_source(fields)
    assert minimum >= 1 - 1e-10
    assert maximum <= 1 + case["rho"] + 1e-5

    validation_x = [x for x, _ in quadrature(config["validation_angular_nodes"])] + [-1.0, 0.0, 1.0]
    validation_p = [[legendre_values(auxiliary, x)[ell] for ell in ells] for x in validation_x]
    angular_defect = 0.0
    full_constraint_residual = 0.0
    for i, r in enumerate(radii):
        if 2 <= i <= intervals - 2:
            laplacians = [fourth_order_laplacian(field, ell, h, i)
                          for field, ell in zip(fields, ells)]
        else:
            laplacians = None
        for x, values in zip(validation_x, validation_p):
            psi = 1 + fsum(fields[j][i] * values[j] for j in range(len(ells)))
            unprojected = source_point(r, x) * psi ** -7
            projected = fsum(coefficients[j][i] * values[j] for j in range(len(ells)))
            angular_defect = max(angular_defect, epsilon_squared * abs(unprojected - projected) / 8)
            if laplacians is not None:
                laplacian = fsum(value * basis for value, basis in zip(laplacians, values))
                full_constraint_residual = max(full_constraint_residual,
                                               abs(laplacian + epsilon_squared * unprojected / 8))

    mass_volume = epsilon_squared * 4 * pi * simpson(
        [r * r * value for r, value in zip(radii, coefficients[0])], h)
    mass_boundary = 32 * pi * radius * fields[0][-1]
    mass2 = epsilon_squared * energy["quadratic"]
    mass4 = mass2 + epsilon_squared ** 2 * energy["quartic"]
    mass6 = mass4 + epsilon_squared ** 3 * energy["sixth"]
    assert mass_volume > 0
    assert mass_volume < mass2
    return {
        **case, "angular_rank_N": 3, "epsilon": sqrt(epsilon_squared),
        "epsilon_squared": epsilon_squared, "h": h,
        "analytic_contraction_bound": contraction_bound,
        "iterations": len(history), "last_fixed_point_change": history[-1],
        "minimum_psi_on_solver_nodes": minimum, "maximum_psi_on_solver_nodes": maximum,
        "sampled_angular_constraint_defect": angular_defect,
        "sampled_fourth_order_constraint_residual": full_constraint_residual,
        "ADM_mass_volume": mass_volume, "ADM_mass_boundary": mass_boundary,
        "boundary_volume_difference": abs(mass_boundary - mass_volume),
        "mass_quadratic": mass2, "mass_through_quartic": mass4, "mass_through_sixth": mass6,
        "scope": "Fixed physical free data, with auxiliary geometry and radial resolution varied.",
    }


def convergence_summary(rows, linear, config):
    linear_ratios = [linear[i]["maximum_linear_green_error"] / linear[i + 1]["maximum_linear_green_error"]
                     for i in range(len(linear) - 1)]
    assert all(3.5 < ratio < 4.5 for ratio in linear_ratios), linear_ratios
    by_name = {row["name"]: row for row in rows}
    refinements = []
    for rho in config["amplitude_rhos"]:
        group = sorted((row for row in rows if row["rho"] == rho and row["auxiliary_J"] == 12
                        and row["angular_nodes"] == 32), key=lambda row: row["radial_intervals"])
        assert len(group) == 3
        coarse, medium, fine = group
        richardson_coarse = (4 * medium["ADM_mass_volume"] - coarse["ADM_mass_volume"]) / 3
        richardson_fine = (4 * fine["ADM_mass_volume"] - medium["ADM_mass_volume"]) / 3
        difference = abs(richardson_fine - richardson_coarse)
        error4 = abs(richardson_fine - fine["mass_through_quartic"])
        error6 = abs(richardson_fine - fine["mass_through_sixth"])
        assert error6 < error4, (rho, error4, error6)
        residual_ratios = [
            coarse["sampled_fourth_order_constraint_residual"] /
            medium["sampled_fourth_order_constraint_residual"],
            medium["sampled_fourth_order_constraint_residual"] /
            fine["sampled_fourth_order_constraint_residual"],
        ]
        assert all(3 < ratio < 5 for ratio in residual_ratios), (rho, residual_ratios)
        assert (fine["boundary_volume_difference"] < medium["boundary_volume_difference"] <
                coarse["boundary_volume_difference"])
        refinements.append({
            "rho": rho, "richardson_mass": richardson_fine,
            "difference_between_two_extrapolations": difference,
            "quartic_series_difference": error4, "sixth_order_series_difference": error6,
            "fourth_order_residual_refinement_ratios": residual_ratios,
            "caution": "Extrapolation differences are diagnostics, not rigorous error bounds.",
        })
    largest = max(config["amplitude_rhos"])
    angular = sorted((row for row in rows if row["rho"] == largest and
                      row["radial_intervals"] == max(config["radial_intervals"]) and
                      row["angular_nodes"] == 32), key=lambda row: row["auxiliary_J"])
    assert [row["auxiliary_J"] for row in angular] == [2, 4, 8, 12]
    defects = [row["sampled_angular_constraint_defect"] for row in angular]
    assert defects[0] > 0
    assert defects[1] < defects[0] / 10
    assert defects[2] < defects[1] / 10
    assert defects[3] <= max(defects[2], 1e-13)
    fine = angular[-1]
    checked = by_name["angular_quadrature_check"]
    quadrature_mass_difference = abs(fine["ADM_mass_volume"] - checked["ADM_mass_volume"])
    assert quadrature_mass_difference < config["angular_quadrature_mass_tolerance"]
    return {"linear_green_refinement_ratios": linear_ratios,
            "amplitude_series_and_radial_refinement": refinements,
            "angular_auxiliary_refinement": [
                {"J": row["auxiliary_J"], "defect": row["sampled_angular_constraint_defect"]}
                for row in angular],
            "independent_angular_quadrature_mass_difference": quadrature_mass_difference}


def main(config):
    start = time.monotonic()
    print(json.dumps({"configuration": config, "python": sys.version, "platform": platform.platform(),
                      "commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID"),
                      "experiment": os.getenv("ORX_EXPERIMENT")}), flush=True)
    report(phase="exact_tt_basis", step=0, of=10, unit="modes")
    modes, mode_data = exact_mode_basis(
        lambda count: report(step=count, m={"exact_TT_residual": 0}))
    source, bound, source_data = axisymmetric_source(modes[1])
    print("EXACT TT BASIS " + json.dumps(mode_data), flush=True)
    print("EXACT GAUNT SOURCE " + json.dumps(source_data), flush=True)

    report(phase="hamiltonian_coefficients", step=0, of=2, unit="quadratures")
    norm = float(F(source_data["source_norm_divided_by_pi_to_3_over_2"])) * pi ** 1.5
    coulomb = float(F(source_data["fourier_coulomb_integral_divided_by_pi_to_3_over_2_over_sqrt2"])) * pi ** 1.5 / sqrt(2)
    integrals = []
    for index, (radial_order, angular_order) in enumerate(config["energy_quadratures"], 1):
        row = energy_integrals(radial_order, angular_order, config["outer_radius"], source)
        row["relative_error_against_exact_fourier_coulomb"] = abs(row["source_times_u2"] - coulomb) / coulomb
        assert row["relative_error_against_exact_fourier_coulomb"] < config["energy_relative_tolerance"]
        assert row["source_times_u2_squared"] > 0
        integrals.append(row)
        print("ENERGY QUADRATURE " + json.dumps(row), flush=True)
        report(step=index, m={"quartic_energy_verification_error": row["relative_error_against_exact_fourier_coulomb"]})
    sixth_difference = abs(integrals[0]["source_times_u2_squared"] -
                           integrals[-1]["source_times_u2_squared"])
    assert sixth_difference / integrals[-1]["source_times_u2_squared"] < config["energy_relative_tolerance"]
    energy = {"quadratic": norm, "quartic": -7 * coulomb,
              "sixth": 77 * integrals[-1]["source_times_u2_squared"]}
    tails = radial_majorant_tails(config["outer_radius"], source)
    assert tails["uniform_linear_potential_tail_bound"] < 1e-30

    report(phase="free_corner_sectors", step=0, of=len(config["area_budgets"]), unit="budgets")
    corner_counts = []
    for index, budget in enumerate(config["area_budgets"], 1):
        row = free_corner_count(budget)
        corner_counts.append(row)
        print("LEADING AREA SECTOR " + json.dumps(row), flush=True)
        report(row="area_sectors", m={"area_budget_B": float(F(row["dimensionless_area_budget_B"])),
                                     "sector_dimension": row["dimension"], "count_residual": 0})
        report(step=index, m={"count_residual": 0})

    report(phase="linear_green_refinement", step=0, of=len(config["radial_intervals"]), unit="grids")
    linear = []
    for index, intervals in enumerate(config["radial_intervals"], 1):
        row = linear_grid_check(intervals, config["outer_radius"], source)
        linear.append(row)
        print("LINEAR GREEN CHECK " + json.dumps(row), flush=True)
        report(row="linear_green", m={"radial_intervals": intervals, "linear_green_error": row["maximum_linear_green_error"]})
        report(step=index, m={"linear_green_error": row["maximum_linear_green_error"]})

    cases = []
    for rho in config["amplitude_rhos"]:
        for intervals in config["radial_intervals"]:
            cases.append({"name": f"rho={rho}, radial={intervals}", "rho": rho,
                          "radial_intervals": intervals, "auxiliary_J": 12, "angular_nodes": 32})
    for auxiliary in (2, 4, 8):
        cases.append({"name": f"auxiliary_J={auxiliary}", "rho": max(config["amplitude_rhos"]),
                      "radial_intervals": max(config["radial_intervals"]),
                      "auxiliary_J": auxiliary, "angular_nodes": 32})
    cases.append({"name": "angular_quadrature_check", "rho": max(config["amplitude_rhos"]),
                  "radial_intervals": max(config["radial_intervals"]),
                  "auxiliary_J": 12, "angular_nodes": 48})
    report(phase="nonlinear_constraints", step=0, of=len(cases), unit="cases")
    nonlinear = []
    for index, case in enumerate(cases, 1):
        row = nonlinear_case(case, source, bound, energy, config)
        nonlinear.append(row)
        print("NONLINEAR CONSTRAINT " + json.dumps(row), flush=True)
        report(row="nonlinear", m={
            "case": case["name"], "rho": case["rho"], "auxiliary_J": case["auxiliary_J"],
            "radial_intervals": case["radial_intervals"],
            "angular_constraint_defect": row["sampled_angular_constraint_defect"],
            "full_constraint_residual": row["sampled_fourth_order_constraint_residual"],
            "ADM_mass": row["ADM_mass_volume"], "mass_flux_difference": row["boundary_volume_difference"]})
        report(step=index, m={"angular_constraint_defect": row["sampled_angular_constraint_defect"],
                             "full_constraint_residual": row["sampled_fourth_order_constraint_residual"]})

    report(phase="export", step=0, of=1, unit="reports")
    convergence = convergence_summary(nonlinear, linear, config)
    for row in convergence["angular_auxiliary_refinement"]:
        report(row="angular_refinement", m={"auxiliary_J": row["J"],
                                           "angular_constraint_defect": row["defect"]})
    result = {
        "config": config, "source_commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID"),
        "model": "Declared reduction-before-cutoff Einstein candidate; asymptotically flat maximal data.",
        "free_mode_basis": mode_data, "exact_source": source_data,
        "energy_coefficients_in_epsilon": energy, "energy_quadratures": integrals,
        "leading_uniform_area_sectors": corner_counts,
        "sixth_integral_quadrature_difference": sixth_difference, "outer_domain_bounds": tails,
        "linear_green_checks": linear, "nonlinear_checks": nonlinear, "convergence": convergence,
        "limitations": [
            "This is a revision of the regulator; equivalence to either historical finite-rank action is not asserted.",
            "The nonlinear numerical family has flat conformal metric and one TT momentum direction. The general canonical map is specified separately with local existence hypotheses.",
            "Dependent geometry contains arbitrarily high harmonics. Auxiliary J is a computation cutoff, not the physical free-data rank N.",
            "The ten profiles do not exhaust radial degrees of freedom. Their oscillator frequencies are basis-dependent, not a new dispersion relation.",
            "Uniform-area counts use the quadratic far-sphere area/ADM relation, a specified vacuum subtraction and a chosen area window. They do not derive fundamental microscopic corner multiplicities.",
            "No complete cubic-through-sixth-order interaction Hamiltonian, interacting one-loop amplitude, physical large-N limit, microscopic corner count or interacting Euclidean reconstruction is established.",
        ],
        "status": "all specified reduced-data constraint and Hamiltonian ingredient checks passed",
        "wall_seconds": time.monotonic() - start,
    }
    output = ROOT / "output"
    output.mkdir(exist_ok=True)
    path = output / "reduced_gaunt_completion_results.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    report(step=1, m={"exact_TT_residual": 0,
                     "full_constraint_residual": nonlinear[-1]["sampled_fourth_order_constraint_residual"],
                     "angular_constraint_defect": nonlinear[-1]["sampled_angular_constraint_defect"],
                     "validated_modes": 10,
                     "validated_nonlinear_cases": len(nonlinear)})
    print("FINAL " + json.dumps({
        "result": str(path), "source_commit": result["source_commit"], "free_angular_rank_N": 3,
        "validated_modes": 10, "nonlinear_cases": len(nonlinear),
        "linear_green_refinement_ratios": convergence["linear_green_refinement_ratios"],
        "auxiliary_angular_defects": convergence["angular_auxiliary_refinement"],
        "energy_coefficients_in_epsilon": energy,
        "status": result["status"], "wall_seconds": result["wall_seconds"],
        "scope": "constraint reconstruction and Hamiltonian ingredients; full quantum goal remains open",
    }, allow_nan=False), flush=True)
