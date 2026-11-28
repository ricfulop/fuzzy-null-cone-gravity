"""Tracked higher-order local geometry and first Green-response validation.

This checks inputs to the gravitational H5/H6 spatial oracle. It does not
evaluate that oracle, a four-point loop, or an angular/radial continuum limit.
All scientific calculations occur only when main(config) is called by run.py.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import fsum, pi, sqrt
from pathlib import Path
import json
import os
import platform
import signal
import sys
import time

from run import report
from higher_local_geometry import local_geometry, compare_legacy
from reduced_gaunt_exact import tensor_square, scale
from reduced_gaunt_vertices import (
    test_cases, tensor_from_labels, integrate, matrix_add, matrix_multiply,
    trace, solve_rational,
)
from reduced_finite_metric import (
    field_values, geometry_at_amplitude, polyval, decimal_fraction,
)

ROOT = Path(__file__).resolve().parent


def coefficient_weights(nodes, degree):
    """Exact coefficient extraction for polynomials up to len(nodes)-1.

    The finite-metric evaluator has no access to these Taylor coefficients.
    No expected low-order term is subtracted from its values.
    """
    matrix = [[F(x)**k for x in nodes] for k in range(len(nodes))]
    rhs = [F(int(k == degree)) for k in range(len(nodes))]
    weights = solve_rational(matrix, rhs)
    assert all(sum((w*F(x)**k for w, x in zip(weights, nodes)), F(0))
               == rhs[k] for k in range(len(nodes)))
    return weights


def point_check(q, s, geometry, point, config, progress=None):
    """Compare high Taylor orders with independent finite exp-metric jets.

    Nine centered samples exactly resolve degrees 0..8. At degree n the
    leading truncation order is 10-n (n even), or 9-n (n odd). Richardson
    eliminates it using two consecutive half-step grids. Decimal precision
    resolves the cancellations in the sixth-order coefficient.
    """
    precision = int(config["geometry_precision"])
    nodes = tuple(range(-4, 5))
    degrees = (2, 3, 4, 5, 6)
    weights = {n: coefficient_weights(nodes, n) for n in degrees}
    with localcontext() as context:
        context.prec = precision
        x = [D(v) for v in point]
        hs = [D(v) for v in config["geometry_steps"]]
        assert len(hs) == 3 and hs[0] == 2*hs[1] == 4*hs[2]
        qv, sv = field_values(q, x), field_values(s, x)
        amplitudes = sorted({D(k)*h for h in hs for k in nodes})
        finite = {}
        for number, h in enumerate(amplitudes, 1):
            finite[h] = geometry_at_amplitude(qv, sv, h, precision)
            if progress:
                progress(number, len(amplitudes))
        envelope = (-sum((a*a for a in x), D(0))/2).exp()
        targets = {}
        for n in (2, 3, 4):
            targets["R"+str(n)] = ("R", None, n,
                polyval(geometry["R"][n], x)*envelope**n)
        for n in (2, 3, 4, 5, 6):
            targets["L"+str(n)] = ("GG", None, n,
                polyval(geometry["L"][n], x)*envelope**n)
        for n in (2, 3):
            for i in range(3):
                targets["b"+str(n)+"_"+str(i)] = ("divergence", i, n,
                    polyval(geometry["b"][n][i], x)*envelope**(n+1))
        commutator = matrix_add(matrix_multiply(q, s), matrix_multiply(s, q), -1)
        twice = matrix_add(matrix_multiply(q, commutator),
                           matrix_multiply(commutator, q), -1)
        raw2 = scale(trace(matrix_multiply(commutator, commutator)), -F(1, 12))
        raw4 = scale(trace(matrix_multiply(twice, twice)), F(1, 360))
        targets["raw_gram2"] = ("raw_gram", None, 2, polyval(raw2, x)*envelope**4)
        targets["raw_gram4"] = ("raw_gram", None, 4, polyval(raw4, x)*envelope**6)
        estimates = []
        for h in hs:
            row = {}
            for name, (key, component, degree, _) in targets.items():
                values = [finite[D(k)*h][key] for k in nodes]
                if component is not None:
                    values = [v[component] for v in values]
                row[name] = sum((decimal_fraction(w)*v
                                 for w, v in zip(weights[degree], values)), D(0))/h**degree
            estimates.append(row)
        refinements = []
        for index in (1, 2):
            row = {}
            for name, (_, _, degree, _) in targets.items():
                exponent = (10 if degree % 2 == 0 else 9)-degree
                ratio = D(2)**exponent
                row[name] = (ratio*estimates[index][name]-estimates[index-1][name])/(ratio-1)
            refinements.append(row)
        errors = {name: abs(refinements[-1][name]-target[3])/(1+abs(target[3]))
                  for name, target in targets.items()}
        differences = {name: abs(refinements[1][name]-refinements[0][name])/(1+abs(target[3]))
                       for name, target in targets.items()}
        maximum = max(errors.values())
        convergence = max(differences.values())
        tolerance = D(str(config["geometry_scaled_tolerance"]))
        determinant_error = max(v["determinant_error"] for v in finite.values())
        connection_error = max(v["unimodular_connection_error"] for v in finite.values())
        jet_tail = max(v["last_jet_term"] for v in finite.values())
        evaluator_tolerance = D(str(config["finite_metric_evaluator_tolerance"]))
        return {
            "point": point, "steps": [str(h) for h in hs], "precision": precision,
            "stencil_nodes": nodes, "finite_amplitudes": len(finite),
            "expected_coefficients": {k: str(v[3]) for k, v in targets.items()},
            "centered_coefficients": [{k: str(v) for k, v in row.items()} for row in estimates],
            "richardson_coefficients": [{k: str(v) for k, v in row.items()} for row in refinements],
            "scaled_errors": {k: str(v) for k, v in errors.items()},
            "refinement_differences": {k: str(v) for k, v in differences.items()},
            "maximum_scaled_error": str(maximum),
            "maximum_refinement_difference": str(convergence),
            "maximum_determinant_error": str(determinant_error),
            "maximum_connection_trace_error": str(connection_error),
            "maximum_last_jet_term": str(jet_tail),
            "maximum_exponential_terms": max(v["exponential_terms"] for v in finite.values()),
            "tolerance": str(tolerance), "evaluator_tolerance": str(evaluator_tolerance),
            "passed": (maximum < tolerance and convergence < tolerance
                       and max(determinant_error, connection_error, jet_tail) < evaluator_tolerance),
        }


def axis_sixth_check(s, config, progress=None):
    """Independent radial-multipole comparison and pure-p sixth-order control.

    The new Cartesian HEAT field and old radial Legendre Green function
    share the physical source, but not their inverse/evaluation algorithms.
    Their full-space potentials are evaluated before the final Gaussian-
    weighted integral is cut off; its omitted tail has a positive bound.
    """
    from gaussian_green_fields import scalar_field, evaluate_source
    from reduced_gaunt_exact import axisymmetric_source
    from reduced_gaunt_completion import (
        quadrature, linear_green, legendre_values, energy_integrals, radial_majorant_tails,
        source_point,
    )
    source, potential_bound, metadata = axisymmetric_source(s)
    poly = tensor_square(s)
    field = scalar_field(poly, F(1))
    values = []
    for r in (0.0, 0.5, 2.0, 6.0, 12.0):
        for mu in ((0.0,) if r == 0 else (0.0, 0.4, 0.9)):
            point = (r*sqrt(1-mu*mu), 0.0, r*mu)
            legendre = legendre_values(4, mu)
            reference = fsum(linear_green(ell, r, source)*legendre[ell] for ell in source)
            current = field.evaluate(point)
            error = abs(current-reference)/(1+abs(reference))
            source_value = evaluate_source(poly, point, F(1))
            source_reference = source_point(r, mu)
            source_error = abs(source_value-source_reference)/(1+abs(source_reference))
            values.append({"r": r, "cos_theta": mu, "heat_value": current,
                           "radial_value": reference, "scaled_error": error,
                           "cartesian_source": source_value, "positive_frame_source": source_reference,
                           "source_scaled_error": source_error})
    radius = float(config["axis_outer_radius"])
    angular_order = int(config["axis_angular_order"])
    rows = []
    for order in config["axis_radial_orders"]:
        quartic, sixth = [], []
        angles = quadrature(angular_order)
        for node, weight in quadrature(order):
            r = radius*(node+1)/2
            for mu, angular_weight in angles:
                point = (r*sqrt(1-mu*mu), 0.0, r*mu)
                density = evaluate_source(poly, point, F(1))
                u = field.evaluate(point)
                # quadrature() weights are normalized to one on [-1,1].
                measure = 4*pi*radius*weight*angular_weight*r*r
                quartic.append(measure*density*u)
                sixth.append(measure*density*u*u)
        direct = {"source_times_u2": fsum(quartic), "source_times_u2_squared": fsum(sixth)}
        reference = energy_integrals(order, angular_order, radius, source)
        errors = {key: abs(value-reference[key])/(1+abs(reference[key])) for key, value in direct.items()}
        exact_coulomb = float(F(metadata["fourier_coulomb_integral_divided_by_pi_to_3_over_2_over_sqrt2"]))*pi**1.5/sqrt(2)
        errors["exact_fourier_coulomb"] = abs(direct["source_times_u2"]-exact_coulomb)/(1+abs(exact_coulomb))
        row = {"radial_order": order, "angular_order": angular_order, "outer_radius": radius,
               "heat_integrals": direct, "radial_integrals": reference, "scaled_errors": errors,
               "pure_p_H6": 77*direct["source_times_u2_squared"]}
        rows.append(row)
        if progress:
            progress(row)
    refinement = max(abs(rows[-1]["heat_integrals"][key]-rows[-2]["heat_integrals"][key])/
                     (1+abs(rows[-1]["heat_integrals"][key])) for key in rows[-1]["heat_integrals"])
    tails = radial_majorant_tails(radius, source)
    sixth_tail = 77*float(potential_bound)**2*tails["source_integral_tail_bound"]
    maximum = max(max(v["scaled_error"] for v in values),
                  max(v["source_scaled_error"] for v in values), refinement,
                  max(error for row in rows for error in row["scaled_errors"].values()))
    tolerance = float(config["axis_scaled_tolerance"])
    return {"points": values, "quadratures": rows, "refinement_scaled_error": refinement,
            "maximum_scaled_error": maximum, "pure_p_H6_tail_bound": sixth_tail,
            "tolerance": tolerance, "tail_tolerance": config["axis_tail_tolerance"],
            "passed": maximum < tolerance and sixth_tail < float(config["axis_tail_tolerance"]),
            "scope": "Independent check of the already derived pure-momentum H6 sector only"}


def main(config):
    from green_field_check import run_green_checks
    started = time.monotonic()
    wall_limit = int(config["wall_limit_seconds"])
    def timeout(signum, frame):
        raise TimeoutError("Declared local wall-clock limit reached")
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(wall_limit)
    print(json.dumps({"configuration": config, "python": sys.version,
                      "platform": platform.platform(), "source_commit": os.getenv("ORX_COMMIT"),
                      "run": os.getenv("ORX_RUN_ID"), "experiment": os.getenv("ORX_EXPERIMENT")}), flush=True)
    chosen = dict(test_cases())
    names = config["geometry_directions"]
    exact, independent = [], []
    sources, vectors = [], []
    pending = []
    # Per direction: 5 connection orders + 4 pointwise curvatures +
    # 5 integrated densities + 4 York source orders, and one legacy check.
    report(phase="higher_geometry", step=0, of=19*len(names), unit="checks")
    counter = 0
    def geometry_progress(row):
        nonlocal counter
        counter += 1
        print("HIGHER GEOMETRY PROGRESS "+json.dumps(row), flush=True)
        report(step=counter)
    for name in names:
        labels = chosen[name]
        q = tensor_from_labels(*labels[:2])
        s = tensor_from_labels(*labels[2:])
        geometry = local_geometry(q, s, progress=geometry_progress)
        legacy = compare_legacy(q, s, geometry)
        row = {"direction": name, "legacy": legacy, "checks": geometry["checks"],
               "curvature_integrals_divided_by_gaussian_unit": {
                   str(n): str(integrate(p, n)) for n, p in geometry["L"].items()},
               "pointwise_curvature_monomials": {str(n): len(p) for n, p in geometry["R"].items()},
               "density_monomials": {str(n): len(p) for n, p in geometry["L"].items()},
               "york_monomials": {str(n): [len(p) for p in vector]
                                   for n, vector in enumerate(geometry["b"])},
               "gaussian_unit": "(2*pi/n)^(3/2), separately for each integrated order n"}
        exact.append(row)
        print("EXACT HIGHER GEOMETRY "+json.dumps(row), flush=True)
        counter += 1
        report(step=counter, m={"exact_geometry_residual": max(legacy["residual_monomials"].values())})
        assert legacy["passed"], row
        # Geometry metadata counts half-Gaussian rates n; Green APIs take
        # exp(-a*r^2). Convert explicitly rather than sharing the wrong rate.
        curvature_rate = F(geometry["gaussian_rates"]["R"][2], 2)
        york_rate = F(geometry["gaussian_rates"]["b"][1], 2)
        assert curvature_rate == york_rate == 1
        sources.extend(((name+":R2", geometry["R"][2], curvature_rate),
                        (name+":S2", tensor_square(s), F(1))))
        vectors.append((name+":b1", tuple(geometry["b"][1]), york_rate))
        pending.append((name, q, s, geometry))
    assert counter == 19*len(names)
    report(phase="finite_metric", step=0, of=len(names)*len(config["geometry_points"]), unit="points")
    for name, q, s, geometry in pending:
        for point in config["geometry_points"]:
            def finite_progress(done, total):
                if done == 1 or done % 8 == 0 or done == total:
                    print("FINITE METRIC PROGRESS "+json.dumps({"direction": name, "point": point,
                                                               "completed_amplitudes": done, "of": total}), flush=True)
            row = {"direction": name, **point_check(q, s, geometry, point, config, finite_progress)}
            independent.append(row)
            print("INDEPENDENT HIGHER COEFFICIENTS "+json.dumps(row), flush=True)
            report(step=len(independent), m={"higher_coefficient_error": float(row["maximum_scaled_error"]),
                                              "higher_coefficient_refinement": float(row["maximum_refinement_difference"])})
            assert row["passed"], row
    report(phase="green_responses", step=0, unit="checks")
    green_count = 0
    def green_progress(row):
        nonlocal green_count
        green_count += 1
        print("GREEN RESPONSE CHECK "+json.dumps(row), flush=True)
        report(step=green_count, m={"green_checks": green_count})
    green = run_green_checks(config, sources, vectors, progress=green_progress)
    print("GREEN RESPONSE SUMMARY "+json.dumps(green), flush=True)
    assert green["passed"], green
    report(phase="axis_sixth_control", step=0, of=len(config["axis_radial_orders"]), unit="quadratures")
    axis_count = 0
    def axis_progress(row):
        nonlocal axis_count
        axis_count += 1
        print("AXIS SIXTH ORDER CONTROL "+json.dumps(row), flush=True)
        report(step=axis_count, m={"axis_sixth_control_error": max(row["scaled_errors"].values())})
    axis_s = next(s for name, q, s, geometry in pending if name == "electric_axis")
    axis = axis_sixth_check(axis_s, config, progress=axis_progress)
    print("AXIS SIXTH ORDER SUMMARY "+json.dumps(axis), flush=True)
    assert axis["passed"], axis
    report(phase="export", step=0, of=1, unit="reports")
    result = {
        "configuration": config, "source_commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID"),
        "physical_free_rank_N": 3, "retained_free_modes": 10,
        "exact_geometry": exact, "independent_finite_metric": independent, "green_responses": green,
        "axis_pure_momentum_sixth_control": axis,
        "maximum_higher_coefficient_error": max(float(r["maximum_scaled_error"]) for r in independent),
        "status": "all specified higher local geometry and first Green-response checks passed",
        "scope": "H5/H6 ingredients; complete higher reduced Hamiltonian and four-point loop remain unevaluated",
        "limitations": ["Two selected directions in the fixed ten-mode basis, not a full tensor reconstruction.",
                        "Green responses invert only polynomial-Gaussian sources; nested nonlinear inverses remain pending.",
                        "No physical large-N/radial/asymptotic scattering limit or microscopic corner multiplicities."],
        "wall_seconds": time.monotonic()-started,
    }
    output = ROOT/"output"
    output.mkdir(exist_ok=True)
    path = output/"higher_geometry_results.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    report(step=1, m={"maximum_higher_coefficient_error": result["maximum_higher_coefficient_error"],
                      "validated_geometry_directions": len(exact)})
    print("FINAL "+json.dumps({"result": str(path), "source_commit": result["source_commit"],
                               "exact_directions": len(exact), "independent_points": len(independent),
                               "green_checks": green_count,
                               "maximum_higher_coefficient_error": result["maximum_higher_coefficient_error"],
                               "status": result["status"], "scope": result["scope"],
                               "wall_seconds": result["wall_seconds"]}), flush=True)
    signal.alarm(0)
