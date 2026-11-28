"""Tracked full H5 spatial assembly with independent geometric/integral checks.

The Hamiltonian functional includes every fifth-order sector of the declared
fixed-mode candidate. Two directional evaluations are not a reconstructed
rank-five tensor, a four-point loop, or a physical large-rank limit.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import atan, pi, sqrt
from pathlib import Path
import json
import os
import platform
import signal
import sys
import time

from run import report
from gaussian_green_fields import HeatField, scalar_field
from gaussian_heat_pairings import PairingAccumulator, collapse_pairing, tail_bound
from reduced_fifth_order import prepare_fields, third_fourth_families, fifth_families, merge_families
from reduced_gaunt_vertices import test_cases, tensor_from_labels, hamiltonian_coefficients, pg
from reduced_gaunt_exact import derivative
from reduced_finite_metric import field_values, polyval, zeros, mscale, exp_jet
from higher_geometry_check import coefficient_weights
from york_operator_fields import york_apply, york_series

ROOT = Path(__file__).resolve().parent


def pairing_controls(options):
    # Exact anisotropic contractions test parity rescue, rate/slot ordering
    # and signed merge behavior before any numerical quadrature is used.
    def signature(accumulator):
        return tuple((g.beta, g.rates, dict(g.coefficients)) for g in accumulator.collapsed())
    x0 = HeatField(F(1), {((1, 0, 0), 0): F(1)})
    y1 = HeatField(F(2), {((0, 1, 0), 1): F(1)})
    xy = {(1, 1, 0): F(1)}
    direct = collapse_pairing(xy, F(3, 2), (x0, y1))
    assert dict(direct.coefficients) == {(4, 0, 1): F(1, 15)}
    reverse = collapse_pairing(xy, F(3, 2), (y1, x0))
    assert direct == reverse
    larger = HeatField(F(1), {((1, 0, 0), 0): F(1), ((0, 0, 1), 0): F(1)})
    assert dict(collapse_pairing(xy, F(3, 2), (larger, y1)).coefficients) == dict(direct.coefficients)
    equal = HeatField(F(1), {((1, 0, 0), 1): F(1)})
    equal_pair = collapse_pairing({(0, 0, 0): F(1)}, F(1), (equal, x0))
    assert dict(equal_pair.coefficients) == {(2, 0, 1): F(1, 3)}
    first = PairingAccumulator().merge(direct)
    second = PairingAccumulator().add(xy, F(3, 2), (x0, y1), -1)
    assert not PairingAccumulator().merge(first).merge(second).collapsed()
    split = PairingAccumulator().merge(first).merge(first)
    together = PairingAccumulator().add(xy, F(3, 2), (x0, y1), 2)
    assert signature(split) == signature(together)
    distinct = PairingAccumulator().merge(first).add(xy, F(2), (x0, y1))
    distinct.add(xy, F(3, 2), (x0, HeatField(F(3), dict(y1.terms))))
    assert len(distinct.collapsed()) == 3
    print("EXACT PAIRING CONTROLS "+json.dumps({"passed": True, "anisotropic_coefficient": "1/15",
        "equal_rate_coefficient": "1/3", "negative_merge_empty": True,
        "split_merge_equal": True, "distinct_rate_groups": 3}), flush=True)
    one = {(0, 0, 0): F(1)}
    beta, a, b = F(3, 2), F(2, 3), F(4, 3)
    beta_f, a_f, b_f = map(float, (beta, a, b))
    f, g = scalar_field(one, a), scalar_field(one, b)
    cases = [
        ("anisotropic_gaussian", {(2, 2, 0): F(1)}, (), pi**1.5/(4*beta_f**3.5)),
        ("single_green", one, (f,), pi**1.5/(16*beta_f*a_f*sqrt(beta_f+a_f))),
        ("two_green", one, (f, g), pi**1.5/(256*sqrt(beta_f)*(a_f*b_f)**1.5)
         *atan(sqrt(a_f*b_f/(beta_f*(beta_f+a_f+b_f))))),
        ("odd_zero", {(1, 0, 0): F(1)}, (f, g), 0.0),
    ]
    results = []
    for name, polynomial, fields, expected in cases:
        accumulator = PairingAccumulator().add(polynomial, beta, fields=fields)
        result = accumulator.integrate(**options)
        error = abs(result["value"]-expected)/(1+abs(expected))
        row = {"name": name, "expected": expected, "scaled_error": error, "integral": result,
               "passed": result["converged"] and error < 1e-10}
        results.append(row)
        print("GAUSSIAN PAIRING CONTROL "+json.dumps(row), flush=True)
        report(step=len(results), m={"pairing_control_error": error})
        assert row["passed"], row
    return results


def york_point_check(data, point, precision=90):
    """Independent finite metric versus operator coefficient convolution.

    A fixed contravariant polynomial test-vector jet is used, so its own
    dependence on gravitational sources cannot contaminate amplitude orders.
    """
    weights = {n: coefficient_weights((-2, -1, 0, 1, 2), n) for n in (1, 2)}
    with localcontext() as context:
        context.prec = precision
        x = [D(a) for a in point]
        vpolys = ({(2, 0, 0): F(1), (0, 1, 0): F(1)},
                  {(1, 0, 1): F(1), (0, 1, 0): F(-2)},
                  {(0, 2, 0): F(1), (1, 0, 0): F(1, 3)})
        w = [polyval(p, x) for p in vpolys]
        dw = [[polyval(derivative(p, i), x) for p in vpolys] for i in range(3)]
        ddw = [[[polyval(derivative(derivative(p, i), j), x) for p in vpolys]
                for j in range(3)] for i in range(3)]
        geometry = data["geometry"]
        r2 = sum((a*a for a in x), D(0))
        def physical(poly, n):
            return polyval(poly, x)*(-D(n)*r2/2).exp()
        metric, dm, ddm, inverse, gamma = [], [], [], [], [None]
        for n in range(3):
            gm = geometry["metric"][n]
            metric.append([[physical(p, n) for p in row] for row in gm])
            dm.append([[[physical(pg(gm[i][j], axis, n), n) for j in range(3)]
                        for i in range(3)] for axis in range(3)])
            ddm.append([[[[physical(pg(pg(gm[i][j], axis, n), other, n), n)
                           for j in range(3)] for i in range(3)]
                         for other in range(3)] for axis in range(3)])
            inverse.append([[physical(p, n) for p in row] for row in geometry["inverse"][n]])
            if n:
                gamma.append([[[physical(p, n) for p in row] for row in tensor]
                              for tensor in geometry["gamma"][n]])
        target = york_series(metric, dm, ddm, inverse, gamma, w, dw, ddw)
        flat = [sum((ddw[k][k][i]+ddw[i][k][k]/3 for k in range(3)), D(0)) for i in range(3)]
        flat_error = max(abs(a-b) for a, b in zip(flat, target["A"][0]))
        trace_error = abs(sum((target["L"][1][i][i] for i in range(3)), D(0))
                          -sum((metric[1][i][j]*target["L"][0][j][i]
                                for i in range(3) for j in range(3)), D(0)))
        q, dq, ddq = field_values(data["q"], x)
        steps = (D("0.0001"), D("0.00005"), D("0.000025"))
        finite = {}
        for h in sorted({D(k)*step for step in steps for k in (-2, -1, 0, 1, 2)}):
            first = [mscale(g, h) for g in dq]+[zeros()]
            second = [[zeros() for _ in range(4)] for _ in range(4)]
            for i in range(3):
                for j in range(3):
                    second[i][j] = mscale(ddq[i][j], h)
            g, dg, ddg, _, jet_tail = exp_jet(mscale(q, h), first, second, precision)
            finite[h] = york_apply(g, dg[:3], [row[:3] for row in ddg[:3]], w, dw, ddw)
            finite[h]["jet_tail"] = jet_tail
        estimates = []
        for h in steps:
            estimates.append({n: [sum((D(weight.numerator)/D(weight.denominator)*finite[D(k)*h]["A"][i]
                                       for k, weight in zip((-2, -1, 0, 1, 2), weights[n])), D(0))/h**n
                                  for i in range(3)] for n in (1, 2)})
        refined = [{n: [(16*estimates[j][n][i]-estimates[j-1][n][i])/15 for i in range(3)]
                    for n in (1, 2)} for j in (1, 2)]
        errors = [abs(refined[-1][n][i]-target["A"][n][i])/(1+abs(target["A"][n][i]))
                  for n in (1, 2) for i in range(3)]
        refinement = max(abs(refined[-1][n][i]-refined[-2][n][i])/(1+abs(target["A"][n][i]))
                         for n in (1, 2) for i in range(3))
        maximum = max(errors)
        determinant_error = max(abs(v["determinant"]-1) for v in finite.values())
        connection_error = max(abs(sum((v["gamma"][k][i][k] for k in range(3)), D(0)))
                               for v in finite.values() for i in range(3))
        jet_tail = max(v["jet_tail"] for v in finite.values())
        return {"point": point, "precision": precision, "finite_amplitudes": len(finite),
                "A1": [str(v) for v in target["A"][1]], "A2": [str(v) for v in target["A"][2]],
                "richardson": {str(n): [str(v) for v in refined[-1][n]] for n in (1, 2)},
                "maximum_scaled_error": str(maximum), "refinement_scaled_error": str(refinement),
                "flat_operator_error": str(flat_error), "linear_tracelessness_error": str(trace_error),
                "determinant_error": str(determinant_error), "connection_trace_error": str(connection_error),
                "maximum_last_jet_term": str(jet_tail),
                "passed": max(maximum, refinement) < D("1e-18")
                          and max(flat_error, trace_error, determinant_error, connection_error, jet_tail) < D("1e-70")}


def calibrate(data, options):
    coefficients, _ = hamiltonian_coefficients(data["q"], data["s"])
    u3, u4 = (2*pi/3)**1.5, pi**1.5/sqrt(2)
    targets = {
        "H3_curvature": u3*float(coefficients["H3_configuration"]),
        "H4_local": u4*float(coefficients["quartic_local_curvature"]+coefficients["quartic_raw_commutator"]),
        "H4_york": u4*float(coefficients["quartic_longitudinal"]),
        "H4_scalar": u4*float(coefficients["quartic_curvature_coulomb"]+coefficients["quartic_scalar_mixed"]+coefficients["H4_momentum"]),
    }
    results = {}
    for name, family in third_fourth_families(data).items():
        integral = family.integrate(**options)
        expected = targets[name]
        error = abs(integral["value"]-expected)/(1+abs(expected))
        row = {"integral": integral, "expected": expected, "scaled_error": error,
               "passed": integral["converged"] and error < 1e-8}
        results[name] = row
        print("H3 H4 PAIRING CALIBRATION "+json.dumps({"family": name, **row}), flush=True)
        assert row["passed"], (name, row)
    return results


def main(config):
    from h5_spatial_reference import integrate_h5_reference
    start = time.monotonic()
    def timeout(signum, frame):
        raise TimeoutError("Declared local H5 wall-clock limit reached")
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(int(config["wall_limit_seconds"]))
    print(json.dumps({"configuration": config, "python": sys.version, "platform": platform.platform(),
                      "source_commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID")}), flush=True)
    options = config["pairing_options"]
    report(phase="pairing_controls", step=0, of=4, unit="cases")
    controls = pairing_controls({**options, "absolute_tolerance": 1e-12, "relative_tolerance": 1e-10})
    names = config["h5_directions"]
    labels = dict(test_cases())
    prepared = []
    counter = 0
    report(phase="field_preparation", step=0, of=15*len(names), unit="checks")
    for name in names:
        def field_progress(row):
            nonlocal counter
            counter += 1
            print("H5 FIELD PROGRESS "+json.dumps({"direction": name, **row}), flush=True)
            report(step=counter)
        eq, bq, es, bs = labels[name]
        data = prepare_fields(tensor_from_labels(eq, bq), tensor_from_labels(es, bs), progress=field_progress)
        prepared.append((name, data))
        counter += 1
        report(step=counter)
    assert counter == 15*len(names)
    calibration = []
    report(phase="quartic_calibration", step=0, of=len(names), unit="directions")
    for name, data in prepared:
        row = {"direction": name, "families": calibrate(data, options),
               "york_operator": york_point_check(data, config["york_test_point"])}
        print("H5 CALIBRATION SUMMARY "+json.dumps(row), flush=True)
        assert row["york_operator"]["passed"], row["york_operator"]
        calibration.append(row)
        report(step=len(calibration), m={"calibration_error": max(r["scaled_error"] for r in row["families"].values()),
                                         "york_operator_error": float(row["york_operator"]["maximum_scaled_error"])})
    assembled = []
    counter = 0
    report(phase="h5_assembly", step=0, of=6*len(names), unit="families")
    for name, data in prepared:
        def family_progress(row):
            nonlocal counter
            counter += 1
            print("H5 ASSEMBLY "+json.dumps({"direction": name, **row}), flush=True)
            report(step=counter)
        families = fifth_families(data, config["max_pair_products"], progress=family_progress)
        assembled.append((name, data, families, merge_families(families)))
    radial = []
    report(phase="h5_radial", step=0, of=len(names), unit="directions")
    for name, data, families, total in assembled:
        result = total.integrate(**options)
        family_results = {key: value.integrate(**options) for key, value in families.items()}
        row = {"direction": name, "integral": result, "families": family_results}
        radial.append(row)
        print("FULL H5 RADIAL "+json.dumps(row), flush=True)
        assert result["converged"] and all(r["converged"] for r in family_results.values()), row
        report(step=len(radial), m={"H5_value": result["value"], "H5_refinement": result["termwise_refinement_change"]})
    references = []
    levels = len(config["spatial_radial_orders"])
    report(phase="spatial_reference", step=0, of=levels*len(names), unit="quadratures")
    for case_index, (name, data, families, total) in enumerate(assembled):
        # Each retained input is j=2. Rotation-equivariant inverses and
        # contractions make each homogeneous degree-five scalar j<=10.
        # This is specific to the declared ten-profile basis, not a general
        # bound for arbitrary tensor inputs accepted by the assembly module.
        assert all(nphi > 10 and 2*nmu-1 >= 10 for nmu, nphi in config["spatial_angular_orders"])
        def spatial_progress(row):
            print("H5 SPATIAL PROGRESS "+json.dumps({"direction": name, **row}), flush=True)
            if row["panel"] == row["panels"]:
                report(step=case_index*levels+row["level"])
        reference = integrate_h5_reference(data, radial_orders=config["spatial_radial_orders"],
            angular_orders=config["spatial_angular_orders"], radius=config["spatial_radius"],
            panel_width=config["spatial_panel_width"], progress=spatial_progress)
        target = radial[case_index]
        latest = reference["latest"]
        error = abs(latest["total_strong"]-target["integral"]["value"])/(1+abs(target["integral"]["value"]))
        weak_error = abs(latest["york_forms"]["difference"])/(1+abs(latest["york_forms"]["weak"]))
        scalar_ibp_error = abs(latest["scalar_gradient_forms"]["difference"])/(1+abs(latest["scalar_gradient_forms"]["weak"]))
        refinement = max(latest["refinement_from_previous"]["total_strong"]["scaled_delta"],
                         *(r["scaled_delta"] for r in latest["refinement_from_previous"]["families"].values()))
        family_errors = {key: abs(latest["families"][key]-value["value"])/(1+abs(value["value"]))
                         for key, value in target["families"].items()}
        tails = {"total": tail_bound(total, config["spatial_radius"])}
        tails.update({key: tail_bound(value, config["spatial_radius"]) for key, value in families.items()})
        tail_scaled = []
        margins = {}
        for key in tails:
            integral = target["integral"] if key == "total" else target["families"][key]
            bound = tails[key]["absolute_bound"]
            assert bound is not None, tails[key]
            tail_scaled.append(bound/(1+abs(integral["value"])))
            margins[key] = config["roundoff_safety_factor"]*integral["machine_epsilon_times_bin_l1"]/(1+abs(integral["value"]))
        roundoff_gate = max(margins.values()) < config["h5_scaled_tolerance"]/10
        maximum = max(error, weak_error, scalar_ibp_error, refinement, *family_errors.values(), *tail_scaled)
        row = {"direction": name, "reference": reference, "scaled_total_error": error,
               "york_integration_by_parts_scaled_error": weak_error, "family_scaled_errors": family_errors,
               "scalar_integration_by_parts_scaled_error": scalar_ibp_error,
               "reference_radius_tail_bounds": tails, "scaled_roundoff_indicators_with_safety_factor": margins,
               "roundoff_margin_passed": roundoff_gate,
               "angular_exactness_premise": "Five total-j=2 inputs and rotation-equivariant operations imply scalar j<=10; nphi>10 and 2*nmu-1>=10",
               "maximum_scaled_error": maximum, "passed": maximum < config["h5_scaled_tolerance"] and roundoff_gate}
        references.append(row)
        print("H5 INDEPENDENT SPATIAL COMPARISON "+json.dumps(row), flush=True)
        assert row["passed"], row
        report(step=(case_index+1)*levels, m={"H5_spatial_error": maximum})
    report(phase="export", step=0, of=1, unit="reports")
    result = {"configuration": config, "source_commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID"),
              "pairing_controls": controls, "calibration": calibration, "fifth_order_integrals": radial,
              "independent_spatial_references": references,
              "status": "all specified full-H5 directional and independent spatial checks passed",
              "scope": "Two directional evaluations of the full classical H5 functional; H6 and four-point loop remain pending",
              "units": "raw TT field directions, width1 and16*pi*G=1; S=G0^-1 p, not normalized canonical p",
              "wall_seconds": time.monotonic()-start}
    output = ROOT/"output"
    output.mkdir(exist_ok=True)
    path = output/"fifth_order_results.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    report(step=1, m={"H5_directions": len(radial), "H5_maximum_spatial_error": max(r["maximum_scaled_error"] for r in references)})
    print("FINAL "+json.dumps({"result": str(path), "source_commit": result["source_commit"],
                               "H5_directions": len(radial), "status": result["status"],
                               "scope": result["scope"], "wall_seconds": result["wall_seconds"]}), flush=True)
    signal.alarm(0)
