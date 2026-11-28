"""Tracked verification of the full reduced H3/H4 functional ingredients.

These vertices belong to the declared reduction-first regulator. The finite
ten-mode example is not an asymptotic scattering amplitude. No loop integral,
interacting corner spectrum or positive interacting transfer kernel is claimed.
"""

from decimal import Decimal as D
from fractions import Fraction as F
from pathlib import Path
import json
import os
import platform
import sys
import time

from run import report
from reduced_gaunt_vertices import (
    ONE, coulomb, hamiltonian_coefficients, invariants, solve_rational,
    tensor_from_labels, test_cases,
)
from reduced_finite_metric import independent_point_check

ROOT = Path(__file__).resolve().parent


def serial_labels(labels):
    return [[[str(x) for x in row] for row in a] for a in labels]


def main(config):
    start = time.monotonic()
    print(json.dumps({"configuration": config, "python": sys.version,
        "platform": platform.platform(), "commit": os.getenv("ORX_COMMIT"),
        "run": os.getenv("ORX_RUN_ID"), "experiment": os.getenv("ORX_EXPERIMENT")}), flush=True)
    report(phase="exact_vertices", step=0, of=7, unit="directions")
    # Analytic Fourier normalization controls, including an odd field.
    assert coulomb(ONE, ONE) == F(1, 16)
    assert coulomb({(1, 0, 0): F(1)}, {(1, 0, 0): F(1)}) == F(1, 192)
    cases, work = [], {}
    for count, (name, labels) in enumerate(test_cases(), 1):
        q = tensor_from_labels(labels[0], labels[1])
        s = tensor_from_labels(labels[2], labels[3])
        print("EXACT DIRECTION START " + name, flush=True)
        coefficients, polynomials = hamiltonian_coefficients(q, s,
            progress=lambda ingredient: print("EXACT INGREDIENT " + json.dumps(
                {"case": name, "ingredient": ingredient}), flush=True))
        iv = invariants(labels[0], labels[1])
        row = {"name": name, "labels_Eq_Bq_Es_Bs": serial_labels(labels),
               "coefficients": {k: str(v) for k, v in coefficients.items()},
               "invariants": {k: [str(v) for v in x] for k, x in iv.items()},
               "exact_identity_residual": "0"}
        cases.append(row)
        work[name] = (q, s, polynomials)
        print("EXACT VERTICES " + json.dumps(row), flush=True)
        report(row="vertices", m={"case": name,
            "cubic_potential": float(coefficients["H3_configuration"]),
            "quartic_potential": float(coefficients["H4_configuration"]),
            "quartic_mixed": float(coefficients["H4_mixed"]),
            "quartic_momentum": float(coefficients["H4_momentum"]),
            "exact_vertex_residual": 0})
        report(step=count, m={"exact_vertex_residual": 0, "vertex_directions": count})

    report(phase="invariant_reconstruction", step=0, of=1, unit="checks")
    # Two cubic and five quartic rotational invariants span the full
    # parity-even potential on the two selected j=2 multiplets. Fit on
    # dedicated directions; verify on two held-out mixed/rotated choices.
    cubic_indices = (0, 2)
    cubic = solve_rational([cases[i]["invariants"]["cubic"] for i in cubic_indices],
                          [cases[i]["coefficients"]["H3_configuration"] for i in cubic_indices])
    quartic = solve_rational([x["invariants"]["quartic"] for x in cases[:5]],
                            [x["coefficients"]["H4_configuration"] for x in cases[:5]])
    for row in cases:
        for degree, coefficients, key in (("cubic", cubic, "H3_configuration"),
                                          ("quartic", quartic, "H4_configuration")):
            predicted = sum((a*F(b) for a, b in zip(coefficients, row["invariants"][degree])), F(0))
            assert predicted == F(row["coefficients"][key]), (row["name"], degree)
    assert cases[5]["coefficients"] == cases[6]["coefficients"], "Full energy failed a rational SO(3) rotation"
    assert F(cases[0]["coefficients"]["H4_momentum"]) == F(-1305190881, 2621440), "Prior independent momentum result changed"
    invariant_result = {
        "cubic_basis": ["tr(E^3)", "tr(E B^2)"],
        "cubic_coefficients": [str(x) for x in cubic],
        "quartic_basis": ["tr(E^2)^2", "tr(B^2)^2", "tr(E^2) tr(B^2)",
                           "tr(E B)^2", "tr(E^2 B^2)"],
        "quartic_coefficients": [str(x) for x in quartic],
        "fit_directions": [row["name"] for row in cases[:5]],
        "held_out_directions": [row["name"] for row in cases[5:]],
        "invariant_and_full_rotational_residual": "0",
        "H2_unit": "pi^(3/2)", "H3_unit": "(2*pi/3)^(3/2)",
        "H4_unit": "pi^(3/2)/sqrt(2)",
        "canonical_momentum_map": "S=sum (G0^-1 p)^alpha e_alpha; the listed Es,Bs are coefficients of S, not p",
    }
    print("INVARIANT RECONSTRUCTION " + json.dumps(invariant_result), flush=True)
    report(step=1, m={"invariant_residual": 0, "rotational_residual": 0,
                      "nonzero_cubic_invariants": sum(bool(x) for x in cubic)})

    checks = []
    report(phase="independent_finite_metric", step=0, of=len(config["finite_metric_checks"]), unit="points")
    for count, spec in enumerate(config["finite_metric_checks"], 1):
        q, s, polynomials = work[spec["case"]]
        result = independent_point_check(q, s, polynomials, spec["point"],
                                         config["finite_difference_steps"], config["decimal_precision"])
        result["case"] = spec["case"]
        checks.append(result)
        print("INDEPENDENT FINITE METRIC " + json.dumps(result), flush=True)
        error = D(result["maximum_scaled_error"])
        assert error < D(config["coefficient_error_tolerance"]), (spec, str(error))
        assert D(result["maximum_determinant_error"]) < D(config["determinant_tolerance"])
        assert D(result["maximum_connection_trace_error"]) < D(config["determinant_tolerance"])
        report(row="finite_metric", m={"case": spec["case"],
            "finite_metric_error": float(error),
            "metric_determinant_error": float(D(result["maximum_determinant_error"])),
            "exponential_terms": result["maximum_exponential_terms"]})
        report(step=count, m={"finite_metric_error": float(error)})

    report(phase="export", step=0, of=1, unit="reports")
    maximum_error = max(D(row["maximum_scaled_error"]) for row in checks)
    result = {
        "config": config, "source_commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID"),
        "model": "Reduction-first canonical Gaunt candidate; maximal asymptotically flat gravity; 16*pi*G=1 and packet width=1",
        "free_angular_rank_N": 3, "free_modes": 10,
        "scope": "Full H3/H4 functional and evaluated couplings; finite-N loop and full interacting quantum goal remain open",
        "exact_directions": cases, "potential_invariants": invariant_result,
        "independent_finite_metric": checks,
        "max_finite_metric_scaled_error": str(maximum_error),
        "limitations": [
            "This is a declared change of regulator, not equivalence to either historical finite-rank action.",
            "The dependent scalar and longitudinal responses use the untruncated continuum Green operators.",
            "The ten localized modes have one radial profile per tensor parity; they are not asymptotic plane-wave scattering states.",
            "H3/H4 follow from a local canonical constraint chart. Global nonlinear existence and evolution convergence are separate.",
            "H5/H6, ordering and measure terms can contribute to a finite-rank one-loop four-point function and have not been evaluated here.",
            "A negative p^4 Taylor coefficient must not be exponentiated as a globally defined positive quantum Hamiltonian.",
            "No full local corner state derivation, microscopic entropy law, interacting reflection positivity or physical large-N limit is proved."
        ],
        "status": "all specified exact vertex and independent finite-metric checks passed",
        "wall_seconds": time.monotonic()-start,
    }
    output = ROOT/"output"
    output.mkdir(exist_ok=True)
    path = output/"reduced_vertices_results.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    report(step=1, m={"maximum_finite_metric_error": float(maximum_error),
                      "exact_vertex_residual": 0, "verified_hamiltonian_order": 4})
    print("FINAL " + json.dumps({"result": str(path), "source_commit": result["source_commit"],
        "exact_directions": len(cases), "independent_points": len(checks),
        "finite_metric_evaluations": sum(row["finite_amplitudes"] for row in checks),
        "max_finite_metric_scaled_error": str(maximum_error),
        "cubic_invariant_coefficients": invariant_result["cubic_coefficients"],
        "quartic_potential_coefficients": invariant_result["quartic_coefficients"],
        "status": result["status"], "scope": result["scope"], "wall_seconds": result["wall_seconds"]}), flush=True)
