"""PROPOSED local extraction/control run. Do not execute without approval.

This exports the exact stationary radial equations and the discarded response.
It does not solve the radial DAE, provide an admissible non-gauge seed, or
evaluate a nonlinear obstruction.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import time

import numpy as np
import sympy as s

import stationary_system as st

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/oblate-stationary-system-check-2026-11-27"


def exact_zero(values, label):
    if any(s.cancel(s.together(value)) != 0 for value in values):
        raise AssertionError(label)


def check_cross_formula(coords, epsilon):
    _, r, x, _ = coords
    f, w = s.Function("f")(r), s.Function("w")(r)
    even = {label: s.Symbol("u"+str(index), real=True) for index, label in enumerate(st.COLUMNS)}
    gauge = st.odd_coordinate_leg(r, epsilon, f)
    defect = st.n.expected_discarded_generator(gauge, w)
    target = {label: st.response_prefactor(r, epsilon, f, w)*value
              for label, value in st.RESPONSE_TEST.items()}
    exact_zero([value-target.get(label, 0) for label, value in defect.items()], "two-channel map")
    entries = s.symbols("Gtt Grr Gxx Gphiphi Grx", real=True)
    G = s.diag(*entries[:4])
    G[1, 2] = G[2, 1] = entries[4]
    direct = st.n.candidate_source_density(even | gauge, G, coords, epsilon, w)
    factored = st.response_prefactor(r, epsilon, f, w)*st.n.euler_pair_density(
        G, st.n.coefficient_fields(st.RESPONSE_TEST, x), coords, epsilon)
    exact_zero([direct-factored], "mixed parity contraction")


def quadrature_controls(data, r, x, epsilon):
    """Independent numerical angular integration of every exported density."""
    maxima = {}
    count = 0
    for epsilon_value in [s.Rational(1, 5), s.Rational(1, 100)]:
        point = {r: s.Rational(3, 2), epsilon: epsilon_value}
        for node_count in [48, 96]:
            nodes, weights = np.polynomial.legendre.leggauss(node_count)
            largest = 0.0
            for exact_matrix, density_matrix in zip(data["operators"], data["densities"]):
                for expected, density in zip(exact_matrix, density_matrix):
                    if density == 0:
                        continue
                    reference = float(expected.subs(point).evalf(60))
                    function = s.lambdify(x, density.subs(point), "numpy")
                    values = np.broadcast_to(np.asarray(function(nodes), dtype=float), nodes.shape)
                    measured = float(weights @ values / 2)
                    error = abs(measured-reference)/max(1.0, abs(reference))
                    if not np.isfinite(error) or error >= 1e-11:
                        raise AssertionError(("independent angular quadrature", epsilon_value,
                                              node_count, reference, measured))
                    largest = max(largest, error)
                    count += 1
            maxima[f"epsilon={epsilon_value};nodes={node_count}"] = largest
    return {"nonzero_density_comparisons": count, "maximum_scaled_errors": maxima,
            "tolerance": 1e-11, "r": "3/2",
            "purpose": "angular coefficient control; radius is not frozen in the exported operators"}


def verify():
    start = time.monotonic()
    t, x, phi = s.symbols("t x phi", real=True)
    r, epsilon = s.symbols("r epsilon", positive=True)
    coords = [t, r, x, phi]
    check_cross_formula(coords, epsilon)
    data = st.assemble(coords, epsilon)
    operators = data["operators"]
    f = s.Function("f")(r)
    exact_zero(st.apply_operator(operators, st.radial_gauge_profiles(r, epsilon, f), r),
               "radial gauge including response")
    time_rescaling = [0]*8
    time_rescaling[6] = 1
    exact_zero(st.apply_operator(operators, time_rescaling, r), "constant lapse control")
    parameter = st.oblate_parameter_profiles(r, epsilon)
    fields = st.n.coefficient_fields(dict(zip(st.COLUMNS, parameter)), x)
    exact_zero(st.n.linear_metric(fields, coords, epsilon)
               - st.n.background(coords, epsilon).diff(epsilon), "oblate parameter tangent")
    exact_zero(st.apply_operator(operators, parameter, r), "oblate family tangent")
    # This deliberately off-shell profile must not pass the retained rows.
    invalid = [0]*8
    invalid[6] = r*r
    invalid_rows = st.apply_operator(operators, invalid, r)[:-1, :]
    point = {r: s.Rational(3, 2), epsilon: s.Rational(1, 5)}
    magnitudes = [abs(float(value.subs(point).evalf(60))) for value in invalid_rows]
    if not all(np.isfinite(value) for value in magnitudes) or max(magnitudes) <= 1e-8:
        raise AssertionError("off-shell lapse control was not rejected")
    quadrature = quadrature_controls(data, r, x, epsilon)
    output = {
        "status": "stationary differential expressions and controls passed; radial DAE unsolved",
        "column_labels": st.COLUMNS, "retained_row_labels": data["rows"],
        "selected_row_indices": data["selected_row_indices"],
        "automatic_zero_row_indices": data["automatic_row_indices"],
        "operator_convention": "sum_j M_j(r,epsilon) u^(j)(r); no derivative acts on M_j here",
        "last_row": "discarded response D=6 E_X4-E_Cpolar4; not a retained equation",
        "angular_measure": "dx/2; Sigma already occurs once in Euler density",
        "parameter_domain": "r>0, epsilon>0; no explicit regular-Legendre epsilon radius certified",
        "quadrature_controls": quadrature,
        "radial_projection_performed": False, "radial_DAE_reduced": False,
        "smooth_solution_germ_certified": False, "admissible_non_gauge_seed_supplied": False,
        "response_evaluated_on_admissible_non_gauge_seed": False,
        "full_nonlinear_Dirac_classification_established": False,
    }
    # Write arrays of expression strings without SymPy coercing strings back to Expr.
    output["matrices_by_derivative_order"] = [
        [[str(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]
        for matrix in operators]
    sources = [Path(__file__), Path(st.__file__), Path(st.n.__file__),
               ROOT/"sections/nonlinear_action.tex"]
    manifest = {
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "python": platform.python_version(),
        "dependencies": {name: importlib.metadata.version(name) for name in ["numpy", "sympy"]},
        "seconds": round(time.monotonic()-start, 3),
        "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "scientific_invocations": 1, "GPU_jobs": 0, "remote_jobs": 0, "tracked_experiment_launches": 0,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/"stationary-operators.json").write_text(json.dumps(output, indent=2)+"\n")
    (OUT/"stationary-reproducibility.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print("PASS: all 22 retained rows checked; eight selected stationary equations exported, not a rank claim.")
    print("PASS: static dimensional split, two-channel response, exact gauge/family controls.")
    print("PASS: angular moment coefficients agree with independent quadrature; off-shell control rejected.")
    print("No radial DAE solution, admissible non-gauge seed, or nonlinear obstruction established.")


if __name__ == "__main__":
    verify()
