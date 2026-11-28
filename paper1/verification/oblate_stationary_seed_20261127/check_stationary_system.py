"""PROPOSED local extraction/control run. Do not execute without approval.

This exports the exact stationary radial equations and the discarded response.
It does not solve the radial DAE, provide an admissible non-gauge seed, or
evaluate a nonlinear obstruction.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import os
import platform
import subprocess
import sys
import time

import numpy as np
import sympy as s

import stationary_system as st
import extraction_progress as progress

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/oblate-stationary-system-recheck-2026-11-27/outputs"


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


def quadrature_controls(data, r, x, epsilon, recorder):
    """Independent numerical angular integration of every exported density."""
    maxima = {}
    count = 0
    for epsilon_value in [s.Rational(1, 5), s.Rational(1, 100)]:
        point = {r: s.Rational(3, 2), epsilon: epsilon_value}
        for node_count in [48, 96]:
            recorder.event("quadrature_case_begin", epsilon=str(epsilon_value), nodes=node_count)
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
            recorder.event("quadrature_case_end", epsilon=str(epsilon_value), nodes=node_count,
                           maximum_scaled_error=largest)
    return {"nonzero_density_comparisons": count, "maximum_scaled_errors": maxima,
            "tolerance": 1e-11, "r": "3/2",
            "purpose": "angular coefficient control; radius is not frozen in the exported operators"}


def verify(output_directory=OUT):
    if not __debug__:
        raise RuntimeError("Scientific verification requires Python assertions enabled")
    start = time.monotonic()
    t, x, phi = s.symbols("t x phi", real=True)
    r, epsilon = s.symbols("r epsilon", positive=True)
    coords = [t, r, x, phi]
    sources = [Path(__file__), Path(st.__file__), Path(st.n.__file__),
               Path(progress.__file__), ROOT/"sections/nonlinear_action.tex"]
    identity = {
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "python": platform.python_version(), "python_implementation": sys.version,
        "dependencies": {name: importlib.metadata.version(name) for name in ["numpy", "sympy"]},
        "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "coordinates": [str(value) for value in coords],
        "symbol_assumptions": {str(value): value.assumptions0 for value in [*coords, epsilon]},
        "parameter_domain": "r>0, epsilon>0; x=cos(theta)",
        "normalization": "(16 pi G)/(4 pi) I; dx/2 angular average; Sigma occurs once",
        "column_labels": st.COLUMNS,
        "retained_row_labels": sum(st.n.coefficient_labels(), []),
        "source_scope": "original Paper I action; stationary N=3,m=0,even/equatorially even columns",
        "command_arguments": sys.argv,
        "thread_environment": {name: os.environ.get(name) for name in (
            "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")},
    }
    with progress.ExtractionRecorder(output_directory, identity) as recorder:
        with recorder.stage("source_snapshots"):
            for path in sources:
                snapshot = recorder.directory / "sources" / path.relative_to(ROOT)
                snapshot.parent.mkdir(parents=True, exist_ok=True)
                content = path.read_bytes()
                expected = identity["sha256"][str(path.relative_to(ROOT))]
                if hashlib.sha256(content).hexdigest() != expected:
                    raise RuntimeError("Source changed before snapshot was saved")
                with snapshot.open("xb") as stream:
                    stream.write(content)
                if hashlib.sha256(snapshot.read_bytes()).hexdigest() != expected:
                    raise RuntimeError("Saved source snapshot does not match identity")
        return verify_calculation(coords, epsilon, recorder, start, identity, sources)


def verify_calculation(coords, epsilon, recorder, start, identity, sources):
    _, r, x, _ = coords
    with recorder.stage("two_channel_response_control"):
        check_cross_formula(coords, epsilon)
    with recorder.stage("all_columns_assembly"):
        data = st.assemble(coords, epsilon, recorder=recorder)
    operators = data["operators"]
    f = s.Function("f")(r)
    with recorder.stage("radial_gauge_control"):
        exact_zero(st.apply_operator(operators, st.radial_gauge_profiles(r, epsilon, f), r),
                   "radial gauge including response")
    time_rescaling = [0]*8
    time_rescaling[6] = 1
    with recorder.stage("constant_lapse_control"):
        exact_zero(st.apply_operator(operators, time_rescaling, r), "constant lapse control")
    parameter = st.oblate_parameter_profiles(r, epsilon)
    with recorder.stage("oblate_parameter_control"):
        fields = st.n.coefficient_fields(dict(zip(st.COLUMNS, parameter)), x)
        exact_zero(st.n.linear_metric(fields, coords, epsilon)
                   - st.n.background(coords, epsilon).diff(epsilon), "oblate parameter tangent")
        exact_zero(st.apply_operator(operators, parameter, r), "oblate family tangent")
    # This deliberately off-shell profile must not pass the retained rows.
    invalid = [0]*8
    invalid[6] = r*r
    with recorder.stage("off_shell_rejection_control"):
        invalid_rows = st.apply_operator(operators, invalid, r)[:-1, :]
        point = {r: s.Rational(3, 2), epsilon: s.Rational(1, 5)}
        magnitudes = [abs(float(value.subs(point).evalf(60))) for value in invalid_rows]
        if not all(np.isfinite(value) for value in magnitudes) or max(magnitudes) <= 1e-8:
            raise AssertionError("off-shell lapse control was not rejected")
    with recorder.stage("independent_quadrature_controls"):
        quadrature = quadrature_controls(data, r, x, epsilon, recorder)
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
    with recorder.stage("final_expression_serialization"):
        output["matrices_by_derivative_order"] = [
            [[str(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]
            for matrix in operators]
    after = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    if after != identity["sha256"]:
        raise RuntimeError("Source changed during scientific invocation")
    manifest = {
        **identity,
        "seconds_before_final_json_writes": round(time.monotonic()-start, 3),
        "sha256_after": after,
        "checkpoints": recorder.checkpoints,
        "scientific_invocations": 1, "GPU_jobs": 0, "remote_jobs": 0, "tracked_experiment_launches": 0,
    }
    with recorder.stage("final_export"):
        recorder.write_json("stationary-operators.json", output)
        recorder.write_json("stationary-reproducibility.json", manifest)
        recorder.write_json("completion.json", {
            "state": "all checker controls passed; radial DAE unsolved",
            "identity_sha256": recorder.identity_sha256,
            "column_count": len(st.COLUMNS), "retained_row_count": len(data["rows"]),
            "all_checker_controls_passed": True,
            "source_sha256_before": identity["sha256"], "source_sha256_after": after,
            "checkpoints": recorder.checkpoints,
            "export_sha256": {name: hashlib.sha256((recorder.directory/name).read_bytes()).hexdigest()
                              for name in ("stationary-operators.json", "stationary-reproducibility.json")},
            "admissible_seed_certified": False, "nonlinear_consistency_established": False,
        })
    print("PASS: all 22 retained rows checked; eight selected stationary equations exported, not a rank claim.")
    print("PASS: static dimensional split, two-channel response, exact gauge/family controls.")
    print("PASS: angular moment coefficients agree with independent quadrature; off-shell control rejected.")
    print("No radial DAE solution, admissible non-gauge seed, or nonlinear obstruction established.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT,
                        help="Fresh invocation directory; existing directories are refused")
    verify(parser.parse_args().output_dir)
