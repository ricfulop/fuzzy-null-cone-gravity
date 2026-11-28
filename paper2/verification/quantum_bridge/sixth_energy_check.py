"""Tracked pilot of all twenty terms in the reduced sixth-order energy.

The pilot tests new components and records production costs. A completed
number remains an assembled fixed-mode coefficient until its independent
mixed-functional reference has also been evaluated. It is not a loop.
"""

from datetime import datetime, timezone
from decimal import Decimal as D
from fractions import Fraction as F
from math import nextafter
from pathlib import Path
import json
import os
import platform
import signal
import sys
import time

from run import report
from reduced_gaunt_exact import five_polarizations
from reduced_gaunt_vertices import tensor_from_labels, test_cases
from reduced_sixth_sources import prepare_sixth_sources
from reduced_sixth_local import local_families, merge_families, prepare_local_fields, triple_tail
from reduced_sixth_energy import prepare_sixth_data, integrate_sixth, physical_indicator
from gaussian_fourier_bilinears import compile_scalar_bilinear, compile_york_bilinear
from sixth_canonical_terms import canonical_terms, canonical_controls
from sixth_local_controls import run_local_checks
from sixth_source_controls import run_source_checks
from triple_integration_controls import run_triple_integration_controls

ROOT = Path(__file__).resolve().parent


def _directions():
    cases = dict(test_cases())
    axis, zero, _, _ = cases["electric_axis"]
    _, xz, _, _ = cases["noncommuting_labels"]
    return {"pure_momentum_axis": (zero, zero, axis, zero),
            "monopole_mixed": (axis, zero, zero, xz)}


def _basis():
    zero = [[F(0)]*3 for _ in range(3)]
    labels = five_polarizations()
    return [tensor_from_labels(a, zero) for a in labels]+[tensor_from_labels(zero, a) for a in labels]


def _indicator_checks():
    rows = []
    for name, values, factor, target in (
        ("scale_before_conversion", ["1e-400"], "1e100", D("1e-300")),
        ("retain_positive_underflow", ["1e-1000"], 1, D("1e-1000")),
        ("sum_before_conversion", ["1e-400", "2e-400"], "1e100", D("3e-300")),
        ("exact_zero", ["0"], 1, D(0)),
    ):
        actual = physical_indicator(values, factor)
        downward = D.from_float(nextafter(actual, 0.0)) if actual else D(0)
        passed = D.from_float(actual) >= target and (downward <= target if target else actual == 0)
        rows.append({"name": name, "value": actual, "target_decimal": str(target), "passed": passed})
    try:
        physical_indicator(["-1e-1000"])
    except ArithmeticError:
        rejected = True
    else:
        rejected = False
    rows.append({"name": "negative_indicator_rejected", "passed": rejected})
    return {"rows": rows, "passed": all(row["passed"] for row in rows)}


def main(config):
    start = time.monotonic()
    names = config["energy_directions"]
    if names != ["pure_momentum_axis", "monopole_mixed"]:
        raise ValueError("The pilot requires the pure-momentum control before the mixed production direction")
    deadline = datetime.fromisoformat(config["authorization_deadline_utc"].replace("Z", "+00:00"))
    remaining = int((deadline-datetime.now(timezone.utc)).total_seconds())
    wall = min(int(config["wall_limit_seconds"]), remaining)
    if wall <= 0:
        raise ValueError("No positive authorized window remains")
    def timeout(signum, frame):
        raise TimeoutError("Declared local H6 pilot wall-clock limit reached")
    output = ROOT/"output/sixth_energy_results.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {"configuration": config, "source_commit": os.getenv("ORX_COMMIT"),
              "run": os.getenv("ORX_RUN_ID"), "python": sys.version, "platform": platform.platform(),
              "status": "in progress", "directions": [], "triple_control_checkpoints": [],
              "scope": "Complete reduced H6 numerical pilot; independent mixed-functional validation and quantum goal remain open"}
    def checkpoint():
        result["wall_seconds"] = time.monotonic()-start
        temp = output.with_suffix(".tmp")
        temp.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        temp.replace(output)
    def gate(passed, label):
        checkpoint()
        if not passed:
            raise AssertionError(label+" failed; recorded measurements retained in result JSON")
    def progress(row):
        print("SIXTH PILOT PROGRESS "+json.dumps(row), flush=True)
        if "radius_result" in row or "control_result" in row:
            result["triple_control_checkpoints"].append(row)
            checkpoint()
        if "control_result" in row:
            completed = sum("control_result" in item for item in result["triple_control_checkpoints"])
            report(step=completed, m={"H6_controls_completed": completed})
    print("CONFIGURATION "+json.dumps(result), flush=True)
    try:
        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(wall)
        quartic = json.loads((ROOT.parents[1]/"reports/quantum-bridge/data/quartic-loop-results.json").read_text())
        higher = json.loads((ROOT.parents[1]/"reports/quantum-bridge/data/higher-geometry-results.json").read_text())
        if quartic["source_commit"] != "a078aca003f0bb3878c0edf992489424866b0532":
            raise ValueError("Unexpected archived quartic source")
        if higher["source_commit"] != "2a59f1b60f630bd2b4e4d922dc693a7810ab2608":
            raise ValueError("Unexpected archived independent momentum reference")
        report(phase="triple_integral_controls", step=0, of=9, unit="controls")
        result["indicator_controls"] = _indicator_checks()
        gate(result["indicator_controls"]["passed"], "Positive error-indicator conversion")
        result["triple_controls"] = run_triple_integration_controls(config, progress=progress)
        print("TRIPLE INTEGRAL CONTROLS "+json.dumps(result["triple_controls"]), flush=True)
        gate(result["triple_controls"]["passed"], "Triple integral controls")
        report(step=result["triple_controls"]["control_count"])
        basis, direction_labels = _basis(), _directions()
        for name in names:
            current = {"name": name, "accepted_numerical_gates": False}
            result["directions"].append(current)
            labels = direction_labels[name]
            current["labels_Eq_Bq_Es_Bs"] = [[[str(x) for x in row] for row in a] for a in labels]
            def case_progress(row):
                progress({"direction": name, **row})
            report(phase="sixth_preparation", step=0, of=1, unit="directions")
            data = prepare_sixth_data(tensor_from_labels(*labels[:2]), tensor_from_labels(*labels[2:]), progress=case_progress)
            sources = prepare_sixth_sources(data)
            current["source_statistics"] = sources["statistics"]
            prepared_local = prepare_local_fields(data, sources)
            current["sources"] = run_source_checks(data, sources, config["source_points"], config["source_scaled_tolerance"])
            current["local_point_checks"] = run_local_checks(data, sources, prepared_local,
                config["source_points"], scaled_tolerance=config["source_scaled_tolerance"])
            print("SIXTH POINT CHECKS "+json.dumps({"name": name, "sources": current["sources"],
                  "local": current["local_point_checks"]}), flush=True)
            gate(current["sources"]["passed"] and current["local_point_checks"]["passed"], "Pointwise "+name)
            report(step=1)

            report(phase="canonical_terms", step=0, of=10, unit="components")
            def canonical_progress(row):
                case_progress(row)
                if row["stage"] == "canonical_component_complete":
                    report(step=row["completed"])
            canonical = canonical_terms(data, basis, labels, precision=96,
                progress=canonical_progress, **config["canonical_compile_caps"])
            current["canonical"] = canonical
            checkpoint()
            report(phase="canonical_checks", step=0, of=10, unit="components")
            def canonical_check_progress(row):
                case_progress(row)
                if row["stage"] == "canonical_control_component_complete":
                    report(step=row["completed"])
            current["canonical_checks"] = canonical_controls(data, basis, canonical, labels,
                quartic_coefficients=quartic["invariant_coefficients"], progress=canonical_check_progress,
                **config["canonical_compile_caps"])
            print("CANONICAL TERMS "+json.dumps({"name": name, "result": canonical,
                  "checks": current["canonical_checks"]}), flush=True)
            gate(current["canonical_checks"]["passed"], "Canonical "+name)

            report(phase="nested_compilation", step=0, of=2, unit="kernels")
            nested = {}
            current["nested_kernel_statistics"] = {}
            for index, key in enumerate(("scalar_nested", "york_nested"), 1):
                t0 = time.monotonic()
                if key == "scalar_nested":
                    kernel = compile_scalar_bilinear(sources["K"], sources["F"],
                        progress=case_progress, **config["nested_compile_caps"]).scale(-1)
                else:
                    kernel = compile_york_bilinear(sources["J_L"], sources["J_R"],
                        progress=case_progress, **config["nested_compile_caps"]).scale(-2)
                nested[key] = kernel
                current["nested_kernel_statistics"][key] = {**kernel.statistics(), "wall_seconds": time.monotonic()-t0}
                print("SIXTH COMPILED KERNEL "+json.dumps({"name": name, "family": key,
                      **current["nested_kernel_statistics"][key]}), flush=True)
                checkpoint()
                report(step=index)

            report(phase="sixth_local_assembly", step=0, of=15, unit="families")
            def local_progress(row):
                case_progress(row)
                report(step=row["completed"])
            local = local_families(data, sources, prepared=prepared_local,
                max_pair_products=config["max_pair_products"], progress=local_progress)
            if name == "pure_momentum_axis":
                def signature(accumulator, factor):
                    return {(group.beta, group.rates, key): factor*coefficient
                            for group in accumulator.collapsed()
                            for key, coefficient in group.coefficients.items() if coefficient}
                lhs, rhs = signature(local["scalar_B2"], 4), signature(local["scalar_s2_u2"], 7)
                current["pure_momentum_bin_identity"] = {"identity": "4*scalar_B2 = 7*scalar_s2_u2",
                    "left_bin_count": len(lhs), "right_bin_count": len(rhs), "passed": lhs == rhs and bool(lhs)}
                gate(current["pure_momentum_bin_identity"]["passed"], "Exact pure-momentum 49:28 ratio")
            report(phase="sixth_triple_assembly", step=0, of=6, unit="components")
            triple, triple_stats = triple_tail(data, max_pair_products=config["max_pair_products"], progress=local_progress)
            current["assembly_statistics"] = {"local": {key: value.statistics for key, value in local.items()},
                                               "triple_components": triple_stats}
            checkpoint()
            assembly = {"local": local, "local_merged": merge_families(local), "nested": nested,
                        "triple": triple, "canonical": canonical}
            report(phase="sixth_integration", step=0, of=20, unit="families")
            current["completed_ledger"] = []
            def integration_progress(row):
                case_progress(row)
                if row["stage"] == "sixth_energy_family":
                    current["completed_ledger"].append(row)
                    checkpoint()
                    report(row="energy_terms", m={"direction": name, "family": row["family"],
                           "energy": row["physical_value"], "error_indicator": row["absolute_error_indicator"],
                           "outcome": "correct" if row["converged"] else "wrong"})
                    report(step=row["completed"], m={"H6_families_completed": row["completed"]})
            current["energy"] = integrate_sixth(assembly, local_options=config["local_options"],
                fourier_options=config["fourier_options"], triple_options=config["triple_options"],
                absolute_tolerance=config["total_absolute_tolerance"], relative_tolerance=config["total_relative_tolerance"],
                roundoff_safety_factor=config["roundoff_safety_factor"], progress=integration_progress)
            gate(current["energy"]["converged"], "Energy integration "+name)
            if name == "pure_momentum_axis":
                reference_level = higher["axis_pure_momentum_sixth_control"]["quadratures"][-1]
                target = 77*reference_level["radial_integrals"]["source_times_u2_squared"]
                actual = current["energy"]["value"]
                error = abs(actual-target)/max(1,abs(target))
                nonzero = {row["family"] for row in current["energy"]["ledger"] if row["physical_value"] != 0}
                values = {row["family"]: row["physical_value"] for row in current["energy"]["ledger"]}
                family_checks = []
                for family, coefficient in (("scalar_B2", 49), ("scalar_s2_u2", 28)):
                    expected = coefficient*reference_level["radial_integrals"]["source_times_u2_squared"]
                    family_error = abs(values[family]-expected)/max(1, abs(expected))
                    family_checks.append({"family": family, "coefficient": coefficient,
                        "actual": values[family], "reference": expected, "scaled_error": family_error,
                        "passed": family_error <= config["pure_momentum_tolerance"]})
                current["independent_momentum_reference"] = {"source_commit": higher["source_commit"],
                    "reference": target, "actual": actual, "scaled_error": error,
                    "reference_method": "77 times independently evaluated radial source_times_u2_squared",
                    "archived_heat_reference": reference_level["pure_p_H6"],
                    "archived_scaled_refinement": higher["axis_pure_momentum_sixth_control"]["refinement_scaled_error"],
                    "archived_tail_bound": higher["axis_pure_momentum_sixth_control"]["pure_p_H6_tail_bound"],
                    "family_checks": family_checks,
                    "nonzero_families": sorted(nonzero),
                    "passed": error <= config["pure_momentum_tolerance"] and all(row["passed"] for row in family_checks)
                        and nonzero == {"scalar_B2", "scalar_s2_u2"}}
                gate(current["independent_momentum_reference"]["passed"], "Independent momentum reference")
            current["accepted_numerical_gates"] = True
            print("SIXTH ENERGY PILOT "+json.dumps(current), flush=True)
            checkpoint()
        report(phase="export", step=0, of=1, unit="reports")
        result["status"] = "all specified H6 pilot numerical gates passed; independent mixed-functional validation remains pending"
        result["passed"] = True
        result["full_four_point_loop_evaluated"] = False
        result["independent_mixed_functional_validation"] = False
        checkpoint()
        report(step=1, m={"H6_mixed_value": result["directions"][-1]["energy"]["value"],
                         "H6_total_error_indicator": result["directions"][-1]["energy"]["total_absolute_error_indicator"]})
        print("FINAL "+json.dumps({"result": str(output), "status": result["status"],
              "source_commit": result["source_commit"], "wall_seconds": result["wall_seconds"]}), flush=True)
    except Exception as exc:
        result["status"], result["passed"] = "failed", False
        result["failure"] = {"type": type(exc).__name__, "message": str(exc)[:2000]}
        checkpoint()
        raise
    finally:
        signal.alarm(0)
