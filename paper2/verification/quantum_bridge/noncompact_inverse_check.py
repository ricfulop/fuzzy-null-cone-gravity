"""Tracked bounded validation of noncompact first inverses and H6 algebra.

This run does not integrate a mixed H6 functional. The fixed experiment
recipe dispatches here after extracting a committed snapshot. Every source,
level, failed control and adapter diagnostic is checkpointed before advancing.
"""

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import platform
import resource
import signal
import sys
import time

from run import report
from sixth_unreduced_controls import run_unreduced_algebra_controls
from gaussian_radial_controls import run_gaussian_radial_controls
from noncompact_york_controls import run_noncompact_york_controls
from noncompact_first_response import run_first_response_checks
from spherical_york_reference import WorkBudget


ROOT = Path(__file__).resolve().parent


def main(config):
    started = time.monotonic()
    output = ROOT/"output/noncompact_inverse_results.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {"source_commit": os.environ.get("ORX_COMMIT"), "run": os.environ.get("ORX_RUN_ID"),
              "configuration": config, "python": sys.version, "platform": platform.platform(),
              "status": "preparing", "suites": {}, "checkpoints": [],
              "accuracy_certified": False, "independent_mixed_functional_validation": False,
              "four_point_evaluated": False,
              "scope": "Noncompact inverse adapters, manufactured unreduced H6 algebra, and actual first gravitational responses only"}
    counts = {"required_completed": 0, "required_failed": 0, "diagnostic_failed": 0,
              "algebra_groups": 0, "gaussian_groups": 0, "synthetic_levels": 0, "gravitational_levels": 0}
    phase = "preparing"
    memory_cap = None
    old_handler = signal.getsignal(signal.SIGALRM)

    def checkpoint():
        result["wall_seconds"] = time.monotonic()-started
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result["peak_rss_bytes"] = int(rss if sys.platform == "darwin" else 1024*rss)
        result["counts"] = dict(counts)
        temporary = output.with_suffix(".tmp")
        temporary.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        temporary.replace(output)

    def progress(event):
        # Snapshot before retaining; later caller mutations cannot rewrite
        # the meaning of a checkpoint that has already appeared in the log.
        event = json.loads(json.dumps(event, allow_nan=False))
        stage = event.get("stage", "unknown")
        nested = event.get("adapter_event", {})
        full = ("control_result" in event or "row" in event or "level" in event
                or "level_result" in event or "result" in event
                or "source" in event or "sources" in event or "diagnostics" in event
                or stage.endswith(("preflight", "targets", "result"))
                or "result" in nested or "source" in nested or "base" in nested
                or nested.get("stage", "").endswith(("gate", "result")))
        if full:
            result["checkpoints"].append({"suite": phase, **event})
        row = event.get("control_result", event.get("row"))
        if row is not None and "passed" in row:
            required = row.get("acceptance_required", True)
            if required:
                counts["required_completed"] += 1
                counts["required_failed"] += int(not row["passed"])
            else:
                counts["diagnostic_failed"] += int(not row["passed"])
            report(m={"noncompact_required_completed": counts["required_completed"],
                      "noncompact_required_failed": counts["required_failed"],
                      "noncompact_diagnostic_failed": counts["diagnostic_failed"]})
            report(row="noncompact_controls", m={"suite": phase,
                "control": row.get("name", row.get("check", event.get("control", "comparison"))),
                "outcome": "pass" if row["passed"] else "fail",
                "required": "yes" if required else "diagnostic"})
        if stage == "unreduced_algebra_control_complete":
            counts["algebra_groups"] += 1
            report(step=counts["algebra_groups"])
        elif phase == "gaussian_exterior_controls" and "control_result" in event and "exception" not in stage:
            counts["gaussian_groups"] += 1
            report(step=counts["gaussian_groups"])
        elif phase == "noncompact_synthetic" and stage == "level_result":
            counts["synthetic_levels"] += 1
            report(step=counts["synthetic_levels"])
        elif stage == "first_response_level_result":
            counts["gravitational_levels"] += 1
            report(step=counts["gravitational_levels"])
        print("NONCOMPACT PROGRESS "+json.dumps({"suite": phase, **event}, allow_nan=False), flush=True)
        checkpoint()
        if memory_cap is not None and result["peak_rss_bytes"] > memory_cap:
            raise MemoryError("Declared sampled memory ceiling exceeded")

    def gate(passed, message):
        checkpoint()
        if not passed:
            raise AssertionError(message+"; completed evidence retained")

    def timeout(signum, frame):
        raise TimeoutError("Declared local noncompact inverse wall-clock limit reached")

    try:
        if config.get("calculation") != "noncompact_inverse_validation":
            raise ValueError("Unexpected calculation selector")
        wall, memory_cap = config["wall_limit_seconds"], config["memory_limit_bytes"]
        if type(wall) is not int or wall <= 0 or type(memory_cap) is not int or memory_cap <= 0:
            raise ValueError("Wall and memory caps must be positive integers")
        deadline = datetime.fromisoformat(config["authorization_deadline_utc"].replace("Z", "+00:00"))
        remaining = int((deadline-datetime.now(timezone.utc)).total_seconds())
        wall = min(wall, remaining)
        if wall <= 0:
            raise ValueError("No positive authorized execution window remains")
        caps = config["shared_work_caps"]
        expected_caps = {"max_fd_unknowns", "max_point_evaluations", "max_exact_products", "max_stored_values",
                         "max_cartesian_calls", "max_projection_products"}
        if set(caps) != expected_caps or any(type(value) is not int or value <= 0 for value in caps.values()):
            raise ValueError("Every cumulative adapter work cap must be declared")
        budget = WorkBudget(**caps)
        result["source_hashes"] = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in (
            "noncompact_inverse_check.py", "noncompact_first_response.py", "noncompact_york_controls.py",
            "gaussian_radial_controls.py", "gaussian_radial_reference.py", "spherical_york_reference.py",
            "spherical_green_reference.py", "solid_harmonic_reference.py", "kelvin_exterior_reference.py",
            "sixth_unreduced_algebra.py", "sixth_unreduced_controls.py", "higher_local_geometry.py",
            "gaussian_green_fields.py")}
        result["effective_wall_limit_seconds"] = wall
        checkpoint()
        print("CONFIGURATION "+json.dumps(result, allow_nan=False), flush=True)
        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(wall)

        phase = "unreduced_algebra"
        report(phase=phase, step=0, of=18, unit="controls")
        result["suites"][phase] = run_unreduced_algebra_controls(config, progress)
        gate(result["suites"][phase]["passed"] and counts["algebra_groups"] == 18,
             "Unreduced algebra controls failed")

        phase = "gaussian_exterior_controls"
        report(phase=phase, step=0, of=10, unit="controls")
        result["suites"][phase] = run_gaussian_radial_controls(config, progress)
        report(step=counts["gaussian_groups"], of=counts["gaussian_groups"])
        gate(result["suites"][phase]["passed"] and counts["gaussian_groups"] == 10, "Gaussian exterior controls failed")

        phase = "noncompact_synthetic"
        report(phase=phase, step=0, of=9, unit="levels")
        result["suites"][phase] = run_noncompact_york_controls(config, budget, progress)
        gate(result["suites"][phase]["passed"] and counts["synthetic_levels"] == 9,
             "Corrected synthetic inverse controls failed")

        phase = "mixed_first_responses"
        report(phase=phase, step=0, of=4, unit="levels")
        result["suites"][phase] = run_first_response_checks(config, budget, progress)
        gate(result["suites"][phase]["passed"] and counts["gravitational_levels"] == 4,
             "Actual gravitational first-response controls failed")

        phase = "final_checks"
        report(phase=phase, step=0, of=1, unit="reports")
        result["work"] = budget.snapshot()
        gate(counts["required_failed"] == 0 and all(suite["passed"] for suite in result["suites"].values()),
             "Final required comparison aggregation failed")
        result["passed"] = True
        result["status"] = "All declared noncompact adapter and actual first-response gates passed; mixed H6 remains pending"
        checkpoint()
        report(step=1, m={"noncompact_all_required_passed": 1,
                          "noncompact_required_completed": counts["required_completed"],
                          "noncompact_required_failed": 0})
        print("FINAL "+json.dumps({"result": str(output), "source_commit": result["source_commit"],
            "status": result["status"], "counts": counts, "wall_seconds": result["wall_seconds"],
            "peak_rss_bytes": result["peak_rss_bytes"]}), flush=True)
    except Exception as exc:
        result["passed"] = False
        result["status"] = "failed: "+type(exc).__name__+": "+str(exc)
        result["failed_phase"] = phase
        checkpoint()
        report(m={"noncompact_all_required_passed": 0})
        print("FAILED "+json.dumps({"result": str(output), "status": result["status"]}), flush=True)
        raise
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
