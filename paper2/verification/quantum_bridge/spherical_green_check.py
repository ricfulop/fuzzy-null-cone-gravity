"""Tracked independent inverse-primitive validation, with durable checkpoints.

No gravitational source or H6 energy is evaluated here. All modes, joins,
origin patches, exterior data, and resource limits are declared before the
first numerical control. The fixed project entry point dispatches here.
"""

from datetime import datetime, timezone
from math import isfinite
from pathlib import Path
import json
import os
import platform
import resource
import signal
import sys
import time

from run import report
from spherical_angular_controls import run_angular_controls
from spherical_scalar_controls import run_scalar_controls
from spherical_york_controls import run_york_controls


ROOT = Path(__file__).resolve().parent


def _flat(values):
    if isinstance(values, (list, tuple)):
        return [item for value in values for item in _flat(value)]
    value = float(values)
    if not isfinite(value):
        raise ArithmeticError("Nonfinite cross-slice value")
    return [value]


def _compare_slices(reference, changed, settings, kind):
    base = {row["case"]: row for row in reference["levels"]
            if row["intervals"] == reference["settings"]["intervals"][-1]}
    rows = []
    for level in changed["levels"]:
        old = base[level["case"]]
        if "points" not in level or "points" not in old:
            rows.append({"case": level["case"], "slice": kind,
                         "passed": False, "reason": "incomplete source or field solve"})
            continue
        common = sorted(set(old["points"]) & set(level["points"]))
        if not common:
            raise ValueError("No shared physical points for cross-slice comparison")
        for key in ("W", "dW", "ddW", "chi", "div_V", "V"):
            error = max(abs(a-b)/max(1.0, abs(a), abs(b)) for label in common
                        for a, b in zip(_flat(level["points"][label][key]), _flat(old["points"][label][key])))
            tolerance = settings["york_component_tolerance"] if key in ("dW", "ddW", "div_V") else settings["york_value_tolerance"]
            rows.append({"case": level["case"], "slice": kind, "quantity": key,
                         "shared_points": common, "scaled_change": error,
                         "tolerance": tolerance, "passed": error <= tolerance})
        # Exterior norm and flux have R-dependent targets. Compare them across
        # sphere quadrature only; their own exact targets gate every R slice.
        keys = tuple(level["integrals"]) if kind == "angular" else ("source_pairing", "norm_total")
        for key in keys:
            actual, previous = level["integrals"][key], old["integrals"][key]
            error = abs(actual-previous)/max(1.0, abs(actual), abs(previous))
            tolerance = settings["york_integral_tolerance"]
            rows.append({"case": level["case"], "slice": kind, "quantity": key,
                         "actual": actual, "reference_slice": previous, "scaled_change": error,
                         "tolerance": tolerance, "passed": error <= tolerance})
    return {"rows": rows, "passed": bool(rows) and all(row["passed"] for row in rows)}


def main(config):
    started = time.monotonic()
    output = ROOT/"output/spherical_green_results.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {"source_commit": os.environ.get("ORX_COMMIT"), "run": os.environ.get("ORX_RUN_ID"),
              "configuration": config, "python": sys.version, "platform": platform.platform(),
              "status": "preparing", "checkpoints": [], "suites": {}, "cross_slices": [],
              "scope": "Independent scalar/angular/York inverse primitives; full mixed H6 and four-point loop remain pending",
              "independent_mixed_functional_validation": False, "four_point_evaluated": False,
              "accuracy_certified": False}
    counts = {"required_completed": 0, "required_failed": 0, "diagnostic_failed": 0,
              "york_levels": 0, "scalar_controls": 0, "angular_controls": 0}
    phase = None
    memory_cap = None
    old_handler = signal.getsignal(signal.SIGALRM)

    def checkpoint():
        result["wall_seconds"] = time.monotonic()-started
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result["peak_rss_bytes"] = int(rss if sys.platform == "darwin" else 1024*rss)
        result["counts"] = dict(counts)
        temp = output.with_suffix(".tmp")
        temp.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        temp.replace(output)

    def progress(event):
        nonlocal phase
        stage = event.get("stage")
        full = ("control_result" in event or "level_result" in event or "row" in event
                or stage in ("source_result", "level_result"))
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
            report(m={"green_controls_completed": counts["required_completed"],
                      "green_required_failures": counts["required_failed"],
                      "green_coarse_diagnostic_failures": counts["diagnostic_failed"]})
            report(row="inverse_controls", m={"suite": phase,
                "control": row.get("name", row.get("check", event.get("control", "comparison"))),
                "outcome": "correct" if row["passed"] else "wrong",
                "required": "yes" if required else "diagnostic"})
        if stage == "scalar_control_complete":
            counts["scalar_controls"] += 1
            report(step=counts["scalar_controls"])
        if stage == "angular_control_completed":
            counts["angular_controls"] += 1
            report(step=counts["angular_controls"])
        if stage == "level_result":
            counts["york_levels"] += 1
            report(step=counts["york_levels"])
        print("SPHERICAL PROGRESS "+json.dumps(event, allow_nan=False), flush=True)
        if full:
            checkpoint()
        if memory_cap is not None:
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            rss_bytes = rss if sys.platform == "darwin" else rss*1024
            if rss_bytes > memory_cap:
                raise MemoryError("Independent inverse controls exceeded sampled peak RSS budget")

    def gate(condition, message):
        checkpoint()
        if not condition:
            raise AssertionError(message+"; complete measurements retained")

    def timeout(signum, frame):
        raise TimeoutError("Declared independent-inverse wall-clock budget reached")

    try:
        if config["calculation"] != "independent_spherical_green":
            raise ValueError("Wrong calculation dispatch")
        deadline = datetime.fromisoformat(config["authorization_deadline_utc"].replace("Z", "+00:00"))
        remaining = int((deadline-datetime.now(timezone.utc)).total_seconds())
        wall = min(int(config["wall_limit_seconds"]), remaining)
        memory_cap = int(config["memory_limit_bytes"])
        if wall <= 0 or memory_cap <= 0:
            raise ValueError("No positive authorized wall/memory budget")
        scalar = config["scalar_intervals"]
        york = config["york_intervals"]
        if scalar != [256, 512, 1024] or york != [128, 256, 512]:
            raise ValueError("This initial protocol pins the preregistered radial levels")
        if config["angular_ell_max"] != 10:
            raise ValueError("This initial protocol declares 26 angular controls through ell=10")
        point_cap = config["york_global_point_cap"]
        if type(point_cap) is not int or point_cap <= 0:
            raise ValueError("york_global_point_cap must be a positive integer")
        if (config["york_radius"] != 4 or config["york_origin_radius"] != .5
                or config["york_n_mu"] != 8 or config["york_n_phi"] != 16
                or config["york_tail_order"] != 12):
            raise ValueError("Unexpected baseline York grid or origin/tail convention")
        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(wall)
        result["effective_wall_seconds"] = wall
        result["status"] = "in progress"
        print("CONFIGURATION "+json.dumps(result, allow_nan=False), flush=True)
        checkpoint()
        phase = "angular_controls"
        report(phase=phase, step=0, of=26, unit="controls")
        result["suites"][phase] = run_angular_controls(config, progress=progress)
        gate(result["suites"][phase]["passed"], "Angular/exterior algebra controls failed")
        phase = "scalar_controls"
        report(phase=phase, step=0, of=9, unit="controls")
        result["suites"][phase] = run_scalar_controls(config, progress=progress)
        gate(result["suites"][phase]["passed"], "Independent scalar inverse controls failed")
        phase = "york_controls"
        report(phase=phase, step=0, of=12, unit="cases")
        base_config = dict(config, york_max_point_evaluations=min(config["york_max_point_evaluations"], point_cap))
        result["suites"]["york_radial"] = run_york_controls(base_config, progress=progress)
        gate(result["suites"]["york_radial"]["passed"], "Two-stage York radial controls failed")
        point_total = result["suites"]["york_radial"]["work"]["point_evaluations"]
        slices = (("york_angular", {"york_intervals": [512], "york_n_mu": 12, "york_n_phi": 24}, "angular"),
                  ("york_R2", {"york_intervals": [256], "york_radius": 2}, "exterior_radius"),
                  ("york_R6", {"york_intervals": [768], "york_radius": 6}, "exterior_radius"))
        for name, change, kind in slices:
            options = dict(config, **change)
            remaining_points = point_cap-point_total
            if remaining_points <= 0:
                raise RuntimeError("Cumulative York point-evaluation cap exhausted")
            options["york_max_point_evaluations"] = min(options["york_max_point_evaluations"], remaining_points)
            def slice_progress(event):
                progress({"slice_name": name, **event})
            trial = run_york_controls(options, progress=slice_progress)
            result["suites"][name] = trial
            point_total += trial["work"]["point_evaluations"]
            result["total_york_point_evaluations"] = point_total
            comparison = _compare_slices(result["suites"]["york_radial"], trial, config, kind)
            result["cross_slices"].append({"name": name, **comparison})
            for row in comparison["rows"]:
                progress({"stage": "cross_slice_comparison", "slice_name": name, "row": row})
            gate(trial["passed"] and comparison["passed"], name+" comparison failed")
        gate(point_total <= point_cap, "Cumulative York work limit failed")
        gate(counts["required_failed"] == 0, "A mandatory callback comparison failed")
        phase = "export"
        report(phase=phase, step=0, of=1, unit="reports")
        result["status"] = "all specified independent inverse-primitive gates passed; mixed H6 comparison remains pending"
        result["passed"] = True
        checkpoint()
        report(step=1, m={"green_all_required_passed": 1,
                         "green_controls_completed": counts["required_completed"],
                         "green_required_failures": counts["required_failed"]})
        print("FINAL "+json.dumps({"result": str(output), "source_commit": result["source_commit"],
              "status": result["status"], "counts": counts, "wall_seconds": result["wall_seconds"],
              "peak_rss_bytes": result["peak_rss_bytes"]}, allow_nan=False), flush=True)
    except BaseException as exc:
        result["status"] = "failed: "+type(exc).__name__+": "+str(exc)
        result["passed"] = False
        checkpoint()
        print("FAILED "+json.dumps({"result": str(output), "status": result["status"]}), flush=True)
        raise
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


if __name__ == "__main__":
    main(json.loads((ROOT/"config.json").read_text()))
