"""Tracked checks of localized H6 sources, inverse bilinears and exterior tails.

This prepares and checks ingredients of the fixed ten-profile reduction.
It does not assemble the full H6 Hamiltonian or a four-point loop. Every
scientific calculation is inside main or a function it explicitly calls.
"""

from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import json
import os
import platform
import signal
import sys
import time

from run import report
from fourier_bilinear_controls import run_fourier_controls
from heat_tail_controls import run_tail_controls, check_gravitational_tail
from reduced_fifth_order import prepare_fields
from reduced_gaunt_exact import scale
from reduced_gaunt_vertices import test_cases, tensor_from_labels
from reduced_sixth_sources import prepare_sixth_sources
from sixth_source_controls import run_source_checks, evaluate_sources

ROOT = Path(__file__).resolve().parent


def _labels():
    cases = dict(test_cases())
    axis, zero, _, _ = cases["electric_axis"]
    _, xz, _, _ = cases["noncommuting_labels"]
    cases["monopole_mixed"] = (axis, zero, zero, xz)
    cases["pure_configuration"] = (axis, zero, zero, zero)
    cases["pure_momentum"] = (zero, zero, zero, xz)
    return cases


def _flat_sources(sources):
    entries = {"F": sources["F"], "K": sources["K"]}
    for name in ("r", "t", "J_L", "J_R", "A1W"):
        entries.update({name+"_"+str(i): source for i, source in enumerate(sources[name])})
    entries.update({"T1_"+str(i)+str(j): sources["T1"][i][j]
                    for i in range(3) for j in range(3)})
    return entries


def _sector_check(name, sources, points):
    checks = {}
    if name == "pure_configuration":
        checks["F_equals_minus_K"] = not sources["F"].add(sources["K"]).terms
        for key, source in _flat_sources(sources).items():
            if key not in ("F", "K"):
                checks[key+"_zero"] = not source.terms
    elif name == "pure_momentum":
        for key, source in _flat_sources(sources).items():
            if key in ("F", "K") or key.startswith(("r_", "J_R_", "A1W_", "T1_")):
                checks[key+"_zero"] = not source.terms
        for i in range(3):
            checks["J_L_equals_2t_"+str(i)] = not sources["J_L"][i].add(sources["t"][i], -2).terms
        # A nonempty representation alone would not establish a nonzero
        # function because different heat moments can obey recurrences.
        values = [evaluate_sources(sources["t"], point) for point in points]
        checks["t_nonzero_at_a_test_point"] = max(abs(x) for row in values for x in row) > 1e-12
        return {"checks": checks, "passed": all(checks.values()), "t_values": values,
                "interpretation": "The pure-momentum York term vanishes through J_R=0, not through J_L=0"}
    return {"checks": checks, "passed": all(checks.values()),
            "scope": "Specified exact sector identities; empty dictionaries are used only as sufficient zero tests"}


def _momentum_parity(data, sources, progress):
    negative_s = [[scale(poly, -1) for poly in row] for row in data["s"]]
    reversed_data = prepare_fields(data["q"], negative_s, geometry=data["geometry"], progress=progress)
    reversed_sources = prepare_sixth_sources(reversed_data)
    original, reversed_entries = _flat_sources(sources), _flat_sources(reversed_sources)
    checks = {}
    for name, source in original.items():
        parity = 1 if name in ("F", "K") else -1
        checks[name] = not reversed_entries[name].add(source, -parity).terms
    return {"passed": all(checks.values()), "checks": checks,
            "method": "Recompute York sources and first inverses from (Q,-S), reusing only Q geometry",
            "expected": "F,K even; r,t,J_L,J_R,A1W,T1 odd"}


def main(config):
    start = time.monotonic()
    names, tail_names = config["source_directions"], config["tail_directions"]
    required_sources = {"monopole_mixed", "noncommuting_labels", "pure_configuration", "pure_momentum"}
    required_tails = {"monopole_mixed", "noncommuting_labels"}
    if (len(set(names)) != len(names) or not required_sources.issubset(names)
            or len(set(tail_names)) != len(tail_names) or not required_tails.issubset(tail_names)
            or not set(tail_names).issubset(names)):
        raise ValueError("Required distinct source/tail directions must be retained, with tails a subset of sources")
    def timeout(signum, frame):
        raise TimeoutError("Declared local sixth-order backend wall limit reached")
    signal.signal(signal.SIGALRM, timeout)
    wall_limit = int(config["wall_limit_seconds"])
    deadline = datetime.fromisoformat(config["authorization_deadline_utc"].replace("Z", "+00:00"))
    remaining = int((deadline-datetime.now(timezone.utc)).total_seconds())
    if wall_limit <= 0 or remaining <= 0:
        raise ValueError("No positive authorized run window remains")
    signal.alarm(min(wall_limit, remaining))
    result = {"configuration": config, "python": sys.version, "platform": platform.platform(),
              "source_commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID"),
              "status": "in progress", "scope": "H6 ingredient validation; full Hamiltonian and quantum goal remain open"}
    output = ROOT / "output" / "sixth_backend_results.json"
    output.parent.mkdir(parents=True, exist_ok=True)

    def checkpoint():
        result["wall_seconds"] = time.monotonic()-start
        temporary = output.with_suffix(".tmp")
        temporary.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        temporary.replace(output)

    def gate(passed, name):
        checkpoint()
        if not passed:
            raise AssertionError(name+" failed; complete preceding measurements are in the result JSON")

    print("CONFIGURATION "+json.dumps(result), flush=True)
    try:
        report(phase="fourier_controls", step=0, of=9, unit="cases")
        result["fourier_control_rows"] = []
        def fourier_progress(row):
            print("FOURIER CONTROL "+json.dumps(row), flush=True)
            result["fourier_control_rows"].append(row)
            checkpoint()
            metrics = {"fourier_cases": row["completed"]}
            if "scaled_error" in row:
                metrics["fourier_control_error"] = float(row["scaled_error"])
            report(step=row["completed"], m=metrics)
        result["fourier_controls"] = run_fourier_controls(config, progress=fourier_progress)
        gate(result["fourier_controls"]["passed"], "Fourier controls")

        report(phase="tail_controls", step=0, of=1, unit="suites")
        result["tail_controls"] = run_tail_controls(config)
        print("TAIL CONTROLS "+json.dumps(result["tail_controls"]), flush=True)
        report(step=1, m={"tail_control_error": result["tail_controls"]["maximum_scaled_error"]})
        gate(result["tail_controls"]["passed"], "Tail controls")

        labels = _labels()
        prepared = {}
        result["prepared_sources"] = []
        report(phase="field_preparation", step=0, of=len(names), unit="directions")
        for index, name in enumerate(names, 1):
            eq, bq, es, bs = labels[name]
            def preparation_progress(row):
                print("SIXTH SOURCE PREPARATION "+json.dumps({"direction": name, **row}), flush=True)
            data = prepare_fields(tensor_from_labels(eq, bq), tensor_from_labels(es, bs), progress=preparation_progress)
            sources = prepare_sixth_sources(data)
            prepared[name] = (data, sources)
            record = {"name": name, "labels_Eq_Bq_Es_Bs":
                      [[[str(x) for x in row] for row in matrix] for matrix in labels[name]],
                      "source_statistics": sources["statistics"], "conventions": sources["conventions"]}
            result["prepared_sources"].append(record)
            print("PREPARED SIXTH SOURCES "+json.dumps(record), flush=True)
            report(step=index)
            checkpoint()

        points = config["source_points"]
        result["source_checks"] = []
        report(phase="source_checks", step=0, of=len(names), unit="directions")
        for index, name in enumerate(names, 1):
            data, sources = prepared[name]
            checks = run_source_checks(data, sources, points, config["source_scaled_tolerance"])
            sector = _sector_check(name, sources, points)
            record = {"name": name, "pointwise": checks, "sector": sector,
                      "passed": checks["passed"] and sector["passed"]}
            result["source_checks"].append(record)
            print("SIXTH SOURCE CHECK "+json.dumps(record), flush=True)
            report(step=index, m={"source_error": checks["maximum"], "source_directions": index})
            gate(record["passed"], "Source "+name)

        report(phase="momentum_parity", step=0, of=1, unit="directions")
        data, sources = prepared["monopole_mixed"]
        result["momentum_parity"] = _momentum_parity(data, sources,
            lambda row: print("PARITY PREPARATION "+json.dumps(row), flush=True))
        print("SOURCE MOMENTUM PARITY "+json.dumps(result["momentum_parity"]), flush=True)
        report(step=1, m={"parity_failures": sum(not x for x in result["momentum_parity"]["checks"].values())})
        gate(result["momentum_parity"]["passed"], "Momentum parity")

        result["gravitational_tail_checks"] = []
        report(phase="gravitational_tail", step=0, of=len(tail_names), unit="directions")
        for index, name in enumerate(tail_names, 1):
            data, _ = prepared[name]
            record = check_gravitational_tail(data, name, config,
                progress=lambda row: print("GRAVITATIONAL TAIL PROGRESS "+json.dumps(row), flush=True))
            result["gravitational_tail_checks"].append(record)
            print("GRAVITATIONAL TAIL CHECK "+json.dumps(record), flush=True)
            report(step=index, m={"gravitational_tail_error": record["maximum_scaled_error"]})
            gate(record["passed"], "Gravitational tail "+name)

        report(phase="export", step=0, of=1, unit="reports")
        result["status"] = "all specified source, Fourier and tail ingredient checks passed"
        result["passed"] = True
        result["summary"] = {
            "fourier_cases": result["fourier_controls"]["cases"],
            "tail_synthetic_checks": result["tail_controls"]["check_count"],
            "source_directions": len(names),
            "source_pointwise_checks": sum(r["pointwise"]["check_count"] for r in result["source_checks"]),
            "maximum_source_error": max(r["pointwise"]["maximum"] for r in result["source_checks"]),
            "maximum_gravitational_tail_error": max(r["maximum_scaled_error"] for r in result["gravitational_tail_checks"]),
            "accuracy_certified": False,
            "full_mixed_H6_evaluated": False,
            "four_point_loop_evaluated": False,
        }
        checkpoint()
        report(step=1, m={"maximum_source_error": result["summary"]["maximum_source_error"],
                         "maximum_gravitational_tail_error": result["summary"]["maximum_gravitational_tail_error"]})
        print("FINAL "+json.dumps({"result": str(output), "source_commit": result["source_commit"],
              **result["summary"], "status": result["status"], "wall_seconds": result["wall_seconds"]}), flush=True)
    except Exception as exc:
        result["status"] = "failed"
        result["passed"] = False
        result["failure"] = {"type": type(exc).__name__, "message": str(exc)[:2000]}
        checkpoint()
        raise
    finally:
        signal.alarm(0)
