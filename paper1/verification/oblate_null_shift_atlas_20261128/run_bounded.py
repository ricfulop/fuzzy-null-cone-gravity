"""Bounded N=3 angular coefficient atlas from authenticated Paper I source."""

from pathlib import Path
from datetime import datetime, timezone, timedelta
from fractions import Fraction
import math
import hashlib
import json
import os
import platform
import re
import signal
import subprocess
import sys
import threading
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_REPOSITORY = Path("/Users/ricfulop/OpenResearch/fuzzy-null-cone-gravity")
WALL_SECONDS = 900
MEMORY_BYTES = 2 * 1024**3
CUTOFF = datetime.fromisoformat("2026-11-28T16:29:00+00:00")
CERTIFICATE = Path("outputs/null-shift-atlas.json")
SCIENCE_NAMES = ("atlas_core.py", "check_atlas.py", "independent_controls.py")
SOURCE_NAMES = SCIENCE_NAMES + ("run_bounded.py", "test_run_bounded.py", "run.sh",
                                "run-config.json", "view-config.json", "view.xml", "README.md")
ARCHIVE_COMMIT = "d9add22b3509700107601b94c2da9f629050dfb3"
ARCHIVE_ROOT = "reports/oblate-null-map-check-2026-11-27/"
MAP_PATH = ARCHIVE_ROOT + "map-checks.json"
REFERENCE_PATH = ARCHIVE_ROOT + "sources/verification/oblate_nonlinear_null_20261127/null_compatibility.py"
ARCHIVE_PINS = {
    MAP_PATH: "75fc5c1036777716ef80bd1911adc9a90d2ef0582110f0ffe3ccc9ed973305cc",
    REFERENCE_PATH: "539448e6b0fa1cee99da51ef2d5519f324e8a868961531e7de201078c8afa6e5",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


_CONSOLE_LOCK = threading.Lock()


def write_console(data):
    """Best-effort, nonblocking forwarding; the archived child log is authoritative."""
    if isinstance(data, str):
        data = data.encode("utf-8", errors="replace")
    if not _CONSOLE_LOCK.acquire(blocking=False):
        return False
    try:
        fd = sys.stdout.fileno()
        blocking = os.get_blocking(fd)
        try:
            os.set_blocking(fd, False)
            pending = memoryview(data)
            for _ in range(64):
                if not pending:
                    return True
                count = os.write(fd, pending[:65536])
                if count <= 0:
                    return False
                pending = pending[count:]
            return not pending
        finally:
            os.set_blocking(fd, blocking)
    except (OSError, ValueError, TypeError, AttributeError):
        return False
    finally:
        _CONSOLE_LOCK.release()


def report(**fields):
    write_console("::orx:: " + json.dumps(fields) + "\n")


def group_members(pgid, timeout=1):
    raw = subprocess.check_output(["ps", "-axo", "pid=,pgid=,rss="], text=True, timeout=timeout)
    if not raw.strip():
        raise ValueError("Empty process table")
    members = []
    for line in raw.splitlines():
        pid, group, rss = map(int, line.split())
        if group == pgid:
            members.append({"pid": pid, "rss_bytes": rss * 1024})
    return members


def group_state(pgid):
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return {"status": "empty", "probe": "ESRCH"}
    except PermissionError as exc:
        try:
            members = group_members(pgid)
        except Exception as error:
            return {"status": "unknown", "probe": "EPERM", "error": repr(exc),
                    "ps_error": repr(error)}
        return {"status": "live" if members else "empty", "probe": "EPERM",
                "members": members, "error": repr(exc)}
    except OSError as error:
        return {"status": "unknown", "probe": "killpg", "error": repr(error)}
    return {"status": "live", "probe": "killpg"}


def stop_group(process, events):
    for sig, grace in [(signal.SIGTERM, 2.0), (signal.SIGKILL, 1.0)]:
        try:
            os.killpg(process.pid, sig)
            events.append({"signal": int(sig), "result": "sent"})
        except OSError as error:
            events.append({"signal": int(sig), "error": repr(error)})
        until = time.monotonic() + grace
        while True:
            process.poll()
            state = group_state(process.pid)
            if state["status"] == "empty" or time.monotonic() >= until:
                events.append({"after_signal": int(sig), **state})
                break
            time.sleep(0.05)
        if state["status"] == "empty":
            break


def save_json(path, data):
    temporary = path.with_suffix(".json.tmp")
    with temporary.open("x", encoding="utf-8") as stream:
        try:
            stream.write(json.dumps(data, indent=2, allow_nan=False) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
            os.link(temporary, path)
        finally:
            temporary.unlink()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def snapshot(data, destination, expected):
    require(hashlib.sha256(data).hexdigest() == expected, f"Input or source pin mismatch: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as stream:
        stream.write(data)
    require(sha(destination) == expected, f"Snapshot hash mismatch: {destination}")


def prepare(out, record):
    commit = record["source_commit"]
    resolved = subprocess.check_output(
        ["git", "-C", str(SOURCE_REPOSITORY), "rev-parse", commit + "^{commit}"], text=True).strip()
    require(resolved == commit, "Source commit did not resolve to its full immutable revision")
    sources = [HERE / name for name in SOURCE_NAMES] + [ROOT / key for key in ARCHIVE_PINS]
    record["source_sha256_before"] = {}
    for path in sources:
        relative = path.relative_to(ROOT).as_posix()
        committed = subprocess.check_output(
            ["git", "-C", str(SOURCE_REPOSITORY), "show", f"{commit}:paper1/{relative}"])
        digest = hashlib.sha256(committed).hexdigest()
        snapshot(path.read_bytes(), out / "sources" / relative, digest)
        record["source_sha256_before"][relative] = digest

    frozen = out / "sources" / HERE.relative_to(ROOT)
    config_path = frozen / "run-config.json"
    config_bytes = config_path.read_bytes()
    config_sha = record["source_sha256_before"][(HERE / "run-config.json").relative_to(ROOT).as_posix()]
    require(hashlib.sha256(config_bytes).hexdigest() == config_sha, "Committed config snapshot changed")
    config = json.loads(config_bytes)
    for key, fixed in (("wall_seconds", WALL_SECONDS), ("memory_bytes", MEMORY_BYTES),
                       ("approval_cutoff_utc", CUTOFF.isoformat()), ("threads", 1),
                       ("python_version", "3.12"), ("sympy_version", "1.14.0")):
        require(type(config.get(key)) is type(fixed) and config[key] == fixed,
                f"Unexpected {key} in committed config")
    require(config["regression_archive"] == {"commit": ARCHIVE_COMMIT, "paper1_paths": ARCHIVE_PINS},
            "Unexpected archived regression binding")
    for path, expected in ARCHIVE_PINS.items():
        prior = subprocess.check_output(
            ["git", "-C", str(SOURCE_REPOSITORY), "show", f"{ARCHIVE_COMMIT}:paper1/{path}"])
        require(hashlib.sha256(prior).hexdigest() == expected and
                record["source_sha256_before"][path] == expected, "Wrong historical regression reference")
    checker = frozen / "check_atlas.py"
    record.update(configuration=config, checker_input_sha256={
        "config_sha256": config_sha,
        "source_sha256": {name: record["source_sha256_before"][
            (HERE/name).relative_to(ROOT).as_posix()] for name in SCIENCE_NAMES},
    })
    record["command"] = ["uv", "run", "--isolated", "--no-project", "--python", "3.12",
                         "--with", "sympy==1.14.0",
                         "python", str(checker), "--config", str(config_path),
                         "--output", str(out / CERTIFICATE)]


def check_hashes_after(out, record, errors):
    for base, field in ((ROOT, "source_sha256_after"),
                        (out / "sources", "source_snapshot_sha256_after")):
        results = {}
        for key, expected in record["source_sha256_before"].items():
            path = base / key
            try:
                results[key] = sha(path)
                require(results[key] == expected, f"Hash changed: {path}")
            except Exception as error:
                errors.append({"action": "verify_hash", "path": str(path), "error": repr(error)})
        record[field] = results


AXES = ("x", "y", "z")
SCALARS = ("a", "b", "lambda", "X")
VECTORS = ("V_polar", "V_axial", "W_polar", "W_axial")
TENSORS = ("C_polar", "C_axial")
CHECKS = ("affine_retained_with_independent_slope", "rotation_discarded_zero",
          "full_coefficient_reconstruction", "metric_chain_rule", "fourier_reality",
          "axisymmetric_eight_rows", "transverse_scalar_phases", "odd_coordinate_leg_normalization")
SCOPE = ("admissible_seed_supplied", "euler_residual_evaluated",
         "nonlinear_compatibility_evaluated", "dirac_classification_established")
TIMERS = {f"{kind}_{axis}" for kind in ("rotation", "null_shift") for axis in AXES}
MUTATIONS = ("round_trace_from_axial_tensor_defect_rejected", "missing_temporal_W_row_rejected",
             "transverse_y_phase_defect_rejected")


def labels(lower, upper):
    return [[kind, ell, m] for ell in range(lower, upper+1) for m in range(-ell, ell+1)
            for kind in SCALARS + (VECTORS if ell >= 1 else ()) + (TENSORS if ell >= 2 else ())]


def same(actual, expected, what):
    # JSON preserves bool/int/float distinctions that Python equality collapses.
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), what)


def unique_keys(pairs):
    obj = {}
    for key, value in pairs:
        require(key not in obj, f"Duplicate JSON field: {key}")
        obj[key] = value
    return obj


def finite_json_number(token):
    value = float(token)
    require(math.isfinite(value), "Nonfinite atlas JSON number")
    return value


def read_payload(path):
    require(path.is_file() and not path.is_symlink(), f"Expected a regular file: {path}")
    data = path.read_bytes()
    obj = json.loads(data, object_pairs_hook=unique_keys,
                     parse_float=finite_json_number, parse_constant=finite_json_number)
    return obj, hashlib.sha256(data).hexdigest()


def checked_member(directory, filename, reference):
    payload, digest = read_payload(directory / filename)
    same(reference, {"file": "atlas-columns/"+filename, "sha256": digest},
         "Invalid atlas member reference: "+filename)
    return payload, digest


def rational(value):
    require(type(value) is str and re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value),
            "Expected rational string")
    parsed = Fraction(value)
    require(str(parsed) == value, "Expected reduced canonical rational")
    return parsed


def validate_column(column, index, label, rows, catalog_sha):
    same(column["label"], label, "Wrong input label")
    same(column["column"], index, "Wrong column index")
    require(column["catalog_sha256"] == catalog_sha, "Column belongs to a different catalog")
    require(column["atlas_complete"] is False and all(column[k] is True for k in
            ("rotation_zero", "round_reconstruction_checked", "metric_chain_rule_checked")),
            "Invalid column verification flags")
    timers = column["map_seconds"]
    require(timers.keys() == TIMERS and all(type(t) in (float, int) and math.isfinite(t) and t >= 0
                                           for t in timers.values()), "Invalid map durations")
    maps = column["maps"]
    require(maps.keys() == set(AXES), "Missing or unexpected generator")
    counts = {}
    for axis, entries in maps.items():
        require(type(entries) is list, "Expected sparse row list")
        previous = -1
        for entry in entries:
            require(entry.keys() == {"row", "real", "imag"}, "Unexpected sparse entry fields")
            row = entry["row"]
            require(type(row) is int and previous < row < len(rows), "Duplicate, unordered or invalid row")
            previous = row
            delta_m = rows[row][2]-label[2]
            require(delta_m == 0 if axis == "z" else abs(delta_m) == 1, "Invalid Fourier selection")
            real, imaginary = rational(entry["real"]), rational(entry["imag"])
            require(bool(real or imaginary), "Explicit zero in sparse matrix")
        counts[axis] = len(entries)
    return counts


def validate_controls(result, totals, axisymmetric_rows):
    total_entries = sum(totals.values())
    require(totals["x"] == totals["y"] and totals["x"] > 0 and totals["z"] >= 9,
            "Vacuous or inconsistent control counts")
    require(len(axisymmetric_rows) == 9 and len(set(axisymmetric_rows)) == 8,
            "Serialized axisymmetric regression has wrong support counts")
    evidence = {
        "explicit_low_harmonic_basis": {"fixture_count": 6, "signed_modes": [-2,-1,0,1,2]},
        "catalog_and_fourier_selection": {"input_columns": 78, "discarded_labels": 160,
                                         "generator_columns": 234, "nonzero_entries": total_entries},
        "archived_axisymmetric_eight_rows": {"input_columns": 22, "nonzero_output_rows": 8,
                                            "nonzero_entries": 9, "all_other_entries_zero": True,
                                            "temporal_W_rows_included": True,
                                            "round_trace_X_from_C_axial_included": True},
        "all_column_complex_conjugacy": {"column_pairs_checked": 234, "coefficient_comparisons": total_entries,
                                        "negative_m_convention": "conjugate positive mode; no factorial reweighting"},
        "all_column_cartesian_transverse_phases": {"input_columns": 78, "coefficient_comparisons": totals["x"],
                                                  "raising_factor": "-I", "lowering_factor": "I"},
        "all_column_rotation_zero": {"rotation_columns_checked": 234, "all_discarded_zero": True},
        "deliberate_coefficient_defects": {"case_count": 3, "cases": [
            {"case": name, "passed": True} for name in MUTATIONS]},
    }
    same(result, {"status": "independent atlas coefficient controls passed", "case_count": 7,
                  "cases": [{"case": n, "passed": True, "evidence": v} for n,v in evidence.items()],
                  "archived_source": {"path": "paper1/"+REFERENCE_PATH, "function": "expected_discarded_generator",
                                      "convention": "unit w; eight output rows, nine nonzero matrix entries"},
                  "scope": "exact coefficient maps and conventions; no Euler pairing or seed"},
         "Independent controls or counts disagree with serialized maps")


def read_certificate(path, expected_hashes):
    require(not path.parent.is_symlink(), "Unexpected linked output directory")
    final, digest = read_payload(path)
    require(final["state"] == "complete N=3 discarded coefficient maps" and final["atlas_complete"] is True,
            "Missing complete atlas")
    same([final["coefficient_column_count"], final["discarded_row_count"], final["generator_count"]],
         [78,160,3], "Wrong final dimensions")
    same(final["checks"], dict.fromkeys(CHECKS, True), "Incomplete final checks")
    same(final["scope"], dict.fromkeys(SCOPE, False), "Unexpected scientific scope")
    members = path.parent/"atlas-columns"
    require(members.is_dir() and not members.is_symlink(), "Expected checkpoint directory")
    names = [f"column-{i:03d}.json" for i in range(78)]
    require({f.name for f in members.iterdir()} == set(names) | {"catalog.json"},
            "Checkpoint inventory differs from the complete atlas")
    catalog, cat_sha = checked_member(members, "catalog.json", final["catalog"])
    for obj in (final, catalog):
        for key in ("config_sha256", "source_sha256"):
            same(obj[key], expected_hashes[key], "Atlas source/config binding changed")
        same(obj["archived_regression_reference"], {"path": "paper1/"+MAP_PATH, "sha256": ARCHIVE_PINS[MAP_PATH]},
             "Wrong regression reference")
    same(final["runtime"], catalog["runtime"], "Catalog and final runtime differ")
    require(re.fullmatch(r"3\.12\.[0-9]+", final["runtime"]["python"]) and
            final["runtime"]["sympy"] == "1.14.0", "Wrong scientific runtime")
    columns, rows = labels(0, 2), labels(3, 4)
    same(catalog["columns"], columns, "Invalid ordered input basis")
    same(catalog["rows"], rows, "Invalid ordered discarded basis")
    same(catalog["labels"], ["field","ell","m"], "Invalid label convention")
    same(catalog["axes"], list(AXES), "Invalid axes")
    require(catalog["atlas_complete"] is False and catalog["weights"] ==
            "Cartesian generator, unit w, exact real/imaginary rational strings", "Invalid catalog semantics")
    refs = final["columns"]
    require(type(refs) is list and len(refs) == len(names), "Incomplete final column references")
    totals = dict.fromkeys(AXES, 0)
    axisymmetric_rows = []
    for i, filename in enumerate(names):
        member, _ = checked_member(members, filename, refs[i])
        counts = validate_column(member, i, columns[i], rows, cat_sha)
        if columns[i][2] == 0:
            axisymmetric_rows.extend(entry["row"] for entry in member["maps"]["z"])
        for axis in AXES:
            totals[axis] += counts[axis]
    validate_controls(final["independent_controls"], totals, axisymmetric_rows)
    metrics = {"atlas_bundle_passed": 1, "checkpoint_column_count": len(names),
               "independent_control_count": 7, "nonzero_entry_count": sum(totals.values())}
    return digest, final["state"], metrics


def fallback_metrics(record):
    return {"atlas_bundle_passed": 0, "scientific_invocations": record["scientific_invocations"]}


def publish_metrics(metrics):
    numeric = {key: value for key, value in metrics.items() if value is not None}
    if metrics["atlas_bundle_passed"]:
        report(phase="complete", step=1, of=1, unit="certificate", m=numeric)
    else:
        report(m=numeric)


def main():
    run_id = os.environ.get("ORX_RUN_ID")
    commit = os.environ.get("ORX_COMMIT")
    require(run_id and Path(run_id).name == run_id and run_id not in (".", ".."), "Missing or invalid run ID")
    require(commit and re.fullmatch("[0-9a-f]{40}", commit), "Missing or invalid full source commit")
    out = ROOT / "run-output" / run_id
    out.mkdir(parents=True, exist_ok=False)
    overrides = {"PYTHONUNBUFFERED": "1", "PYTHONOPTIMIZE": "0", "PYTHONDONTWRITEBYTECODE": "1",
                 "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                 "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1",
                 "VECLIB_MAXIMUM_THREADS": "1"}
    record = {"run_id": run_id, "source_commit": commit, "source_repository": str(SOURCE_REPOSITORY),
              "wall_seconds": WALL_SECONDS, "memory_bytes": MEMORY_BYTES, "approval_cutoff_utc": CUTOFF.isoformat(),
              "environment_overrides": overrides, "supervisor_python": sys.version,
              "scientific_invocations": 0, "remote_jobs": 0, "GPU_jobs": 0,
              "rss_measurement": "sampled child process-group RSS; not a hard OS cap",
              "start_utc": datetime.now(timezone.utc).isoformat()}
    try:
        require(not sys.flags.optimize and os.environ.get("PYTHONOPTIMIZE") in (None, "", "0"),
                "Optimized Python would disable scientific assertions")
        require(datetime.now(timezone.utc) + timedelta(seconds=WALL_SECONDS + 10) < CUTOFF,
                "The bounded invocation would exceed the approval cutoff")
        record["platform"] = platform.platform()
        prepare(out, record)
        save_json(out / "launch.json", record)
    except Exception as error:
        metrics = fallback_metrics(record)
        record.update(status="failed", stop_reason="preflight failed", supervisor_error=repr(error),
                      exit_code=None, group_disposition={"status": "not_spawned"},
                      end_utc=datetime.now(timezone.utc).isoformat(), metrics=metrics)
        save_json(out / "execution.json", record)
        save_json(out / "resource-samples.json", [])
        save_json(out / "metrics.json", metrics)
        publish_metrics(metrics)
        return 1
    report(phase="setup", step=1, of=1, unit="environment")
    return supervise(out, record, {**os.environ, **overrides}, CUTOFF)


def supervise(out, record, env, cutoff):
    start = time.monotonic()
    samples = []
    process = None
    thread = None
    relay_started = False
    relay_active = False
    log = None
    stop_reason = None
    relay_errors = []
    received_signals = []
    cleanup_events = []
    cleanup_errors = []
    supervisor_error = None
    disposition = {"status": "not_spawned"}
    old_handlers = {}
    accepting_cancellation = True

    def interrupted(signum, frame):
        if accepting_cancellation:
            received_signals.append(int(signum))

    def relay():
        try:
            with log, process.stdout:
                for line in iter(process.stdout.readline, b""):
                    log.write(line)
                    log.flush()
                    write_console(line)  # A stalled live console cannot stop archival drain.
        except Exception as error:
            relay_errors.append(repr(error))

    try:
        for sig in (signal.SIGTERM, signal.SIGINT):
            old_handlers[sig] = signal.signal(sig, interrupted)
        log = (out / "verification.log").open("xb")
        if datetime.now(timezone.utc) + timedelta(seconds=record["wall_seconds"] + 10) >= cutoff:
            raise RuntimeError("Setup consumed the approved invocation window")
        if received_signals:
            raise InterruptedError("Supervisor signalled before spawning")
        process = subprocess.Popen(record["command"], cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        record.update(scientific_invocations=1, pid=process.pid)
        if received_signals:
            stop_reason = "supervisor received signal"
        thread = threading.Thread(target=relay, daemon=True)
        thread.start()
        relay_started = True
        while not stop_reason and process.poll() is None:
            elapsed = time.monotonic() - start
            remaining = record["wall_seconds"] - elapsed
            if remaining <= 0:
                stop_reason = "wall budget reached"
                break
            rss = sum(member["rss_bytes"] for member in
                      group_members(process.pid, timeout=min(1, remaining)))
            samples.append({"seconds": round(elapsed, 3), "rss_bytes": rss})
            if time.monotonic() - start >= record["wall_seconds"]:
                stop_reason = "wall budget reached"
            elif rss > record["memory_bytes"]:
                stop_reason = "sampled process-group RSS exceeded threshold"
            elif relay_errors:
                stop_reason = "log relay failed"
            elif received_signals:
                stop_reason = "supervisor received signal"
            else:
                time.sleep(0.25)
        if stop_reason is None and time.monotonic() - start >= record["wall_seconds"]:
            stop_reason = "wall budget reached"
    except Exception as error:
        supervisor_error = repr(error)
        stop_reason = stop_reason or "supervisor exception"
    finally:
        if process is not None:
            disposition = {"status": "unknown"}
            try:
                initial = group_state(process.pid)
            except Exception as error:
                initial = {"status": "unknown"}
                cleanup_errors.append({"action": "initial_probe", "error": repr(error)})
            cleanup_events.append({"action": "initial_probe", **initial})
            if initial["status"] != "empty":
                stop_reason = stop_reason or "child group required cleanup"
                try:
                    stop_group(process, cleanup_events)
                except Exception as error:
                    cleanup_errors.append({"action": "terminate_group", "error": repr(error)})
            try:
                process.wait(timeout=0.5)
            except Exception as error:
                cleanup_errors.append({"action": "wait_child", "error": repr(error)})
            try:
                disposition = group_state(process.pid)
            except Exception as error:
                cleanup_errors.append({"action": "final_probe", "error": repr(error)})
            if disposition["status"] != "empty":
                stop_reason = stop_reason or "child group " + disposition["status"]

        try:
            if relay_started:
                thread.join(timeout=2)
                relay_active = thread.is_alive()
                if relay_active:
                    stop_reason = stop_reason or "log relay remained active after group cleanup"
            else:
                if log is not None:
                    log.close()
                if process is not None:
                    process.stdout.close()
        except Exception as error:
            cleanup_errors.append({"action": "finish_relay", "error": repr(error)})

        try:
            check_hashes_after(out, record, cleanup_errors)
        except Exception as error:
            cleanup_errors.append({"action": "verify_hashes", "error": repr(error)})
        if cleanup_errors:
            stop_reason = stop_reason or "supervisor cleanup failed"
        if relay_errors:
            stop_reason = stop_reason or "log relay failed"
        try:
            save_json(out / "resource-samples.json", samples)
        except Exception as error:
            cleanup_errors.append({"action": "save_samples", "error": repr(error)})
            stop_reason = stop_reason or "resource evidence write failed"
        finally:
            try:
                if received_signals:
                    stop_reason = stop_reason or "supervisor received signal"
                exit_code = None if process is None else process.returncode
                success = (exit_code == 0 and stop_reason is None and disposition["status"] == "empty")
                metrics = fallback_metrics(record)
                if success:
                    try:
                        digest, status, flags = read_certificate(out / CERTIFICATE, record["checker_input_sha256"])
                        record.update(certificate_sha256=digest, atlas_status=status)
                        metrics.update(flags)
                    except Exception as error:
                        success = False
                        stop_reason = "terminal certificate validation failed"
                        record["certificate_error"] = repr(error)
                # Freeze cancellation before publishing any terminal evidence.
                # A later signal is outside this documented terminal decision.
                accepting_cancellation = False
                terminal_signals = list(received_signals)
                if terminal_signals:
                    success = False
                    stop_reason = stop_reason or "supervisor received signal"
                    metrics = fallback_metrics(record)
                try:
                    save_json(out / "metrics.json", metrics)
                except Exception as error:
                    success = False
                    stop_reason = stop_reason or "metrics evidence write failed"
                    cleanup_errors.append({"action": "save_metrics", "error": repr(error)})
                    metrics = fallback_metrics(record)
                record.update({"status": "passed" if success else "failed", "exit_code": exit_code,
                               "elapsed_seconds": time.monotonic() - start,
                               "end_utc": datetime.now(timezone.utc).isoformat(),
                               "stop_reason": stop_reason, "supervisor_error": supervisor_error,
                               "group_disposition": disposition, "cleanup_events": cleanup_events,
                               "cleanup_errors": cleanup_errors, "relay_errors": relay_errors,
                               "relay_active_after_join": relay_active,
                               "received_signals": terminal_signals, "metrics": metrics,
                               "cancellation_boundary": "before terminal evidence publication",
                               "terminal_commit_marker": "terminal-evidence.json",
                               "rss_sample_count": len(samples),
                               "sampled_max_rss_bytes": max((p["rss_bytes"] for p in samples), default=None)})
                save_json(out / "execution.json", record)
                if (out / "metrics.json").exists():
                    # Standalone metrics are provisional until this last marker
                    # authenticates both terminal files. No pass is emitted sooner.
                    save_json(out / "terminal-evidence.json", {
                        "status": record["status"],
                        "execution_sha256": sha(out / "execution.json"),
                        "metrics_sha256": sha(out / "metrics.json"),
                    })
            finally:
                for sig, handler in old_handlers.items():
                    signal.signal(sig, handler)
    if not relay_active:
        publish_metrics(metrics)
        write_console("Saved execution evidence: " + str(out) + "\n")
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
