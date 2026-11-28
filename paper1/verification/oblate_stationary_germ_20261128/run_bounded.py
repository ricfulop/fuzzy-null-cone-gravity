"""One supervised terminal-germ check using committed source and two pinned JSON inputs."""

from pathlib import Path
from datetime import datetime, timezone, timedelta
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
CERTIFICATE = Path("outputs/germ-certificate.json")
INPUT_NAMES = {"stationary_operators": "stationary-operators.json", "dae_closure": "dae-closure.json"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def report(**fields):
    print("::orx:: " + json.dumps(fields), flush=True)


def group_members(pgid):
    raw = subprocess.check_output(["ps", "-axo", "pid=,pgid=,rss="], text=True, timeout=1)
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
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


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
    sources = [HERE / name for name in ("run_bounded.py", "test_run_bounded.py", "run.sh",
                                         "certify_germ.py", "run-config.json", "view-config.json", "view.xml")]
    sources.append(ROOT / "verification/oblate_stationary_reduction_20261128/run.sh")
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
                       ("approval_cutoff_utc", CUTOFF.isoformat())):
        require(config.get(key, fixed) == fixed, f"Unexpected {key} in committed config")
    terms = config.get("max_atan_terms", 256)
    require(type(terms) is int and terms > 0, "max_atan_terms must be a positive integer")
    record.update(configuration=config, input_paths={}, input_sha256_before={})
    for key, filename in INPUT_NAMES.items():
        entry = config["inputs"][key]
        require(isinstance(entry["path"], str) and entry["path"], f"Missing {key} input path")
        source = (ROOT / entry["path"]).resolve(strict=True)
        digest = entry["sha256"]
        require(isinstance(digest, str) and re.fullmatch("[0-9a-f]{64}", digest),
                f"Missing or invalid {key} SHA-256")
        snapshot(source.read_bytes(), out / "inputs" / filename, digest)
        record["input_paths"][key] = str(source)
        record["input_sha256_before"][key] = digest

    record["checker_input_sha256"] = {"config": config_sha, **record["input_sha256_before"]}
    record["command"] = [
        "uv", "run", "--isolated", "--no-project", "--with", "sympy==1.14.0",
        "python", str(frozen / "certify_germ.py"), "--config", str(config_path),
        "--operators", str(out / "inputs" / INPUT_NAMES["stationary_operators"]),
        "--closure", str(out / "inputs" / INPUT_NAMES["dae_closure"]),
        "--output", str(out / CERTIFICATE), "--max-atan-terms", str(terms),
    ]


def check_hashes_after(out, record, errors):
    # Check both the original bytes and the copies actually passed to the child.
    for group in ("source", "input"):
        originals, copies = {}, {}
        for key, expected in record[f"{group}_sha256_before"].items():
            if group == "source":
                original, saved = ROOT / key, out / "sources" / key
            else:
                original = Path(record["input_paths"][key])
                saved = out / "inputs" / INPUT_NAMES[key]
            for path, results in ((original, originals), (saved, copies)):
                try:
                    results[key] = sha(path)
                    require(results[key] == expected, f"Hash changed: {path}")
                except Exception as error:
                    errors.append({"action": "verify_hash", "path": str(path), "error": repr(error)})
        record[f"{group}_sha256_after"] = originals
        record[f"{group}_snapshot_sha256_after"] = copies


def read_certificate(path, expected_hashes):
    data = path.read_bytes()
    cert = json.loads(data)
    require(cert["input_sha256"] == expected_hashes, "Certificate input identity mismatch")
    require(cert["state"] == "terminal germ certificate passed", "Missing terminal pass")
    require(cert["manufactured_controls"]["state"] == "terminal germ manufactured controls passed",
            "Missing manufactured-control pass")
    cases = cert["manufactured_controls"]["cases"]
    require(isinstance(cases, list) and cases, "Missing manufactured-control cases")
    for key in ("smooth_solution_germ_certified", "all_22_retained_rows_zero",
                "x0_equation_included_in_22_retained_rows", "discarded_response_matches_reducer_serialization",
                "constant_control_augmented_state_supported"):
        require(cert[key] is True, f"Missing terminal flag: {key}")
    require(cert["epsilon0"] == "1/100" and cert["nonlinear_consistency_established"] is False
            and cert["instantaneous_configuration_only_claim"] is False, "Unexpected certificate scope")
    response = cert["response_certificate"]
    status = response["status"]
    states = {
        "silent_on_terminal_subfamily": (False, False, False),
        "response_nonzero_at_certified_r0": (True, True, True),
        "formal_response_nonzero_but_unresolved_at_sampled_r0": (True, None, False),
    }
    require(status in states, "Unknown response certificate status")
    formal, measured, certified = states[status]
    for subkey, key, expected in (
        ("formal_response_nonzero", "formal_discarded_response_nonzero", formal),
        ("response_nonzero", "discarded_response_nonzero", measured),
        ("certified_response_nonzero", "certified_discarded_response_nonzero", certified),
    ):
        require(response[subkey] is expected and cert[key] is expected, f"Inconsistent {key}")
    require(cert["nonintegrability_conclusion_supported"] is certified
            and response["admissible_non_gauge_seed_supplied"] is certified
            and response["universal_zero_claim"] is False, "Inconsistent response scope")
    metrics = {"terminal_germ_certificate_passed": 1, "manufactured_case_count": len(cases),
               "formal_response_nonzero": int(formal), "certified_response_nonzero": int(certified),
               "response_determined": int(measured is not None),
               "discarded_response_nonzero": None if measured is None else int(measured),
               "all_22_retained_rows_zero": 1, "nonintegrability_conclusion_supported": int(certified)}
    for key in ("free_state_count", "constant_control_count"):
        require(type(cert[key]) is int and cert[key] >= 0, f"Invalid {key}")
        metrics[key] = cert[key]
    return hashlib.sha256(data).hexdigest(), status, metrics


def fallback_metrics(record):
    return {"terminal_germ_certificate_passed": 0, "formal_response_nonzero": None,
            "certified_response_nonzero": None, "discarded_response_nonzero": None,
            "response_determined": 0, "scientific_invocations": record["scientific_invocations"]}


def publish_metrics(metrics):
    numeric = {key: value for key, value in metrics.items() if value is not None}
    if metrics["terminal_germ_certificate_passed"]:
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
                 "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
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

    def interrupted(signum, frame):
        received_signals.append(int(signum))

    def relay():
        try:
            with log, process.stdout:
                for line in iter(process.stdout.readline, b""):
                    log.write(line)
                    log.flush()
                    # A blocked console must not hold Python's shared stdout lock
                    # while the supervising thread saves evidence or exits.
                    remaining = memoryview(line)
                    while remaining:
                        count = os.write(sys.stdout.fileno(), remaining)
                        if count <= 0:
                            raise OSError("Console relay made no progress")
                        remaining = remaining[count:]
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
            rss = sum(member["rss_bytes"] for member in group_members(process.pid))
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
                        record.update(certificate_sha256=digest, response_status=status)
                        metrics.update(flags)
                    except Exception as error:
                        success = False
                        stop_reason = "terminal certificate validation failed"
                        record["certificate_error"] = repr(error)
                if received_signals:
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
                               "received_signals": received_signals, "metrics": metrics,
                               "rss_sample_count": len(samples),
                               "sampled_max_rss_bytes": max((p["rss_bytes"] for p in samples), default=0)})
                save_json(out / "execution.json", record)
            finally:
                for sig, handler in old_handlers.items():
                    signal.signal(sig, handler)
    if not relay_active:
        publish_metrics(metrics)
        print("Saved execution evidence: " + str(out), flush=True)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
