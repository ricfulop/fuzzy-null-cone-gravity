"""One supervised JSON-only reduction from an immutable local run snapshot."""

from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib
import json
import os
import platform
import signal
import subprocess
import sys
import threading
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_REPOSITORY = Path("/Users/ricfulop/OpenResearch/fuzzy-null-cone-gravity")
UPSTREAM_COMPLETION_SHA256 = "640ae2f7c9e1cd3702f8330b074df361ee3fdff701ee2a566ef2b84c4e3a6b4d"


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
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def main():
    config = json.loads((HERE / "run-config.json").read_text())
    run_id = os.environ.get("ORX_RUN_ID")
    commit = os.environ.get("ORX_COMMIT")
    if not run_id or not commit:
        raise RuntimeError("Launch through the recorded experiment; run identity is required")
    if sys.flags.optimize or os.environ.get("PYTHONOPTIMIZE") not in (None, "", "0"):
        raise RuntimeError("Optimized Python would disable scientific assertions")
    cutoff = datetime.fromisoformat(config["approval_cutoff_utc"])
    if datetime.now(timezone.utc) + timedelta(seconds=config["wall_seconds"] + 10) >= cutoff:
        raise RuntimeError("The bounded invocation would exceed the approved launch window")
    out = ROOT / "run-output" / run_id
    out.mkdir(parents=True, exist_ok=False)
    sources = sorted(HERE.glob("*.py")) + [HERE / "run.sh", HERE / "run-config.json"]
    source_hashes = {str(p.relative_to(ROOT)): sha(p) for p in sources}
    resolved = subprocess.check_output(
        ["git", "-C", str(SOURCE_REPOSITORY), "rev-parse", commit + "^{commit}"], text=True).strip()
    if resolved != commit:
        raise RuntimeError("Run commit must identify the full immutable revision")
    for p in sources:
        relative = str(p.relative_to(ROOT))
        committed = subprocess.check_output(
            ["git", "-C", str(SOURCE_REPOSITORY), "show", commit + ":paper1/" + relative])
        if hashlib.sha256(committed).hexdigest() != source_hashes[relative]:
            raise RuntimeError("Materialized source differs from recorded commit: " + relative)
        copied = out / "sources" / p.relative_to(ROOT)
        copied.parent.mkdir(parents=True, exist_ok=True)
        copied.write_bytes(p.read_bytes())
        if sha(copied) != source_hashes[relative]:
            raise RuntimeError("Source snapshot mismatch")
    completion = ROOT / config["input_directory"] / "completion.json"
    if sha(completion) != UPSTREAM_COMPLETION_SHA256:
        raise RuntimeError("Input completion marker differs from the reviewed export")
    overrides = {"PYTHONUNBUFFERED": "1", "PYTHONOPTIMIZE": "0",
                 "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                 "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
    env = os.environ.copy()
    env.update(overrides)
    command = ["uv", "run", "--isolated", "--no-project", "--with", "sympy==1.14.0",
               "python", str(out / "sources" / HERE.relative_to(ROOT) / "reduce_stationary.py"), "--input-dir",
               str(ROOT / config["input_directory"]), "--output-dir", str(out / "outputs"),
               "--stage", config["stage"], "--max-steps", str(config["max_steps"])]
    record = {"run_id": run_id, "source_commit": commit, "command": command,
              "configuration": config, "source_sha256_before": source_hashes,
              "source_snapshot_sha256_before": dict(source_hashes),
              "source_repository": str(SOURCE_REPOSITORY),
              "input_completion_sha256": UPSTREAM_COMPLETION_SHA256,
              "supervisor_python": sys.version, "platform": platform.platform(),
              "environment_overrides": overrides, "scientific_invocations": 0,
              "remote_jobs": 0, "GPU_jobs": 0,
              "rss_measurement": "sampled sum for child process group; not a hard OS cap",
              "start_utc": datetime.now(timezone.utc).isoformat()}
    (out / "launch.json").write_text(json.dumps(record, indent=2) + "\n")
    report(phase="setup", step=0, of=1, unit="environment")
    return supervise(out, record, env, sources, cutoff)


def supervise(out, record, env, sources, cutoff):
    start = time.monotonic()
    config = record["configuration"]
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
        if datetime.now(timezone.utc) + timedelta(seconds=config["wall_seconds"] + 10) >= cutoff:
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
            remaining = config["wall_seconds"] - elapsed
            if remaining <= 0:
                stop_reason = "wall budget reached"
                break
            rss = sum(member["rss_bytes"] for member in
                      group_members(process.pid, timeout=min(1, remaining)))
            samples.append({"seconds": round(elapsed, 3), "rss_bytes": rss})
            if time.monotonic() - start >= config["wall_seconds"]:
                stop_reason = "wall budget reached"
            elif rss > config["memory_bytes"]:
                stop_reason = "sampled process-group RSS exceeded threshold"
            elif relay_errors:
                stop_reason = "log relay failed"
            elif received_signals:
                stop_reason = "supervisor received signal"
            else:
                time.sleep(0.25)
        if stop_reason is None and time.monotonic() - start >= config["wall_seconds"]:
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

        after = {}
        snapshot_after = {}
        try:
            after = {str(p.relative_to(ROOT)): sha(p) for p in sources}
            if after != record["source_sha256_before"]:
                stop_reason = stop_reason or "source files changed during invocation"
            # The shared diagnostic wrapper labels committed sources separately
            # from its original/copied input files in the broader after-check set.
            snapshot_expected = record.get("source_snapshot_sha256_before",
                                           record.get("committed_source_sha256",
                                                      record["source_sha256_before"]))
            snapshot_after = {relative: sha(out / "sources" / relative)
                              for relative in snapshot_expected}
            if snapshot_after != snapshot_expected:
                stop_reason = stop_reason or "executed source snapshot changed during invocation"
        except Exception as error:
            cleanup_errors.append({"action": "source_hashes", "error": repr(error)})
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
                # The terminal decision stops accepting cancellation here. Later
                # signals cannot produce contradictory status/metrics/exit evidence.
                accepting_cancellation = False
                terminal_signals = list(received_signals)
                if terminal_signals:
                    stop_reason = stop_reason or "supervisor received signal"
                exit_code = None if process is None else process.returncode
                success = (exit_code == 0 and stop_reason is None and disposition["status"] == "empty")
                record.update({"status": "passed" if success else "failed", "exit_code": exit_code,
                               "elapsed_seconds": time.monotonic() - start,
                               "end_utc": datetime.now(timezone.utc).isoformat(),
                               "stop_reason": stop_reason, "supervisor_error": supervisor_error,
                               "group_disposition": disposition, "cleanup_events": cleanup_events,
                               "cleanup_errors": cleanup_errors, "relay_errors": relay_errors,
                               "relay_active_after_join": relay_active,
                               "received_signals": terminal_signals, "source_sha256_after": after,
                               "source_snapshot_sha256_after": snapshot_after,
                               "cancellation_boundary": "before terminal evidence publication",
                               "rss_sample_count": len(samples),
                               "sampled_max_rss_bytes": max((p["rss_bytes"] for p in samples), default=None)})
                save_json(out / "execution.json", record)
            finally:
                for sig, handler in old_handlers.items():
                    signal.signal(sig, handler)
    if not relay_active:
        write_console("Saved execution evidence: " + str(out) + "\n")
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
