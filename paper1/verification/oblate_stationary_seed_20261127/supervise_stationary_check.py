"""Proposed single local invocation; requires a separately approved source/budget."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/oblate-stationary-system-recheck-2026-11-27"
PREP = ROOT / "reports/stationary-extraction-diagnostics-2026-11-27"
COMMAND = ["uv", "run", "--isolated", "--no-project", "--with",
           "numpy==2.3.5", "--with", "sympy==1.14.0", "python",
           "verification/oblate_stationary_seed_20261127/check_stationary_system.py",
           "--output-dir", str(OUT / "outputs")]
LIMIT_SECONDS = 900
LIMIT_RSS_BYTES = 2 * 1024**3
SAMPLE_SECONDS = 0.25


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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def group_exists(pgid):
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    return True


def signal_group(pgid, sig):
    try:
        os.killpg(pgid, sig)
    except ProcessLookupError:
        pass


def stop_group(process, grace_seconds=2):
    """Stop descendants even if the launcher exits before the grace period."""
    pgid = process.pid
    signal_group(pgid, signal.SIGTERM)
    deadline = time.monotonic()+grace_seconds
    while group_exists(pgid) and time.monotonic() < deadline:
        process.poll()  # Reap the launcher if possible; do not equate it with the group.
        time.sleep(0.05)
    if group_exists(pgid):
        signal_group(pgid, signal.SIGKILL)
    process.wait(timeout=0.5)
    return {"group_visible_after_signals": group_exists(pgid)}


def main(pin, approval_message):
    # An existing log is an invocation sentinel, even if a previous run failed.
    if (OUT / "verification.log").exists():
        raise RuntimeError("Invocation already recorded; no automatic rerun permitted")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    if head != pin:
        raise RuntimeError(f"Pinned commit mismatch: {head}")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True):
        raise RuntimeError("Worktree must be clean")
    if sys.flags.optimize or os.environ.get("PYTHONOPTIMIZE") not in (None, "", "0"):
        raise RuntimeError("Python optimization must not disable scientific assertions")
    expected = json.loads((PREP / "preparation-provenance.json").read_text())["sha256"]
    for name, digest in expected.items():
        if sha(ROOT / name) != digest:
            raise RuntimeError(f"Source hash mismatch: {name}")
    OUT.mkdir(parents=True, exist_ok=False)
    for name, digest in expected.items():
        source = ROOT / name
        frozen = OUT / "sources" / name
        frozen.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, frozen)
    for name in ["outputs/stationary-operators.json", "outputs/stationary-reproducibility.json"]:
        if (OUT / name).exists():
            raise RuntimeError(f"Existing result must not be overwritten: {name}")
    env = os.environ.copy()
    overrides = {"PYTHONUNBUFFERED": "1", "PYTHONOPTIMIZE": "0", "OPENBLAS_NUM_THREADS": "1",
                 "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                 "VECLIB_MAXIMUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}
    env.update(overrides)
    record = {
        "source_commit": pin,
        "approval_messages": [approval_message],
        "command": COMMAND, "cwd": str(ROOT),
        "uv_path": shutil.which("uv"),
        "uv_version": subprocess.check_output(["uv", "--version"], text=True).strip(),
        "platform": platform.platform(), "supervisor_python": sys.version,
        "environment_overrides": overrides,
        "source_sha256_before": expected,
        "stdout_stderr": "verification.log (merged, unmodified bytes)",
        "wall_budget_seconds": LIMIT_SECONDS,
        "memory_budget_bytes": LIMIT_RSS_BYTES,
        "rss_measurement": "sum of ps RSS for the launched process group; sampled, not a hard OS memory cap",
        "sampling_period_seconds": SAMPLE_SECONDS,
        "scientific_invocations": 0,
        "tracked_experiment_launches": 0, "GPU_jobs": 0, "remote_jobs": 0,
        "start_utc": datetime.now(timezone.utc).isoformat(),
    }
    start = time.monotonic()
    peak = 0
    samples = []
    stop_reason = None
    process = None
    log = None
    cleanup_errors = []
    supervisor_error = None
    received_signals = []
    old_handlers = {}
    accepting_cancellation = True

    def interrupted(signum, frame):
        if accepting_cancellation:
            received_signals.append(int(signum))

    try:
        for sig in (signal.SIGTERM, signal.SIGINT):
            old_handlers[sig] = signal.signal(sig, interrupted)
        log = (OUT / "verification.log").open("xb")
        if received_signals:
            raise InterruptedError("Supervisor signalled before spawning")
        process = subprocess.Popen(COMMAND, cwd=ROOT, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        record["pid"] = process.pid
        record["scientific_invocations"] = 1
        write_console(f"Started single checker pid={process.pid}\n")
        last_notice = 0
        while process.poll() is None:
            elapsed = time.monotonic() - start
            if received_signals:
                stop_reason = "supervisor received signal"
                break
            if elapsed >= LIMIT_SECONDS:
                stop_reason = "wall budget reached"
                break
            raw = subprocess.check_output(["ps", "-axo", "pgid=,rss="], text=True,
                                          timeout=min(1, LIMIT_SECONDS - elapsed))
            rss = sum(int(parts[1]) * 1024 for line in raw.splitlines()
                      if len(parts := line.split()) == 2 and int(parts[0]) == process.pid)
            peak = max(peak, rss)
            samples.append({"seconds": round(elapsed, 3), "process_group_rss_bytes": rss})
            if time.monotonic() - start >= LIMIT_SECONDS:
                stop_reason = "wall budget reached"
            elif rss > LIMIT_RSS_BYTES:
                stop_reason = "sampled process-group RSS exceeded memory budget"
            if stop_reason:
                break
            if elapsed - last_notice >= 30:
                write_console(f"Running: {elapsed:.1f}s; sampled peak RSS {peak/1024**2:.1f} MiB\n")
                last_notice = elapsed
            time.sleep(SAMPLE_SECONDS)
        if stop_reason is None and time.monotonic() - start >= LIMIT_SECONDS:
            stop_reason = "wall budget reached"
    except Exception as error:
        supervisor_error = repr(error)
        stop_reason = stop_reason or "supervisor exception"
    finally:
        try:
            if received_signals:
                stop_reason = stop_reason or "supervisor received signal"
            if process is not None:
                try:
                    live_group = group_exists(process.pid)
                except Exception as error:
                    live_group = True  # An unknown group still requires cleanup.
                    cleanup_errors.append({"action": "probe_group", "error": repr(error)})
                if live_group:
                    stop_reason = stop_reason or "process-group cleanup required after launcher exit or supervisor exception"
                    try:
                        record["final_group_cleanup"] = stop_group(process)
                    except Exception as error:
                        cleanup_errors.append({"action": "stop_group", "error": repr(error)})
                try:
                    process.wait(timeout=0.5)
                except Exception as error:
                    cleanup_errors.append({"action": "wait_child", "error": repr(error)})
            if log is not None:
                try:
                    log.close()
                except Exception as error:
                    cleanup_errors.append({"action": "close_log", "error": repr(error)})
            after = {}
            for name in expected:
                try:
                    after[name] = sha(ROOT / name)
                except Exception as error:
                    cleanup_errors.append({"action": "source_hash", "path": name, "error": repr(error)})
            if after != expected:
                stop_reason = stop_reason or "source files changed during invocation"
            supervisor_hash = None
            try:
                supervisor_hash = sha(Path(__file__))
            except Exception as error:
                cleanup_errors.append({"action": "supervisor_hash", "error": repr(error)})
            if cleanup_errors:
                stop_reason = stop_reason or "supervisor cleanup failed"
            accepting_cancellation = False
            terminal_signals = list(received_signals)
            if terminal_signals:
                stop_reason = stop_reason or "supervisor received signal"
            record.update({
                "status": "passed" if process is not None and process.returncode == 0 and stop_reason is None else "failed",
                "exit_code": None if process is None else process.returncode,
                "elapsed_seconds": round(time.monotonic() - start, 3),
                "end_utc": datetime.now(timezone.utc).isoformat(),
                "stop_reason": stop_reason,
                "sampled_peak_process_group_rss_bytes": peak if samples else None,
                "received_signals": terminal_signals,
                "cancellation_boundary": "before terminal evidence publication",
                "source_sha256_after": after,
                "supervisor_sha256": supervisor_hash,
                "supervisor_error": supervisor_error,
                "cleanup_errors": cleanup_errors,
            })
            (OUT / "execution.json").write_text(json.dumps(record, indent=2) + "\n")
            (OUT / "resource-samples.json").write_text(json.dumps(samples, indent=2) + "\n")
        finally:
            for sig, handler in old_handlers.items():
                signal.signal(sig, handler)
    write_console(json.dumps({key: record[key] for key in ["exit_code", "elapsed_seconds",
          "sampled_peak_process_group_rss_bytes", "stop_reason"]}) + "\n")
    sys.exit(0 if record["exit_code"] == 0 and not stop_reason else 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--approval-message", required=True)
    arguments = parser.parse_args()
    main(arguments.source_commit, arguments.approval_message)
