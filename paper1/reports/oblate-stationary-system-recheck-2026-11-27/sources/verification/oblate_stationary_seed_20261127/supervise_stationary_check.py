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
    process.wait()
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
    with (OUT / "verification.log").open("xb") as log:
        try:
            process = subprocess.Popen(COMMAND, cwd=ROOT, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            record["pid"] = process.pid
            record["scientific_invocations"] = 1
            print(f"Started single checker pid={process.pid}", flush=True)
            last_notice = 0
            while process.poll() is None:
                elapsed = time.monotonic() - start
                raw = subprocess.check_output(["ps", "-axo", "pgid=,rss="], text=True)
                rss = sum(int(parts[1]) * 1024 for line in raw.splitlines()
                          if len(parts := line.split()) == 2 and int(parts[0]) == process.pid)
                peak = max(peak, rss)
                samples.append({"seconds": round(elapsed, 3), "process_group_rss_bytes": rss})
                if elapsed >= LIMIT_SECONDS:
                    stop_reason = "wall budget reached"
                elif rss > LIMIT_RSS_BYTES:
                    stop_reason = "sampled process-group RSS exceeded memory budget"
                if stop_reason:
                    record["group_cleanup"] = stop_group(process)
                    break
                if elapsed - last_notice >= 30:
                    print(f"Running: {elapsed:.1f}s; sampled peak RSS {peak/1024**2:.1f} MiB", flush=True)
                    last_notice = elapsed
                time.sleep(SAMPLE_SECONDS)
            process.wait()
        finally:
            if process is not None and group_exists(process.pid):
                if stop_reason is None:
                    stop_reason = "process-group cleanup required after launcher exit or supervisor exception"
                record["final_group_cleanup"] = stop_group(process)
            if process is not None:
                process.wait()
            record.update({
                "exit_code": None if process is None else process.returncode,
                "elapsed_seconds": round(time.monotonic() - start, 3),
                "end_utc": datetime.now(timezone.utc).isoformat(),
                "stop_reason": stop_reason,
                "sampled_peak_process_group_rss_bytes": peak,
                "source_sha256_after": {name: sha(ROOT / name) for name in expected},
                "supervisor_sha256": sha(Path(__file__)),
            })
            (OUT / "execution.json").write_text(json.dumps(record, indent=2) + "\n")
            (OUT / "resource-samples.json").write_text(json.dumps(samples, indent=2) + "\n")
    print(json.dumps({key: record[key] for key in ["exit_code", "elapsed_seconds",
          "sampled_peak_process_group_rss_bytes", "stop_reason"]}), flush=True)
    sys.exit(0 if record["exit_code"] == 0 and not stop_reason else 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--approval-message", required=True)
    arguments = parser.parse_args()
    main(arguments.source_commit, arguments.approval_message)
