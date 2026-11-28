"""Single authorized local invocation; stdlib supervision, no scientific imports."""
from pathlib import Path
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

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
PIN = "de24fa2bb5c21c2b027498581872e64cb65da445"
PREP = ROOT / "reports/oblate-stationary-seed-preparation-2026-11-27"
COMMAND = ["uv", "run", "--isolated", "--no-project", "--with",
           "numpy==2.3.5", "--with", "sympy==1.14.0", "python",
           "verification/oblate_stationary_seed_20261127/check_stationary_system.py"]
LIMIT_SECONDS = 300
LIMIT_RSS_BYTES = 2 * 1024**3
SAMPLE_SECONDS = 0.25


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    # An existing log is an invocation sentinel, even if a previous run failed.
    if (OUT / "verification.log").exists():
        raise RuntimeError("Invocation already recorded; no automatic rerun permitted")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    assert head == PIN, head
    expected = json.loads((PREP / "preparation-provenance.json").read_text())["sha256"]
    for name, digest in expected.items():
        source = ROOT / name
        assert sha(source) == digest, name
        frozen = OUT / "sources" / name
        frozen.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, frozen)
    for name in ["stationary-operators.json", "stationary-reproducibility.json"]:
        assert not (OUT / name).exists(), name
    env = os.environ.copy()
    overrides = {"PYTHONUNBUFFERED": "1", "OPENBLAS_NUM_THREADS": "1",
                 "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                 "VECLIB_MAXIMUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}
    env.update(overrides)
    record = {
        "source_commit": PIN,
        "approval_messages": ["msg_fa2e97d4"],
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
        "scientific_invocations": 1,
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
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                    break
                if elapsed - last_notice >= 30:
                    print(f"Running: {elapsed:.1f}s; sampled peak RSS {peak/1024**2:.1f} MiB", flush=True)
                    last_notice = elapsed
                time.sleep(SAMPLE_SECONDS)
            process.wait()
        finally:
            if process is not None and process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)
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
    main()
