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


def report(**fields):
    print("::orx:: " + json.dumps(fields), flush=True)


def group_exists(pgid):
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    return True


def stop_group(process):
    for sig, grace in [(signal.SIGTERM, 2.0), (signal.SIGKILL, 1.0)]:
        try:
            os.killpg(process.pid, sig)
        except ProcessLookupError:
            break
        until = time.monotonic() + grace
        while group_exists(process.pid) and time.monotonic() < until:
            process.poll()
            time.sleep(0.05)
    process.wait()


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
               "python", str(HERE / "reduce_stationary.py"), "--input-dir",
               str(ROOT / config["input_directory"]), "--output-dir", str(out / "outputs"),
               "--stage", config["stage"], "--max-steps", str(config["max_steps"])]
    record = {"run_id": run_id, "source_commit": commit, "command": command,
              "configuration": config, "source_sha256_before": source_hashes,
              "source_repository": str(SOURCE_REPOSITORY),
              "input_completion_sha256": UPSTREAM_COMPLETION_SHA256,
              "supervisor_python": sys.version, "platform": platform.platform(),
              "environment_overrides": overrides, "scientific_invocations": 0,
              "remote_jobs": 0, "GPU_jobs": 0,
              "rss_measurement": "sampled sum for child process group; not a hard OS cap",
              "start_utc": datetime.now(timezone.utc).isoformat()}
    (out / "launch.json").write_text(json.dumps(record, indent=2) + "\n")
    report(phase="setup", step=0, of=1, unit="environment")
    start = time.monotonic()
    samples = []
    process = None
    thread = None
    stop_reason = None
    relay_errors = []
    received_signals = []

    def interrupted(signum, frame):
        received_signals.append(signum)

    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    with (out / "verification.log").open("xb") as log:
        try:
            if datetime.now(timezone.utc) + timedelta(seconds=config["wall_seconds"] + 10) >= cutoff:
                raise RuntimeError("Setup consumed the approved invocation window")
            if received_signals:
                raise InterruptedError("Supervisor signalled before spawning")
            process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            record["scientific_invocations"] = 1
            record["pid"] = process.pid
            if received_signals:
                stop_reason = "supervisor received signal"
                stop_group(process)

            def relay():
                try:
                    for line in iter(process.stdout.readline, b""):
                        log.write(line)
                        log.flush()
                        sys.stdout.buffer.write(line)
                        sys.stdout.buffer.flush()
                except BaseException as error:
                    relay_errors.append(repr(error))

            thread = threading.Thread(target=relay, daemon=True)
            thread.start()
            while process.poll() is None:
                elapsed = time.monotonic() - start
                raw = subprocess.check_output(["ps", "-axo", "pgid=,rss="], text=True)
                rss = sum(int(parts[1]) * 1024 for line in raw.splitlines()
                          if len(parts := line.split()) == 2 and int(parts[0]) == process.pid)
                samples.append({"seconds": round(elapsed, 3), "rss_bytes": rss})
                if elapsed >= config["wall_seconds"]:
                    stop_reason = "wall budget reached"
                elif rss > config["memory_bytes"]:
                    stop_reason = "sampled process-group RSS exceeded threshold"
                elif relay_errors:
                    stop_reason = "log relay failed"
                if received_signals:
                    stop_reason = "supervisor received signal"
                if stop_reason:
                    stop_group(process)
                    break
                time.sleep(0.25)
            process.wait()
            thread.join(timeout=3)
            if thread.is_alive():
                stop_reason = stop_reason or "log relay did not finish"
        finally:
            if received_signals:
                stop_reason = stop_reason or "supervisor received signal"
            if process is not None and group_exists(process.pid):
                stop_reason = stop_reason or "surviving child group required cleanup"
                stop_group(process)
            if thread is not None:
                thread.join(timeout=3)
                if thread.is_alive():
                    stop_reason = stop_reason or "log relay remained active after group cleanup"
            after = {str(p.relative_to(ROOT)): sha(p) for p in sources}
            if after != source_hashes:
                stop_reason = stop_reason or "source files changed during invocation"
            record.update({"exit_code": None if process is None else process.returncode,
                           "elapsed_seconds": time.monotonic() - start,
                           "end_utc": datetime.now(timezone.utc).isoformat(),
                           "stop_reason": stop_reason, "relay_errors": relay_errors,
                           "received_signals": received_signals,
                           "source_sha256_after": after,
                           "sampled_max_rss_bytes": max((p["rss_bytes"] for p in samples), default=0)})
            (out / "execution.json").write_text(json.dumps(record, indent=2) + "\n")
            (out / "resource-samples.json").write_text(json.dumps(samples, indent=2) + "\n")
    print("Saved execution evidence: " + str(out), flush=True)
    return 0 if record["exit_code"] == 0 and not stop_reason and not relay_errors else 1


if __name__ == "__main__":
    sys.exit(main())
