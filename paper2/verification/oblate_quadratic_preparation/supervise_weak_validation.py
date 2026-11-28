"""Bound one local weak-curvature check and preserve its execution record.

Resource accounting is sampled process-group RSS, not a hard memory limit.
The scientific process runs from the framework's immutable source snapshot.
"""

import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time


MAX_SECONDS = 900
MAX_RSS_BYTES = 2 * 1024**3
APPROVAL_END = dt.datetime(2026, 11, 28, 16, 29, tzinfo=dt.timezone.utc)
SOURCE_NAMES = (
    "verify_weak_forms.py", "weak_forms.py", "potential_hessian.py",
    "kinematic_blocks.py", "supervise_weak_validation.py", "run_weak_validation.sh",
)


def utc_now():
    return dt.datetime.now(dt.timezone.utc)


def save_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def process_group_rss(group):
    raw = subprocess.check_output(
        ["ps", "-axo", "pid=,pgid=,rss="], text=True, timeout=5
    )
    selected = []
    for line in raw.splitlines():
        fields = line.split()
        if len(fields) == 3 and int(fields[1]) == group:
            selected.append({"pid": int(fields[0]), "rss_bytes": int(fields[2]) * 1024})
    return sum(row["rss_bytes"] for row in selected), selected


def stop_group(group, sig):
    try:
        os.killpg(group, sig)
    except ProcessLookupError:
        pass


def group_exists(group):
    try:
        os.killpg(group, 0)
    except ProcessLookupError:
        return False
    return True


def main():
    started = utc_now()
    # Reserve cleanup time; never start a run beyond the approval window.
    remaining = (APPROVAL_END - started).total_seconds() - 10
    if remaining <= 0:
        raise RuntimeError("The overnight authorization window has expired")
    budget = min(MAX_SECONDS, remaining)
    if not os.environ.get("ORX_RUN_ID") or not os.environ.get("ORX_COMMIT"):
        raise RuntimeError("Launch through the recorded local experiment contract")
    script_dir = Path(__file__).resolve().parent
    paper_dir = script_dir.parent.parent
    output = paper_dir / "reports" / "oblate-weak-curvature-validation" / "runtime"
    output.mkdir(parents=True, exist_ok=False)
    uv = shutil.which("uv")
    if not uv:
        raise RuntimeError("The prepared isolated uv environment requires uv")
    command = [uv, "run", "--isolated", "--no-project", "--python", "3.12",
               "--with", "numpy==2.3.5", "--with", "sympy==1.14.0", "python", "-u",
               str(script_dir / "verify_weak_forms.py"), "--output-dir", str(output)]
    hashes = {name: hashlib.sha256((script_dir / name).read_bytes()).hexdigest()
              for name in SOURCE_NAMES}
    manifest = {
        "started_utc": started.isoformat(), "approval_message": "msg_ad636446",
        "approval_end_utc": APPROVAL_END.isoformat(), "wall_budget_seconds": budget,
        "sampled_group_rss_budget_bytes": MAX_RSS_BYTES,
        "memory_budget_is_hard_limit": False, "command": command,
        "cwd": str(paper_dir), "source_sha256": hashes,
        "source_commit": os.environ.get("ORX_COMMIT"),
        "run_id": os.environ.get("ORX_RUN_ID"), "supervisor_python": sys.version,
        "source_equations_changed_by_supervisor": False,
        "status": "prepared_to_launch",
    }
    save_json(output / "execution.json", manifest)
    print("::orx:: " + json.dumps({"phase": "setup"}), flush=True)
    print("EXECUTION " + json.dumps(manifest), flush=True)
    samples = []
    child = None
    group = None
    termination_requested = [None]

    def on_signal(signum, frame):
        termination_requested[0] = "supervisor_signal_%d" % signum

    old_handlers = {sig: signal.signal(sig, on_signal) for sig in (signal.SIGINT, signal.SIGTERM)}
    monotonic_start = time.monotonic()
    stop_reason = None
    return_code = None
    try:
        child = subprocess.Popen(command, cwd=paper_dir, start_new_session=True)
        group = child.pid
        manifest.update({"status": "running", "process_group": group})
        save_json(output / "execution.json", manifest)
        while child.poll() is None:
            elapsed = time.monotonic() - monotonic_start
            rss, members = process_group_rss(group)
            samples.append({"elapsed_seconds": elapsed, "rss_bytes": rss, "members": members})
            if termination_requested[0]:
                stop_reason = termination_requested[0]
            elif elapsed >= budget or utc_now() >= APPROVAL_END:
                stop_reason = "wall_budget"
            elif rss > MAX_RSS_BYTES:
                stop_reason = "sampled_group_rss_budget"
            if stop_reason:
                break
            time.sleep(0.25)
        if stop_reason:
            stop_group(group, signal.SIGTERM)
            grace_end = time.monotonic() + 5
            while group_exists(group) and time.monotonic() < grace_end:
                child.poll()
                time.sleep(0.05)
            if group_exists(group):
                stop_group(group, signal.SIGKILL)
        return_code = child.wait(timeout=10)
        # A child that exited must not leave a worker in its process group.
        if group_exists(group):
            stop_group(group, signal.SIGKILL)
            stop_reason = stop_reason or "orphaned_process_group"
    except BaseException as exc:
        stop_reason = stop_reason or "supervisor_exception: " + repr(exc)
        if group is not None:
            stop_group(group, signal.SIGKILL)
        if child is not None:
            try:
                return_code = child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
        manifest.update({
            "ended_utc": utc_now().isoformat(),
            "elapsed_seconds": time.monotonic() - monotonic_start,
            "child_return_code": return_code, "stop_reason": stop_reason,
            "peak_sampled_group_rss_bytes": max((row["rss_bytes"] for row in samples), default=0),
            "sample_count": len(samples),
            "status": "child_exit_zero" if return_code == 0 and stop_reason is None else "failed",
        })
        save_json(output / "resource-samples.json", samples)
        save_json(output / "execution.json", manifest)
        print("EXECUTION_FINAL " + json.dumps(manifest), flush=True)
    return 0 if return_code == 0 and stop_reason is None else 1


if __name__ == "__main__":
    raise SystemExit(main())
