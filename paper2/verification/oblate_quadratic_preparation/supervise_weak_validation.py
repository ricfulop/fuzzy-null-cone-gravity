"""Bound one local angular weak-slot diagnostic and preserve its execution record.

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


MAX_SECONDS = 300
MAX_RSS_BYTES = 2 * 1024**3
CHECK_MODE = "profile-arithmetic-probe"
CHECKER_NAME = "export_angular_weak_slots.py"
# Renewed user approval; see profile-arithmetic-probe-preparation.md.
LAUNCH_APPROVAL_MESSAGE = "msg_c7ac8e6c"
APPROVAL_END = dt.datetime(2026, 11, 28, 16, 29, tzinfo=dt.timezone.utc)
SOURCE_NAMES = (
    CHECKER_NAME, "weak_forms.py", "potential_hessian.py",
    "kinematic_blocks.py", "supervise_weak_validation.py", "run_weak_validation.sh",
)


def utc_now():
    return dt.datetime.now(dt.timezone.utc)


def save_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def process_group_members(group):
    raw = subprocess.check_output(
        ["ps", "-axo", "pid=,pgid=,rss="], text=True, timeout=5
    )
    selected = []
    for line in raw.splitlines():
        fields = line.split()
        if len(fields) == 3 and int(fields[1]) == group:
            selected.append({"pid": int(fields[0]), "rss_bytes": int(fields[2]) * 1024})
    return selected


def process_group_rss(group):
    selected = process_group_members(group)
    return sum(row["rss_bytes"] for row in selected), selected


def signal_group(group, sig, cleanup_events):
    try:
        os.killpg(group, sig)
        cleanup_events.append({"action": "signal", "signal": int(sig), "status": "sent"})
        return True
    except ProcessLookupError:
        cleanup_events.append({
            "action": "signal", "signal": int(sig), "status": "no_such_group"
        })
        return False
    except PermissionError as exc:
        members = None
        status = "permission_error_unknown"
        try:
            members = process_group_members(group)
            status = "permission_error_live" if members else "permission_error_empty"
        except Exception as member_exc:
            cleanup_events.append({
                "action": "list_members", "status": "failed",
                "type": type(member_exc).__name__, "message": str(member_exc),
            })
        cleanup_events.append({
            "action": "signal", "signal": int(sig), "status": status,
            "errno": exc.errno, "message": str(exc), "members": members,
        })
        return False


def probe_group(group, cleanup_events):
    try:
        members = process_group_members(group)
    except Exception as exc:
        cleanup_events.append({
            "action": "probe", "status": "unknown",
            "type": type(exc).__name__, "message": str(exc), "members": None,
        })
        return True, None
    try:
        os.killpg(group, 0)
        status = "live"
        exists = True
    except ProcessLookupError:
        status = "empty"
        exists = False
    except PermissionError as exc:
        exists = bool(members)
        status = "permission_error_live" if exists else "permission_error_empty"
        cleanup_events.append({
            "action": "probe", "status": status, "errno": exc.errno,
            "message": str(exc), "members": members,
        })
    else:
        cleanup_events.append({"action": "probe", "status": status, "members": members})
    return exists, members


def terminate_group(group, child, cleanup_events, grace_seconds=5):
    signal_group(group, signal.SIGTERM, cleanup_events)
    grace_end = time.monotonic() + grace_seconds
    live = True
    members = []
    while time.monotonic() < grace_end:
        live, members = probe_group(group, cleanup_events)
        if not live:
            return False, members
        child.poll()
        time.sleep(0.05)
    if live:
        signal_group(group, signal.SIGKILL, cleanup_events)
        live, members = probe_group(group, cleanup_events)
    return live, members


def final_disposition(group, child, cleanup_events):
    live, members = probe_group(group, cleanup_events)
    if live:
        signal_group(group, signal.SIGKILL, cleanup_events)
        live, members = probe_group(group, cleanup_events)
        for _ in range(3):
            if not live:
                break
            child.poll()
            time.sleep(0.05)
            live, members = probe_group(group, cleanup_events)

    status = "unknown" if members is None else "live" if live else "empty"
    cleanup_events.append({
        "action": "final_disposition", "status": status, "members": members,
    })
    return live, members, status


def main():
    if not LAUNCH_APPROVAL_MESSAGE:
        raise RuntimeError("This angular diagnostic still requires per-run user approval")
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
    output = paper_dir / "reports" / "oblate-angular-weak-slots" / "profile-arithmetic-probe-runtime"
    output.mkdir(parents=True, exist_ok=False)
    uv = shutil.which("uv")
    if not uv:
        raise RuntimeError("The prepared isolated uv environment requires uv")
    command = [uv, "run", "--isolated", "--no-project", "--python", "3.12",
               "--with", "numpy==2.3.5", "--with", "sympy==1.14.0", "python", "-u",
               str(script_dir / CHECKER_NAME), "--output-dir", str(output),
               "--mode", CHECK_MODE]
    hashes = {name: hashlib.sha256((script_dir / name).read_bytes()).hexdigest()
              for name in SOURCE_NAMES}
    manifest = {
        "started_utc": started.isoformat(), "approval_message": LAUNCH_APPROVAL_MESSAGE,
        "approval_chat_id": "chat_5351fc4a-9bc3-4992-9e39-0a9505c63137",
        "approval_quote": "ok I am going back to sleep, you have my approval for any of the upcoming experiments that come up on the chat or other requests to keep going till I wake up",
        "standing_window_message": "msg_ad636446",
        "approval_end_utc": APPROVAL_END.isoformat(), "wall_budget_seconds": budget,
        "sampled_group_rss_budget_bytes": MAX_RSS_BYTES,
        "memory_budget_is_hard_limit": False, "command": command,
        "cwd": str(paper_dir), "source_sha256": hashes,
        "source_commit": os.environ.get("ORX_COMMIT"),
        "run_id": os.environ.get("ORX_RUN_ID"), "supervisor_python": sys.version,
        "source_equations_changed_by_supervisor": False,
        "status": "prepared_to_launch", "checker_mode": CHECK_MODE,
    }
    save_json(output / "execution.json", manifest)
    print("::orx:: " + json.dumps({
        "phase": "setup", "step": 0, "of": 1, "unit": "environment"
    }), flush=True)
    print("EXECUTION " + json.dumps(manifest), flush=True)
    samples = []
    cleanup_events = []
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
            live, members = terminate_group(group, child, cleanup_events)
            if live:
                cleanup_events.append({
                    "action": "cleanup", "status": "unresolved_live_group",
                    "members": members,
                })
        try:
            return_code = child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            cleanup_events.append({"action": "child_wait", "status": "timeout"})
        # A child that exited must not leave a worker in its process group.
        if group is not None:
            live, members, disposition = final_disposition(group, child, cleanup_events)
        else:
            live, members, disposition = False, [], "empty"
        if live:
            stop_reason = stop_reason or "orphaned_process_group"
            cleanup_events.append({
                "action": "cleanup",
                "status": "post_wait_unknown_group"
                if disposition == "unknown" else "post_wait_live_group",
                "members": members,
            })
    except BaseException as exc:
        cleanup_events.append({
            "action": "supervisor_exception", "status": "raised",
            "type": type(exc).__name__, "message": str(exc),
        })
        stop_reason = stop_reason or "supervisor_exception: " + repr(exc)
        if group is not None:
            signal_group(group, signal.SIGKILL, cleanup_events)
        if child is not None:
            try:
                return_code = child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                cleanup_events.append({"action": "exception_child_wait", "status": "timeout"})
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
        manifest.update({
            "ended_utc": utc_now().isoformat(),
            "elapsed_seconds": time.monotonic() - monotonic_start,
            "child_return_code": return_code, "stop_reason": stop_reason,
            "peak_sampled_group_rss_bytes": max((row["rss_bytes"] for row in samples), default=0),
            "sample_count": len(samples), "cleanup_events": cleanup_events,
            "status": "child_exit_zero" if return_code == 0 and stop_reason is None else "failed",
        })
        save_json(output / "resource-samples.json", samples)
        save_json(output / "execution.json", manifest)
        print("EXECUTION_FINAL " + json.dumps(manifest), flush=True)
    return 0 if return_code == 0 and stop_reason is None else 1


if __name__ == "__main__":
    raise SystemExit(main())
