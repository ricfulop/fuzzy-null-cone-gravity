"""Bounded selected endpoint checks from the authenticated Paper I source snapshot."""

from pathlib import Path
from datetime import datetime, timezone, timedelta
from fractions import Fraction
import ast
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
WALL_SECONDS = 300
MEMORY_BYTES = 2 * 1024**3
CUTOFF = datetime.fromisoformat("2026-11-28T16:29:00+00:00")
CERTIFICATE = Path("outputs/lapse-endpoint-checks.json")
SOURCE_NAMES = ("run_bounded.py", "test_run_bounded.py", "run.sh", "check_lapse_endpoints.py",
                "run-config.json", "view-config.json", "view.xml", "README.md")
REPORT = Path("reports/oblate-common-domain-preparation-2026-11-28/README.md")
REPORT_COMMIT = "ae32314d27e20ac682c5ed1b71b4d7f11bee7133"
REPORT_SHA256 = "fa659f17f3ba9a8dce1ff0d2201cb5418fa2a10faeb64251ebcc20f690d3bc7d"
MATH_REVISION = "773ccc4e3e8a1a89e06086163042fbf45e171c76"


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
    sources = [HERE / name for name in SOURCE_NAMES] + [ROOT / REPORT]
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
    require(config["report"] == {"path": "paper1/" + REPORT.as_posix(), "commit": REPORT_COMMIT,
                                 "sha256": REPORT_SHA256}, "Unexpected preparation report binding")
    report_bytes = subprocess.check_output(
        ["git", "-C", str(SOURCE_REPOSITORY), "show", f"{REPORT_COMMIT}:paper1/{REPORT}"])
    require(hashlib.sha256(report_bytes).hexdigest() == REPORT_SHA256
            and record["source_sha256_before"][REPORT.as_posix()] == REPORT_SHA256,
            "Preparation report differs from its published revision")
    require(config["derivation_commit"] == MATH_REVISION, "Unexpected derivation revision")
    checker = frozen / "check_lapse_endpoints.py"
    record.update(configuration=config, checker_input_sha256={"config": config_sha, "checker": sha(checker)})
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


ZERO_FIELDS = {
    "metric_determinant_sigma_squared": ("residual",),
    "background_ricci_zero": ("residual",),
    "scalar_curvature_green_density": ("residual",),
    "two_hessian_density": ("residual",),
    "scalar_dirichlet_radial_flux": ("residual",),
    "hessian_dirichlet_radial_flux": ("residual",),
    "scalar_angular_pole_fluxes": ("south_residual", "north_residual"),
    "hessian_angular_pole_fluxes": ("south_residual", "north_residual"),
    "clamped_fluxes": ("scalar_residual", "hessian_residual"),
}
CONTROL_NAMES = {
    "linearized Ricci scalar constructed from varied Christoffels",
    "strong-minus-adjoint density equals radial plus angular flux divergence",
    "two-Hessian density equals radial plus angular flux divergence",
    "Dirichlet radial fluxes match reviewed endpoint predictions",
    "angular pole fluxes vanish",
    "clamped endpoint slopes kill both fluxes",
    "sign and factor defects stay nonzero",
    "outer-only and inner-only endpoint orientations have expected signs",
}
IDENTITIES = ("metric_determinant_matches_sigma_squared", "linearized_ricci_from_varied_christoffels",
              "scalar_curvature_density_divergence", "two_hessian_density_divergence",
              "angular_pole_fluxes_vanish", "dirichlet_radial_fluxes_match")
SCOPE = ("full_operator_domain_established", "evolution_invariance_established",
         "physical_null_boundary_selected", "final_Dirac_classification_established",
         "physical_DOF_count_established")


def passed_cases(rows, names):
    require(isinstance(rows, list) and len(rows) == len(names), "Unexpected case count")
    by_name = {row["name"]: row for row in rows}
    require(by_name.keys() == set(names) and all(row["passed"] is True for row in rows),
            "Missing, duplicate or failed check")
    return by_name


def exact_formula(value):
    """Check serialized field syntax only. Exact identities are tested by the pinned child."""
    require(isinstance(value, str) and bool(value.strip()), "Empty expression")
    nodes = list(ast.walk(ast.parse(value, mode="eval")))
    allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Load,
               ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd)
    for node in nodes:
        require(isinstance(node, allowed), "Unexpected expression syntax")
        if isinstance(node, ast.Constant):
            require(type(node.value) is int, "Expected integer literal in exact expression")
        if isinstance(node, ast.Name):
            require(node.id in {"r", "x", "F1", "G1"}, "Unexpected or nonfinite expression symbol")
    require(any(isinstance(node, ast.Name) for node in nodes), "Expected variable flux expression")


def exact_sample(value):
    require(type(value) is str and re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value),
            "Expected an exact rational sample")
    return Fraction(value)


def reject_inexact(literal):
    raise ValueError(f"Unexpected inexact JSON number: {literal}")


def read_certificate(path, expected_hashes):
    data = path.read_bytes()
    cert = json.loads(data, parse_constant=reject_inexact, parse_float=reject_inexact)
    require(cert["input_sha256"] == expected_hashes, "Checker or config digest mismatch")
    require(cert["state"] == "prepared lapse endpoint identities verified", "Unexpected endpoint state")
    provenance = {"report": "paper1/" + REPORT.as_posix(), "report_commit": REPORT_COMMIT,
                  "report_sha256": REPORT_SHA256, "source_of_derivation_commit": MATH_REVISION}
    require(cert["reviewed_mathematical_definitions"] == provenance
            and cert["reviewed_report"] == provenance["report"]
            and cert["reviewed_report_commit"] == REPORT_COMMIT
            and cert["reviewed_report_sha256"] == REPORT_SHA256
            and cert["source_of_derivation_commit"] == MATH_REVISION, "Unexpected source links")
    runtime = cert["runtime"]
    version = runtime["python_version_info"]
    require(set(version) == {"major", "minor", "micro"} and
            all(type(n) is int and n >= 0 for n in version.values()), "Invalid child version info")
    child_version = f"3.12.{version['micro']}"
    require(version["major"] == 3 and version["minor"] == 12 and
            runtime["python_version"] == child_version and runtime["sys_version"].split()[:1] == [child_version]
            and runtime["sympy_version"] == "1.14.0" and
            isinstance(runtime["python_executable"], str) and Path(runtime["python_executable"]).is_absolute(),
            "Unexpected child interpreter or SymPy version")
    params = cert["fixed_parameters"]
    require(params == {"N": 3, "m": 0, "ell": 2, "epsilon": "1/100",
                       "radial_shell": ["1", "2"], "Y": "3*x**2/2 - 1/2"}
            and all(type(params[k]) is int for k in ("N", "m", "ell")), "Wrong selected sector")
    require(cert["normalization"] == {"scalar_curvature_identity_prefactor": "1/(8*pi)",
            "hessian_identity_prefactor": "1/(4*pi)"}, "Wrong normalized prefactors")
    for key, names, flag in (("identities", IDENTITIES, True), ("scope", SCOPE, False)):
        require(cert[key].keys() == set(names) and all(cert[key][name] is flag for name in names),
                f"Unexpected {key} claims")

    cases = passed_cases(cert["case_results"], set(ZERO_FIELDS) | {
        "sign_and_factor_defect_controls", "endpoint_orientation_signs"})
    controls = passed_cases(cert["controls"], CONTROL_NAMES)
    for name, fields in ZERO_FIELDS.items():
        require(all(cases[name][field] == "0" for field in fields), f"Nonzero residual in {name}")
    for key in ("scalar_dirichlet_radial", "hessian_dirichlet_radial"):
        exact_formula(cert["endpoint_fluxes_before_prefactors"][key])
    for key in ("scalar_wrong_sign", "hessian_wrong_factor"):
        exact_formula(cases["sign_and_factor_defect_controls"][key])

    samples = cert["sign_samples"]
    require((samples["x"], samples["inner_r"], samples["outer_r"]) == ("0", "1", "2")
            and samples["orientation_reverses_inner_endpoint"] is True, "Wrong sign-control points")
    oriented = cases["endpoint_orientation_signs"]
    for identity, sign in (("scalar", 1), ("hessian", -1)):
        inner = exact_sample(samples[identity + "_inner_coordinate"])
        outer = exact_sample(samples[identity + "_outer_coordinate"])
        require(sign * inner > 0 and sign * outer > 0, f"Wrong {identity} flux signs")
        require(exact_sample(oriented[identity + "_outer"]) == outer and
                exact_sample(oriented[identity + "_oriented_inner"]) == -inner,
                f"Inconsistent oriented {identity} flux")
    metrics = {"endpoint_checks_passed": 1, "endpoint_case_count": len(cases),
               "control_count": len(controls)}
    return hashlib.sha256(data).hexdigest(), cert["state"], metrics


def fallback_metrics(record):
    return {"endpoint_checks_passed": 0, "scientific_invocations": record["scientific_invocations"]}


def publish_metrics(metrics):
    numeric = {key: value for key, value in metrics.items() if value is not None}
    if metrics["endpoint_checks_passed"]:
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
                        record.update(certificate_sha256=digest, endpoint_status=status)
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
