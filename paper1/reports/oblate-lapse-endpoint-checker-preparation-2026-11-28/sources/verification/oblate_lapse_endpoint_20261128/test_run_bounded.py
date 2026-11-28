"""Mocked supervisor regressions. All child, ps, and signal operations are replaced."""

from contextlib import ExitStack
from datetime import datetime, timedelta, timezone
from itertools import chain, count
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import io
import json
import unittest
from unittest.mock import Mock, call, patch

import run_bounded as runner


def manufactured_certificate(hashes):
    """Synthetic wire-format example. No scientific functions or data are invoked."""
    rows = [{"name": name, "passed": True, **dict.fromkeys(fields, "0")}
            for name, fields in runner.ZERO_FIELDS.items()]
    rows.extend([{"name": "sign_and_factor_defect_controls", "passed": True,
                  "scalar_wrong_sign": "2*x", "hessian_wrong_factor": "x/2"},
                 {"name": "endpoint_orientation_signs", "passed": True, "scalar_outer": "2/3",
                  "scalar_oriented_inner": "-1", "hessian_outer": "-5/3", "hessian_oriented_inner": "2"}])
    report = "paper1/" + runner.REPORT.as_posix()
    return {
        "state": "prepared lapse endpoint identities verified", "input_sha256": hashes,
        "reviewed_report": report, "reviewed_report_commit": runner.REPORT_COMMIT,
        "reviewed_report_sha256": runner.REPORT_SHA256, "source_of_derivation_commit": runner.MATH_REVISION,
        "reviewed_mathematical_definitions": {"report": report, "report_commit": runner.REPORT_COMMIT,
            "report_sha256": runner.REPORT_SHA256, "source_of_derivation_commit": runner.MATH_REVISION},
        "runtime": {"sys_version": "3.12.11 (mock build)", "python_version": "3.12.11", "python_executable": "/mock/bin/python",
                    "python_version_info": {"major": 3, "minor": 12, "micro": 11},
                    "sympy_version": "1.14.0"},
        "fixed_parameters": {"N": 3, "m": 0, "ell": 2, "epsilon": "1/100",
                             "radial_shell": ["1", "2"], "Y": "3*x**2/2 - 1/2"},
        "normalization": {"scalar_curvature_identity_prefactor": "1/(8*pi)",
                          "hessian_identity_prefactor": "1/(4*pi)"},
        "scope": dict.fromkeys(runner.SCOPE, False), "identities": dict.fromkeys(runner.IDENTITIES, True),
        "controls": [{"name": name, "passed": True} for name in sorted(runner.CONTROL_NAMES)],
        "case_results": rows,
        "endpoint_fluxes_before_prefactors": {"scalar_dirichlet_radial": "F1*r*x",
                                             "hessian_dirichlet_radial": "-F1*G1/r"},
        "sign_samples": {"x": "0", "outer_r": "2", "inner_r": "1", "scalar_outer_coordinate": "2/3",
                         "scalar_inner_coordinate": "1", "hessian_outer_coordinate": "-5/3",
                         "hessian_inner_coordinate": "-2", "orientation_reverses_inner_endpoint": True},
    }


class GroupProbeTests(unittest.TestCase):
    def test_permission_denied_requires_a_readable_process_table(self):
        cases = [
            ("1 1 50\n2 2 80\n", "empty"),
            ("1 1 50\n111 77 80\n", "live"),
            ("", "unknown"),
            ("1 1 50\ntruncated\n", "unknown"),
            (runner.subprocess.TimeoutExpired("ps", 1), "unknown"),
        ]
        for table, expected in cases:
            with self.subTest(table=table), \
                 patch.object(runner.os, "killpg", side_effect=PermissionError()), \
                 patch.object(runner.subprocess, "check_output") as ps:
                if isinstance(table, str):
                    ps.return_value = table
                else:
                    ps.side_effect = table
                state = runner.group_state(77)
                self.assertEqual(state["status"], expected)
                self.assertEqual(state["probe"], "EPERM")
                if expected == "live":
                    self.assertEqual(state["members"], [{"pid": 111, "rss_bytes": 81920}])
                self.assertGreater(ps.call_args.kwargs["timeout"], 0)

    def test_stubborn_descendant_escalated_after_launcher_exit(self):
        process = Mock(pid=77)
        process.poll.return_value = 0
        clock = count()
        events = []
        with patch.object(runner.os, "killpg", side_effect=PermissionError()) as signals, \
             patch.object(runner, "group_state", return_value={"status": "live"}), \
             patch.object(runner.time, "monotonic", side_effect=lambda: next(clock)), \
             patch.object(runner.time, "sleep"):
            runner.stop_group(process, events)
        self.assertEqual(signals.call_args_list, [call(77, runner.signal.SIGTERM),
                                                 call(77, runner.signal.SIGKILL)])
        self.assertEqual(events[-1]["status"], "live")


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.patches = ExitStack()
        self.addCleanup(self.patches.close)
        self.out = Path(self.patches.enter_context(TemporaryDirectory()))
        self.handlers = {}
        self.received_during_install = False

        def install(sig, handler):
            self.handlers[sig] = handler
            if self.received_during_install and callable(handler):
                handler(sig, None)
            return runner.signal.SIG_DFL

        self.replace(runner.signal, "signal", side_effect=install)
        self.killpg = self.replace(runner.os, "killpg", side_effect=ProcessLookupError())
        self.ps = self.replace(runner.subprocess, "check_output", return_value="55 55 50\n")
        self.child = Mock(pid=77, returncode=0, stdout=io.BytesIO())
        self.child.poll.return_value = 0
        self.launch = self.replace(runner.subprocess, "Popen", return_value=self.child)
        self.stop_implementation = runner.stop_group
        self.stop = self.replace(runner, "stop_group")
        self.console = self.replace(runner.os, "write", side_effect=lambda fd, data: len(data))
        console_stream = self.replace(runner.sys, "stdout", new=Mock())
        console_stream.fileno.return_value = 1
        self.replace(runner.os, "get_blocking", return_value=True)
        self.replace(runner.os, "set_blocking")
        self.record = {"wall_seconds": 10, "memory_bytes": 2048,
                       "scientific_invocations": 0, "command": ["dummy"],
                       "source_sha256_before": {}, "checker_input_sha256": {}}
        self.certificate = self.out / runner.CERTIFICATE
        self.certificate.parent.mkdir()
        self.certificate.write_text(json.dumps(manufactured_certificate({})))

    def replace(self, obj, name, **kwargs):
        return self.patches.enter_context(patch.object(obj, name, **kwargs))

    def supervise(self):
        rc = runner.supervise(self.out, self.record, {},
                              datetime.now(timezone.utc) + timedelta(minutes=10))
        execution = json.loads((self.out / "execution.json").read_bytes())
        samples = json.loads((self.out / "resource-samples.json").read_bytes())
        self.assertEqual(execution["metrics"], json.loads((self.out / "metrics.json").read_bytes()))
        self.assertEqual(execution["rss_sample_count"], len(samples))
        self.assertEqual(rc == 0, execution["status"] == "passed")
        self.assertEqual(set(self.handlers.values()), {runner.signal.SIG_DFL})
        if self.launch.called and execution["scientific_invocations"]:
            self.assertGreater(self.child.wait.call_args.kwargs["timeout"], 0)
        return execution, samples

    def test_normal_exit_eperm_empty_and_exact_raw_log(self):
        self.killpg.side_effect = PermissionError()
        self.child.stdout = io.BytesIO(b"stage\nlast")
        self.console.side_effect = lambda fd, data: 1   # partial writes
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["group_disposition"]["status"], "empty")
        self.assertEqual((self.out / "verification.log").read_bytes(), b"stage\nlast")
        self.assertGreaterEqual(self.console.call_count, len(b"stage\nlast"))
        self.stop.assert_not_called()

    def test_timeout_preserves_measured_sample_when_cleanup_raises(self):
        self.record["wall_seconds"] = 1
        clock = chain([0, 0.5], count(2))
        self.replace(runner.time, "monotonic", side_effect=lambda: next(clock))
        self.child.poll.return_value = None
        self.child.returncode = None
        self.ps.return_value = "77 77 12\n"
        self.killpg.side_effect = PermissionError()
        self.stop.side_effect = PermissionError("manufactured termination error")
        self.child.wait.side_effect = runner.subprocess.TimeoutExpired("dummy", 0.5)
        evidence, samples = self.supervise()
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["stop_reason"], "wall budget reached")
        self.assertEqual(evidence["group_disposition"]["status"], "live")
        self.assertEqual([e["action"] for e in evidence["cleanup_errors"]],
                         ["terminate_group", "wait_child"])
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["rss_bytes"], 12 * 1024)

    def test_timeout_with_eperm_after_term_keeps_both_records(self):
        self.record["wall_seconds"] = 1
        clock = chain([0, 0.5], count(2))
        self.replace(runner.time, "monotonic", side_effect=lambda: next(clock))
        self.child.poll.side_effect = [None, -15]
        self.child.returncode = -15
        self.stop.side_effect = self.stop_implementation
        self.killpg.side_effect = [None, None, PermissionError(), PermissionError()]
        self.ps.side_effect = ["77 77 12\n", "1 1 50\n", "1 1 50\n"]
        evidence, samples = self.supervise()
        self.assertEqual(evidence["group_disposition"]["status"], "empty")
        self.assertEqual(evidence["group_disposition"]["probe"], "EPERM")
        self.assertEqual(evidence["cleanup_errors"], [])
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["stop_reason"], "wall budget reached")
        self.assertEqual(evidence["sampled_max_rss_bytes"], samples[0]["rss_bytes"])
        self.assertEqual(self.killpg.call_args_list,
                         [call(77, 0), call(77, runner.signal.SIGTERM), call(77, 0), call(77, 0)])

    def test_zero_child_exit_with_unknown_group_is_failure(self):
        self.killpg.side_effect = PermissionError()
        self.ps.side_effect = runner.subprocess.TimeoutExpired("ps", 1)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["exit_code"], 0)
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["group_disposition"]["status"], "unknown")
        self.stop.assert_called_once()

    def test_signal_before_spawn(self):
        self.received_during_install = True
        evidence, _ = self.supervise()
        self.launch.assert_not_called()
        self.assertTrue(evidence["received_signals"])
        self.assertEqual(evidence["scientific_invocations"], 0)
        self.assertEqual(evidence["group_disposition"]["status"], "not_spawned")

    def test_signal_during_spawn_still_cleans_child_and_drains_log(self):
        def launch(*args, **kwargs):
            self.handlers[runner.signal.SIGTERM](runner.signal.SIGTERM, None)
            return self.child

        self.launch.side_effect = launch
        self.killpg.side_effect = [None, ProcessLookupError()]
        self.child.stdout = io.BytesIO(b"before stop\n")
        evidence, samples = self.supervise()
        self.assertEqual(evidence["received_signals"], [int(runner.signal.SIGTERM)])
        self.assertEqual(evidence["scientific_invocations"], 1)
        self.assertEqual(evidence["stop_reason"], "supervisor received signal")
        self.assertEqual(evidence["group_disposition"]["status"], "empty")
        self.assertFalse(evidence["relay_active_after_join"])
        self.assertEqual((self.out / "verification.log").read_bytes(), b"before stop\n")
        self.assertEqual(samples, [])
        self.stop.assert_called_once()

    def test_stalled_relay_join_is_bounded_and_fails_closed(self):
        relay = self.replace(runner.threading, "Thread").return_value
        relay.is_alive.return_value = True
        evidence, _ = self.supervise()
        self.assertGreater(relay.join.call_args.kwargs["timeout"], 0)
        self.assertTrue(evidence["relay_active_after_join"])
        self.assertEqual(evidence["status"], "failed")
        # The mock did not start the worker that ordinarily owns these streams.
        relay_target = runner.threading.Thread.call_args.kwargs["target"]
        relay_target()


    def test_zero_exit_without_certificate_fails(self):
        self.certificate.unlink()
        evidence, _ = self.supervise()
        self.assertEqual(evidence["stop_reason"], "terminal certificate validation failed")
        self.assertEqual(evidence["metrics"]["endpoint_checks_passed"], 0)
        self.assertNotIn("endpoint_case_count", evidence["metrics"])

    def test_zero_exit_observed_after_deadline_fails(self):
        with patch.object(runner.time, "monotonic", side_effect=[0, 11, 11]):
            evidence, _ = self.supervise()
        self.assertEqual(evidence["exit_code"], 0)
        self.assertEqual(evidence["stop_reason"], "wall budget reached")
        self.assertEqual(evidence["metrics"]["endpoint_checks_passed"], 0)

    def test_endpoint_payload_rejections(self):
        mutations = [
            (("input_sha256",), {"checker": "changed"}),
            (("state",), "failed"),
            (("runtime", "python_version"), "3.13.0"),
            (("runtime", "sys_version"), "3.11.9 (mock build)"),
            (("runtime", "irrelevant_value"), float("nan")),
            (("runtime", "sympy_version"), "1.13.0"),
            (("runtime", "python_version_info", "micro"), True),
            (("reviewed_report_commit",), runner.MATH_REVISION),
            (("reviewed_mathematical_definitions", "report_sha256"), ""),
            (("fixed_parameters", "epsilon"), 0.01),
            (("fixed_parameters", "m"), False),
            (("normalization", "scalar_curvature_identity_prefactor"), "1/(4*pi)"),
            (("identities", runner.IDENTITIES[0]), 1), (("scope", runner.SCOPE[0]), True),
            (("case_results", 0, "passed"), False), (("case_results", 0, "residual"), "1"),
            (("case_results", 7, "south_residual"), 0),
            (("case_results", -1, "scalar_oriented_inner"), "1"),
            (("case_results", -2, "scalar_wrong_sign"), "0"),
            (("case_results", 0, "name"), "background_ricci_zero"),
            (("controls",), []), (("controls", 0, "passed"), 1),
            (("sign_samples", "hessian_outer_coordinate"), "1"),
            (("sign_samples", "scalar_inner_coordinate"), "nan"),
            (("sign_samples", "hessian_inner_coordinate"), "-2.0"),
            (("endpoint_fluxes_before_prefactors", "scalar_dirichlet_radial"), "oo*x"),
            (("endpoint_fluxes_before_prefactors", "scalar_dirichlet_radial"), "0.5*x"),
        ]
        for keys, value in mutations:
            with self.subTest(keys=keys):
                cert = manufactured_certificate({})
                parent = cert
                for key in keys[:-1]:
                    parent = parent[key]
                parent[keys[-1]] = value
                self.certificate.write_text(json.dumps(cert))
                with self.assertRaises(ValueError):
                    runner.read_certificate(self.certificate, {})

    def test_pass_metrics_follow_authenticated_terminal_records(self):
        def inspect(metrics):
            self.assertEqual(metrics["endpoint_checks_passed"], 1)
            self.assertEqual(metrics["endpoint_case_count"], 11)
            self.assertEqual(metrics["control_count"], 8)
            marker = json.loads((self.out / "terminal-evidence.json").read_bytes())
            self.assertEqual(marker["status"], "passed")
            self.assertEqual(marker["execution_sha256"], runner.sha(self.out / "execution.json"))
            execution = json.loads((self.out / "execution.json").read_bytes())
            self.assertEqual(execution["certificate_sha256"], runner.sha(self.certificate))
        emit = self.replace(runner, "publish_metrics", side_effect=inspect)
        self.supervise()
        emit.assert_called_once()

    def test_resource_stop_keeps_sample_and_fails(self):
        self.child.poll.side_effect = [None, 0]
        self.ps.return_value = "77 77 9\n"
        evidence, samples = self.supervise()
        self.assertEqual(evidence["stop_reason"], "sampled process-group RSS exceeded threshold")
        self.assertEqual(samples[0]["rss_bytes"], 9 * 1024)
        self.assertEqual(evidence["metrics"]["endpoint_checks_passed"], 0)

    def test_existing_terminal_json_is_preserved(self):
        path = self.out / "standalone.json"
        runner.save_json(path, {"value": 1})
        with self.assertRaises(FileExistsError):
            runner.save_json(path, {"value": 2})
        self.assertEqual(json.loads(path.read_text()), {"value": 1})
        self.assertFalse(path.with_suffix(".json.tmp").exists())

    def test_console_backpressure_preserves_archived_log(self):
        self.console.side_effect = BlockingIOError("console reader stalled")
        self.child.stdout = io.BytesIO(b"archived regardless of console\n" * 1024)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["relay_errors"], [])
        self.assertEqual((self.out / "verification.log").read_bytes(),
                         b"archived regardless of console\n" * 1024)

    def test_empty_rss_sample_set_is_unknown(self):
        evidence, samples = self.supervise()
        self.assertEqual(samples, [])
        self.assertIsNone(evidence["sampled_max_rss_bytes"])

    def test_sampling_timeout_is_capped_by_remaining_budget(self):
        self.record["wall_seconds"] = 1
        self.child.poll.side_effect = [None, 0]
        self.replace(runner.time, "sleep")
        with patch.object(runner.time, "monotonic", side_effect=[0, 0.75, 0.8, 0.9, 1]):
            evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(self.ps.call_args.kwargs["timeout"], 0.25)

    def test_signal_after_terminal_boundary_is_outside_the_decision(self):
        original_save = runner.save_json
        def inject(path, payload):
            if path.name == "metrics.json":
                self.handlers[runner.signal.SIGTERM](runner.signal.SIGTERM, None)
            return original_save(path, payload)
        self.replace(runner, "save_json", side_effect=inject)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["received_signals"], [])
        self.assertEqual(evidence["cancellation_boundary"], "before terminal evidence publication")

    def test_terminal_marker_authenticates_both_records(self):
        evidence, _ = self.supervise()
        marker = json.loads((self.out / "terminal-evidence.json").read_bytes())
        self.assertEqual(marker["status"], evidence["status"])
        for key, filename in [("execution_sha256", "execution.json"), ("metrics_sha256", "metrics.json")]:
            self.assertEqual(marker[key], hashlib.sha256((self.out / filename).read_bytes()).hexdigest())

    def test_failed_execution_write_cannot_commit_pass_metrics(self):
        original_save = runner.save_json
        def fail_execution(path, payload):
            if path.name == "execution.json":
                raise OSError("synthetic execution-write failure")
            return original_save(path, payload)
        self.replace(runner, "save_json", side_effect=fail_execution)
        with self.assertRaisesRegex(OSError, "synthetic execution-write"):
            self.supervise()
        self.assertFalse((self.out / "terminal-evidence.json").exists())

    def test_signal_during_certificate_validation_still_fails(self):
        original_read = runner.read_certificate
        def cancel(*args, **kwargs):
            self.handlers[runner.signal.SIGTERM](runner.signal.SIGTERM, None)
            return original_read(*args, **kwargs)
        self.replace(runner, "read_certificate", side_effect=cancel)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["received_signals"], [int(runner.signal.SIGTERM)])


class DispatchTests(unittest.TestCase):
    """Only synthetic files enter the supervisor. Git, process and clock are mocked."""

    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(TemporaryDirectory()))
        self.here = self.root / "verification/endpoint"
        self.here.mkdir(parents=True)
        self.commit = "a" * 40
        self.out = self.root / "run-output/mock-run"
        self.config = json.loads((runner.HERE / "run-config.json").read_bytes())
        report = self.root / runner.REPORT
        report.parent.mkdir(parents=True)
        report.write_bytes(b"mock report")
        report_hash = runner.sha(report)
        self.stack.enter_context(patch.object(runner, "REPORT_SHA256", report_hash))
        self.config["report"]["sha256"] = report_hash
        for name in runner.SOURCE_NAMES:
            (self.here / name).write_text("mock source " + name)
        (self.here / "run-config.json").write_text(json.dumps(self.config))
        self.committed = {path.relative_to(self.root).as_posix(): path.read_bytes()
                          for path in self.root.rglob("*") if path.is_file()}

        def git(args, **kwargs):
            self.assertEqual(args[0], "git")
            if args[3] == "rev-parse":
                return self.commit
            self.assertEqual(args[3], "show")
            full_id, path = args[4].split(":paper1/")
            self.assertIn(full_id, (self.commit, runner.REPORT_COMMIT))
            return self.committed[path]

        self.stack.enter_context(patch.object(runner, "ROOT", self.root))
        self.stack.enter_context(patch.object(runner, "HERE", self.here))
        self.stack.enter_context(patch.object(runner.subprocess, "check_output", side_effect=git))
        self.popen = self.stack.enter_context(patch.object(runner.subprocess, "Popen", side_effect=AssertionError("unexpected child")))
        self.dispatch = self.stack.enter_context(patch.object(runner, "supervise", return_value=0))
        self.stack.enter_context(patch.object(runner, "report"))
        self.stack.enter_context(patch.object(runner.platform, "platform", return_value="mock host"))
        self.now = self.stack.enter_context(patch.object(runner, "datetime"))
        self.now.now.return_value = datetime(2026, 11, 28, tzinfo=timezone.utc)
        self.stack.enter_context(patch.dict(runner.os.environ, {"ORX_RUN_ID": "mock-run",
                            "ORX_COMMIT": self.commit, "PYTHONOPTIMIZE": "0"}))

    def test_only_snapshot_bytes_dispatched_and_fresh_directory_required(self):
        self.assertEqual(runner.main(), 0)
        out, record, env, cutoff = self.dispatch.call_args.args
        self.assertEqual((out, cutoff), (self.out, runner.CUTOFF))
        self.assertEqual((record["wall_seconds"], record["memory_bytes"]), (300, 2 * 1024**3))
        argv = record["command"]
        self.assertEqual(argv[:8], ["uv", "run", "--isolated", "--no-project", "--python", "3.12",
                                    "--with", "sympy==1.14.0"])
        frozen = out / "sources/verification/endpoint"
        self.assertEqual(argv[argv.index("python") + 1:],
                         [str(frozen / "check_lapse_endpoints.py"), "--config",
                          str(frozen / "run-config.json"), "--output", str(out / runner.CERTIFICATE)])
        self.assertEqual(record["checker_input_sha256"],
                         {"checker": runner.sha(frozen / "check_lapse_endpoints.py"),
                          "config": runner.sha(frozen / "run-config.json")})
        self.assertEqual(env["PYTHONOPTIMIZE"], "0")
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "MKL_NUM_THREADS"):
            self.assertEqual(env[key], "1")
        errors = []
        runner.check_hashes_after(out, record, errors)
        self.assertEqual(errors, [])
        self.assertEqual(record["source_sha256_before"], record["source_sha256_after"])
        self.assertEqual(record["source_sha256_before"], record["source_snapshot_sha256_after"])
        (out / "sources" / runner.REPORT).write_bytes(b"changed copy")
        (self.here / "check_lapse_endpoints.py").write_bytes(b"changed original")
        runner.check_hashes_after(out, record, errors)
        self.assertEqual(len(errors), 2)
        with self.assertRaises(FileExistsError):
            runner.main()
        self.assertEqual(self.dispatch.call_count, 1)

    def check_preflight_failure(self):
        self.assertEqual(runner.main(), 1)
        self.dispatch.assert_not_called()
        self.popen.assert_not_called()
        evidence = json.loads((self.out / "execution.json").read_bytes())
        self.assertEqual(evidence["scientific_invocations"], 0)
        self.assertEqual(evidence["group_disposition"]["status"], "not_spawned")
        self.assertNotIn("endpoint_case_count", evidence["metrics"])
        return evidence["supervisor_error"]

    def test_report_mismatch_cannot_dispatch(self):
        (self.root / runner.REPORT).write_bytes(b"wrong report")
        self.assertIn("pin mismatch", self.check_preflight_failure())

    def test_modified_source_cannot_dispatch(self):
        (self.here / "check_lapse_endpoints.py").write_bytes(b"wrong source")
        self.assertIn("pin mismatch", self.check_preflight_failure())

    def test_expired_window_cannot_dispatch(self):
        self.now.now.return_value = runner.CUTOFF
        self.assertIn("cutoff", self.check_preflight_failure())

    def test_missing_identity_cannot_dispatch(self):
        for key in ("ORX_RUN_ID", "ORX_COMMIT"):
            with self.subTest(key=key), patch.dict(runner.os.environ, {key: ""}):
                with self.assertRaises(ValueError):
                    runner.main()
        self.dispatch.assert_not_called()
        self.assertFalse(self.out.exists())


if __name__ == "__main__":
    unittest.main()
