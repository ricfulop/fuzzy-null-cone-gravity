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


def manufactured_certificate(hashes, status="silent_on_forward_state_germ_subfamily"):
    formal, measured, certified = {
        "silent_on_forward_state_germ_subfamily": (False, False, False),
        "response_nonzero_at_certified_r0": (True, True, True),
        "formal_response_nonzero_unresolved": (True, None, False),
    }[status]
    states = {"silent_on_forward_state_germ_subfamily": "state_germ_forward_subfamily_silent",
              "response_nonzero_at_certified_r0": "state_germ_with_nonzero_discarded_response",
              "formal_response_nonzero_unresolved": "state_germ_response_unresolved"}
    regular = {"r0": "1", "terms": 3, "atan_interval": ["1/12", "1/10"],
               "denominator_count": 1, "certified_denominators": [{"sign": 1}]}
    response = {"status": status, "formal_response_nonzero": formal, "response_nonzero": measured,
                "certified_response_nonzero": certified, "stationary_universal_zero_claim": False}
    if certified:
        response.update(responding_initial_vector_supplied=True, r0="1", certified_column=2,
                        unit_vector_kind="free_state", unit_vector_index=2,
                        free_state_initial=["0", "0", "1", "0", "0", "0", "0", "0"],
                        certified_numerator={"column": 2, "r0": "1", "sign": 1, "interval": ["1", "1"]},
                        atan_interval=regular["atan_interval"], regularity_certificate=regular.copy())
    return {"state": states[status], "input_sha256": hashes, "state_germ_certificate_passed": True,
            "epsilon0": "1/100", "free_state_count": 8,
            "manufactured_controls": {"cases": ["toy"], "case_count": 1, "state": "direct state-germ manufactured controls passed"}, "smooth_solution_germ_certified": True,
            "stationary_universal_zero_claim": False, "nonlinear_consistency_established": False,
            "physical_boundary_domain_certified": False, "final_dirac_classification_established": False,
            "formal_discarded_response_nonzero": formal, "discarded_response_nonzero": measured,
            "certified_discarded_response_nonzero": certified, "admissible_non_gauge_seed_supplied": certified,
            "response_certificate": response, "B_shape": [14, 0], "remaining_control_count": 0,
            "proof_matrix_flags": dict.fromkeys(("CM_zero", "JM_identity", "SM_zero", "profile_lift_zero",
                                                 "all_22_retained_rows_zero"), True),
            "retained_euler_certificate": {"all_22_retained_rows_zero": True, "retained_row_count": 22, "x0_equation_included_in_22_retained_rows": True,
                                           "operator_identity": "(E0 Pv + E1 Pw + E2 Pw K) M"},
            "regularity_certificate": regular, "same_radius_witness_required": True,
            "constraint_rank_certificate": {"rank": 6, "minor_columns": [0, 1, 2, 7, 8, 9],
                                            "r0": "1", "terms": 3, "atan_interval": regular["atan_interval"],
                                            "row_swap_count": 0, "determinant_interval": ["-2", "-1"],
                                            "determinant_sign": -1}}



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
                       "source_sha256_before": {}, "input_sha256_before": {}, "checker_input_sha256": {}}
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
        self.assertEqual(evidence["metrics"]["state_germ_certificate_passed"], 0)
        self.assertIsNone(evidence["metrics"]["discarded_response_nonzero"])

    def test_zero_exit_observed_after_deadline_fails(self):
        with patch.object(runner.time, "monotonic", side_effect=[0, 11, 11]):
            evidence, _ = self.supervise()
        self.assertEqual(evidence["exit_code"], 0)
        self.assertEqual(evidence["stop_reason"], "wall budget reached")
        self.assertEqual(evidence["metrics"]["state_germ_certificate_passed"], 0)

    def test_unresolved_certificate_is_a_pass_with_unknown_response(self):
        status = "formal_response_nonzero_unresolved"
        self.certificate.write_text(json.dumps(manufactured_certificate({}, status)))
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["response_status"], status)
        self.assertEqual(evidence["certificate_sha256"], runner.sha(self.certificate))
        self.assertIsNone(evidence["metrics"]["discarded_response_nonzero"])
        self.assertEqual(evidence["metrics"]["certified_response_nonzero"], 0)

    def test_certificate_requires_pins_and_consistent_boolean_flags(self):
        for field, value in (("input_sha256", {"config": "different"}), ("state_germ_certificate_passed", 1),
                             ("discarded_response_nonzero", None), ("stationary_universal_zero_claim", True)):
            with self.subTest(field=field):
                cert = manufactured_certificate({})
                cert[field] = value
                self.certificate.write_text(json.dumps(cert))
                with self.assertRaises(ValueError):
                    runner.read_certificate(self.certificate, {})

    def test_rejects_missing_proof_and_incompatible_witness(self):
        for group, key, value in (
            ("proof_matrix_flags", "SM_zero", False),
            ("retained_euler_certificate", "x0_equation_included_in_22_retained_rows", False),
            ("constraint_rank_certificate", "minor_columns", [0]*6),
            ("constraint_rank_certificate", "determinant_interval", ["-1", "0"]),
            ("constraint_rank_certificate", "r0", "2"),
            ("response_certificate", "r0", "2"),
            ("response_certificate", "certified_numerator", {"r0": "3"}),
            ("response_certificate", "free_state_initial", [0]*8),
            ("response_certificate", "unit_vector_index", 1),
            ("response_certificate", "unit_vector_kind", "control"),
        ):
            with self.subTest(group=group, key=key):
                cert = manufactured_certificate({}, "response_nonzero_at_certified_r0")
                cert[group][key] = value
                self.certificate.write_text(json.dumps(cert))
                with self.assertRaises(ValueError):
                    runner.read_certificate(self.certificate, {})

    def test_silent_unresolved_and_responding_metrics(self):
        for status, measured, certified in (
            ("silent_on_forward_state_germ_subfamily", 0, 0),
            ("formal_response_nonzero_unresolved", None, 0),
            ("response_nonzero_at_certified_r0", 1, 1),
        ):
            with self.subTest(status=status):
                self.certificate.write_text(json.dumps(manufactured_certificate({}, status)))
                digest, actual, metrics = runner.read_certificate(self.certificate, {})
                self.assertEqual(digest, runner.sha(self.certificate))
                self.assertEqual(actual, status)
                self.assertEqual(metrics["discarded_response_nonzero"], measured)
                self.assertEqual(metrics["certified_response_nonzero"], certified)
                self.assertEqual(metrics["state_germ_certificate_passed"], 1)
                self.assertEqual(metrics["response_determined"], int(measured is not None))
                with patch.object(runner, "report") as emit:
                    runner.publish_metrics(metrics)
                numeric = emit.call_args.kwargs["m"]
                self.assertEqual("discarded_response_nonzero" in numeric, measured is not None)


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
        self.here = self.root / "verification/germ"
        self.here.mkdir(parents=True)
        self.commit = "a" * 40
        self.out = self.root / "run-output/mock-run"
        self.config = {"inputs": {}}
        for key, filename in runner.INPUT_NAMES.items():
            path = self.root / filename
            data = json.dumps({"manufactured": key}).encode()
            path.write_bytes(data)
            self.config["inputs"][key] = {"path": filename, "sha256": hashlib.sha256(data).hexdigest()}
        for name in ("run.sh", "run_bounded.py", "test_run_bounded.py", "certify_state_germ.py",
                     "exact_decimal.py", "test_exact_decimal.py", "view.xml", "view-config.json", "README.md"):
            (self.here / name).write_text("mock source " + name)
        (self.here / "run-config.json").write_text(json.dumps(self.config))
        self.committed = {path.relative_to(self.root).as_posix(): path.read_bytes()
                          for path in (self.root / "verification").rglob("*") if path.is_file()}

        def git(args, **kwargs):
            self.assertEqual(args[0], "git")
            if args[3] == "rev-parse":
                return self.commit
            self.assertEqual(args[3], "show")
            full_id, path = args[4].split(":paper1/")
            self.assertEqual(full_id, self.commit)
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
        self.assertEqual((record["wall_seconds"], record["memory_bytes"]), (900, 2 * 1024**3))
        argv = record["command"]
        self.assertEqual(argv[argv.index("python") + 1], str(out / "sources/verification/germ/certify_state_germ.py"))
        self.assertEqual(argv[argv.index("--max-atan-terms") + 1], "256")
        self.assertEqual(argv[argv.index("--output") + 1], str(out / runner.CERTIFICATE))
        for key, flag in (("stationary_operators", "--operators"), ("terminal_state", "--state")):
            saved = Path(argv[argv.index(flag) + 1])
            self.assertEqual(saved, out / "inputs" / runner.INPUT_NAMES[key])
            self.assertEqual(runner.sha(saved), self.config["inputs"][key]["sha256"])
        self.assertEqual(runner.sha(Path(argv[argv.index("--config") + 1])), record["checker_input_sha256"]["config"])
        self.assertEqual(env["PYTHONOPTIMIZE"], "0")
        errors = []
        runner.check_hashes_after(out, record, errors)
        self.assertEqual(errors, [])
        for group in ("source", "input"):
            self.assertEqual(record[group + "_sha256_before"], record[group + "_sha256_after"])
            self.assertEqual(record[group + "_sha256_before"], record[group + "_snapshot_sha256_after"])
        (out / "inputs" / runner.INPUT_NAMES["terminal_state"]).write_bytes(b"changed copy")
        (self.root / runner.INPUT_NAMES["stationary_operators"]).write_bytes(b"changed original")
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
        self.assertIsNone(evidence["metrics"]["discarded_response_nonzero"])
        return evidence["supervisor_error"]

    def test_input_pin_mismatch_cannot_dispatch(self):
        (self.root / runner.INPUT_NAMES["terminal_state"]).write_bytes(b"wrong terminal state")
        self.assertIn("pin mismatch", self.check_preflight_failure())

    def test_modified_source_cannot_dispatch(self):
        (self.here / "certify_state_germ.py").write_bytes(b"wrong source")
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
