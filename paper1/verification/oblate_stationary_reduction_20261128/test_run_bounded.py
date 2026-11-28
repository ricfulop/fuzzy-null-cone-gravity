"""Mocked supervisor regressions. All child, ps, and signal operations are replaced."""

from contextlib import ExitStack
from datetime import datetime, timedelta, timezone
from itertools import chain, count
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib
import io
import json
import unittest
from unittest.mock import Mock, call, patch

import run_bounded as runner


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


class FrozenSourceTests(unittest.TestCase):
    def test_main_executes_the_authenticated_copy(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            here = root / "verification" / "reducer"
            here.mkdir(parents=True)
            config = {"input_directory": "inputs", "wall_seconds": 10,
                      "memory_bytes": 2048, "stage": "principal", "max_steps": 1,
                      "approval_cutoff_utc": (datetime.now(timezone.utc) +
                                              timedelta(minutes=5)).isoformat()}
            (here / "run-config.json").write_text(json.dumps(config))
            (here / "run.sh").write_text("dummy entrypoint")
            (here / "reduce_stationary.py").write_text("dummy checker; never executed")
            (root / "inputs").mkdir()
            (root / "inputs" / "completion.json").write_bytes(b"{}")
            committed = {"paper1/" + str(path.relative_to(root)): path.read_bytes()
                         for path in here.iterdir()}
            def git(command, **kwargs):
                return "a" * 40 if "rev-parse" in command else committed[command[-1].split(":", 1)[1]]
            with patch.object(runner, "ROOT", root), patch.object(runner, "HERE", here), \
                 patch.object(runner, "UPSTREAM_COMPLETION_SHA256", hashlib.sha256(b"{}").hexdigest()), \
                 patch.object(runner.subprocess, "check_output", side_effect=git), \
                 patch.object(runner.platform, "platform", return_value="synthetic"), \
                 patch.object(runner, "report"), patch.object(runner, "supervise", return_value=0) as launch, \
                 patch.dict(runner.os.environ, {"ORX_RUN_ID": "synthetic", "ORX_COMMIT": "a" * 40,
                                                "PYTHONOPTIMIZE": "0"}):
                self.assertEqual(runner.main(), 0)
            record = launch.call_args.args[1]
            copied = root / "run-output/synthetic/sources/verification/reducer/reduce_stationary.py"
            self.assertIn(str(copied), record["command"])
            self.assertNotIn(str(here / "reduce_stationary.py"), record["command"])
            self.assertEqual(copied.read_bytes(), committed["paper1/verification/reducer/reduce_stationary.py"])


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.patches = ExitStack()
        self.addCleanup(self.patches.close)
        self.out = Path(self.patches.enter_context(TemporaryDirectory()))
        self.sources = []
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
        self.record = {"configuration": {"wall_seconds": 10, "memory_bytes": 2048},
                       "scientific_invocations": 0, "command": ["dummy"],
                       "source_sha256_before": {}}

    def replace(self, obj, name, **kwargs):
        return self.patches.enter_context(patch.object(obj, name, **kwargs))

    def supervise(self):
        rc = runner.supervise(self.out, self.record, {}, self.sources,
                              datetime.now(timezone.utc) + timedelta(minutes=10))
        execution = json.loads((self.out / "execution.json").read_bytes())
        samples = json.loads((self.out / "resource-samples.json").read_bytes())
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
        self.record["configuration"]["wall_seconds"] = 1
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
        self.record["configuration"]["wall_seconds"] = 1
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
        self.record["configuration"]["wall_seconds"] = 1
        self.child.poll.side_effect = [None, 0]
        self.replace(runner.time, "sleep")
        with patch.object(runner.time, "monotonic", side_effect=[0, 0.75, 0.8, 0.9, 1]):
            evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(self.ps.call_args.kwargs["timeout"], 0.25)

    def test_signal_after_terminal_boundary_is_outside_the_decision(self):
        original_save = runner.save_json
        def inject(path, payload):
            if path.name == "execution.json":
                self.handlers[runner.signal.SIGTERM](runner.signal.SIGTERM, None)
            return original_save(path, payload)
        self.replace(runner, "save_json", side_effect=inject)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["received_signals"], [])
        self.assertEqual(evidence["cancellation_boundary"], "before terminal evidence publication")

    def test_zero_exit_observed_after_deadline_fails(self):
        with patch.object(runner.time, "monotonic", side_effect=[0, 11, 11]):
            evidence, _ = self.supervise()
        self.assertEqual(evidence["exit_code"], 0)
        self.assertEqual(evidence["stop_reason"], "wall budget reached")

    def test_mutated_executed_snapshot_fails_even_when_original_is_unchanged(self):
        root = self.out / "original"
        source = root / "verification" / "synthetic.py"
        source.parent.mkdir(parents=True)
        source.write_bytes(b"reviewed source")
        self.sources = [source]
        relative = str(source.relative_to(root))
        self.record["source_sha256_before"] = {relative: hashlib.sha256(source.read_bytes()).hexdigest()}
        copied = self.out / "sources" / relative
        copied.parent.mkdir(parents=True)
        copied.write_bytes(b"mutated executed source")
        self.replace(runner, "ROOT", new=root)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["stop_reason"], "executed source snapshot changed during invocation")
        self.assertEqual(evidence["source_sha256_before"], evidence["source_sha256_after"])
        self.assertNotEqual(evidence["source_sha256_before"], evidence["source_snapshot_sha256_after"])

    def test_shared_wrapper_input_paths_are_not_treated_as_source_copies(self):
        root = self.out / "original"
        source = root / "verification" / "synthetic.py"
        source.parent.mkdir(parents=True)
        source.write_bytes(b"reviewed source")
        input_path = root / "input.json"
        input_path.write_bytes(b"reviewed input")
        self.sources = [source, input_path]
        before = {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in self.sources}
        self.record["source_sha256_before"] = before
        relative = str(source.relative_to(root))
        self.record["committed_source_sha256"] = {relative: before[relative]}
        copied = self.out / "sources" / relative
        copied.parent.mkdir(parents=True)
        copied.write_bytes(source.read_bytes())
        self.replace(runner, "ROOT", new=root)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["source_sha256_after"], before)
        self.assertEqual(evidence["source_snapshot_sha256_after"], self.record["committed_source_sha256"])


if __name__ == "__main__":
    unittest.main()
