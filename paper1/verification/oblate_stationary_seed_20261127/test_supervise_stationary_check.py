"""Mocked stdlib supervisor regressions; no child or scientific module runs."""

from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch
import hashlib
import json
import signal
import unittest

import supervise_stationary_check as supervisor


class SupervisorTests(unittest.TestCase):
    def exercise(self, *, signal_on_spawn=False, mutate_source=False,
                 delete_source=False, cleanup_failure=False, supervisor_hash_failure=False):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "verification" / "synthetic_checker.py"
            source.parent.mkdir()
            original_source = b"# synthetic source; never executed\n"
            source.write_bytes(original_source)
            relative_source = str(source.relative_to(root))
            preparation = root / "preparation"
            preparation.mkdir()
            (preparation / "preparation-provenance.json").write_text(json.dumps({
                "sha256": {relative_source: hashlib.sha256(original_source).hexdigest()}
            }))
            output = root / "attempt"
            process = Mock(pid=43210, returncode=0)
            process.poll.side_effect = [0] if signal_on_spawn else [None, 0]
            process.wait.return_value = 0
            initial_handlers = {signal.SIGTERM: object(), signal.SIGINT: object()}
            handlers = initial_handlers.copy()
            log_bytes = b"synthetic child log\x00\n"

            def register(sig, handler):
                previous = handlers[sig]
                handlers[sig] = handler
                return previous

            def fake_spawn(command, **kwargs):
                self.assertEqual(command, ["synthetic-never-execute"])
                self.assertIs(kwargs["start_new_session"], True)
                self.assertEqual(kwargs["cwd"], root)
                self.assertTrue(all(callable(handlers[sig]) for sig in initial_handlers))
                kwargs["stdout"].write(log_bytes)
                if mutate_source:
                    source.write_bytes(b"# persistent synthetic source mutation\n")
                if delete_source:
                    source.unlink()
                if signal_on_spawn:
                    handlers[signal.SIGTERM](signal.SIGTERM, None)
                return process

            def fake_command(command, **kwargs):
                if command == ["git", "rev-parse", "HEAD"]:
                    return "synthetic-pin\n"
                if command == ["git", "status", "--porcelain"]:
                    return ""
                if command == ["uv", "--version"]:
                    return "uv synthetic\n"
                if command == ["ps", "-axo", "pgid=,rss="]:
                    return "43210 64\n99999 1024\n"
                raise AssertionError(f"Unexpected external command: {command}")

            with ExitStack() as stack:
                for name, value in (("ROOT", root), ("PREP", preparation), ("OUT", output),
                                    ("COMMAND", ["synthetic-never-execute"])):
                    stack.enter_context(patch.object(supervisor, name, value))
                launch = stack.enter_context(patch.object(
                    supervisor.subprocess, "Popen", side_effect=fake_spawn))
                commands = stack.enter_context(patch.object(
                    supervisor.subprocess, "check_output", side_effect=fake_command))
                stack.enter_context(patch.object(supervisor.signal, "signal", side_effect=register))
                stack.enter_context(patch.object(supervisor.os, "killpg", side_effect=AssertionError(
                    "Tests must never signal a real process group")))
                stack.enter_context(patch.object(supervisor, "group_exists",
                                                return_value=signal_on_spawn or cleanup_failure))
                cleanup = stack.enter_context(patch.object(supervisor, "stop_group", return_value={
                    "group_visible_after_signals": False}))
                if cleanup_failure:
                    cleanup.side_effect = OSError("synthetic cleanup failure")
                if supervisor_hash_failure:
                    original_sha = supervisor.sha

                    def failing_hash(path):
                        if path == Path(supervisor.__file__):
                            raise FileNotFoundError("synthetic missing supervisor")
                        return original_sha(path)

                    stack.enter_context(patch.object(supervisor, "sha", side_effect=failing_hash))
                stack.enter_context(patch.object(supervisor, "write_console", return_value=False))
                stack.enter_context(patch.object(supervisor.shutil, "which", return_value="/synthetic/uv"))
                stack.enter_context(patch.object(supervisor.platform, "platform", return_value="synthetic"))
                stack.enter_context(patch.object(supervisor.time, "monotonic", return_value=0.0))
                stack.enter_context(patch.object(supervisor.time, "sleep"))
                stack.enter_context(patch.dict(supervisor.os.environ, {"PYTHONOPTIMIZE": "0"}))
                with self.assertRaises(SystemExit) as exited:
                    supervisor.main("synthetic-pin", "synthetic-test-only")
                launch.assert_called_once()
                process.wait.assert_called_once()
                self.assertEqual(handlers, initial_handlers)
                if signal_on_spawn or cleanup_failure:
                    cleanup.assert_called_once_with(process)
                else:
                    cleanup.assert_not_called()
                self.assertEqual((output / "verification.log").read_bytes(), log_bytes)
                self.assertEqual((output / "sources" / relative_source).read_bytes(), original_source)
                record = json.loads((output / "execution.json").read_text())
                self.assertEqual(record["status"], "passed" if exited.exception.code == 0 else "failed")
                samples = json.loads((output / "resource-samples.json").read_text())
                self.assertEqual(record["scientific_invocations"], 1)  # Mocked invocation only.
                self.assertEqual(record["source_commit"], "synthetic-pin")
                ps_calls = [call for call in commands.call_args_list if call.args[0][0] == "ps"]
            return exited.exception.code, record, samples, ps_calls

    def test_sigterm_during_spawn_fails_and_preserves_cleanup_evidence(self):
        code, record, samples, ps_calls = self.exercise(signal_on_spawn=True)
        self.assertEqual(code, 1)
        self.assertEqual(record["exit_code"], 0)  # Cancellation still overrides a clean child exit.
        self.assertEqual(record["stop_reason"], "supervisor received signal")
        self.assertEqual(record["received_signals"], [int(signal.SIGTERM)])
        self.assertEqual(record["final_group_cleanup"], {"group_visible_after_signals": False})
        self.assertEqual(samples, [])
        self.assertEqual(ps_calls, [])

    def test_persistent_source_mutation_prevents_successful_exit(self):
        code, record, samples, _ = self.exercise(mutate_source=True)
        self.assertEqual(code, 1)
        self.assertEqual(record["exit_code"], 0)
        self.assertEqual(record["stop_reason"], "source files changed during invocation")
        self.assertNotEqual(record["source_sha256_before"], record["source_sha256_after"])
        self.assertEqual(record["received_signals"], [])
        self.assertEqual(samples[0]["process_group_rss_bytes"], 64 * 1024)

    def test_clean_completion_succeeds_with_matching_source_and_resource_evidence(self):
        code, record, samples, ps_calls = self.exercise()
        self.assertEqual(code, 0)
        self.assertEqual(record["exit_code"], 0)
        self.assertIsNone(record["stop_reason"])
        self.assertEqual(record["received_signals"], [])
        self.assertEqual(record["source_sha256_before"], record["source_sha256_after"])
        self.assertEqual(record["sampled_peak_process_group_rss_bytes"], 64 * 1024)
        self.assertEqual(len(samples), 1)
        self.assertEqual(len(ps_calls), 1)
        self.assertLessEqual(ps_calls[0].kwargs["timeout"], 1)

    def test_deleted_source_preserves_failed_execution_and_restores_handlers(self):
        code, record, samples, _ = self.exercise(delete_source=True)
        self.assertEqual(code, 1)
        self.assertEqual(record["stop_reason"], "source files changed during invocation")
        self.assertEqual(record["source_sha256_after"], {})
        self.assertEqual(record["cleanup_errors"][0]["action"], "source_hash")
        self.assertIn("FileNotFoundError", record["cleanup_errors"][0]["error"])
        self.assertEqual(len(samples), 1)

    def test_cleanup_exception_preserves_failed_execution_and_restores_handlers(self):
        code, record, samples, _ = self.exercise(cleanup_failure=True)
        self.assertEqual(code, 1)
        self.assertEqual(record["cleanup_errors"][0]["action"], "stop_group")
        self.assertIn("synthetic cleanup failure", record["cleanup_errors"][0]["error"])
        self.assertEqual(record["source_sha256_before"], record["source_sha256_after"])
        self.assertEqual(len(samples), 1)

    def test_supervisor_hash_failure_preserves_failed_execution_and_restores_handlers(self):
        code, record, _, _ = self.exercise(supervisor_hash_failure=True)
        self.assertEqual(code, 1)
        self.assertEqual(record["stop_reason"], "supervisor cleanup failed")
        self.assertIsNone(record["supervisor_sha256"])
        self.assertEqual(record["cleanup_errors"][0]["action"], "supervisor_hash")


if __name__ == "__main__":
    unittest.main(verbosity=2)
