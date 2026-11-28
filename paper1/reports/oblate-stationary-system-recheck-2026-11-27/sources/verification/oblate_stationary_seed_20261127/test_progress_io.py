"""Filesystem-only regression checks; no scientific modules or expressions run."""
from contextlib import redirect_stderr
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, call, patch
import hashlib
import io
import json
import unittest

from extraction_progress import ExtractionRecorder, atomic_json_new, canonical_bytes
import supervise_stationary_check as supervisor


class ProgressIOTests(unittest.TestCase):
    def test_existing_evidence_is_not_overwritten(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence.json"
            atomic_json_new(path, {"original": True})
            before = path.read_bytes()
            with self.assertRaises(FileExistsError):
                atomic_json_new(path, {"original": False})
            self.assertEqual(path.read_bytes(), before)
            self.assertEqual(list(path.parent.glob(".partial-*")), [])
            with self.assertRaises(FileExistsError):
                ExtractionRecorder(path.parent, {})

    def test_interrupted_publication_does_not_create_a_checkpoint(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "column.json"
            with patch("extraction_progress.os.link", side_effect=OSError("simulated interruption")):
                with self.assertRaises(OSError):
                    atomic_json_new(path, {"exact_expression": "1/(r**2+epsilon)"})
            self.assertFalse(path.exists())
            self.assertEqual(list(path.parent.glob(".partial-*")), [])

    def test_failure_preserves_partial_evidence_without_completion(self):
        with TemporaryDirectory() as temporary, redirect_stderr(io.StringIO()):
            directory = Path(temporary) / "attempt"
            with self.assertRaisesRegex(RuntimeError, "manufactured failure"):
                with ExtractionRecorder(directory, {"source": "toy"}) as recorder:
                    with recorder.stage("toy_column", column=0):
                        recorder.save_column(0, ("toy", 0), {
                            "operators_by_derivative_order": [["1/3"], ["0"], ["sqrt(epsilon)"]],
                            "densities_by_derivative_order": [["x**2"], ["0"], ["sqrt(epsilon)"]],
                            "all_global_controls_passed": False,
                        })
                    with recorder.stage("toy_control"):
                        raise RuntimeError("manufactured failure")
            events = [json.loads(line) for line in (directory/"progress.jsonl").read_text().splitlines()]
            self.assertTrue(any(event["event"] == "error" and event["stage"] == "toy_control"
                                for event in events))
            self.assertTrue((directory/"failure.json").exists())
            self.assertFalse((directory/"completion.json").exists())
            self.assertFalse((directory/"stationary-operators.json").exists())
            checkpoint = directory/"checkpoints/column-00-toy-0.json"
            envelope = json.loads(checkpoint.read_text())
            self.assertEqual(envelope["payload_sha256"],
                             hashlib.sha256(canonical_bytes(envelope["payload"])).hexdigest())
            self.assertEqual(envelope["payload"]["operators_by_derivative_order"][2], ["sqrt(epsilon)"])
            self.assertFalse(envelope["payload"]["all_global_controls_passed"])
            saved = next(event for event in events if event["event"] == "checkpoint")
            self.assertEqual(saved["file_sha256"], hashlib.sha256(checkpoint.read_bytes()).hexdigest())

    def test_invalid_json_or_escape_cannot_create_result(self):
        with TemporaryDirectory() as temporary, redirect_stderr(io.StringIO()):
            directory = Path(temporary) / "attempt"
            with ExtractionRecorder(directory, {"source": "toy"}) as recorder:
                with self.assertRaises(ValueError):
                    recorder.write_json("bad.json", {"value": float("nan")})
                with self.assertRaises(ValueError):
                    recorder.write_json("../outside.json", {})
            self.assertFalse((directory/"bad.json").exists())
            self.assertFalse((directory.parent/"outside.json").exists())
            self.assertFalse((directory/"completion.json").exists())

    def test_optimized_python_preflight_cannot_launch_child(self):
        with TemporaryDirectory() as temporary:
            with patch.object(supervisor, "OUT", Path(temporary)/"new"), \
                 patch.object(supervisor.subprocess, "check_output", side_effect=["pin\n", ""]), \
                 patch.object(supervisor.subprocess, "Popen") as launch, \
                 patch.dict(supervisor.os.environ, {"PYTHONOPTIMIZE": "1"}):
                with self.assertRaisesRegex(RuntimeError, "disable scientific assertions"):
                    supervisor.main("pin", "test-only")
                launch.assert_not_called()

    def test_descendant_cleanup_does_not_trust_launcher_exit(self):
        process = Mock(pid=12345)
        process.poll.return_value = 0
        with patch.object(supervisor, "group_exists", side_effect=[True, True, True, False]), \
             patch.object(supervisor, "signal_group") as signals, \
             patch.object(supervisor.time, "monotonic", side_effect=[0, 0.01, 3]), \
             patch.object(supervisor.time, "sleep"):
            result = supervisor.stop_group(process)
        self.assertEqual(signals.call_args_list,
                         [call(12345, supervisor.signal.SIGTERM), call(12345, supervisor.signal.SIGKILL)])
        process.wait.assert_called_once()
        self.assertFalse(result["group_visible_after_signals"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
