"""Stdlib-only dispatch tests; the scientific supervisor is always mocked."""
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import run_bounded as runner


class DispatchTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(TemporaryDirectory())) / "paper1"
        self.here = self.root / "verification/probe"
        self.shared = self.root / "verification/reduction/run_bounded.py"
        self.here.mkdir(parents=True)
        self.shared.parent.mkdir(parents=True)
        for name in ("run_bounded.py", "probe_terminal_embedding.py", "run.sh", "view-config.json", "view.xml"):
            (self.here / name).write_text("fixture")
        self.shared.write_text("fixture")
        (self.shared.parent / "run.sh").write_text("fixture")
        self.input = self.root / "checkpoint.json"
        self.input.write_bytes(b'{"fixture":true}\n')
        self.digest = hashlib.sha256(self.input.read_bytes()).hexdigest()
        self.config = {"wall_seconds": 300, "memory_bytes": 2147483648,
            "approval_cutoff_utc": "2099-01-01T00:00:00+00:00", "checkpoint_path": "checkpoint.json",
            "checkpoint_sha256": self.digest, "stage": "diagnostic"}
        self.write_config()
        for name, value in (("HERE", self.here), ("ROOT", self.root), ("SHARED", self.shared)):
            self.stack.enter_context(patch.object(runner, name, value))
        self.stack.enter_context(patch.dict(runner.os.environ,
            {"ORX_RUN_ID": "fixture-run", "ORX_COMMIT": "a" * 40, "PYTHONOPTIMIZE": "0"}))
        self.launch = self.stack.enter_context(patch.object(runner.common, "supervise", return_value=0))
        self.stack.enter_context(patch.object(runner.common, "report"))
        self.git = self.stack.enter_context(patch.object(runner.subprocess, "check_output", side_effect=self.git_bytes))

    def write_config(self):
        (self.here / "run-config.json").write_text(json.dumps(self.config))

    def git_bytes(self, command, **kwargs):
        if "rev-parse" in command:
            return "a" * 40 + "\n"
        return (self.root / command[-1].split(":paper1/", 1)[1]).read_bytes()

    def test_dispatches_verified_copy_and_tracks_both_input_paths(self):
        self.assertEqual(runner.main(), 0)
        out, record, env, paths, cutoff = self.launch.call_args.args
        saved = out / "inputs/dae-step-02.json"
        self.assertEqual(saved.read_bytes(), self.input.read_bytes())
        self.assertIn(str(saved), record["command"])
        self.assertNotIn(str(self.input), record["command"])
        self.assertEqual(record["scientific_invocations"], 0)
        self.assertEqual(len(record["input_sha256_before"]), 2)
        self.assertIn(saved, paths)
        self.assertIn(self.input, paths)
        self.assertEqual(env["OMP_NUM_THREADS"], "1")

    def test_bad_checkpoint_pin_never_dispatches(self):
        self.config["checkpoint_sha256"] = "0" * 64
        self.write_config()
        with self.assertRaisesRegex(ValueError, "digest mismatch"):
            runner.main()
        self.launch.assert_not_called()
        evidence = json.loads((self.root / "run-output/fixture-run/execution.json").read_bytes())
        self.assertEqual(evidence["scientific_invocations"], 0)
        self.assertEqual(evidence["stop_reason"], "preflight failed")

    def test_wrong_committed_source_never_dispatches(self):
        original = self.git_bytes
        self.git.side_effect = lambda command, **kwargs: (b"different" if "show" in command
                                                        else original(command, **kwargs))
        with self.assertRaisesRegex(ValueError, "digest mismatch"):
            runner.main()
        self.launch.assert_not_called()

    def test_expired_window_never_dispatches(self):
        self.config["approval_cutoff_utc"] = "2000-01-01T00:00:00+00:00"
        self.write_config()
        with self.assertRaisesRegex(ValueError, "approved window"):
            runner.main()
        self.launch.assert_not_called()


if __name__ == "__main__":
    unittest.main()
