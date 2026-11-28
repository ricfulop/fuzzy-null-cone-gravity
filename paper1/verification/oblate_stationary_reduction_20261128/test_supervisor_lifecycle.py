"""Synthetic POSIX process/pipe checks; no scientific child or module is used."""

from datetime import datetime, timedelta, timezone
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch


VERIFY = Path(__file__).resolve().parent.parent
RUNNERS = ["oblate_stationary_reduction_20261128", "oblate_stationary_germ_20261128",
           "oblate_stationary_state_germ_20261128"]


def load_runner(folder):
    spec = importlib.util.spec_from_file_location(folder + "_synthetic_runner", VERIFY / folder / "run_bounded.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@unittest.skipUnless(os.name == "posix", "process groups and nonblocking pipes require POSIX")
class SyntheticLifecycleTests(unittest.TestCase):
    def test_full_console_pipe_never_blocks_forwarding(self):
        script = '''
import importlib.util, os, sys
spec = importlib.util.spec_from_file_location("runner", sys.argv[1])
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
read_fd, write_fd = os.pipe()
try:
    os.set_blocking(write_fd, False)
    try:
        while True:
            os.write(write_fd, b"x" * 65536)
    except BlockingIOError:
        pass
    os.set_blocking(write_fd, True)
    original = sys.stdout
    with os.fdopen(write_fd, "wb", buffering=0, closefd=False) as stream:
        sys.stdout = stream
        try:
            assert runner.write_console(b"not allowed to block\\n") is False
            assert os.get_blocking(write_fd) is True
        finally:
            sys.stdout = original
finally:
    os.close(write_fd)
    os.close(read_fd)
'''
        for folder in RUNNERS:
            with self.subTest(runner=folder):
                result = subprocess.run([sys.executable, "-B", "-c", script,
                                         str(VERIFY / folder / "run_bounded.py")],
                                        capture_output=True, timeout=3)
                self.assertEqual(result.returncode, 0, result.stderr.decode())

    def test_timeout_observes_actual_descendant_exit(self):
        # The launcher reaps its child on TERM, so an exited descendant cannot
        # remain a zombie and confuse the independent kill-zero observation.
        launcher = '''
import os, signal, subprocess, sys
code = "import os,signal,sys,time; from pathlib import Path; signal.signal(signal.SIGTERM,signal.SIG_IGN); Path(sys.argv[1]).write_text(str(os.getpid())); time.sleep(60)"
child = subprocess.Popen([sys.executable, "-B", "-c", code, sys.argv[1]])
def stop(*_):
    child.kill()
    child.wait()
    raise SystemExit(0)
signal.signal(signal.SIGTERM, stop)
child.wait()
'''
        for folder in RUNNERS:
            with self.subTest(runner=folder), TemporaryDirectory() as temporary:
                runner = load_runner(folder)
                output = Path(temporary)
                pid_file = output / "descendant.pid"
                record = {"wall_seconds": 1, "memory_bytes": 2 * 1024**3,
                          "configuration": {"wall_seconds": 1, "memory_bytes": 2 * 1024**3},
                          "command": [sys.executable, "-B", "-c", launcher, str(pid_file)],
                          "scientific_invocations": 0, "source_sha256_before": {},
                          "input_sha256_before": {}, "checker_input_sha256": {}}
                processes = []
                original_popen = subprocess.Popen
                def launch(*args, **kwargs):
                    process = original_popen(*args, **kwargs)
                    if kwargs.get("start_new_session"):
                        processes.append(process)
                    return process
                try:
                    cutoff = datetime.now(timezone.utc) + timedelta(minutes=2)
                    with patch.object(runner.subprocess, "Popen", side_effect=launch), \
                         patch.object(runner, "write_console", return_value=True):
                        args = [output, record, os.environ.copy()]
                        if folder == RUNNERS[0]:
                            args.append([])
                        rc = runner.supervise(*args, cutoff)
                    evidence = json.loads((output / "execution.json").read_bytes())
                    self.assertEqual(rc, 1)
                    self.assertEqual(evidence["stop_reason"], "wall budget reached")
                    self.assertEqual(evidence["group_disposition"]["status"], "empty")
                    descendant = int(pid_file.read_text())
                    with self.assertRaises(ProcessLookupError):
                        os.kill(descendant, 0)
                finally:
                    for process in processes:
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        process.wait(timeout=3)
                    if pid_file.exists():
                        try:
                            os.kill(int(pid_file.read_text()), signal.SIGKILL)
                        except ProcessLookupError:
                            pass


if __name__ == "__main__":
    unittest.main()
