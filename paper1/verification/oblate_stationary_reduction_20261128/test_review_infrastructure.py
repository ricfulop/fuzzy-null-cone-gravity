"""Stdlib-only buffer and evidence regressions; no scientific modules imported."""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import authenticated_export
from reduction_progress import ReductionRecorder


class AuthenticatedBufferTests(unittest.TestCase):
    def fixture(self):
        operators = b'{"coefficient": "reviewed"}'
        reproducibility = b'{"source": "reviewed"}'
        payloads = {"stationary-operators.json": operators,
                    "stationary-reproducibility.json": reproducibility}
        marker = json.dumps({"export_sha256": {
            name: hashlib.sha256(data).hexdigest() for name, data in payloads.items()
        }}).encode()
        return {"completion.json": marker, **payloads}, hashlib.sha256(marker).hexdigest()

    def test_parse_and_identity_use_each_authenticated_buffer_once(self):
        payloads, pin = self.fixture()
        reads = []
        def read(path):
            reads.append(path.name)
            if reads.count(path.name) > 1:
                raise AssertionError("authenticated input was reread")
            return payloads[path.name]
        with patch.object(Path, "read_bytes", read):
            _, decoded, hashes = authenticated_export.read_completed_buffers("unused", pin)
        self.assertEqual(decoded["stationary-operators.json"], {"coefficient": "reviewed"})
        self.assertEqual(reads, list(payloads))
        self.assertEqual(hashes, {name: hashlib.sha256(data).hexdigest()
                                  for name, data in payloads.items()})

    def test_bad_marker_is_rejected_before_json_or_payload_reads(self):
        with patch.object(Path, "read_bytes", return_value=b"not JSON") as read, \
             patch.object(authenticated_export.json, "loads") as decode:
            with self.assertRaisesRegex(ValueError, "completion marker"):
                authenticated_export.read_completed_buffers("unused", "0" * 64)
        self.assertEqual(read.call_count, 1)
        decode.assert_not_called()

    def test_unreviewed_payload_never_reaches_json_parser(self):
        payloads, pin = self.fixture()
        payloads["stationary-operators.json"] = b"unreviewed malformed JSON"
        original = json.loads
        with patch.object(Path, "read_bytes", lambda p: payloads[p.name]), \
             patch.object(authenticated_export.json, "loads", wraps=original) as decode:
            with self.assertRaisesRegex(ValueError, "stationary-operators.json hash"):
                authenticated_export.read_completed_buffers("unused", pin)
        self.assertEqual(decode.call_count, 1)
        self.assertEqual(decode.call_args.args, (payloads["completion.json"],))

    def test_reproducibility_digest_is_required_before_payload_decoding(self):
        payloads, pin = self.fixture()
        payloads["stationary-reproducibility.json"] = b"{}"
        with patch.object(Path, "read_bytes", lambda p: payloads[p.name]):
            with self.assertRaisesRegex(ValueError, "stationary-reproducibility.json hash"):
                authenticated_export.read_completed_buffers("unused", pin)


class RecorderEvidenceTests(unittest.TestCase):
    def test_closed_stderr_preserves_original_failure_record(self):
        with TemporaryDirectory() as temporary:
            output = Path(temporary) / "attempt"
            with patch("reduction_progress.sys.stderr") as console:
                console.write.side_effect = BrokenPipeError("reader disconnected")
                with self.assertRaisesRegex(ValueError, "original failure"):
                    with ReductionRecorder(output, {}) as recorder:
                        with recorder.stage("synthetic"):
                            raise ValueError("original failure")
            failure = json.loads((output / "failure.json").read_bytes())
            self.assertEqual(failure["error_type"], "ValueError")
            self.assertEqual(failure["error"], "original failure")
            events = [json.loads(line)["event"] for line in (output / "progress.jsonl").read_text().splitlines()]
            self.assertEqual(events, ["invocation_begin", "begin", "error", "invocation_error"])

    def test_journal_failure_is_not_suppressed_as_console_failure(self):
        with TemporaryDirectory() as temporary:
            recorder = ReductionRecorder(Path(temporary) / "attempt", {})
            try:
                with patch.object(recorder.stream, "write", side_effect=OSError("journal failed")):
                    with self.assertRaisesRegex(OSError, "journal failed"):
                        recorder.event("synthetic")
            finally:
                recorder.stream.close()


if __name__ == "__main__":
    unittest.main()
