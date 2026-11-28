"""Stdlib-only JSON persistence for stationary reduction attempts."""

from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import sys
import tempfile
import time
import uuid


def report(**fields):
    print("::orx:: " + json.dumps(fields, allow_nan=False), flush=True)


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    return sha256_bytes(Path(path).read_bytes())


def atomic_json_new(path, value):
    """Publish one complete JSON file without replacing existing evidence."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=".partial-",
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


class ReductionRecorder:
    """Fresh-directory progress recorder for exact reduction probes."""

    def __init__(self, directory, identity):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=False)
        self.identity = json.loads(canonical_bytes(identity))
        self.identity_sha256 = sha256_bytes(canonical_bytes(self.identity))
        self.started = time.monotonic()
        self.sequence = 0
        self.stage_count = 0
        self.attempt_id = uuid.uuid4().hex
        self.stream = (self.directory / "progress.jsonl").open("x", encoding="utf-8")
        self.write_json(
            "invocation.json",
            {
                "state": "started; no reduction completion asserted",
                "identity": self.identity,
                "identity_sha256": self.identity_sha256,
                "attempt_id": self.attempt_id,
                "start_utc": datetime.now(timezone.utc).isoformat(),
                "checkpoint_policy": "write only; no automatic reuse or resume",
            },
        )

    def write_json(self, name, value):
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Output must stay within the fresh attempt directory")
        atomic_json_new(self.directory / relative, value)

    def event(self, event, **details):
        self.sequence += 1
        record = {
            "sequence": self.sequence,
            "event": event,
            "seconds": round(time.monotonic() - self.started, 6),
            "utc": datetime.now(timezone.utc).isoformat(),
            **details,
        }
        line = json.dumps(record, sort_keys=True, allow_nan=False)
        self.stream.write(line + "\n")
        self.stream.flush()
        print("PROGRESS " + line, file=sys.stderr, flush=True)

    @contextmanager
    def stage(self, name, **details):
        if name in {"manufactured_controls", "fixed_epsilon_field_conversion",
                    "radial_x0_gauge_right_kernel", "principal_block_probe", "dae_closure"}:
            report(phase=name, step=0, of=1, unit="check")
        self.stage_count += 1
        stage_id = self.stage_count
        started = time.monotonic()
        self.event("begin", stage=name, stage_id=stage_id, **details)
        try:
            yield
        except BaseException as error:
            self.event(
                "error",
                stage=name,
                stage_id=stage_id,
                duration_seconds=round(time.monotonic() - started, 6),
                error_type=type(error).__name__,
                error=str(error),
                **details,
            )
            raise
        else:
            self.event(
                "end",
                stage=name,
                stage_id=stage_id,
                duration_seconds=round(time.monotonic() - started, 6),
                **details,
            )
            if name in {"manufactured_controls", "fixed_epsilon_field_conversion",
                        "radial_x0_gauge_right_kernel", "principal_block_probe", "dae_closure"}:
                report(phase=name, step=1, of=1, unit="check")

    def __enter__(self):
        self.event("invocation_begin")
        return self

    def __exit__(self, error_type, error, traceback):
        try:
            if error_type is not None:
                self.event(
                    "invocation_error",
                    error_type=error_type.__name__,
                    error=str(error),
                )
                self.write_json(
                    "failure.json",
                    {
                        "state": "exception; no successful reduction",
                        "error_type": error_type.__name__,
                        "error": str(error),
                        "identity_sha256": self.identity_sha256,
                    },
                )
            else:
                self.event("invocation_return")
        finally:
            self.stream.close()
        return False
