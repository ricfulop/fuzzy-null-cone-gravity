"""Stdlib-only timing and write-only checkpoints, not a scientific verifier.

Every invocation requires a new output directory. Column checkpoints preserve
exact expression strings but are explicitly pending the global controls. There
is no automatic resume, expression parser, or reuse of earlier results here.
"""
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


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def atomic_json_new(path, value):
    """Publish a complete JSON file without replacing any existing evidence."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                         dir=path.parent, prefix=".partial-",
                                         delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        # Same-filesystem link publishes atomically and fails if path exists.
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


class ExtractionRecorder:
    def __init__(self, directory, identity):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=False)
        self.identity = json.loads(canonical_bytes(identity))
        self.identity_sha256 = hashlib.sha256(canonical_bytes(self.identity)).hexdigest()
        self.started = time.monotonic()
        self.sequence = 0
        self.stage_count = 0
        self.attempt_id = uuid.uuid4().hex
        self.checkpoints = []
        self.stream = (self.directory / "progress.jsonl").open("x", encoding="utf-8")
        self.write_json("invocation.json", {
            "state": "started; no scientific completion asserted",
            "identity": self.identity, "identity_sha256": self.identity_sha256,
            "attempt_id": self.attempt_id,
            "start_utc": datetime.now(timezone.utc).isoformat(),
            "checkpoint_policy": "write only; no automatic reuse or resume",
        })

    def write_json(self, name, value):
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Output must stay within the fresh invocation directory")
        atomic_json_new(self.directory / relative, value)

    def event(self, event, **details):
        self.sequence += 1
        record = {"sequence": self.sequence, "event": event,
                  "seconds": round(time.monotonic()-self.started, 6),
                  "utc": datetime.now(timezone.utc).isoformat(), **details}
        line = json.dumps(record, sort_keys=True, allow_nan=False)
        self.stream.write(line+"\n")
        self.stream.flush()
        print("PROGRESS "+line, file=sys.stderr, flush=True)

    @contextmanager
    def stage(self, name, **details):
        self.stage_count += 1
        stage_id = self.stage_count
        started = time.monotonic()
        self.event("begin", stage=name, stage_id=stage_id, **details)
        try:
            yield
        except BaseException as error:
            self.event("error", stage=name, stage_id=stage_id,
                       duration_seconds=round(time.monotonic()-started, 6),
                       error_type=type(error).__name__, error=str(error), **details)
            raise
        else:
            self.event("end", stage=name, stage_id=stage_id,
                       duration_seconds=round(time.monotonic()-started, 6), **details)

    def save_column(self, column, label, payload):
        body = {"schema_version": 1,
                "kind": "stationary_column_checkpoint",
                "state": "column assembled; global controls pending; not an operator export",
                "identity_sha256": self.identity_sha256,
                "attempt_id": self.attempt_id,
                "column_index": column, "column_label": label,
                "radial_DAE_reduced": False, "admissible_seed_certified": False,
                "nonlinear_consistency_established": False, **payload}
        envelope = {"payload": body,
                    "payload_sha256": hashlib.sha256(canonical_bytes(body)).hexdigest()}
        name = f"checkpoints/column-{column:02d}-{label[0]}-{label[1]}.json"
        self.write_json(name, envelope)
        file_hash = hashlib.sha256((self.directory / name).read_bytes()).hexdigest()
        self.checkpoints.append({"file": name, "column": column, "sha256": file_hash})
        self.event("checkpoint", column=column, label=label, file=name,
                   payload_sha256=envelope["payload_sha256"], file_sha256=file_hash)

    def __enter__(self):
        self.event("invocation_begin")
        return self

    def __exit__(self, error_type, error, traceback):
        try:
            if error_type is not None:
                self.event("invocation_error", error_type=error_type.__name__, error=str(error))
                self.write_json("failure.json", {
                    "state": "exception; no successful completion",
                    "error_type": error_type.__name__, "error": str(error),
                    "identity_sha256": self.identity_sha256,
                })
            else:
                self.event("invocation_return")
        finally:
            self.stream.close()
        return False
