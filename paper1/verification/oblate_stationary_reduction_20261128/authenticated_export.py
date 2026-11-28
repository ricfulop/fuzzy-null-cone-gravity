"""Authenticate completed JSON export buffers before decoding any coefficients."""

import hashlib
import json
from pathlib import Path


def read_completed_buffers(directory, completion_sha256):
    directory = Path(directory)
    marker_bytes = (directory / "completion.json").read_bytes()
    digests = {"completion.json": hashlib.sha256(marker_bytes).hexdigest()}
    if digests["completion.json"] != completion_sha256:
        raise ValueError("completion marker differs from the independently reviewed upstream export")
    marker = json.loads(marker_bytes)
    payloads = {}
    for name in ("stationary-operators.json", "stationary-reproducibility.json"):
        data = (directory / name).read_bytes()
        digests[name] = hashlib.sha256(data).hexdigest()
        if marker.get("export_sha256", {}).get(name) != digests[name]:
            raise ValueError(name + " hash does not match completion marker")
        payloads[name] = data
    # Only authenticated buffers reach the JSON parser. No file is reread.
    return marker, {name: json.loads(data) for name, data in payloads.items()}, digests
