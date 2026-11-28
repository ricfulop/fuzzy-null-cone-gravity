"""Bounded diagnostic dispatch, reusing the reviewed reduction supervisor."""
from datetime import datetime, timezone, timedelta
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SHARED = ROOT / "verification/oblate_stationary_reduction_20261128/run_bounded.py"
spec = importlib.util.spec_from_file_location("stationary_supervisor", SHARED)
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)


def require(value, message):
    if not value:
        raise ValueError(message)


def pinned_copy(source, destination, expected):
    data = source.read_bytes()
    require(hashlib.sha256(data).hexdigest() == expected, "Checkpoint digest mismatch")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as stream:
        stream.write(data)
    require(common.sha(destination) == expected, "Checkpoint snapshot mismatch")


def main():
    run_id, commit = os.environ.get("ORX_RUN_ID"), os.environ.get("ORX_COMMIT")
    require(run_id and Path(run_id).name == run_id and run_id not in (".", ".."), "Invalid run ID")
    require(commit and re.fullmatch("[0-9a-f]{40}", commit), "Missing immutable source commit")
    require(not sys.flags.optimize and os.environ.get("PYTHONOPTIMIZE") in (None, "", "0"),
            "Scientific assertions must remain enabled")
    config = json.loads((HERE / "run-config.json").read_bytes())
    cutoff = datetime.fromisoformat(config["approval_cutoff_utc"])
    require(config["wall_seconds"] == 300 and config["memory_bytes"] == 2147483648,
            "Unexpected diagnostic resource envelope")
    require(datetime.now(timezone.utc) + timedelta(seconds=310) < cutoff,
            "Diagnostic would exceed approved window")
    out = ROOT / "run-output" / run_id
    out.mkdir(parents=True, exist_ok=False)
    sources = sorted(HERE.glob("*.py")) + [HERE / n for n in
        ("run-config.json", "view-config.json", "view.xml", "run.sh")]
    sources += [SHARED, SHARED.parent / "run.sh"]
    record = {"run_id": run_id, "source_commit": commit, "configuration": config,
              "scientific_invocations": 0, "remote_jobs": 0, "GPU_jobs": 0,
              "diagnostic_only": True, "dae_closure_certified": False,
              "supervisor_python": sys.version, "start_utc": datetime.now(timezone.utc).isoformat(),
              "source_repository": str(common.SOURCE_REPOSITORY),
              "rss_measurement": "sampled process-group RSS, not a hard OS cap"}
    try:
        resolved = subprocess.check_output(
            ["git", "-C", str(common.SOURCE_REPOSITORY), "rev-parse", commit + "^{commit}"], text=True).strip()
        require(resolved == commit, "Source revision mismatch")
        record["source_sha256_before"] = {}
        for p in sources:
            relative = p.relative_to(ROOT).as_posix()
            committed = subprocess.check_output(
                ["git", "-C", str(common.SOURCE_REPOSITORY), "show", commit + ":paper1/" + relative])
            digest = hashlib.sha256(committed).hexdigest()
            target = out / "sources" / relative
            pinned_copy(p, target, digest)
            record["source_sha256_before"][relative] = digest
        record["committed_source_sha256"] = dict(record["source_sha256_before"])
        require(config == json.loads((out / "sources" / HERE.relative_to(ROOT) / "run-config.json").read_bytes()),
                "Dispatched configuration differs from committed snapshot")
        checkpoint = ROOT / config["checkpoint_path"]
        digest = config["checkpoint_sha256"]
        require(isinstance(digest, str) and re.fullmatch("[0-9a-f]{64}", digest), "Invalid checkpoint pin")
        saved = out / "inputs/dae-step-02.json"
        pinned_copy(checkpoint, saved, digest)
        # The shared supervisor verifies both original and dispatched input bytes after cleanup.
        sources += [checkpoint, saved]
        record["input_sha256_before"] = {str(checkpoint.relative_to(ROOT)): digest,
                                         str(saved.relative_to(ROOT)): digest}
        record["source_sha256_before"].update(record["input_sha256_before"])
        record["source_hash_scope"] = "committed sources plus the two separately labeled input paths"
        record["command"] = ["uv", "run", "--isolated", "--no-project", "--with", "sympy==1.14.0",
            "python", str(HERE / "probe_terminal_embedding.py"), "--checkpoint", str(saved),
            "--checkpoint-sha256", digest, "--output", str(out / "outputs/embedding-probe.json")]
        overrides = {"PYTHONUNBUFFERED": "1", "PYTHONOPTIMIZE": "0", "PYTHONDONTWRITEBYTECODE": "1",
                     "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                     "VECLIB_MAXIMUM_THREADS": "1"}
        record["environment_overrides"] = overrides
        common.save_json(out / "launch.json", record)
    except Exception as error:
        record.update(status="failed", stop_reason="preflight failed", supervisor_error=repr(error),
                      exit_code=None, group_disposition={"status": "not_spawned"})
        common.save_json(out / "execution.json", record)
        common.save_json(out / "resource-samples.json", [])
        raise
    common.report(phase="setup", step=1, of=1, unit="environment")
    return common.supervise(out, record, {**os.environ, **overrides}, sources, cutoff)


if __name__ == "__main__":
    sys.exit(main())
