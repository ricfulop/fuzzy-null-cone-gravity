"""Full N=3 coefficient-map atlas, prepared for a separately bounded invocation."""
import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import platform
import tempfile
import time

import sympy as s

from atlas_core import Label, affine_check, column_map, labels
from independent_controls import verify_independent_controls


AXES = ("x", "y", "z")
HERE = Path(__file__).resolve().parent
ARCHIVED_REGRESSION = {
    "path": "paper1/reports/oblate-null-map-check-2026-11-27/map-checks.json",
    "sha256": "75fc5c1036777716ef80bd1911adc9a90d2ef0582110f0ffe3ccc9ed973305cc",
}
SCOPE = dict(admissible_seed_supplied=False, euler_residual_evaluated=False,
             nonlinear_compatibility_evaluated=False, dirac_classification_established=False)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def publish(path, payload):
    data = (json.dumps(payload, indent=2, allow_nan=False) + "\n").encode()
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".atlas-", delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
            os.link(temporary, path)
        finally:
            temporary.unlink()
    return digest(data)


def progress(**fields):
    print("::orx:: " + json.dumps(fields, allow_nan=False), flush=True)


def detail(event, **fields):
    print("PROGRESS " + json.dumps(dict(event=event, **fields), allow_nan=False), flush=True)


@contextmanager
def phase(name, unit):
    start = time.monotonic()
    progress(phase=name, step=0, of=1, unit=unit)
    yield
    progress(phase=name, step=1, of=1, unit=unit, m={"phase_seconds": time.monotonic()-start})


def anchor_controls(lookup):
    # These independent scalar highest-mode checks fix the two transverse phases.
    for field in ("a", "b", "lambda"):
        column = Label(field, 2, 2)
        row = Label(field, 3, 3)
        require(lookup(column, "x").get(row) == -s.I/5 and lookup(column, "y").get(row) == -s.S.One/5,
                f"Nonaxisymmetric scalar transport {field}")
    # On the archived coordinate leg, factor rho*f from c1 and w from the map.
    rows = lookup(Label("C_axial", 2, 0), "z")
    mixed = {row: -s.Rational(2, 3)*weight for row, weight in rows.items()}
    # The f' amplitude is independent of rho*f; its mapped coefficient must vanish.
    require(not lookup(Label("V_axial", 2, 0), "z"), "Radial odd leg has discarded entries")
    require(mixed == {Label("X", 4, 0): 6*s.Rational(4, 105),
                      Label("C_polar", 4, 0): -s.Rational(4, 105)}, "odd-leg map normalization")


def encode_map(rows, row_indices):
    entries = []
    for row, weight in sorted(rows.items(), key=lambda item: row_indices[item[0]]):
        real, imag = s.expand_complex(weight).as_real_imag()
        require(isinstance(real, s.Rational) and isinstance(imag, s.Rational),
                "A projected coefficient is outside Q(i)")
        entries.append({"row": row_indices[row], "real": str(real), "imag": str(imag)})
    return entries


def verify(config, output):
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    checkpoints = output.parent/"atlas-columns"
    checkpoints.mkdir(exist_ok=False)
    with phase("catalog", "catalog"):
        config_bytes = config.read_bytes()
        require(isinstance(json.loads(config_bytes), dict), "Expected config object")
        sources = {name: digest((HERE/name).read_bytes()) for name in ("atlas_core.py", "check_atlas.py", "independent_controls.py")}
        identity = {"config_sha256": digest(config_bytes), "source_sha256": sources,
                    "runtime": {"python": platform.python_version(), "sympy": s.__version__}}
        columns = labels(0, 2)
        rows = labels(3, 4)
        require(len(set(columns)) == len(columns) == 78 and len(set(rows)) == len(rows) == 160,
                "Unexpected coefficient basis dimensions")
        row_indices = {row: i for i, row in enumerate(rows)}
        catalog_hash = publish(checkpoints/"catalog.json", {
            **identity, "labels": ["field", "ell", "m"], "columns": columns, "rows": rows,
            "axes": AXES, "weights": "Cartesian generator, unit w, exact real/imaginary rational strings",
            "archived_regression_reference": ARCHIVED_REGRESSION, "atlas_complete": False})
    with phase("affine_background", "check"):
        for axis in AXES:
            for rotation in (False, True):
                detail("affine_begin", axis=axis, rotation=rotation)
                affine_check(axis, rotation)
                detail("affine_end", axis=axis, rotation=rotation)

    records, maps = [], {}
    progress(phase="atlas_columns", step=0, of=len(columns), unit="column")
    for i, column in enumerate(columns):
        started = time.monotonic()
        detail("column_begin", index=i, label=column)
        mapped = {}
        timings = {}
        for axis in AXES:
            for rotation in (False, True):
                mode = "rotation" if rotation else "null_shift"
                timestamp = time.monotonic()
                detail("map_begin", index=i, axis=axis, kind=mode)
                result = column_map(column, axis, rotation)
                maps[column, axis, rotation] = result
                if rotation:
                    require(not result, f"Discarded rotation {axis} {column}")
                else:
                    mapped[axis] = result
                elapsed = time.monotonic()-timestamp
                timings[f"{mode}_{axis}"] = elapsed
                detail("map_end", index=i, axis=axis, kind=mode, seconds=elapsed)
        entries = {axis: encode_map(mapped[axis], row_indices) for axis in AXES}
        filename = f"column-{i:03d}.json"
        checksum = publish(checkpoints/filename, {
            "catalog_sha256": catalog_hash, "column": i, "label": column, "maps": entries,
            "rotation_zero": True, "round_reconstruction_checked": True, "metric_chain_rule_checked": True,
            "map_seconds": timings, "atlas_complete": False})
        records.append({"file": "atlas-columns/"+filename, "sha256": checksum})
        detail("column_published", index=i, sha256=checksum, seconds=time.monotonic()-started)
        progress(phase="atlas_columns", step=i+1, of=len(columns), unit="column",
                 m={"checkpoint_column_count": i+1})
    def lookup(label, axis, rotation=False):
        # Isolate stored results from deliberate mutations in the independent controls.
        return maps[label, axis, rotation].copy()

    with phase("independent_controls", "check"):
        anchor_controls(lookup)
        independent = verify_independent_controls(lookup)
        require(isinstance(independent, dict), "Missing independent controls record")
    with phase("final_export", "bundle"):
        final = {"state": "complete N=3 discarded coefficient maps", "atlas_complete": True, **identity,
                 "catalog": {"file": "atlas-columns/catalog.json", "sha256": catalog_hash},
                 "columns": records, "coefficient_column_count": len(columns), "discarded_row_count": len(rows),
                 "generator_count": 3, "archived_regression_reference": ARCHIVED_REGRESSION,
                 "independent_controls": independent,
                 "checks": {"affine_retained_with_independent_slope": True, "rotation_discarded_zero": True,
                            "full_coefficient_reconstruction": True, "metric_chain_rule": True,
                            "fourier_reality": True, "axisymmetric_eight_rows": True,
                            "transverse_scalar_phases": True, "odd_coordinate_leg_normalization": True},
                 "scope": SCOPE}
        final_hash = publish(output, final)
    detail("atlas_exported", sha256=final_hash, columns=len(columns))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    verify(args.config, args.output)
