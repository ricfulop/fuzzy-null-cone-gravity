"""Standard-library integrity and layout checks for exact-leading exports.

The sidecar is produced alongside the unchanged list-shaped data file. It
binds that file to the producer source and explicit ordered basis conventions;
it is not an additional symbolic verification of the matrix entries.
"""

from pathlib import Path
import hashlib
import json


SCHEMA = "oblate-exact-leading-v1"
SOURCE_FILES = ("all_modes.py", "leading_certificate.py")
COMPONENTS = ("lambda", "X", "V_polar", "V_axial", "C_polar", "C_axial")
SECTORS = tuple((L, m) for L in (2, 3) for m in range(L + 1))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def source_fingerprints(directory):
    directory = Path(directory)
    return {name: sha256_bytes((directory / name).read_bytes())
            for name in SOURCE_FILES}


def shift_labels(L, m):
    labels = []
    for ell in range(m, L + 1):
        labels.append(["b", ell])
        if ell >= 1:
            labels.extend([["W_polar", ell], ["W_axial", ell]])
    return labels


def normal_labels(L):
    return [[name, ell] for ell in (L + 1, L + 2) for name in COMPONENTS]


def basis_layout():
    return [{"L": L, "m": m, "shift_labels": shift_labels(L, m),
             "normal_row_labels": normal_labels(L)} for L, m in SECTORS]


def _same_layout(actual, expected):
    # JSON spelling distinguishes bool/float labels from the required integers.
    return (json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False))


def _validate_sources(sources):
    require(isinstance(sources, dict) and set(sources) == set(SOURCE_FILES),
            "exact-leading certificate source inventory changed")
    require(all(isinstance(value, str) and len(value) == 64 and
                all(char in "0123456789abcdef" for char in value)
                for value in sources.values()),
            "exact-leading certificate has invalid source digests")


def _matrix(value, rows, columns, name):
    require(isinstance(value, list) and len(value) == rows,
            f"{name}: unexpected row count")
    for row in value:
        require(isinstance(row, list) and len(row) == columns,
                f"{name}: unexpected column count")
        require(all(isinstance(entry, str) and entry.strip() for entry in row),
                f"{name}: expected nonempty exact-expression strings")


def validate_rows(rows):
    require(isinstance(rows, list) and len(rows) == len(SECTORS),
            "exact-leading certificate must contain all seven ordered sectors")
    totals = {2: 0, 3: 0}
    for row, (L, m) in zip(rows, SECTORS):
        require(isinstance(row, dict), "exact-leading sector must be an object")
        require(all(type(row.get(key)) is int for key in ("L", "N", "m")) and
                (row["L"], row["N"], row["m"]) == (L, L + 1, m),
                "exact-leading sector order or identity changed")
        labels = shift_labels(L, m)
        require(_same_layout(row.get("labels"), labels), "exact-leading shift label order changed")
        require(_same_layout(row.get("normal_row_labels"), normal_labels(L)),
                "exact-leading normal basis row order changed")
        columns = len(labels)
        _matrix(row.get("leading_DeWitt_Gram"), columns, columns, "leading Gram")
        _matrix(row.get("normal_coefficients"), 12, columns, "normal coefficients")
        expected = 3 if m == L else 4
        require(all(type(row.get(key)) is int and row[key] == expected for key in
                    ("exact_leading_rank", "exact_leading_normal_map_rank")),
                f"unexpected exact-leading sector rank for L={L}, m={m}")
        selected = row.get("nonzero_principal_minor_indices")
        require(isinstance(selected, list) and len(selected) == expected and
                all(type(index) is int and 0 <= index < columns for index in selected) and
                len(set(selected)) == expected,
                "exact-leading principal minor indices are invalid")
        determinant = row.get("nonzero_principal_minor_determinant")
        require(isinstance(determinant, str) and determinant.strip(),
                "exact-leading principal minor determinant is missing")
        minors = row.get("positive_Sylvester_minors")
        require(isinstance(minors, list) and len(minors) == expected and
                all(isinstance(value, str) and value.strip() for value in minors),
                "exact-leading Sylvester minor list is invalid")
        totals[L] += (1 if m == 0 else 2) * row["exact_leading_rank"]
    require(all(totals[L] == 8 * L + 2 for L in totals),
            "exact-leading real rank totals changed")


def manifest_path(path):
    path = Path(path)
    return path.with_name(path.stem + "-manifest.json")


def write_certificate(path, rows, sources):
    """Write data plus metadata only after all structural checks pass."""
    path = Path(path)
    _validate_sources(sources)
    validate_rows(rows)
    data = (json.dumps(rows, indent=2, allow_nan=False) + "\n").encode("utf-8")
    manifest = {"schema": SCHEMA, "producer": "all_modes.exact_leading",
                "source_sha256": sources, "basis_layout": basis_layout(),
                "payload_file": path.name, "payload_sha256": sha256_bytes(data)}
    path.write_bytes(data)
    manifest_path(path).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def load_certificate(path, expected_sources):
    """Read each file once, check metadata, then authenticate before parsing data."""
    path = Path(path)
    _validate_sources(expected_sources)
    sidecar = manifest_path(path)
    try:
        manifest_data = sidecar.read_bytes()
    except FileNotFoundError as error:
        raise ValueError(
            "exact-leading export lacks its certificate manifest; legacy exports are "
            "not accepted: regenerate with the current all_modes.py producer"
        ) from error
    manifest = json.loads(manifest_data)
    require(isinstance(manifest, dict), "exact-leading manifest must be an object")
    require(manifest.get("schema") == SCHEMA and
            manifest.get("producer") == "all_modes.exact_leading",
            "unsupported exact-leading certificate; regenerate with the current producer")
    require(manifest.get("source_sha256") == expected_sources,
            "exact-leading producer source hashes differ; regenerate with the current producer")
    require(_same_layout(manifest.get("basis_layout"), basis_layout()),
            "exact-leading certificate basis layout changed")
    require(manifest.get("payload_file") == path.name,
            "exact-leading certificate names a different data file")
    data = path.read_bytes()
    digest = sha256_bytes(data)
    require(manifest.get("payload_sha256") == digest,
            "exact-leading data hash does not match its certificate")
    rows = json.loads(data)
    validate_rows(rows)
    return rows, {"schema": SCHEMA, "source_sha256": expected_sources.copy(),
                  "payload_sha256": digest,
                  "manifest_sha256": sha256_bytes(manifest_data)}
