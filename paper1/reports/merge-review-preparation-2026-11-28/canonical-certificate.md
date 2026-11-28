# Exact-leading certificate migration

The maintained producer now checks each exact sector rank against its expected
value (three for `m=L`, four otherwise) and the real total against `8L+2`.
It retains the list-shaped `exact-leading.json` format, adds explicit ordered
normal-row labels, and writes an `exact-leading-manifest.json` sidecar binding
the payload digest, producer-source digests and ordered basis layout.

The lapse consumer requires that sidecar, checks its source and layout metadata,
authenticates the payload before parsing it, and validates sector identities,
labels, matrix shapes, ranks and principal-minor indices. Its provenance records
the digests of the exact payload and sidecar bytes consumed. The helper loads by
its path beside `all_modes.py`, including when a caller loads that module with
`importlib.spec_from_file_location` outside the source directory.

All checked-in historical result JSON, report text and source snapshots remain
unchanged. The new consumer intentionally rejects legacy exports without the
sidecar: a future lapse execution first needs newly generated authenticated
producer output from the matching sources, in a separately authorized execution
that preserves the historical files. No historical certificate was retrofitted,
and no scientific calculation was rerun for this preparation.

Validation here is limited to synthetic standard-library tests and syntax
inspection. These checks cover integrity, provenance and data layout; they do
not independently recompute or certify the symbolic matrix entries.
