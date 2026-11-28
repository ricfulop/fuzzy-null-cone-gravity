# Direct-state stationary germ checker and supervisor

This is a prepared alternative certificate for the fixed-$\epsilon=1/100$, zero-control oblate system. It starts from the original stationary operator export and a candidate fourteen-state snapshot. The [mathematical contract](../../reports/oblate-stationary-state-germ-preparation-2026-11-28/README.md) gives the interval-local uniqueness argument for the full ODE. This preparation supplies no physical-radius certificate or responding seed.

## Dispatch and evidence

The entry point is `run.sh` in this directory. A tracked invocation must supply a run ID, a full immutable source commit, and a committed sibling `run-config.json`. Under `inputs`, both `stationary_operators` and `terminal_state` require a `path` and a lowercase SHA-256 pin. Relative input paths resolve from the Paper I root. The supervisor verifies the committed local sources, checks their copies, and dispatches the checker on fresh snapshot bytes:

```text
certify_state_germ.py --config CONFIG --operators OPERATORS_COPY --state STATE_COPY --output OUTPUT --max-atan-terms 256
```

One isolated local child uses SymPy 1.14.0 and single-thread settings. The wall budget is 900 seconds. The sampled process-group RSS limit is 2 GiB. The supervisor reserves cleanup time before the 2026-11-28 16:29 UTC cutoff. The command, raw merged child log, input and source snapshots, hashes, samples, group disposition, final metrics and execution outcome belong to the fresh run directory. The certificate path there is `outputs/state-germ-certificate.json`.

## What a pass means

The independent check must establish the kernel and invariance identities, the physical lift, all 22 original retained rows, and a regular physical interval with a rank-six minor of the saved constraint matrix. The supervisor accepts the exported certificate only after a zero child exit, verified empty process group, and successful provenance checks. It checks the result against the actual three input digests and requires consistent response flags.

A certified nonzero response supplies a witness on the forward family. An exact zero establishes silence on that family. If a formally nonzero coefficient cannot be separated at a regular candidate radius, the response remains JSON `null`; numeric telemetry omits its truth value. The response and the certificate pass are separate outcomes. Global radial existence, a converse covering all stationary solutions, and nonlinear consistency remain open.

## Mocked supervisor validation

From `paper1/`, the following commands exited 0. The unittest command produced the full combined output shown below. The shell syntax check produced no output.

```sh
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_state_germ_20261128 -p test_run_bounded.py -q
sh -n verification/oblate_stationary_state_germ_20261128/run.sh
```

```text
----------------------------------------------------------------------
Ran 20 tests in 0.231s

OK
```

The tests use manufactured certificate data and replace child dispatch and process-group operations. They cover copied-input command paths and digests, fresh output, stop and cleanup evidence, and the three response outcomes. The [validation record](../../reports/oblate-stationary-state-germ-checker-preparation-2026-11-28/runner-validation.json) preserves the merged output. Scientific execution remains pending.

| Validated file | SHA-256 |
|---|---|
| `run_bounded.py` | `819631f2cac2977ddbf3e074425e7cf14ebad4d39c0d9d6c7c8edb01a1653be5` |
| `test_run_bounded.py` | `8839bba6daa272b9420d15ff177024e465e19c232d9cad1e92f5c1f7a2b83231` |
| `run.sh` | `baff5b85437ecf781c8729fb83b33c5f199ecbda26545fdd430a51dd98082893` |
| `view-config.json` | `96683101486e99547dc5b3fa6370a4ef0c90a15f9519ba083613f2f8abb1d1b0` |
| `view.xml` | `0b7c584a47b5bd3d61496513e22975d99adae48b818364ed4e1e93e6ea1dfbeb` |
