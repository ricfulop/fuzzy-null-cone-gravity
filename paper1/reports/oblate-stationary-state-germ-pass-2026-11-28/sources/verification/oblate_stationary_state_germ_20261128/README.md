# Direct-state stationary germ checker and supervisor

This is an alternative certificate for the fixed-$\epsilon=1/100$, zero-control oblate system. It starts from the original stationary operator export and a candidate fourteen-state snapshot. The [mathematical contract](../../reports/oblate-stationary-state-germ-preparation-2026-11-28/README.md) gives the interval-local uniqueness argument for the full ODE. No physical-radius certificate or responding seed has yet been exported.

## First attempt and interval repair

The [first direct-state attempt](../../reports/oblate-stationary-state-germ-interval-timeout-2026-11-28/README.md), at source `0ae1c15e11b2b8c9a2f82ca86b6774b68f7e0e79`, passed 18 manufactured controls and the exact kernel, preservation, physical-lift and all-22 retained-row identities. It found a formally nonzero discarded response and collected 57 denominator polynomials. The 900-second wall limit stopped it during the combined regularity, rank-minor and response search. No physical-radius or solution-germ certificate was produced.

The cache repair stores exact polynomial objects and their per-radius coefficient lists, advances the same alternating arctangent brackets incrementally, and retains the original rational Horner and elimination arithmetic. Radius, refinement and response-column order are unchanged. Detailed progress events identify preparation, denominator batches, pivots and serialization without printing large endpoints. Only a dedicated undecided-enclosure exception triggers refinement. The subsequent run `9afcd73b-3499-4f08-8185-0f162c716cfb` passed all 19 mathematical controls, including the cache/reference and determinant-sign cases. It selected a responding candidate at radius 1 with 25 arctangent terms, then failed while converting a minor interval endpoint to decimal text. No certificate was exported. Inputs, proof identities, 900-second wall limit and sampled 2 GiB memory limit are unchanged.

## Dispatch and evidence

The entry point is `run.sh` in this directory. A tracked invocation must supply a run ID, a full immutable source commit, and a committed sibling `run-config.json`. Under `inputs`, both `stationary_operators` and `terminal_state` require a `path` and a lowercase SHA-256 pin. Relative input paths resolve from the Paper I root. The supervisor verifies the committed local sources, checks their copies, and dispatches the checker on fresh snapshot bytes:

```text
certify_state_germ.py --config CONFIG --operators OPERATORS_COPY --state STATE_COPY --output OUTPUT --max-atan-terms 256
```

One isolated local child uses SymPy 1.14.0 and single-thread settings. The wall budget is 900 seconds. The sampled process-group RSS limit is 2 GiB. The supervisor reserves cleanup time before the 2026-11-28 16:29 UTC cutoff. The command, raw merged child log, input and source snapshots, hashes, samples, group disposition, final metrics and execution outcome belong to the fresh run directory. The certificate path there is `outputs/state-germ-certificate.json`.

## What a pass means

The independent check must establish the kernel and invariance identities, the physical lift, all 22 original retained rows, and a regular physical interval with a rank-six minor of the saved constraint matrix. The supervisor accepts the exported certificate only after a zero child exit, verified empty process group, and successful provenance checks. It checks the result against the actual three input digests and requires consistent response flags.

A certified nonzero response supplies a witness on the forward family. An exact zero establishes silence on that family. If a formally nonzero coefficient cannot be separated at a regular candidate radius, the response remains JSON `null`; numeric telemetry omits its truth value. The response and the certificate pass are separate outcomes. Global radial existence, a converse covering all stationary solutions, and nonlinear consistency remain open.

`terminal-evidence.json` is the terminal commit marker. Its SHA-256 fields must match both `execution.json` and `metrics.json` before a reader accepts a pass; standalone metrics are provisional. Cancellation is accepted through certificate validation and frozen before terminal evidence publication. Later signals cannot change that decision. Console forwarding is nonblocking and best effort, while the child log is archived independently. No RSS samples means an unknown (`null`) maximum.

## Initial preparation validation

At preparation commit `0ae1c15`, the following commands exited 0 from `paper1/`. The unittest command produced the full combined output shown below. The shell syntax check produced no output. These recorded hashes describe that preparation, before the later supervisor review fixes.

```sh
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_state_germ_20261128 -p test_run_bounded.py -q
sh -n verification/oblate_stationary_state_germ_20261128/run.sh
```

```text
----------------------------------------------------------------------
Ran 20 tests in 0.231s

OK
```

The tests use manufactured certificate data and replace child dispatch and process-group operations. They cover copied-input command paths and digests, fresh output, stop and cleanup evidence, and the three response outcomes. The [validation record](../../reports/oblate-stationary-state-germ-checker-preparation-2026-11-28/runner-validation.json) preserves the original merged output. The maintained suite adds console-backpressure, terminal-commit, cancellation and missing-measurement regressions. These infrastructure tests execute no scientific checker.

| File at original supervisor validation | SHA-256 |
|---|---|
| `run_bounded.py` | `819631f2cac2977ddbf3e074425e7cf14ebad4d39c0d9d6c7c8edb01a1653be5` |
| `test_run_bounded.py` | `8839bba6daa272b9420d15ff177024e465e19c232d9cad1e92f5c1f7a2b83231` |
| `run.sh` | `baff5b85437ecf781c8729fb83b33c5f199ecbda26545fdd430a51dd98082893` |
| `view-config.json` | `96683101486e99547dc5b3fa6370a4ef0c90a15f9519ba083613f2f8abb1d1b0` |
| `view.xml` | `0b7c584a47b5bd3d61496513e22975d99adae48b818364ed4e1e93e6ea1dfbeb` |

## Integrated cache retry

The [integrated preparation record](../../reports/oblate-stationary-state-germ-cache-retry-preparation-2026-11-28/README.md) combined the reviewed interval-cache checker with the corrected supervisor from `698d7fcd65c3e7eb9f25f6d311d96af27a73e069`. Its captured 27-test mocked infrastructure validation is separate from the original 20-test record above. At preparation time its 19 mathematical controls were unexecuted; their subsequent pass belongs to run `9afcd73b`, which failed during certificate serialization. The input pins, proof obligations, search order, 900-second wall limit and sampled 2 GiB limit remain unchanged. An invocation must produce an authenticated physical-radius certificate before any germ or response conclusion is accepted.

## Exact endpoint serialization repair

The [serialization failure archive](../../reports/oblate-stationary-state-germ-serialization-failure-2026-11-28/README.md) preserves the raw candidate-selection event, the Python 4,300-digit conversion error, and authenticated failed terminal evidence. That run supplied no certificate. The [new preparation record](../../reports/oblate-stationary-state-germ-serialization-preparation-2026-11-28/README.md) changes rational endpoint output to exact base-$10^9$ chunk formatting. It leaves the field identities, interval arithmetic and search decisions unchanged, and does not alter Python's global digit-limit setting. Both the helper and its tests are included in the supervisor's authenticated source snapshots. Six stdlib serialization tests and the 27 mocked supervisor tests pass; an added large-rational mathematical control brings the next tracked invocation to 20 controls. That new scientific control has not yet run.
