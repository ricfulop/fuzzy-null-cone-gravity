# Stationary radial extraction: incomplete at the wall-time limit

**The single approved invocation was stopped at its five-minute limit.** It did not complete the stationary operator export or final controls. This is an incomplete verification, not a counterexample to the proposed identities or a nonlinear-consistency result.

The source was frozen at `de24fa2bb5c21c2b027498581872e64cb65da445`. User approval `msg_fa2e97d4` covered one local CPU invocation with a 300-second wall budget and a 2-GiB monitored memory budget. The original [preparation](../oblate-stationary-seed-preparation-2026-11-27/README.md) and [approved proposal](approved-run-proposal.json) are preserved.

## Observed outcome

| Item | Recorded outcome |
|---|---|
| Local date | 27 November 2026, America/New_York |
| UTC start | 28 November 2026, 03:58:54.154718 |
| UTC stop | 28 November 2026, 04:03:54.336729 |
| Supervised elapsed time | 300.181 seconds |
| Stop reason | Wall budget reached |
| Child exit code | 143 |
| Supervisor exit code | 1 |
| Largest of 1,031 process-group RSS samples | 277,037,056 bytes, 264.203125 MiB |
| Scientific invocations | 1 |
| GPU, remote, and tracked experiment jobs | 0 |

The [full merged stdout/stderr](verification.log) records these five column-completion messages:

1. `lambda, 0`
2. `lambda, 2`
3. `X, 0`
4. `X, 2`
5. `V_polar, 2`

The log contains no assertion traceback and no final PASS message. Three of the eight planned columns have no completion message. The final global row, gauge/family, rejection, and quadrature controls follow the complete assembly and were not reached. Neither `stationary-operators.json` nor `stationary-reproducibility.json` was produced. The in-memory partial matrices were not checkpointed, so the column messages are not a recoverable partial operator export.

The [execution record](execution.json), [resource samples](resource-samples.json), and [supervisor](supervise.py) preserve the command, approval, process information, selected environment overrides, source hashes, and stopping procedure. Memory was sampled approximately every 0.25 seconds; its maximum is not a measurement of an unobserved inter-sample peak or a hard OS memory cap. The observed memory use did not trigger the stop.

## What can and cannot be inferred

The symbolic work exceeded the initial runtime estimate. The exact expensive operation at termination was not instrumented, and no profiling or scientific residual inspection has been performed after the stop. It would be premature to attribute the delay to a specific simplification, angular integral, or curvature expression.

No assertion failure was reported before termination, but the incomplete invocation does not certify the full stationary extraction. It establishes neither a radial differential-algebraic reduction nor an admissible non-gauge solution germ. The discarded response has not been evaluated on such a solution. Prior completed rank, lapse, map/control, and companion kinematic results are unchanged; the full nonlinear release gate remains open.

## Provenance and next decision

The command ran from `paper1/`:

```sh
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 python verification/oblate_stationary_seed_20261127/check_stationary_system.py
```

The command pins NumPy 2.3.5 and SymPy 1.14.0. The checker did not reach its successful-run dependency/version manifest, so this archive distinguishes the requested dependency pins from a completed checker-generated environment certificate. The supervisor's Python version, platform, launcher version, and thread-limit overrides are recorded separately.

All four scientific source/action hashes matched before and after execution. Their exact snapshots are under [sources](sources/). No scientific source was changed, no automatic retry occurred, and no additional computation was launched after the stop.

The proposed next step is permission to inspect the computational bottleneck and prepare per-column timing/checkpointing or another justified optimization while preserving the mathematical equations and controls. A concrete rerun would receive separate approval. The current archive authorizes neither that correction nor a larger compute budget.
