# Stationary radial operator extraction: passed

**The single approved, instrumented invocation completed all eight columns and every prepared control.** It exported three exact $23\times8$ radial coefficient matrices: 22 retained axisymmetric Euler rows and one separate discarded response row. The radial differential-algebraic system is not yet solved; no admissible non-gauge solution germ or nonlinear obstruction has been established.

The scientific source is `1cd1a69c258810d38c12adc334733aff1fac6b7a`. The earlier 300-second timeout remains separately archived at [its original report](../oblate-stationary-system-check-2026-11-27/README.md). This invocation starts from scratch and does not reuse a checkpoint. It was authorized by user `msg_4c4bcc8d` and the subsequent approval of other work during the same overnight window, conservatively ending 28 November 2026 at 16:29 UTC.

## Measured outcome

| Item | Outcome |
|---|---|
| Supervisor / child exit | 0 / 0; no stop reason |
| Supervised elapsed time | 463.842 seconds |
| Checker time before final JSON writes | 462.145 seconds |
| Largest of 1,593 process-group RSS samples | 288,686,080 bytes = 275.3125 MiB |
| Completed column checkpoints | 8 |
| Raw operator dimensions | 3 derivative orders, each $23\times8$ |
| Independent nonzero-density quadrature comparisons | 772 |
| Largest scaled quadrature error | $1.769991488445588\times10^{-13}$ |
| Quadrature acceptance threshold | $10^{-11}$ |
| Scientific invocations / remote / GPU / tracked experiment jobs | 1 / 0 / 0 / 0 |

The RSS maximum is a sampled aggregate, not a measured instantaneous peak or a hard OS limit. The prepared budget was 900 seconds with a 2-GiB sampled RSS threshold. The [full merged log](verification.log) records the stage events and all three final PASS lines. The [execution record](execution.json), [resource samples](resource-samples.json), [startup identity](outputs/invocation.json), and [completion marker](outputs/completion.json) separate process exit, checkpoint publication, and successful export.

The actual checker environment was Python 3.12.11, NumPy 2.3.5, and SymPy 1.14.0. The supervisor environment is recorded separately. All frozen source hashes match before and after execution; export and checkpoint hashes are recorded in the completion marker.

## What was checked

The stationary columns are
$$
u=(\lambda_0,\lambda_2,X_0,X_2,V_{\mathrm p2},C_{\mathrm p2},a_0,a_2).
$$
The exported raw differential expression is
$$
E(u)=M_0(r,\epsilon)u+M_1(r,\epsilon)u'+M_2(r,\epsilon)u'',
$$
with no derivative acting on $M_j$ in this convention. The last row is
$$
D=6\mathscr E_{X,4}-\mathscr E_{C_{\mathrm p},4};
$$
it is a discarded response, not an extra retained equation. The angular average is $dx/2$, and the background density $\Sigma$ already occurs once in the Euler pairing.

The completed controls cover:

- The mixed-parity two-channel response identity.
- The full four-dimensional/static spatial-lapse split for every column. This is a dimensional cross-check using the same curvature implementation, not an independent curvature implementation.
- All 22 retained tests, including exact vanishing of the 14 non-selected rows across all columns and derivative orders.
- Arbitrary radial gauge profiles, constant lapse rescaling, and the exact oblate-parameter tangent, including the discarded response row.
- Rejection of a deliberately off-shell lapse profile, with finite-value checks.
- Exact angular moments against independent 48- and 96-node quadrature for every nonzero density at $\epsilon=1/5,1/100$ and $r=3/2$. This numerical control does not freeze $r$ in the symbolic export.

The eight selected rows are selected by symmetry, not asserted to be eight independent equations. Their rank, differential consequences, and local solution space require the next reduction.

## Timing diagnosis

[Stage timings](timing-summary.json) locate the dominant work: four-dimensional curvature construction took 244.017 seconds across the columns, and the static-split controls took 171.664 seconds. Together these two non-overlapping stages account for about 90% of the supervised runtime. Euler-density contraction took 4.923 seconds, radial coefficient extraction 2.819 seconds, and angular averaging 9.648 seconds. Column checkpoint publication took 0.463 seconds in total.

The former hypothesis that repeated row contractions might dominate is not supported by these measurements. No optimization was needed to complete this invocation. The timing table also contains enclosing stages; their durations overlap with nested entries and must not all be summed. These measurements concern this single run, not a controlled performance comparison with the earlier timeout.

## Reproducibility and next scientific step

From `paper1/` the supervisor was invoked as:
```sh
python3 verification/oblate_stationary_seed_20261127/supervise_stationary_check.py \
  --source-commit 1cd1a69c258810d38c12adc334733aff1fac6b7a \
  --approval-message msg_4c4bcc8d
```

It launched the pinned local command recorded in `execution.json`. The supervisor and checker preserve source snapshots, complete stdout/stderr, actual environment metadata, exact column checkpoints, and the final [operator export](outputs/stationary-operators.json). The [reproducibility manifest](outputs/stationary-reproducibility.json) and [completion marker](outputs/completion.json) link those objects by hash. A future invocation must use a new output directory; the archive must not be overwritten.

The next step is an exact reduction of the radial equations, retaining all original rows as the certificate. A responding seed must solve those equations throughout an open radial interval. A consistent jet or numerical nullspace at one radius is insufficient. A possible outcome is that the discarded response vanishes identically on the stationary solution space; the export does not prejudge that test.

The prior oblate rank and linearized-lapse results are unchanged. The final nonlinear constraint classification, canonical boundary realization, physical degree-of-freedom count, and full Paper I completion goal remain open.
