# Oblate stationary reduction preparation

This directory prepares a JSON-only reducer for the stationary oblate
operator export.  It deliberately imports no stationary extraction modules:
the only scientific input is a completed
`paper1/reports/oblate-stationary-system-recheck-2026-11-27/outputs`
directory containing:

- `completion.json`
- `stationary-operators.json`
- `stationary-reproducibility.json`

The reducer fixes `epsilon=1/100` and rewrites every exported expression in
the differential field `Q(r, T)`, where

```text
T = atan(1/(10*r)),       dT/dr = -10/(100*r**2 + 1).
```

Its first obligations are structural:

1. Refuse checkpoint-only input.
2. Recompute the hashes recorded in the completion marker.
3. Check the expected 8 columns, 22 retained rows, 3 derivative orders, and
   final discarded response row.
4. Convert expressions at fixed epsilon and verify the conversion is reversible
   before replacing `atan(1/(10*r))` by `T`.
5. Verify the radial right-kernel on all retained rows and on the response row.
6. Remove the `X0` column only after the gauge check, while keeping the old
   `X0` equation among the selected rows.

The prepared DAE closure works with

```text
z = (v, v')
```

for the seven gauge-fixed profiles.  It starts from the eight selected retained
equations

```text
M2 v'' + M1 v' + M0 v = 0,
```

solves their highest-derivative block for as many accelerations as possible,
records residual constraints, and iteratively preserves each new constraint
by reducing

```text
C' z + C K z + C B f = 0.
```

Remaining controls are retained explicitly.  The final response is serialized
as a function of both reduced free coordinates and those controls, so a probe
cannot hide a source by setting all controls to zero.

No radial solution germ, non-gauge seed, nonlinear obstruction, or canonical
Dirac classification follows from preparing this reducer.

## First tracked invocation

`run-config.json` selects only the principal probe. The fixed command is
`bash paper1/verification/oblate_stationary_reduction_20261128/run.sh`.
Run it through the project experiment engine, not directly: the supervisor
requires a run ID and full commit, verifies the materialized source bytes
against that commit, and uses one isolated local SymPy 1.14.0 child.
The reviewed upstream completion marker is pinned independently at SHA-256
`640ae2f7c9e1cd3702f8330b074df361ee3fdff701ee2a566ef2b84c4e3a6b4d`.

The envelope is 900 seconds, a 2 GiB sampled process-group RSS stop threshold,
and no GPU or remote compute. The supervisor also rejects a launch whose
budget and cleanup margin would exceed the temporary approval window ending
2026-11-28 16:29 UTC. The RSS threshold is not a hard operating-system cap.

Two manufactured systems exercise hidden differentiated constraints, retention
of a free acceleration in the response, and variable-coefficient constraint
derivatives. The actual principal probe then checks exact rational-field
membership, reversible conversion, and the full radial gauge identity before
computing the highest-derivative rank. This rank is not a Dirac rank.

The complete raw child output, sampled resources, source bytes and before/after
hashes are saved under the immutable run snapshot's `paper1/run-output/RUN_ID/`.
A successful principal probe still supplies no DAE closure, radial solution,
or non-gauge compatibility source. Future closure is a separate child node
using the same run command with a changed committed configuration.

Static review: the algebraic elimination was reviewed in `msg_98e342dd`;
its rational-field and input-pinning requests were incorporated. Runner review
`msg_879491a2` requested commit binding, a second cutoff guard, and post-cleanup
log draining; those changes are incorporated. No scientific function was
executed during this preparation.

## Closure child

The child configuration now selects `closure`, building on the passed rank-4
principal probe (run e04ffc4a-d9c3-4c8a-ad61-3bf8da079ff3). The unchanged
formal elimination path differentiates constraints and checks the original
22 rows plus all remaining controls. Final JSON and telemetry report only
formal dimensions and response; an actual regular interval is still required.
The fixed command and resource envelope are unchanged.
