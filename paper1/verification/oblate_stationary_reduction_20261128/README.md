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

The first invocation's `run-config.json` selected only the principal probe. The fixed command is
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

## Repair after the incomplete closure attempt

Run `ef17ca1d-2663-4c0d-acc8-2f8607ae7881` did not complete closure. Its
original evidence is preserved separately in
`reports/oblate-stationary-closure-attempt-2026-11-28/`. The scientific log
finished step 0 and began step 1; process-group cleanup then failed with
`PermissionError`, including inside the final evidence-writing path.

The first retry replaced expression-by-expression row reduction by exact
`DomainMatrix` operations over `QQ(r,T)`. It clears row denominators and uses
fraction-free elimination on the whole augmented `[A | C]` matrix. Only pivots
in `A` solve controls. All remaining state-only rows are retained. Reducing the
upper state block can change the displayed control formulas by multiples of
the residual state constraints; the solution set is unchanged on that
constraint surface. The differential derivation and all 22 original-row
terminal checks are unchanged.

Independent-constraint selection uses one transpose pivot calculation and
returns the original selected candidate rows in their original order. Five
additional manufactured row-operation tests cover control/state pivots, a free
control, empty dimensions, dependent candidate rows, and inexact-input
rejection. The two prior differential-reduction controls remain. These tests
run before the actual closure within the bounded scientific invocation.

New timings separate solving controls, applying substitutions, selecting
constraints, checkpoint serialization, and differentiation. An interim
principal-rank headline is available before the final closure signals. Partial
step counts never certify a solution. That invocation passed its controls but
was stopped after its first constraint rows grew substantially and their
differentiation exceeded five minutes. Its complete evidence is preserved in
`reports/oblate-stationary-closure-basis-growth-2026-11-28/`.

## Control-only rational-field arithmetic

The next variant restores the original first-nonzero control pivot order,
paired row swaps, normalization, and elimination order. It converts both input
blocks once to `QQ(r,T)` field elements and performs exact field arithmetic
inside that historical loop. It never pivots the state block. By induction
over the row operations, the resulting entries are the same rational functions
as in the original expression-based implementation. Batched independent-row
selection still returns the original candidate rows in order.

Six row-operation controls now check the historical residual rows entrywise,
including a row swap with distinct state entries, a free control, preservation
of the state rows when there are no control columns, empty dimensions,
dependent candidate rows, and inexact-input rejection. The two original DAE
controls and every final retained-row identity remain. Run `8fa984ea-3c9b-4e85-b984-c23956fc95d8` completed three elimination
steps, then reached the 900-second wall budget during the combined terminal
embedding stage. Its clean shutdown, partial checkpoints and full evidence
are preserved in `reports/oblate-stationary-closure-terminal-timeout-2026-11-28/`.


## Terminal kernel repair

The isolated diagnostic `6128c4ef-4b35-4ddf-9763-1b40369bd6f6` verified the
saved 6-by-14 constraint matrix's rank-six kernel in 117.565 seconds. It used
C-only fraction-free RREF and an exact-field product with the original C.
The diagnostic checked neither dynamics nor the original retained equations.
Its report is `reports/oblate-terminal-embedding-probe-2026-11-28/`.

The prepared full reducer uses that terminal kernel method while preserving
the historical control-only discovery basis and every terminal equation check.
It starts again from the original completed operator export. Each kernel and
subsequent terminal operation is timed separately. No earlier checkpoint is
loaded as a substitute for reconstructing the dynamics. A successful full
closure still requires a separate interval-regularity and responding-germ
certificate before any nonlinear conclusion.


## Terminal preservation repair

The full run `157032da-1468-473e-8cd7-37763ff86822` completed the terminal
kernel in 114.775 seconds, then stopped at the 900-second wall budget during
construction of DC+CK. The log does not distinguish the derivative, product,
and simplification costs. The subsequent multiplication by M was unreached.
The failure is archived in `reports/oblate-stationary-preservation-timeout-2026-11-28/`.

The source prepared for run `5d9b96ad` changed only terminal preservation to sparse exact-field
arithmetic, with D=partial_r-10/(100r^2+1) partial_T. It requires both CB=0 and
(DC+CK)M=0 and omits zero derivative entries from sparse storage. Discovery
uses its original function and basis. All 22 original-row, reduced-dynamics,
reconstruction and response checks remain. New positive and negative controls
exercise a nonzero DC+CK whose projection vanishes, independent CB/SM defects,
empty shapes and exact serialization.

After the successful kernel, a fresh no-overwrite terminal-state JSON records
C,K,B,M and their labels/identity before preservation. It is explicitly
uncertified and write-only: there is no loader, resume path, closure or germ
claim from its presence. The controls were unexecuted at preparation time;
they subsequently passed in `5d9b96ad`, together with the terminal preservation
identities. That run's evidence is preserved in
`reports/oblate-stationary-dynamics-timeout-2026-11-28/`.

## Terminal Euler and reduced-dynamics field arithmetic

The 5d9b96ad attempt passed exact preservation in 75.004 seconds and all original retained-row checks in 317.070 seconds, then reached the wall budget in reduced dynamics. The implementation frozen at `96916c7` retains discovery through preservation and the discarded response unchanged. Only terminal Euler products and reduced-dynamics/reconstruction arithmetic use the exact field.

For the original 22-row blocks E0,E1,E2, seven acceleration rows Ka,Ba, and terminal embedding M, the helper checks ([E0|E1]+E2 Ka)M=0 and E2 Ba=0 independently. Removing the X0 column never removes its equation row. The 23rd discarded response is excluded from the zero requirement.

The helper computes W=KM-DM once, selects ordered free rows J to obtain Kred=JW and Bred=JB, then checks MKred-W=0 and MBred-B=0 on every row. Products, derivatives, extraction and conversion to expressions have separate timers. Controls test retained X0 and outside-discovery row defects, nonzero remaining controls, non-free-row reconstruction defects, empty widths and inexact/nonfinite entries. These controls passed in `c2a1482f`. The original-row stage took 50.614 seconds and KM took 240.467 seconds; the run then reached its wall limit during DM. Reconstruction, response, and formal closure did not complete. See `reports/oblate-stationary-field-dynamics-timeout-2026-11-28/` for the unchanged raw evidence.

That invocation used a 900-second wall limit. The current configuration, frozen at `c763d3a`, raises the limit to 1800 seconds with the same scientific Python files and shell entry point, a sampled 2-GiB limit, and the existing approval cutoff. See `reports/oblate-stationary-extended-budget-preparation-2026-11-28/`. It remains a fresh original-export calculation without partial-state loading. Formal closure and any independent analytic-germ certificate remain separate milestones.

## Maintained supervisor and input validation

The maintained wrapper now executes the verified source copies under its fresh
output directory and rechecks both the originals and executed copies afterward.
The reducer hashes each input byte buffer before decoding it,
and its recorded digests describe those same buffers; input paths are not reread
between authentication and parsing. The standalone manufactured-control entry
point rejects optimized Python, which would disable its assertions.

Child output is archived before best-effort nonblocking console forwarding.
A closed stderr mirror cannot suppress the reducer's original failure record;
errors writing the authoritative journal still fail the attempt. The supervisor
checks the deadline after observing child exit, caps process-table sampling by
the remaining budget, and records an absent RSS sample set as unknown. It stops
accepting cancellation immediately before publishing terminal evidence, so
signals received earlier determine one consistent status and exit decision.
The accompanying mocked and synthetic-process tests exercise infrastructure
only; these changes add no scientific result and do not alter archived sources.
