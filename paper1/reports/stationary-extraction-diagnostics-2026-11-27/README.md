# Stationary extraction: diagnostic preparation

**Preparation complete; the scientific checker has not been rerun.** User message `msg_4ba97146` authorizes bottleneck inspection and preparation of timing, checkpoints, or a justified optimization, and removes the need for future approval for diagnostic/preparation work of this kind. Scientific reruns and experiment launches retain their separate approval requirement.

The [original timeout](../oblate-stationary-system-check-2026-11-27/README.md) remains unchanged at result commit `7901d07a2abe291ce4a4073f3f7cd0593a50882e`, with scientific source `de24fa2bb5c21c2b027498581872e64cb65da445`. All changes here are on a child branch. No manuscript claims or mathematical definitions have been changed.

## What the diagnosis establishes

The old log contains five completed-column messages and no operation-level timestamps. The supervisor stopped the single invocation at 300.181 seconds; the largest sampled process-group RSS was 264.203125 MiB. There is no assertion traceback, completed operator export, or saved partial matrix. These facts identify a wall-budget stop but **do not identify the expensive symbolic operation**. The sixth column had no completion message; curvature construction, its static-split control, density contraction, jet extraction, and angular averaging cannot be distinguished from that log. The old code had no per-column serialization; final serialization was not reached.

Static inspection exposes possible repeated work: each of 23 row tests raises the same four-dimensional Einstein tensor; the fixed test tensors recur across columns; the fixed background metric, inverse, and connection are reconstructed repeatedly. These are optimization candidates, not measured bottlenecks. Quadrature follows full assembly and was not reached. No symbolic algorithm has been optimized in this preparation.

## Prepared changes

The assembler keeps the original eight columns, all 22 retained rows, and the additional discarded response
$$
D=6\mathscr E_{X,4}-\mathscr E_{C_{\mathrm p},4}.
$$
It differentiates the same reconstructed metric before extracting radial coefficients. The moment recurrence, static-split equations, full continuum Euler evaluator, harmonic definitions, and action are unchanged.

The added observation and storage code provides:

- Flushed, timestamped begin/end/error events for coefficient reconstruction, four-dimensional curvature, static-split verification, each Euler-density contraction, radial coefficient extraction, each angular average, checkpoint serialization/publication, and the final controls. Nested stage IDs distinguish a stopped inner operation from its enclosing stage. Quadrature case markers identify deformation and node count.
- A fresh output directory for every invocation. The checker refuses an existing output directory. Startup provenance records the actual Python/dependency versions, symbol assumptions, source hashes, coordinate and density conventions, command arguments, selected thread limits, and an attempt ID before expensive assembly.
- An exact JSON checkpoint after each completed column. Each includes all three derivative orders of all 23 density coefficients and averaged coefficients, including zeros, as strings; the row/column order, selected and automatic row indices, and static-split status are explicit. Moment requests are labeled **cumulative through the column**, not column-local. Exact nonfinite expressions are rejected before publication.
- Atomic, non-overwriting publication of checkpoints, with file and payload hashes. The file is flushed and fsynced before publication. Progress confirms a checkpoint only after publication succeeds. An interrupted current column need not be recoverable; completed published columns remain available.
- A final completion marker only after all original controls, unchanged-source checks, and final export writes succeed. It records the source hashes, export hashes, and checkpoint list. Partial checkpoints explicitly say that global controls are pending and do not establish a seed or nonlinear result.

The checkpoint implementation is **write-only**. There is no automatic resume, expression parser, or scientific cache reuse. A future loader must validate source/dependency identity, symbol assumptions, every row/column and shape, and exact expressions, then rerun the final controls. The earlier timed-out invocation cannot be recovered from these new checkpoints because it never wrote them.

The completion marker certifies export/control completion, not a clean process exit: a later logging error or supervisor stop must still be reported. A hard stop can also leave a published checkpoint before its progress event, or an unpublished `.partial-*` file. Inspect the published envelope and hashes rather than inferring publication solely from the last event. Incomplete source snapshots are not valid provenance; their hashes must match the startup identity. The checker records time through expression serialization before JSON writes; the supervisor supplies end-to-end duration.

The unchanged final controls are the two-channel identity, static dimensional split in every column, all 14 non-selected retained rows, radial gauge including the response, constant lapse, oblate-parameter family, finite nonzero off-shell rejection, and 48/96-node independent quadrature at both deformations. The partial checkpoints do not certify these final controls.

## Preparation validation

[The validation record](preparation-validation.json) distinguishes software checks from scientific execution. The final six stdlib-only infrastructure regression tests passed in 0.016 seconds. They exercise refusal to overwrite evidence, interrupted atomic publication, preservation of partial evidence after a manufactured failure without a completion marker, rejection of nonfinite JSON or paths escaping the output directory, rejection of optimized-Python launch settings, and cleanup of mocked surviving descendants despite launcher exit. The process tests use mocks and launch no child or OS signal. The expression strings are inert toy data; no SymPy or NumPy module was imported by those tests.

AST inspection checks syntax and confirms unchanged analytic helper definitions. The full continuum evaluator and action bytes still match the frozen scientific source. Independent reviews are recorded in [review.json](review.json). No scientific checker, profiling run, seed calculation, remote job, or GPU job has been executed in this stage. The instrumentation's scientific path remains to be checked in the separately approved invocation.

## Proposed single check

The [run proposal](run-proposal.json) requests one supervised local CPU invocation from the committed preparation source, with one process, a **900-second wall budget** and **2 GiB sampled process-group RSS threshold**. This is a proposed budget, not a measured runtime prediction or a hard OS memory cap. The original five-minute limit was insufficient; completion within 15 minutes is still uncertain. Checkpoints and stage records preserve useful evidence if the bound is reached again.

From `paper1/`, the prepared supervisor takes the pinned commit and the new approval message:
```sh
python3 verification/oblate_stationary_seed_20261127/supervise_stationary_check.py \
  --source-commit <preparation-commit> --approval-message <new-run-approval>
```

It uses explicit launch guards for a clean worktree and source hashes, rejects settings that disable Python assertions, launches the same pinned NumPy/SymPy command, retains merged stdout/stderr, samples RSS, and terminates the process group at the wall or memory threshold, including surviving descendants after launcher exit. Termination has a short cleanup grace period beyond the stop trigger. It has no automatic retry. The scientific command is:
```sh
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 \
  python verification/oblate_stationary_seed_20261127/check_stationary_system.py \
  --output-dir reports/oblate-stationary-system-recheck-2026-11-27/outputs
```

A success would validate and export the stationary differential expressions. It would still not solve the radial DAE, certify an admissible non-gauge solution germ, evaluate $B_\xi$ on such a germ, establish a canonical endpoint domain, or complete the nonlinear Dirac/degree-of-freedom analysis.
