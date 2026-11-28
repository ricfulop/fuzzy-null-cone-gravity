# Stationary closure stopped during terminal preservation

One fresh local CPU invocation from source `e48324b06a7edc16bd2bff2d7dd16125762fd96c`, run `157032da-1468-473e-8cd7-37763ff86822`, stopped at its 900-second wall budget. The supervisor recorded **900.093 seconds**, child exit 143, empty final process group, no cleanup or relay error, and no scientific assertion traceback. The largest of 3,064 process-group RSS samples was **443,760,640 bytes (423.2 MiB)**. The 2-GiB setting was a sampled threshold; the instantaneous peak was not measured.

The new kernel method completed successfully inside the full reducer in **114.774936 seconds**, ending at child elapsed time 474.315756 seconds. Its original C*M exact-field check took 88.320807 seconds, and the C-only fraction-free RREF took 23.658989 seconds. All six kernel substages completed. The following `terminal_constraint_preservation` stage began but emitted no completion before the wall stop.

That stage combines differentiation DC, multiplication CK and simplification of their sum. The log does not separate those costs. Multiplication by M belongs to the next verification stage and was not reached. This is a resource-limited incomplete calculation, not a failed preservation identity or a no-go result.

## Completed and missing work

The six row-operation control groups, five new terminal-kernel control groups, and two manufactured DAE controls passed before fixed-field conversion and reduction of the original equations. Fixed-field conversion, the full radial gauge identity, and principal rank four passed. Three elimination checkpoints completed; they are byte-identical to the corresponding checkpoints from the prior control-only attempt `8fa984ea-3c9b-4e85-b984-c23956fc95d8`. The final partial checkpoint records six state constraints and zero remaining acceleration controls.

The kernel was constructed and verified in memory but was not exported as a terminal dynamics checkpoint. The original 22 retained-row checks, reduced dynamics/reconstruction checks, discarded response and actual interval regularity remain uncompleted. There is no `dae-closure.json` or final `summary.json`. No solution germ or nonlinear compatibility conclusion is established.

## Preserved evidence

- [Complete engine log](engine.log), [child log](verification.log), [nested stage timings](timing-summary.json).
- [Execution record](execution.json), [launch record](launch.json), [all RSS samples](resource-samples.json).
- Source snapshots under `sources/`; partial exact matrices and controls under `outputs/`.
- [Archive status and independent process check](archive-status.json).

All eleven monitored source snapshots match their before/after records and the bytes at the executed commit. The independently pinned upstream completion digest is unchanged. The supervisor environment is Python 3.13.7; the child invocation records Python 3.12.11 and SymPy 1.14.0. A later process-table check found no member of the original process group. Archiving used file copying, log/JSON inspection and hashing only; no scientific code was rerun.

## Next bounded repair

The next preparation targets only terminal preservation arithmetic over QQ(r,T), retaining discovery equations and every downstream check. A write-only uncertified C/K/B/M snapshot is also warranted before expensive terminal verification; no checkpoint loader or automatic resume is proposed. Separate derivative, product and projection timings will distinguish the remaining costs. The prior failures and the successful isolated kernel diagnostic remain unchanged.

## Independent review

Review `msg_6c38802b` found no material archive, provenance or scope issue using file, JSON and log inspection plus byte comparisons and hashes. No scientific code was imported or rerun. The later watcher flag `f590a9e37` was settled by retaining this confirmed failed record; this does not endorse a closure result.
