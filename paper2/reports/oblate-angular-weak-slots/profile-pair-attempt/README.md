# Independent profile-pair probe: bounded timeout

Run `5e8d761b-ccb2-440a-8f9d-07098036d1af`, source `8eab4afd1a153331c2252669e3cfde2b6017b13f`, stopped at its 300-second wall limit. This second attempt used pairwise integration of independently profile-weighted tensor fields, restricted to `L=2,m=0,r=3/2,epsilon=1/4`. It produced no completed snapshot or final probe result.

## Verified progress

The [raw log](verification.log), [execution record](execution.json), [progress](angular-weak-slot-progress.jsonl), and [control checkpoint](angular-weak-slot-checkpoint.json) agree: three control batches completed, zero of one snapshot completed. The checkpoint retains its historical running status because it predates the terminal execution record.

Symbolic slots took 0.594622 seconds; density cancellation/quadrature took 121.717378 seconds and flux 13.769565 seconds. All six left/right/lapse profile-family builds completed in a combined 0.521677 seconds. The first potential-density block, family `(0,0)`, completed rows 0 through 9: at least 120 of its 144 pair integrals. Row 10 started and did not finish before termination. No direct pair block completed, so no reconstructed family value was accepted. Individual values from those partial terms were not saved as a result. See [parsed timing events](timing-events.json).

The row timer covers density construction plus evaluator substitution, expression expansion/cancellation, rationality/root checks and quadrature. It does not isolate which of these operations dominates. Symbolic normalization remains a source-based hypothesis; jet construction is demonstrably small in this attempt.

## Resources and integrity

Supervisor elapsed time was 300.389361 seconds, return code 143, stop reason `wall_budget`, peak sampled group RSS 286,277,632 bytes (273.015625 MiB), across 1,025 samples. This is lower sampled peak RSS than the first dense attempt (1003.15625 MiB), but the two runs stopped at different computational work and this is not a completed-runtime comparison.

Final cleanup reports an empty group after a transient permission-error probe; a later independent process-table check also found no members. All six source hashes match the executed snapshot. The [manifest](manifest.json) covers every archive payload except itself. Only named scientific runtime records, the raw scientific output, and the six explicitly pinned source files were copied; no environment-launch record is present.

## Interpretation and next decision

The pairwise representation reduces the measured memory peak and exposes progress, but does not yet complete the first snapshot under the approved limit. Initial refinement zero remains an initialization value, not a measured convergence result. No mathematical assertion failed. Both failed attempts remain preserved.

After two incomplete bounded attempts, another full-case retry is not justified by this record alone. The next proposal should isolate a small pair/row and measure density construction versus normalization, with exact raw-versus-compact algebra checks if componentwise rational compaction or a different exact representation is proposed. Keep the independent actual-profile path, all ultimate controls and tolerances, and existing resource limits. No follow-up scientific run is part of this archive.

This is an original metric-coefficient Galerkin diagnostic. It establishes no six-snapshot angular export, radial operator/projection, physical endpoint domain, finite K/E/D, finite-basis kernel preservation, spectrum, admissible seed, or nonlinear theorem.
