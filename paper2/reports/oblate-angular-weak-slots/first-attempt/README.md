# Angular weak-slot diagnostic: first timeout

Run `b092637e-9aa2-406e-86a8-963593e30f47` from source `ffc5fcd0855932596e4d51fc6232aa754af2c2b4` stopped at its 300-second wall limit. This is an incomplete angular diagnostic of the original coefficient-Galerkin regulator. It is not an angular-export pass or a failure of a mathematical identity.

## Verified outcome

The raw [verification log](verification.log), [execution record](execution.json), and [progress stream](angular-weak-slot-progress.jsonl) agree: all three control batches completed; zero of six snapshots completed. The preserved [checkpoint](angular-weak-slot-checkpoint.json) contains the Legendre/regularity controls, exact round weak anchor, and manufactured realification controls. Its historical `status: running` reflects when it was written; the execution record gives the terminal failure.

For the first case, m=0 and r=3/2, symbolic slot construction completed in 0.602 seconds, density cancellation/quadrature in 117.314 seconds, and flux cancellation/quadrature in 12.950 seconds. The last substage is `direct_profile_oracle start`. Available timers do not distinguish combined-field jet construction, scalar-density construction, normalization, or quadrature inside that oracle. Large symbolic expression expansion is a source-based hypothesis, not a measured inner hotspot.

The supervisor recorded 300.398998 seconds, child return code 143, stop reason `wall_budget`, and peak sampled group RSS 1,051,885,568 bytes (1003.15625 MiB) across 1,028 samples. A transient permission error occurred during termination; final disposition was empty. A later independent process-table check also found no group members. See [archive status](archive-status.json) and [resource samples](resource-samples.json).

No final result JSON or complete snapshot file was written. The reported zero refinement metric was an initialization value, not measured convergence. Earlier weak-form and curvature validation results retain their separate scope.

## Preserved evidence

The archive contains only the four named scientific runtime records, the raw scientific output log, the six source files named by the supervisor, and this report/status/manifest. The copied source bytes match all six execution hashes. No environment-launch record is included. Raw records are unedited. The manifest records SHA-256 for every payload except itself.

## Proposed next step

Keep all mathematical controls and tolerance gates. Add timings inside the independent combined-profile oracle and prepare a bounded first-case diagnostic. Consider algebra-preserving compaction of independently differentiated fields before bilinear expansion, after source review. No source repair or follow-up run is represented by this archive. A new successful first-case diagnostic would still not certify all six snapshots.

The fixed-radius angular scope excludes radial projection, physical endpoint-domain choice, finite K/E/D assembly, finite-basis kernel preservation, spectra, admissible seeds, and nonlinear consistency.
