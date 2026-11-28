# First weak-curvature validation attempt: incomplete

The first tracked attempt stopped during the symbolic pure-gauge control after exceeding its monitored memory budget. Cleanup then raised a separate permission error. This is not a failed gravitational identity or a completed validation.

- Run: `85aa6197-97a5-4a4a-9239-4119a350d7ef`.
- Source: `42200f8c135e9ae5837a25515ba33b3c4055e1fe`.
- Supervisor elapsed time: 667.526 seconds; framework exit code 1.
- Primary recorded stop: `sampled_group_rss_budget`, with peak sampled process-group RSS 2,150,842,368 bytes against 2,147,483,648 bytes allowed. This was sampled enforcement, not an OS hard memory cap.
- Completed: the three representative L=2, m=0,1,2 weak-sector batches. Their exact sampled point, labels, profiles and counts remain in `weak-form-progress.jsonl`.
- Incomplete: m=0 symbolic pure-gauge control, subsequent gauge controls, Cartesian TT normalization, ablation, and final result export.

The complete log contains no scientific assertion failure. It does show `PermissionError` while checking group existence after SIGTERM and again during exception cleanup. The scientific child return code was not recovered and remains unknown. A later process-table inspection found no original PID or process-group members; that observation and its time are in `failure-summary.json`. The later observation does not reconstruct the child's exact exit status.

All six source hashes match the frozen run snapshot. The durable archive `artifacts/paper2-oblate-weak-first-attempt-2026-11-28/` preserves the exact sources and the 2,271 RSS samples. The proposed next step is an algebra-preserving reduction of symbolic expression growth plus a supervisor repair that distinguishes an empty process group from unresolved cleanup, under the existing overnight approval. No mathematical test is to be weakened to obtain a pass.

This archive records one attempt. Read the raw [verification log](verification.log), [execution record](execution.json), [partial checkpoints](weak-form-progress.jsonl) and [failure summary](failure-summary.json). The code's sampled Green checks are not a global identity theorem, radial operator export, admissible-seed certificate, or nonlinear consistency proof.
