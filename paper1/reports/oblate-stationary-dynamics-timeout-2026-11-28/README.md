# Stationary reduction: preservation and original rows passed, dynamics incomplete

Run `5d9b96ad-c712-487d-bc28-678e1a026243` from source `8c693bc2ede0783ada2b872fbd8b0d2b6f0f23ef` reached the 900-second wall budget during `terminal_reduced_dynamics`. The child exited 143; the supervisor recorded 900.223 seconds, no supervisor, relay or cleanup errors, and an empty final process group. A later process-table check also found no members. The 3074 process-group RSS samples reached 522,436,608 bytes (498.234 MiB); this is a sampled maximum, not an instantaneous peak.

## Completed checks

The invocation passed six exact row-operation controls, five kernel controls, eight new terminal-preservation control groups, and two manufactured DAE controls before fixed-field conversion and actual reduction. The radial-gauge identity and rank-four principal block passed again. All three elimination checkpoints are byte-identical to the preceding run.

The terminal kernel completed in 115.051 seconds. The new write-only snapshot was durably published in 2.321 seconds. Exact terminal preservation then passed in 75.004 seconds: DC took 0.235 seconds, CK 5.410 seconds, and (DC+CK)M 62.931 seconds; CB and both exact zero checks also completed. These are nested durations and must not be summed with their enclosing stages.

The unchanged original retained-row stage completed in 317.070 seconds, ending at child time 878.085 seconds. The frozen source requires both the state residual on M and the remaining-control residual to vanish in all 22 original equations, including the X0 equation. The separate discarded response row is not included in that zero requirement. The next stage, reduced dynamics, began but did not finish before the wall stop.

## Saved intermediate state

`outputs/terminal-state-before-preservation.json` contains exact C (6 by 14), K (14 by 14), B (14 by 0), M (14 by 8), free-row indices, fourteen explicit profile/derivative labels, and matching source/input/attempt identity. Its SHA-256 is `f3f899f32ef77b3419fe1c777198db5db5261c3fd6a1e9f3617335f372c0e157`. It was written **before** preservation; its false preservation and original-row flags correctly describe that publication moment and remain unchanged. Later completion of those two checks is evidenced by the raw progress log, not by retrospectively editing the snapshot. The file remains intermediate evidence, with no loader or resume path.

No `dae-closure.json`, `summary.json`, or independent germ certificate was produced. Reduced-dynamics construction, its full reconstruction checks, the discarded response, and interval regularity remain unfinished. Passing the preservation and original-row checks does not supply a responding solution germ, a nonlinear nonintegrability conclusion, a physical radial domain, or a Dirac/degree-of-freedom count.

## Provenance and next step

All 13 supervisor source files match their before/after hashes and Git bytes at the executed commit. All three upstream input digests match the completed stationary export. The archive includes the exact child bundle, engine log and framework exit/wrapper records; no scientific calculation was rerun during archiving. `timing-summary.json` labels nested timing data, and `archive-manifest.json` authenticates the archived files.

The measured 317-second original-row stage and unfinished dynamics motivate a narrow exact-field repair of those terminal products and full reconstruction checks. Discovery, kernel, preservation, original equations and the response definition are retained. Every formal terminal check and subsequent independent interval/germ certificate remains required.

Independent review `msg_40088be0` reconciled all archived payloads and the durable mirror, source/input identity, snapshot metadata, completed-stage boundary and resource records. It found no material issue. Watcher flag `f6cfb92db` was settled by retaining this confirmed incomplete record. No scientific code was imported or rerun by the reviewer.
