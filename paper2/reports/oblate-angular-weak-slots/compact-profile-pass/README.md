# First compact angular snapshot: accepted local diagnostic

Run `47d66179-da34-4f4e-809e-2eb2d7f19c96` from source `7d4b90ced87473ae3124878a75a74cef012e0562` passed `compact-profile-probe`. Independent raw-result review `msg_c0ce11d6` and source/result interpretation `msg_12d4d048` agree. This is one fixed-radius angular snapshot for the original metric-coefficient Galerkin regulator: **L=2, m=0, r=3/2, epsilon=1/4**. All remaining five planned snapshots are still pending. This local checkpoint is held for later integration. Manuscript release/submission is on hold; the user permits reviewed GitHub data backups, coordinated separately.

## Evidence and scope

The output contains all four potential and lapse derivative-slot matrices plus separate radial flux slots: 12 spatial labels, 3 lapse labels, 1,440 scalar matrix entries in the four families at order 96. The quadrature ladder is 24/48/96 with angular measure dx/2, Sigma included once, and slot mapping P=11, Q=01, R=10, S=00. Each final real/imaginary scalar integrand passed exact rational numerator/denominator guards and the endpoint/interior denominator-root check before quadrature.

Actual profile-weighted tensor fields were differentiated before exact radius/epsilon substitution, then compacted componentwise. All **2,064 spatial plus 12 lapse equalities** passed exactly before density/flux use. The direct oracle retained all four profile-family combinations for each of potential density, potential flux, lapse density and lapse flux. Its 576+576+144+144 pair integrals completed. All 16 comparisons against the independently exported raw slot matrices passed; the largest scaled reconstruction residual was 1.9806674970262414e-16 under the unchanged 2e-10 gate. Potential bulk Hermiticity and left-derivative flux structural-zero residuals were zero in the exported arrays. Profile-span determinants were nonzero. The m=0 actual real-field oracle is explicitly unperformed because this sector is already real; no m>0 field oracle or radical compaction was tested.

The maximum 48/96 scaled refinement indicator was **5.1181281435219717e-14**, under the unchanged 2e-10 threshold. Several 24/48-to-48/96 absolute-error ratios were nonmonotone at small floating-point differences; these are retained as metadata and were never an acceptance gate. They do not support a rigorous quadrature error bound or monotone-convergence theorem.

Three control batches passed: Legendre norms with strict final-scalar and component-compaction controls, the exact round weak-density/flux anchor, and manufactured realification with wrong-sign rejection. The final result says status=passed, and the stream reaches setup1/1, controls3/3, angular_export1/1, export1/1 with compact_profile_probe_passed=1. The early notification's zero-snapshot counters were stale. The checkpoint file is a pre-export running checkpoint and is preserved unchanged; the completed result and final execution record establish termination.

## Time, resources and provenance

- Result-construction elapsed: 252.407091332832 seconds.
- Final export progress elapsed: 252.414857874857 seconds (recorded slightly later).
- Supervised elapsed: 253.311114584096 seconds, under the 300-second cap.
- Peak sampled process-group RSS: 230981632 bytes = 220.28125000 MiB, across 867 samples; the 2-GiB threshold is sampled, not an OS hard limit.
- Raw density-slot integration: 118.823275958 seconds; raw flux-slot integration: 12.975257334 seconds.
- Actual-profile oracle: 119.699899750 seconds, including 5.345077166 seconds for compaction/equalities.
- Exit zero, no stop reason, final group 19378 empty; a later independent process-table check also found no members.

The six executed sources are copied from the isolated run and authenticated against both the execution hashes and committed Git blobs. The physics modules remain unchanged. Runtime records and the complete 58,855-byte raw log are copied verbatim. This archive and its durable mirror are byte-identical; manifest.json hashes every payload and intentionally excludes itself. No launch-environment or credential record is included.

The supervisor preserves the actual user standing quote “you dont need to ask me for approval till 10am”, interpreted in America/New_York as November 28, 2026 at 15:00 UTC. The message ref msg_908bb84e can expand in the interface to the following assistant launch reply; the literal user quote and previously read author context supply the authority, not that expanded assistant text.

## Remaining obligations

This establishes the first angular snapshot and its stated controls only. It does not establish the remaining five snapshots, arbitrary radial coefficient functions, a radial operator or projection, a common physical endpoint domain, finite K/E/D assembly, finite-basis kernel preservation, rank or spectrum, an admissible seed, or nonlinear consistency. It does not identify this original regulator with the reduction-first Gaunt candidate. The next proposed check changes only the fixed m=0 radius to r=2 under the same numerical/resource gates.
