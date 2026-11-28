# Complete curvature diagnostic after single-entry pass

The three-batch single-entry diagnostic `9b6ab336-ed3c-460b-9308-b2dc21c470e4` passed from source `3445f140ac8352e7e83a07edefb9ecdcaad83e73`. Its full log, JSON, source hashes, all 125 resource samples and clean process disposition are preserved in [ricci-entry-pass/README.md](ricci-entry-pass/README.md). The observed whole invocation lasted 36.748044 seconds with 156.1875 MiB peak sampled group RSS. The single global gauge Ricci entry (1,2) normalized in 2.029145 seconds and passed the original exact-zero test. This supports returning to the existing complete curvature diagnostic, without changing the checker again.

Only the supervisor's selected mode changes from `ricci-entry-probe` to `curvature-probe`. The verifier remains byte-identical with SHA256 `3bb2a49863cda8b6f5f3b7a1a718ec580cc6916ed6f5cde13d682d2571ac7113`. The tested weak forms, background/tensor map, chart normalizer, exact-zero predicate and raw dense Ricci oracle remain unchanged.

One local CPU invocation is bounded by the same 300-second wall cap and 2-GiB sampled process-group RSS budget, with the same isolated environment and cleanup. The complete workload is not yet timed; the cap is not a promised completion time. The user's overnight authorization ends at 2026-11-28 16:29 UTC.

The five diagnostic batches are:

1. The existing chart arithmetic, exact nonzero and unsupported-expression controls.
2. All nine raw and normalized contracted Ricci entries compared with the dense raw oracle at the recorded nonzero offshell m=1 point.
3. All nine Einstein components and the scalar curvature for the original global symbolic gauge perturbation at m=1.
4. The same complete gauge test at m=0.
5. The same complete gauge test at m=2.

The diagnostic writes `curvature-probe-checks.json` and emits `diagnostic_passed` only after that export. It does not run the sampled weak potential/lapse Green identities, radial endpoint flux, TT or weak-density ablation batches. A pass therefore supports a later decision on the eight-batch full weak-form validator; it is not that validator's pass and establishes no finite radial operators, physical boundary domain, spectrum, admissible seed, or nonlinear consistency theorem. Historical failed attempts and the single-entry pass remain distinct.
