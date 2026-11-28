# Complete weak-form recheck after curvature diagnostics

The five-batch curvature diagnostic `025b8511-9538-4bc7-a189-5c053156fb27` passed from source `1d3c606a8d066c2c850565616734ce9ea5c5132a`: 142.570241 supervisor seconds, 141.581383 checker seconds, 165.515625 MiB peak sampled process-group RSS, all five checkpoints and a clean exit. The raw log contains all 27 Einstein-entry and three scalar exact-zero completions for the three prepared symbolic gauge vectors. Its record is preserved in [curvature-recheck-pass/README.md](curvature-recheck-pass/README.md). This pass is distinct from complete weak-form validation.

The next invocation changes only the supervisor's mode to `full`. The verifier remains byte-identical with SHA256 `3bb2a49863cda8b6f5f3b7a1a718ec580cc6916ed6f5cde13d682d2571ac7113`. The physics modules, strict chart arithmetic, exact-zero predicate and dense Ricci oracle are unchanged. Two source/result reviewers found no blocker in the full-mode wiring or the preceding diagnostic record (`msg_8a6fddec`, `msg_5cad2b38`).

One tracked local CPU invocation keeps the 300-second wall cap and sampled 2-GiB process-group RSS budget, isolated Python 3.12 / SymPy 1.14.0 / NumPy 2.3.5, assertions enabled, provenance and cleanup. The expected workload is comparable to the measured diagnostic plus the sampled weak sectors; the complete current-mode runtime is unmeasured, so 300 seconds is a cap, not a promised duration. The user's standing overnight authorization ends at 2026-11-28 16:29 UTC.

The runtime plan reports setup (one environment), validation (eight batches), then export (one checker result). The eight batches are:

1. Representative L=2 weak-density, potential/lapse Green and radial endpoint-flux checks at m=0.
2. The corresponding representative m=1 checks.
3. The corresponding representative m=2 checks.
4. The prepared symbolic pure-gauge test at m=0.
5. The prepared symbolic pure-gauge test at m=1.
6. The prepared symbolic pure-gauge test at m=2.
7. Cartesian transverse-traceless sign/normalization and wrong-sign rejection.
8. Nonzero rejection of omitted cross/trace terms.

The sampled sectors use the raw shared-Ricci lane and recorded nonconstant profiles; the gauge tests use the independently varied-connection path with combined normalization. The TT and ablation controls remain enabled. The result file is `weak-form-checks.json`; only successful writing of that file permits `validation_passed=1`. Diagnostic filenames and success signals are not used for this invocation.

A pass would validate the specified local derivative-slot densities and flux identities with the stated sampled and symbolic controls. It would not prove the identities for all retained labels and all chart points, assemble finite radial K/E/D, choose physical endpoint spaces, certify preserved finite-basis kernels, construct an admissible seed, compute a spectrum, count physical degrees of freedom, or establish nonlinear gravitational consistency. The original Paper I coefficient-Galerkin action remains distinct from the reduction-before-cutoff Gaunt candidate.
