# First-case independent profile reconstruction probe

Status: independently source-reviewed and approved for one bounded probe; scientific execution not yet performed.

The user renewed authorization in `msg_c7ac8e6c`, in chat `chat_5351fc4a-9bc3-4992-9e39-0a9505c63137`: “ok I am going back to sleep, you have my approval for any of the upcoming experiments that come up on the chat or other requests to keep going till I wake up”. That message was read directly. This preparation uses it for one reviewed local diagnostic with a 300-second wall cap and 2-GiB sampled process-group RSS threshold. The existing 2026-11-28 16:29 UTC cutoff remains an additional limit. No remote/GPU budget is used.

## Evidence motivating the change

The [first-attempt archive](first-attempt/README.md) preserves run `b092637e-9aa2-406e-86a8-963593e30f47`, source `ffc5fcd0855932596e4d51fc6232aa754af2c2b4`. All three control batches completed, but zero angular snapshots completed. The first `m=0,r=3/2` case finished density and flux integration in 117.314 and 12.950 seconds, then stopped within `direct_profile_oracle` at the 300-second limit. The exact operation within that oracle was not timed. Peak sampled RSS was 1003.15625 MiB; shutdown and a later process-table check found the group empty.

The source-based hypothesis is that forming and expanding a single dense profile-pair expression causes costly symbolic cancellation. It is not a demonstrated inner hotspot and does not establish a mathematical failure.

## Algebra-preserving change

For each of the original two profile families, build every actual field `c_I F_I(r) T_I(r,x)` separately and apply `tensor_jet(field, data, 1, 0, gamma)` before substituting the fixed radius and oblate parameter. Thus the actual profile and basis radial derivatives remain in the direct path. The lapse analog independently constructs and differentiates each scalar profile times its harmonic.

For each of the four left/right family combinations, evaluate the prepared potential density and radial flux on every independent spatial column pair, and the lapse density and flux on every independent lapse/spatial pair. Integrate each contribution through the unchanged rationality, pole-certification and quadrature evaluator. Sum the real and imaginary numerical contributions separately with stable summation. The direct values must not read or be reconstructed from exported slot matrices.

This distribution uses sesquilinearity: for a conjugate-left form `B`, `B(sum_I h_I, sum_J k_J) = sum_IJ B(h_I,k_J)`. Coefficients occur once in the actual fields; the existing density/flux formulas conjugate the left field. No extra coefficient or conjugation is applied during accumulation. The independent direct sum is compared to the original slot reconstruction with the unchanged scaled tolerance. Floating sums remain numerical checks, not exact integration certificates.

The probe retains all three existing control batches, all four profile-family combinations, all retained labels of the first case, all four derivative slots, the separate flux slots, exact profile-span determinants, rational integrand and denominator guards, potential Hermiticity, and the existing acceptance tolerances. Full mode retains its actual real-field controls for the nonzero m sectors; this `m=0` probe does not claim those were exercised.

## Probe and reporting contract

- Fixed target: `L=2`, `m=0`, `r=3/2`, `epsilon=1/4`.
- The same pairwise reconstruction helper serves full mode and probe mode.
- The fixed experiment command remains unchanged; the supervisor passes `--mode direct-profile-probe` to the exporter.
- Plan: setup 1 environment; controls 3 batches; angular export 1 snapshot; export 1 result.
- Result: `direct-profile-probe-checks.json` under `direct-profile-probe-runtime`.
- Success: `direct_profile_probe_passed` only after writing the final JSON; no `angular_export_passed` in probe mode.
- Timings identify profile-jet construction and each family/block/row of pair integrations, with expected and completed term counts. They do not advance completed snapshots.
- The full six-case mode and its distinct result remain available. The earlier failed snapshot and source records remain immutable.

Runtime of the pairwise oracle is unmeasured. This is a bounded first-case test, not a promise that the complete six-case workload fits the same cap.

## Scope

A probe pass would establish only the first fixed-radius snapshot and its specified controls. It would not finish all six angular snapshots or perform radial projection, choose a physical endpoint domain, establish finite-basis kernel preservation, assemble finite K/E/D, compute spectra, construct an admissible seed, or prove nonlinear consistency. The original coefficient-Galerkin regulator remains distinct from the reduction-first Gaunt candidate.

## Static review and freeze

Independent reviews `msg_25b87cf7` and `msg_36bbf637` found no launch blocker in the paired direct oracle, conjugation, radial product rule, mode isolation, or final row-timing instrumentation. The failed-attempt archive was independently verified in `msg_da598ec4`. AST/in-memory compilation, JSON/view-fragment parsing, unchanged physics/shim/checker comparisons, whitespace checks, and all nine stdlib-only mocked supervisor tests passed. These tests did not execute scientific modules. See [probe preflight](direct-profile-probe-preflight.json) for the reviewed source hashes.
