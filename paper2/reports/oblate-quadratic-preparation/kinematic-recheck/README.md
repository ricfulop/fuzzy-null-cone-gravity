# Corrected oblate kinematic checker: passed

The single rerun authorized by `msg_da46f19d` passed with exit **0**. The input was correction commit `2986a52960840cc846f88a0665acbc44e4c6bf6e`, following the preserved failed attempt at `5c92074dd6cd230156c04a6481a6a83a9cc6e934`. Only the verifier's exact-zero normalization and regression controls changed; the action, kinematic implementation, and pinned Paper I reference did not.

## Completed checks

The [full log](verification.log) and [structured result](kinematic-checks.json) record:

- Three exact-zero regressions accepted, including the original failed residual; five deliberately nonzero expressions rejected, including an exact `10^-30` offset and a square-root branch control.
- All `L=2` magnetic sectors `m=0,1,2` passed direct geometric and pinned-reference comparisons: **28 complex spatial tangent columns**, **16 complex shift columns**, and **three reference blocks**. These are implementation-check counts, not physical degree-of-freedom counts.
- The actual DeWitt density path passed the deformed rational-point comparison against an independent inverse/determinant contraction. The round `L=2,m=2,r=2` sector passed exact angular integration against the reference pairing.
- Shift-profile derivative columns passed through that integrated density path and had zero normal-Hessian derivative rows and columns after Schur projection.
- The manufactured indefinite-Gram Legendre check passed, and omission of the shift-norm term failed as intended. The synthetic imaginary cross block passed positive-sine realification.

The system timing wrapper records **3.18 seconds** including the runner and **81,969,152 bytes** maximum resident set size. The checker records **2.3439 seconds** internally, Python **3.12.11**, NumPy **2.3.5**, and SymPy **1.14.0**. See [execution provenance](execution.json) for the command, input hashes, approval, and timing boundaries.

## What the correction established

The two radicands in the original failure, `-(x-1)*(x+1)` and `1-x**2`, are now put into the same polynomial representation before exact cancellation. The regression passes without modifying any harmonic or tensor coefficient. Nonzero and branch controls prevent accepting an expression merely because it is small or because an invalid radical identity was assumed.

The independent pre-execution review (`msg_a0174d43`) found no blocker in the concrete predicate/control diff. The original failure and its provenance are unchanged. Exactly one full rerun was executed after the new approval; there was no separate scientific diagnostic run and no remote/GPU compute.

## Limits and next proposed gate

This closes the stated kinematic/reference and selected density/normalization checkpoint. It does **not** validate the potential or lapse weak forms against independent curvature variations, verify their Green fluxes, select a radial projection order or common boundary domain, export a full finite radial `K,E,D` system, construct an admissible non-gauge seed, or establish a physical spectrum or nonlinear constraint classification.

The next proposed scientific task is a separately approved local validation of the potential and lapse weak forms against independent curvature variations and Green identities, retaining boundary fluxes explicitly. No such calculation has been launched by this approval.
