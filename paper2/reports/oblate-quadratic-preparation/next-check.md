# Proposed first operator check — not launched

Preparation approval: `msg_7c4562eb`. This proposal requests a **separate** approval for one local invocation. No remote resources or tracked experiment have been launched by this preparation.

## Question

Do the general-r coefficient tangent and shift maps, canonical Legendre factors, and real harmonic conversion agree with direct geometric definitions and the pinned Paper I specialization? This should be resolved before angular/radial matrix assembly or a physical-spectrum calculation.

## Fixed input and command

Input physics is Paper I commit `157ec0f04866ea36509b325ec65f0a68890f9b17`. The unmodified angular reference is copied with its source hash into `verification/oblate_quadratic_preparation/reference/all_modes_157ec0f.py`. Its scientific verifiers are not called by this check. The new preparation commit must be recorded before launch.

Run from `paper2/` after approval:

```bash
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 \
  python verification/oblate_quadratic_preparation/verify_kinematics.py \
  --output-dir reports/oblate-quadratic-preparation/first-kinematic-check
```

Capture the complete stdout/stderr and exit status alongside the output JSON. Preserve any failed invocation. Do not fix and retry after a failure without the user's decision.

## Coverage and controls

- Use `N=3` (`L=2`) and all three nonnegative complex magnetic sectors. `r` and `epsilon` remain symbolic; this is not a fixed-radius mode calculation.
- Differentiate the actual nonlinear spatial metric map and compare every retained tangent column.
- Evaluate the coordinate definition of the Lie derivative with an arbitrary radial profile. Compare with `U0 s + T J1 s'`, including the radial derivative terms.
- At `r=1`, compare the resulting tensors and labels with the pinned Paper I source. This is a source-consistency control, separate from the direct geometric control.
- Exercise the DeWitt density assembly at the rational deformed point `r=2, x=1/3, epsilon=1/5` against a full-matrix inverse/determinant contraction. Exactly integrate the round `L=2,m=2,r=2` sector and compare with the pinned round pairing after the frame/volume conversion.
- In that integrated sector, include the shift-profile derivative columns before Schur projection and require their normal-Hessian rows and columns to vanish. This checks the actual density path in addition to the pointwise tensor tangency identity.
- Verify the Legendre factors with a manufactured indefinite invertible spatial Gram. Deliberately omitting the shift norm must fail that identity.
- Verify the positive-sine realification by direct angular integration of a synthetic imaginary cross block.
- Require exact symbolic cancellation on the stated coordinate chart; use no floating-point tolerance to declare these identities zero.

The script produces source hashes, environment versions, elapsed time, and check counts only after successful completion. It does not fabricate partial success summaries on an exception. A pass does not validate the spatial-potential implementation, establish the common radial domain, or compute nonlinear obstruction values.

## Resources

One local Python process, CPU only, with no paid external API. There is no measured runtime for this new checker. Planning estimate: seconds to a few minutes and under approximately 1 GB RAM for the symbolic `N=3` maps; this is an estimate, not a benchmark. If the process exceeds a three-minute inspection threshold, inspect progress/resource use and report before considering expansion. No automatic retry or hardware escalation is part of the request.

## Subsequent gates

After this check is reviewed, propose a separate direct curvature-variation/weak-form test. Only after potential coefficients, boundary traces, canonical normalization, and radial projection order are explicit should a small radial operator assembly be proposed. The exact-rational reducer cannot consume rationalized floating quadrature as proof of a structural rank.
