# Oblate angular weak-slot snapshot diagnostics

Status: **the first fixed-radius angular snapshot passed; five of the six planned snapshots remain pending**. Two earlier timeouts and the successful twelve-pair arithmetic diagnostic are preserved unchanged.
This child follows the completed [bounded local weak-form validation](../oblate-weak-curvature-validation/full-recheck-pass/README.md).
It belongs to the original Paper I metric-coefficient Einstein–Hilbert Galerkin regulator. It does not test the distinct reduction-before-cutoff Gaunt candidate.

The question is whether the validated local densities can be integrated into consistently labelled angular matrices at a small, declared set of radii. This is an interface diagnostic, not a finite radial operator calculation.

The [first attempt](first-attempt/README.md) preserves run `b092637e-9aa2-406e-86a8-963593e30f47` unchanged. Its source hashes and cleanup were verified. The [narrow follow-up preparation](direct-profile-probe-preparation.md) targeted only `m=0,r=3/2` under the same resource and numerical limits. Its [preserved attempt](profile-pair-attempt/README.md) also hit the wall limit: profile jets built in 0.522 seconds, and at least 120 pair integrals in the first 144-term potential-density family block completed before row 10 was interrupted. No block reconstruction or snapshot passed. The subsequent [arithmetic comparison passed](arithmetic-pass/README.md): 559 exact component equalities, twelve exact raw-versus-compact density identities and all quadrature gates on row 10. Recorded symbolic normalization/integration work was substantially smaller after component compaction, with the one-process timing caveat retained. The subsequent [compact first-snapshot test passed](compact-profile-pass/README.md) at m=0,r=3/2: all 2,076 component identities, all 16 profile-family reconstruction comparisons, strict angular regularity, Hermiticity and refinement gates passed. The saved matrices cover all twelve spatial and three lapse labels. The six-snapshot target below remains incomplete; m=0,r=2 is the next same-method check.

## Fixed scientific scope

- Cutoff `L=2`; all retained labels in complex sectors `m=0,1,2`.
- Radii `r=3/2,2`, each with `epsilon=1/4`: six angular snapshots.
- Spatial column counts `12,10,6`; lapse row counts `3,2,1`, respectively.
- All four derivative slots `00,01,10,11` for the potential and rectangular lapse bilinear.
- Separate radial flux matrices in the same four-slot layout. Left-derivative flux slots are evaluated and checked; they are not folded into the bulk density.
- Conjugate-left angular measure `dx/2` after equal-`m` azimuthal averaging. The prepared densities already contain `Sigma` once; there is no extra `r^2` or `2*pi` multiplier.

For radial coefficient vectors `v,u`, slot `ab` means derivative order `a` on the conjugated left profile and `b` on the right profile. Thus the later formal weak expression would be

```
v'* P u' + v* Q u' + v'* R u + v* S u,
P = B11, Q = B01, R = B10, S = B00.
```

No symmetry mirroring supplies uncomputed entries. The potential bulk checks are `B00=B00^dagger`, `B11=B11^dagger`, and `B01=B10^dagger`. The complete weak bulk form is already Hermitian; its Hermiticity does not require subtracting the radial Green flux. Rectangular lapse blocks have no analogous Hermiticity requirement.

## Evaluation and controls

Tensor jets are differentiated with symbolic `r` before substituting a radius, retaining the radial dependence of the angular tensor basis. The exporter then substitutes exact rational `r,epsilon` before constructing and cancelling paired integrands. At these fixed radii, accepted integrands must reduce to rational functions of `x`; unsupported radical remnants are rejected. Cancelled real and imaginary denominators are checked exactly for endpoint or interior zeros before numerical quadrature. This regularity check concerns the exported scalar angular integrands, not a physical radial endpoint domain.

The quadrature is float64 Gauss–Legendre at orders `24,48,96`. The acceptance comparison is componentwise `48` versus `96`, scaled by `max(1,abs(left),abs(right))`, with threshold `2e-10`; the `24` versus `48` difference is also retained. Refinement differences are numerical indicators, not certified quadrature error bounds. Nonfinite values fail.

Reconstruction uses the same `2e-10` scaled threshold. Potential Hermiticity uses absolute tolerance `2e-9`, actual-field realification uses `1e-9`, structural flux zeros and `m=0` imaginary leakage use `1e-11`, and the manufactured realification anchor uses `1e-12`. The exact round anchor, Legendre symbolic norms, profile-span determinants and denominator-root tests use exact arithmetic. These thresholds are fixed before execution.

The control surface includes:

1. Exact associated-Legendre normalization under `dx/2`, using `(ell+m)!/((2*ell+1)*(ell-m)!)`, checked through the actual quadrature routine at all three orders, plus an exact round radial weak-density/flux anchor.
2. Fresh polynomial-profile tensor jets, comparing direct integrated density/flux pair sums for the declared field combinations with the raw exported slot contractions. The accepted m=0 first snapshot compacts each differentiated, fixed-radius component only after proving its exact equality to the raw component. Two profile families on each side have an exactly nonzero two-slot jet determinant at each selected radius; all four family combinations are checked. These are assembly controls; their recorded combinations do not constitute a theorem for arbitrary profiles.
3. Positive-sine realification of square and rectangular matrices, checked against direct real-field azimuthal averaging and a manufactured imaginary rectangular block. The convention is `[[Re(H),Im(H)],[-Im(H),Re(H)]]` for `m>0`, corresponding to complex coefficients `a-i*b` for real cosine/sine coefficients `(a,b)`. The actual-field control transforms the already differentiated complex jets and lapse gradients, including the azimuthal derivative, then contracts real fields before angular averaging. A wrong-sign control must fail.
4. Independent potential Hermiticity, exact angular regularity, and quadrature refinement checks for each snapshot.

The exact round anchor has passed independent analytic/source review. For `epsilon=0`, `m=0`, the `ell=0` radial tensor `h=F(r) dr^2`, right tensor `k=G(r) dr^2`, and constant lapse harmonic give the identities

```
potential density = -2*F*G - r*(F*G' + F'*G)
potential radial flux = r*F*G
lapse density = -r*f'*G
lapse radial flux = r*f*G.
```

All formulas here use the prepared action normalization. The source review and executable exact anchor, rather than this prose alone, determine acceptance.

The first control batch also exercises the same rationality/regularity guard with two positive expressions and eight typed rejection cases, including endpoint/interior poles, unsupported functions and coefficients, accidental floating coefficients and nonfinite values.

## Resource and provenance plan

One local scientific process is approved, with a **300-second wall cap** and **2-GiB sampled process-group RSS threshold**. RSS sampling is not an OS hard memory limit. Runtime for this new workload is unknown; the parent’s local weak-validation timing is not a performance estimate for angular integration.

The inherited fixed launch command continues through the existing supervisor. This child's wrapper selects the new exporter and distinct output directory. Its cleanup, source-hash recording, isolated Python 3.12 environment, NumPy 2.3.5 and SymPy 1.14.0 remain pinned. The physics modules and prior scientific verifier are unchanged.

Historical preparation and failed-attempt records retain their original approval and unrun/partial status. The first attempt used the user’s “run it” instruction, the row arithmetic probe used explicit approval msg_37a2dc0f, and the first compact snapshot used the directly read standing instruction “you dont need to ask me for approval till 10am” (msg_908bb84e). The current local cutoff is November 28, 2026 at 10:00 America/New_York, equivalent to 15:00 UTC. Execution records preserve the literal quote because that ref can render as the following assistant reply in the interface. Physics modules and numerical tolerances remain unchanged.

Full six-case success is `angular_export_passed`; the accepted one-case run instead emitted `compact_profile_probe_passed`. Both are emitted only after their distinct completed result JSON is written. Partial case files, stage timings and resource records do not count as a completed diagnostic. The result records the exact labels, radii, slot convention, angular measure, quadrature orders, control coverage and source hashes.

## Completed preparation checks

Independent source reviews `msg_3412be83` and `msg_24db192f` found no remaining mathematical or provenance blocker for this bounded diagnostic. The exporter and supervisor pass AST parsing and in-memory bytecode compilation. All nine stdlib-only mocked supervisor tests pass, including the unapproved-launch guard; no scientific module was imported or executed by these tests. The plan JSON, view XML and whitespace checks pass. See [design review](design-review.json) and [preflight record](preflight.json).

Those historical preparation records are source-only. The later accepted [first snapshot archive](compact-profile-pass/README.md) supplies its own completed runtime evidence: 252.407 seconds to result construction, 253.311 seconds supervised, 220.28125 MiB peak sampled RSS, and maximum scaled 48/96 refinement difference 5.1181281435219717e-14. This is a numerical refinement indicator, not a certified error bound; 12 of 16 block-level absolute-error ratio flags were nonmonotone at small floating-point differences and are retained as metadata. The remaining five snapshots have no accepted result.

## Limits and later decisions

Two fixed-radius matrices are diagnostic snapshots, not coefficient functions over a radial interval. This work does not select radial basis functions, perform radial projection or integration, choose an endpoint domain, establish finite-basis kernel preservation, invert a constraint bracket, assemble finite `K,E,D`, compute rank or spectrum, identify an admissible seed, or establish nonlinear consistency. No change to either paper's action or regulator is proposed.

Later radial work must explicitly choose the projection order and test-space inclusions described in [canonical assembly](../oblate-quadratic-preparation/canonical-assembly.md) and [boundary/domain preparation](../oblate-quadratic-preparation/boundary-domain.md). A passing angular snapshot diagnostic alone will not settle those choices.

Reviewed work may be backed up to GitHub under the user’s explicit clarification; manuscript release/submission remains held pending a supported joint story. These original-regulator checks do not establish equivalence with the reduction-first Gaunt candidate.
