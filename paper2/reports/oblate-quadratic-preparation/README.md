# Oblate operator and boundary preparation

Prepared under user approval `msg_7c4562eb`. Physics is pinned to the original Paper I metric-coefficient Galerkin action at `157ec0f04866ea36509b325ec65f0a68890f9b17`. The action has not been replaced by an alternative canonical or reduction-before-cutoff prescription.

## Deliverables

- [Canonical assembly](canonical-assembly.md): the actual coefficient Legendre transform, full formal `K,E,D`, general-r shift derivative injections, real harmonic convention, and distinct radial projection orders.
- [Spatial-curvature Hessian](potential-hessian.md): the complete covariant flat-background operator, first-derivative weak form, connection coefficients, and Green flux.
- [Boundary/domain analysis](boundary-domain.md): compatible-domain requirements, endpoint obstructions, canonical dual bases, and clearly labelled shell-diagnostic options.
- [Machine-readable interface](interface.json): field/normalization metadata, formal blocks, unresolved choices, and the seed handoff to Paper I's nonlinear analysis.
- [First check proposal](next-check.md): the original preparation-stage command, fixed coverage, source provenance, and local resource estimate.
- [First kinematic check](first-kinematic-check/README.md): the initial invocation failed at an exact-zero reference comparison; its original log and source provenance are preserved unchanged.
- [Exact-zero correction](zero-comparison-correction.md): the later approved predicate repair and zero/nonzero controls, without changing the physics implementation.
- [Corrected kinematic check](kinematic-recheck/README.md): the single approved rerun passed all prescribed sectors and controls.
- [Bounded weak-form validation](../oblate-weak-curvature-validation/full-recheck-pass/README.md): the full eight-batch recheck passed representative sampled potential/lapse Green identities, three fixed symbolic gauge controls, TT normalization, and negative controls. Radial-domain and finite-assembly questions remain open.

Implementation preparation is in `paper2/verification/oblate_quadratic_preparation/`:

- `kinematic_blocks.py`: retained tensor columns, general-r shift split, DeWitt densities, post-angular-integration canonical blocks, and realification.
- `weak_forms.py`: potential and lapse weak derivative-slot densities and separate endpoint fluxes.
- `potential_hessian.py`: formula/connection definitions and a coefficient-export container.
- `verify_kinematics.py`: the prepared first local check; it does not call the archived scientific verifiers.
- `verify_weak_forms.py`: independent varied-connection curvature, representative weak/Green comparisons, symbolic gauge controls, and TT/ablation checks; diagnostic and full modes have separate outputs.
- `reference/all_modes_157ec0f.py`: an unmodified pinned Paper I reference with SHA256 `55ba3d4ba6c44104668362258aa96465bd485fef8f407f65283e12aa0c643ebe`.

## Result status

The derivations define how to assemble the missing stationary operator. The preparation at `ee7d0a3` was unrun. The first separately authorized local check passed `m=0` and failed at the `m=1` reference equality predicate; its original log remains intact. The later predicate correction at `2986a52` passed one approved rerun of all `m=0,1,2` sectors, selected density/integration checks, and the normalization/realification controls. This completed the bounded kinematic checkpoint.

The subsequent full weak-form run `32b93879-a4f2-426f-82ac-2fa455844477`, from source `b59be9535b38e25d59a8437efe9873094b2bc4bd`, passed all eight batches. It covers 108 potential-pair and 36 lapse-pair Green comparisons exactly at recorded rational points, separate radial flux comparisons, three fixed symbolic gauge controls, Cartesian TT normalization, and wrong-sign/omitted-term rejection. The [result archive](../oblate-weak-curvature-validation/full-recheck-pass/README.md) distinguishes sampled Green evidence from symbolic controls and preserves the earlier incomplete attempts separately. No finite radial `K,E,D` matrices, physical spectrum, on-shell seed, nonlinear compatibility value, or degree-of-freedom count have been produced by this preparation.

Two preparation reviews reconciled the potential normalization, realification convention, and adjoint/boundary distinction. The review also exposed and corrected a missing one-half in the draft potential cross terms and errors in a draft finite momentum-conversion description before scientific execution. Source review and the later bounded measured checks are separate evidence.

## What the parallel tracks need from one another

Paper I retains the second-order nonlinear-preservation derivation and the full continuum linear-Euler residual evaluator. This package supplies the operator/boundary interface needed to construct a retained on-shell linear seed. The seed must include the reconstructed spacetime perturbation and sufficient time/radial jets; discarded temporal residuals cannot be inferred from retained Hamiltonian matrices alone. Pure continuum-gauge fields are useful controls but cannot expose a nonzero compatibility pairing when their full linear Einstein residual vanishes.

The two operator preparation tasks can proceed analytically in parallel, but finite assembly needs one explicit radial projection order and one boundary/test-space choice. A common continuum physical evolution domain has not been proved. Equal radial bubble spaces can break exact derivative/product inclusions; those approximation defects must be measured separately from a failure of the angular theory.

## Decision after preparation

The kinematic and weak-form checkpoints are complete at their recorded bounded scopes; all earlier failed or incomplete invocations remain preserved. The next preparation is a raw four-slot angular coefficient/flux interface, followed by finite-basis derivative/product inclusion diagnostics with an explicitly declared radial projection and test space. A compact-support or homogeneous-shell diagnostic does not choose physical endpoint conditions. Physical boundary choices, finite radial assembly, numerical rank policy, and a genuinely on-shell seed remain subsequent decisions; the present result supplies none of those by implication.
