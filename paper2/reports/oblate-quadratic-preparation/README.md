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
- [Corrected kinematic check](kinematic-recheck/README.md): the single approved rerun passed all prescribed sectors and controls; curvature/Green and radial-domain verification remain open.

Implementation preparation is in `paper2/verification/oblate_quadratic_preparation/`:

- `kinematic_blocks.py`: retained tensor columns, general-r shift split, DeWitt densities, post-angular-integration canonical blocks, and realification.
- `weak_forms.py`: potential and lapse weak derivative-slot densities and separate endpoint fluxes.
- `potential_hessian.py`: formula/connection definitions and a coefficient-export container.
- `verify_kinematics.py`: the prepared first local check; it does not call the archived scientific verifiers.
- `reference/all_modes_157ec0f.py`: an unmodified pinned Paper I reference with SHA256 `55ba3d4ba6c44104668362258aa96465bd485fef8f407f65283e12aa0c643ebe`.

## Result status

The derivations define how to assemble the missing stationary operator. The preparation at `ee7d0a3` was unrun. The first separately authorized local check passed `m=0` and failed at the `m=1` reference equality predicate; its original log remains intact. The later predicate correction at `2986a52` passed one approved rerun of all `m=0,1,2` sectors, selected density/integration checks, and the normalization/realification controls. This is a bounded kinematic verification. No finite radial `K,E,D` matrices, physical spectrum, on-shell seed, nonlinear compatibility value, or degree-of-freedom count have been produced.

Two preparation reviews reconciled the potential normalization, realification convention, and adjoint/boundary distinction. The review also exposed and corrected a missing one-half in the draft potential cross terms and errors in a draft finite momentum-conversion description before scientific execution. Review is not a substitute for the proposed independent calculations.

## What the parallel tracks need from one another

Paper I retains the second-order nonlinear-preservation derivation and the full continuum linear-Euler residual evaluator. This package supplies the operator/boundary interface needed to construct a retained on-shell linear seed. The seed must include the reconstructed spacetime perturbation and sufficient time/radial jets; discarded temporal residuals cannot be inferred from retained Hamiltonian matrices alone. Pure continuum-gauge fields are useful controls but cannot expose a nonzero compatibility pairing when their full linear Einstein residual vanishes.

The two operator preparation tasks can proceed analytically in parallel, but finite assembly needs one explicit radial projection order and one boundary/test-space choice. A common continuum physical evolution domain has not been proved. Equal radial bubble spaces can break exact derivative/product inclusions; those approximation defects must be measured separately from a failure of the angular theory.

## Decision after preparation

The first failed invocation and the successful approved rerun are preserved in separate directories. The proposed next scientific action is independent curvature/Green-identity validation of the potential and lapse forms, subject to a separate user decision. Physical boundary choices, finite radial assembly, numerical rank policy, and a genuinely on-shell seed remain subsequent decisions. No spectrum run, SSH work, or B200 allocation is authorized by the completed kinematic check.
