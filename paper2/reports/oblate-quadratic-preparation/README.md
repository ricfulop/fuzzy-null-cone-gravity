# Oblate operator and boundary preparation

Prepared under user approval `msg_7c4562eb`. Physics is pinned to the original Paper I metric-coefficient Galerkin action at `157ec0f04866ea36509b325ec65f0a68890f9b17`. The action has not been replaced by an alternative canonical or reduction-before-cutoff prescription.

## Deliverables

- [Canonical assembly](canonical-assembly.md): the actual coefficient Legendre transform, full formal `K,E,D`, general-r shift derivative injections, real harmonic convention, and distinct radial projection orders.
- [Spatial-curvature Hessian](potential-hessian.md): the complete covariant flat-background operator, first-derivative weak form, connection coefficients, and Green flux.
- [Boundary/domain analysis](boundary-domain.md): compatible-domain requirements, endpoint obstructions, canonical dual bases, and clearly labelled shell-diagnostic options.
- [Machine-readable interface](interface.json): field/normalization metadata, formal blocks, unresolved choices, and the seed handoff to Paper I's nonlinear analysis.
- [First check proposal](next-check.md): an exact command, fixed coverage, source provenance, and local resource estimate. **Not launched.**

Implementation preparation is in `paper2/verification/oblate_quadratic_preparation/`:

- `kinematic_blocks.py`: retained tensor columns, general-r shift split, DeWitt densities, post-angular-integration canonical blocks, and realification.
- `weak_forms.py`: potential and lapse weak derivative-slot densities and separate endpoint fluxes.
- `potential_hessian.py`: formula/connection definitions and a coefficient-export container.
- `verify_kinematics.py`: the prepared first local check; it does not call the archived scientific verifiers.
- `reference/all_modes_157ec0f.py`: an unmodified pinned Paper I reference with SHA256 `55ba3d4ba6c44104668362258aa96465bd485fef8f407f65283e12aa0c643ebe`.

## Result status

The derivations define how to assemble the missing stationary operator. The new implementation has **not** been scientifically executed. No finite radial `K,E,D` matrices, physical spectrum, on-shell seed, nonlinear compatibility value, or degree-of-freedom count have been produced. Syntax/metadata inspection, if recorded in the artifact manifest, is not a scientific verification.

Two preparation reviews reconciled the potential normalization, realification convention, and adjoint/boundary distinction. The review also exposed and corrected a missing one-half in the draft potential cross terms and errors in a draft finite momentum-conversion description before scientific execution. Review is not a substitute for the proposed independent calculations.

## What the parallel tracks need from one another

Paper I retains the second-order nonlinear-preservation derivation and the full continuum linear-Euler residual evaluator. This package supplies the operator/boundary interface needed to construct a retained on-shell linear seed. The seed must include the reconstructed spacetime perturbation and sufficient time/radial jets; discarded temporal residuals cannot be inferred from retained Hamiltonian matrices alone. Pure continuum-gauge fields are useful controls but cannot expose a nonzero compatibility pairing when their full linear Einstein residual vanishes.

The two operator preparation tasks can proceed analytically in parallel, but finite assembly needs one explicit radial projection order and one boundary/test-space choice. A common continuum physical evolution domain has not been proved. Equal radial bubble spaces can break exact derivative/product inclusions; those approximation defects must be measured separately from a failure of the angular theory.

## Decision after preparation

The smallest next action is the separately approved kinematic check in `next-check.md`. A pass would justify the next independent curvature/Green-identity validation. It would not authorize a spectrum run. Physical boundary choices, finite radial assembly, numerical rank policy, and a genuinely on-shell seed remain subsequent decisions. Local exact work is sufficient for the proposed check; no SSH or B200 allocation is requested.
