# Bounded oblate weak-form validation passed

Run `32b93879-a4f2-426f-82ac-2fa455844477` passed all eight planned validation batches from frozen source `b59be9535b38e25d59a8437efe9873094b2bc4bd`. This checks the prepared weak potential/lapse densities and radial fluxes for the original Paper I coefficient-Galerkin action, pinned to `157ec0f04866ea36509b325ec65f0a68890f9b17`. It does not validate the distinct reduction-before-cutoff Gaunt candidate.

Evidence: [structured result](weak-form-checks.json), [raw log](verification.log), [execution manifest](execution.json), [eight checkpoints](weak-form-progress.jsonl), [resource samples](resource-samples.json), and [SHA256 manifest](sha256-manifest.json).

The supervisor recorded exit 0 in 129.526260 seconds, no stop reason, 169,082,880 bytes (161.25 MiB) peak sampled process-group RSS across 439 samples, and an empty final process group. A later independent process-table check also found no members. The checker recorded 128.458574 seconds under Python 3.12.11, SymPy 1.14.0 and NumPy 2.3.5. The 2-GiB memory budget was sampled, not a hard OS limit.

## What passed

Each m=0,1,2 sector at L=2 selects one retained tensor from each of six species: lambda, X, polar/axial V, and polar/axial C. The result lists all selected labels. Independent background Christoffels and tensor jets are compared with the prepared paths, including radial derivatives of the basis and separate radial coefficient/derivative slots. Curvature is computed from the independently varied connection/Ricci path, with an additional comparison to compact linearized formulas.

The sector-dependent strong/weak potential and lapse Green residuals, including the local radial flux expressions, are evaluated exactly at these recorded points:

| m | r | x | epsilon | potential pairs | lapse pairs |
|---|---|---|---|---:|---:|
| 0 | 3/2 | 1/5 | 1/4 | 36 | 18 |
| 1 | 23/14 | 1/5 | 1/4 | 36 | 12 |
| 2 | 25/14 | 1/5 | 1/4 | 36 | 6 |

Thus 18 selected tensor instances cover 108 potential pairs and 36 lapse pairs. These are exact sampled checks, not identities established for every chart point or every retained label. Nonconstant profiles are used consistently: left r^2-2r+5/4, right r^3-r+3/2, and lapse r^2+r/3+2/5. The Green comparison includes the complete divergence of the conjugate-paired flux; the radial component is compared separately without imposing endpoint conditions.

The three prepared continuum gauge test vectors each pass all nine Einstein-component and one scalar-curvature exact-zero checks symbolically in the chart variables. This covers those fixed vectors, not arbitrary gauge functions or a finite projected gauge kernel. The Cartesian transverse-traceless test passes its strong/weak/flux identity, weak f'^2 normalization and strong -f*f'' normalization, while rejecting the wrong sign. The nonzero lambda/X ablation rejects omission of the cross/trace terms.

All eight checkpoints and the final `weak-form-checks.json` exist. The log reports validation progress 0 through 8, then `validation_passed=1` only after the result write. The full source hashes match both the immutable snapshot and pinned Git commit. The checker hash is identical to the two preceding passing curvature diagnostics, and the three physics modules remain unchanged from the preparation.

## What this closes and leaves open

This closes the prepared bounded local potential/lapse weak-form and Green-flux validation checkpoint. Earlier failures remain archived: a memory stop followed by cleanup error, a blanket-normalization time regression, and an incomplete five-batch curvature diagnostic. The later single-entry and complete-curvature passes remain separate supporting records. No failed record has been replaced.

It does not perform angular integration, finite radial projection, endpoint-domain selection, integrated endpoint-boundary testing, finite-basis kernel preservation, finite K/E/D assembly, an admissible non-gauge seed construction, a spectrum, a degree-of-freedom count, or nonlinear consistency classification. No independent finite-deformation curvature oracle is added by this run. Those are subsequent tasks with their own tests and assumptions.

Before finite radial assembly, the natural next preparation is an explicit raw four-slot angular coefficient/flux interface and a finite-basis inclusion check under a declared projection/test-space convention. A diagnostic shell is not a physical boundary model. The result supports that preparation without making those choices implicitly.

## Review and preservation

Independent read-only result and archive audits found no discrepancy (`msg_05c31715`, `msg_6aa9451a`). They checked the raw log, eight checkpoints, result/source hashes, sampled resource maximum, empty process group, and byte agreement between the run snapshot, repository archive, and durable artifact copy. A separate source-scope review (`msg_d88daac5`) confirmed the bounded coverage and remaining operator/domain requirements. See [review metadata](review.json). These agent reviews did not rerun scientific code and are not external peer review.

The executed experiment source is frozen. This archive and current preparation-status updates are recorded separately on `paper2/oblate-weak-validation-results-2026-11-28`; the earlier incomplete and diagnostic-only records remain intact.
