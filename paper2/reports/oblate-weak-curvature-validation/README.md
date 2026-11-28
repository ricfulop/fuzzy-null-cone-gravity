# Bounded oblate weak-curvature validation

Status: both full-validation invocations remain incomplete. The first passed three representative weak-density/Green batches, then reached its monitored memory limit in the symbolic gauge control. The retry timed out in the second sampled sector after completing only the first. A subsequent targeted curvature diagnostic passed its arithmetic and nonzero sampled Ricci comparisons, then reached its 300-second cap inside the symbolic m=1 Ricci entry (1,2). It completed two of five diagnostic batches and produced no final result. The records are preserved in [first-attempt/README.md](first-attempt/README.md), [second-attempt/README.md](second-attempt/README.md), and [curvature-probe-attempt/README.md](curvature-probe-attempt/README.md). Preparation documents are not passing results.

This check tests the potential and lapse densities prepared at `510e7d87d0e8cd1581ffc7233122a7bedafea969`, for the original Paper I coefficient-Galerkin action and its oblate flat background. It does not test the separate reduction-before-cutoff Gaunt candidate. The source action and retained harmonic conventions remain pinned to the prior preparation; no scientific equation in `weak_forms.py`, `potential_hessian.py`, or `kinematic_blocks.py` is changed here.

The hypothesis is that the prepared first-derivative densities and radial Green fluxes agree with an independently assembled metric/Christoffel variation in the recorded representative sectors. Independent here means the strong curvature is obtained from varying the connection and Ricci tensor, while the tested weak density is the previously prepared first-derivative implementation. Both use the same stated metric and retained tensor basis.

## Planned evidence and limits

- At $L=2$, test six representative tensor types in each of $m=0,1,2$: the first available scalar lambda and X columns, polar/axial vectors, and polar/axial tensors. This is not every angular coefficient at arbitrary cutoff.
- Recompute the background connection using zero azimuthal mode. Compare first-derivative tensor jets symbolically, including radial derivatives of the basis.
- Evaluate the curvature-formula and potential/lapse Green residuals exactly at one recorded rational point per $m$, using consistent polynomial radial profiles. These are exact sampled checks, not a global symbolic proof of the Green identities.
- Include the full spatial divergence in the pointwise Green check. The azimuthal mode of a conjugated equal-$m$ bilinear is zero. Exported radial fluxes are checked separately without imposing boundary conditions.
- Check three local continuum pure-gauge perturbations symbolically and a Cartesian transverse-traceless normalization through the actual prepared and independent evaluator paths. The TT control explicitly requires the positive weak density and correct strong normalization.
- Require a nonzero response to omitted trace/cross terms and reject a wrong-sign TT operator. A positive control alone cannot exclude an always-zero evaluator.

The planned coarse batches are three sector checks, three pure-gauge checks, one TT check and one ablation check. The checker records its selected labels, rational points and profiles. No angular integral, radial discretization, finite `K,E,D`, common self-adjoint domain, physical spectrum, admissible non-gauge seed, degree-of-freedom count or nonlinear consistency theorem is produced.

## Execution contract

The fixed command for this independent experiment baseline is:

```sh
bash paper2/verification/oblate_quadratic_preparation/run_weak_validation.sh
```

Run from an immutable local source snapshot. The shim changes to `paper2/`; the supervisor launches isolated Python 3.12 with SymPy 1.14.0 and NumPy 2.3.5. Assertions stay enabled and numerical worker counts are limited to one. Scientific output goes to the fresh snapshot-local directory `paper2/reports/oblate-weak-curvature-validation/runtime/`.

The two full-validation invocations had a 900-second wall cap. The current targeted diagnostic uses a 300-second cap. Both use a 2-GiB sampled process-group RSS budget, monitored rather than enforced as an OS hard cap. Setup/dependency resolution is inside the wall budget. The supervisor saves source hashes and environment provenance before launch, resource samples and terminal status afterward, and cleans up the child's process group on interruption or budget exhaustion. The archived full run log is the source of truth for controls and exceptions.

The user's overnight authorization is `msg_ad636446`, followed by approval for other project decisions. The shared conservative cutoff is **2026-11-28 16:29 UTC**. The supervisor reduces the wall budget near the cutoff and refuses a new launch after it. The 12-hour window is an authorization window, not a request to inflate this check's workload. No paid remote compute is used.

## Preflight review

The draft verifier was corrected during source review before any scientific execution: background versus perturbation azimuthal derivatives, the raised connection index, coordinate Lie derivatives, consistent second-derivative profile data, equal-mode bilinear divergence, genuine TT transversality and trace-sensitive negative controls. These were checker implementation defects caught in preparation, not failed scientific runs or changes to the gravitational action.

The weak-density formulas received a separate analytic review. The original `preflight-review.json` records the first source snapshot and remains historical provenance. The retry has a separate preparation record. Static syntax checks are not scientific verification. Results, including any failure, must be read from the subsequent run before changing the status above.
