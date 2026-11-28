# Bounded oblate weak-curvature validation

Status: **the full eight-batch bounded local validation passed** in run `32b93879-a4f2-426f-82ac-2fa455844477`, from frozen source `b59be9535b38e25d59a8437efe9873094b2bc4bd`. The [full result and provenance](full-recheck-pass/README.md) record 108 potential-pair and 36 lapse-pair exact sampled checks, three fixed symbolic gauge controls, Cartesian TT normalization and wrong-sign rejection, and a nonzero omitted-cross/trace ablation. Supervised time was 129.526260 seconds, peak sampled process-group RSS was 161.25 MiB, and cleanup ended with an empty process group.

The earlier [first full attempt](first-attempt/README.md), [second full attempt](second-attempt/README.md), and [incomplete curvature probe](curvature-probe-attempt/README.md) remain unchanged. They stopped on monitored memory or wall budgets and did not complete their configured validations. The later [single-entry diagnostic](ricci-entry-pass/README.md) and [complete five-batch curvature diagnostic](curvature-recheck-pass/README.md) passed before the full recheck. The checker bytes are identical across all three passing runs: SHA256 `3bb2a49863cda8b6f5f3b7a1a718ec580cc6916ed6f5cde13d682d2571ac7113`, verified against each archived source and Git commit. Their [single-entry](ricci-entry-pass/execution.json), [curvature](curvature-recheck-pass/execution.json), and [full-validation](full-recheck-pass/execution.json) execution manifests record the same checker hash; within the six pinned source files, only the supervisor's fixed mode changed.

This check tests the potential and lapse densities prepared at `510e7d87d0e8cd1581ffc7233122a7bedafea969`, for the original Paper I coefficient-Galerkin action and its oblate flat background. It does not test the separate reduction-before-cutoff Gaunt candidate. The source action and retained harmonic conventions remain pinned to the prior preparation; no scientific equation in `weak_forms.py`, `potential_hessian.py`, or `kinematic_blocks.py` is changed here.

The hypothesis is that the prepared first-derivative densities and radial Green fluxes agree with an independently assembled metric/Christoffel variation in the recorded representative sectors. Independent here means the strong curvature is obtained from varying the connection and Ricci tensor, while the tested weak density is the previously prepared first-derivative implementation. Both use the same stated metric and retained tensor basis.

## Completed evidence and limits

- At $L=2$, six representative tensor types were tested in each of $m=0,1,2$: the first available scalar lambda and X columns, polar/axial vectors, and polar/axial tensors. This is not every angular coefficient at arbitrary cutoff.
- The background connection was recomputed using zero azimuthal mode, and first-derivative tensor jets were compared symbolically, including radial derivatives of the basis.
- Curvature-formula and potential/lapse Green residuals passed exactly at one recorded rational point per $m$, using consistent polynomial radial profiles. These are exact sampled checks, not a global symbolic proof of the Green identities.
- The pointwise Green check included the full spatial divergence. The azimuthal mode of a conjugated equal-$m$ bilinear is zero. Exported radial fluxes were checked separately without imposing boundary conditions.
- Three fixed continuum pure-gauge perturbations passed symbolic checks of all nine Einstein entries and the scalar curvature. The Cartesian transverse-traceless control passed through the actual prepared and independent evaluator paths with the positive weak density and correct strong normalization.
- Omitted trace/cross terms gave a nonzero residual, and the wrong-sign TT operator was rejected.

The eight completed batches are three sector checks, three pure-gauge checks, one TT check and one ablation check. The checker records its selected labels, rational points and profiles. No angular integral, radial discretization, finite `K,E,D`, common self-adjoint domain, physical spectrum, admissible non-gauge seed, degree-of-freedom count or nonlinear consistency theorem was produced. The next preparation is the raw four-slot angular coefficient/flux interface and finite-basis inclusion diagnostics under an explicit projection/test-space convention.

## Execution contract

The fixed command for this independent experiment baseline is:

```sh
bash paper2/verification/oblate_quadratic_preparation/run_weak_validation.sh
```

Run from an immutable local source snapshot. The shim changes to `paper2/`; the supervisor launches isolated Python 3.12 with SymPy 1.14.0 and NumPy 2.3.5. Assertions stay enabled and numerical worker counts are limited to one. Scientific output goes to the fresh snapshot-local directory `paper2/reports/oblate-weak-curvature-validation/runtime/`.

The first two full-validation invocations had a 900-second wall cap. The targeted diagnostics and final full recheck used a 300-second cap. All used a 2-GiB sampled process-group RSS budget, monitored rather than enforced as an OS hard cap. Setup/dependency resolution is inside the wall budget. The supervisor saves source hashes and environment provenance before launch, resource samples and terminal status afterward, and cleans up the child's process group on interruption or budget exhaustion. The archived full run log is the source of truth for controls and exceptions.

The user's overnight authorization is `msg_ad636446`, followed by approval for other project decisions. The shared conservative cutoff is **2026-11-28 16:29 UTC**. The supervisor reduces the wall budget near the cutoff and refuses a new launch after it. The 12-hour window is an authorization window, not a request to inflate this check's workload. No paid remote compute is used.

## Preflight review

The draft verifier was corrected during source review before any scientific execution: background versus perturbation azimuthal derivatives, the raised connection index, coordinate Lie derivatives, consistent second-derivative profile data, equal-mode bilinear divergence, genuine TT transversality and trace-sensitive negative controls. These were checker implementation defects caught in preparation, not failed scientific runs or changes to the gravitational action.

The weak-density formulas received a separate analytic review. The original `preflight-review.json` records the first source snapshot and remains historical provenance. Later preparations and results have separate records. The final result received an independent read-only provenance/scope audit (`msg_05c31715`) and a source-scope review (`msg_d88daac5`). Static syntax and source reviews are distinct from the measured validation.

The answered experiment source remains frozen at `b59be9535b38e25d59a8437efe9873094b2bc4bd`. This result archive and current-status documentation are maintained on the separate `paper2/oblate-weak-validation-results-2026-11-28` branch. Further scientific changes require a new child experiment rather than rewriting the answered node.
