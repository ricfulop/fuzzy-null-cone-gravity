# Resolve the nonlinear consistency of the finite-state gravity model

**Requested by Ric Fulop, 27 November 2026.**

**Status:** active research plan; the initial checkpoint and an actual-action primary–secondary constraint calculation are complete. Tertiary stabilization and the full objective remain open.

## Objective and completion standard

Resolve the finite-rank gravitational constraints, the boundary quantization and gluing, and the connection to the simulation in one mutually consistent formulation. Preserve the citation corrections and all archived simulation results. The starting scientific source is repair PR #6, commit `ab303dd3f756c7f3f8d393ef090d4da35a7fdae2`.

A result can establish an exact construction, a controlled approximation with stated errors, or an impossibility theorem with explicit hypotheses. A numerical residual, a toy model, or an obstruction to one candidate must not be described as completing the gravitational theory. This goal is complete only when the three components below have compatible definitions and their stated completion tests pass, or when a decisive obstruction to the requested formulation is proved and the exact assumptions that must change are documented.

## Completion tests

| Component | Required deliverable | Completion test |
|---|---|---|
| Gravitational constraints | Canonical analysis of the actual proposed bulk action, retaining lapse, shift, and all metric components until variation | Derive every primary and secondary constraint; classify first/second class; compute the brackets on the constraint surface; count physical degrees of freedom; prove the claimed algebra on the specified domain and test representative finite ranks |
| Boundary quantum theory | Specified null boundary data, generator normalization, joint terms, reduced symplectic form, representation, inner product, and gluing map | Show differentiability of the action, compatibility of the constraints and state cutoff, positivity of the physical inner product, and consistent sewing with the correct dual orientation; specify exactly which symmetries survive |
| Simulation connection | Derivation of a sampled or otherwise computationally accessible theory from that same constrained action | Account for gauge fixing, determinants, integration contour, couplings, and observables; match quadratic, cubic, and quartic terms and the corresponding symmetry identities; distinguish exact equivalence from controlled effective-theory matching |

## First checkpoint: completed

See [the research note](reports/nonlinear-consistency-2026-11-27/README.md) and its JSON results.

- An exact counterexample rules out the naive projected angular Poisson bracket at harmonic cutoff L=3.
- A matrix-commutator candidate closes that area-preserving angular subalgebra, with a measured departure from its classical bracket. It does not yet supply the full gravitational constraints.
- Trace arguments exclude an exact canonical pair or a nontrivial unitary sl(2,R) representation on a finite-dimensional space. An exact identity quantifies errors from compressing observables.
- A finite reduced R x SU(2) boundary example implements area matching, dual orientation, a positive projection, and associative sewing. It treats the boost angle as a gauge parameter. Its derivation from the gravitational corner phase space remains open.
- The real Euclidean conformal-sign comparison rules out direct unreduced identification with the positive scalar sector of the archived HMC action. An effective or complex-contour relation remains an independent research question.
- The first canonical calculation proves that angular velocity truncation and the Legendre transform do not commute. In a spin-two kinetic block, a nonconstant lapse acquires nonlinear Hamiltonian dependence through discarded-mode leakage. The first-order constrained formulation must therefore be derived independently.
- The current bulk pullback has a guaranteed exact subgroup of base diffeomorphisms and base-dependent SO(3) rotations. This is a proved symmetry of the retained ansatz, not yet a classification of all canonical constraints or physical degrees of freedom.
- Companion Paper II PR #7 was independently reproduced at commit `5e0ab55afed69f24432bccffd6ac2a781f29cfbf`: the continuum weak-field four-graviton benchmark passes, as do 20 additional nonlinear curvature-identity configurations. Its projected Gaunt product and kinematic Ward projection must not be identified with the unprojected nonlinear action or a full finite-rank Dirac-algebra proof.

Reproduce this checkpoint on CPU:

```bash
python verification/closure_audit_20261127/run_all.py
```

## Next work, in dependency order

The [actual-action canonical note](reports/canonical-analysis-2026-11-27/README.md) now derives the exact ADM variable map, a regular spatial Legendre domain, and the temporal primary and secondary equations. At explicit vacuum data satisfying all primaries and secondaries, the scalar lapse Hessian has rank $4(N-1)$: 8 at $N=3$ and 12 at $N=4$. A general harmonic proof and independent local checks support this result. It does not determine the final Dirac surface or physical degree-of-freedom count. The next proposed calculation is to stabilize the secondaries, including their null directions and radial derivatives, before changing the regulator. It awaits the user's next copilot decision.

Reproduce this calculation from `paper1/` with `python verification/canonical_analysis_20261127/actual_adm.py` (NumPy and SymPy; no GPU jobs).

1. **Derive the canonical constraints first.** Start from the first-order ADM formulation with a stated angular coefficient ansatz. Compare it with the Legendre transform of the current Lagrangian Galerkin restriction; truncation and the Legendre transform need not commute. Use N=3 and N=4 as explicit tests, then formulate the general-rank statement. Do not identify an SU(N) commutator with the hypersurface-deformation algebra without deriving the mixed and Hamiltonian brackets.
2. **Test reduction before the state cutoff.** Preserve the gravitational constraint structure on the kinematical space, reduce it, and only then test whether the clock/occupation restriction defines the claimed finite physical sector. Prove preservation, or quantify leakage, under dynamics and the physical observables. This is a working route, not an established solution.
3. **Derive the boundary reduction.** Choose null boundary conditions and polarization explicitly. Determine whether the tested area-character/compact-spin model follows from that reduction, or which omitted generators and edge modes must be restored. Do not require an operator conjugate to area inside a fixed-area sector merely because the unreduced classical phase space has that pair.
4. **Derive the computational measure.** Reuse the independently checked continuum vertices from Paper II PR #7 as weak-field matching targets, with their stated momentum-domain restrictions. Determine the real or complex integration contour and the reduced effective action. Solving constraints does not by itself guarantee a bounded real Euclidean action. If the existing HMC model fails matching, implement a new evaluator with a separate name and preserve the old runs as benchmarks of their original model.
5. **Integrate supported results.** Update the manuscript and figures only at the strength established by the derivations. Rebuild the PDFs and retain a source-to-claim audit. Keep citation fixes intact.

## Checkpoint rules

- Record assumptions, proofs, counterexamples, source versions, commands, and numerical tolerances.
- Use exact symbolic calculations where feasible and independent numerical checks where they test a distinct failure mode.
- Do not use smaller residuals alone as proof of constraint closure or quantum equivalence.
- Do not count trivial added gauge variables as recovery of Einstein's constraint algebra.
- Begin with bounded local CPU calculations. Any proposed large external computation should have a concrete scientific purpose, cost estimate, and validation plan before execution.
- If a candidate fails, preserve the counterexample and state precisely which formulation it excludes. Continue with the next justified construction; do not silently redefine the goal as a completed theory.
