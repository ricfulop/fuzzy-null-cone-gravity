# Paper II: Toward Interacting Einstein Gravity from Finite Causal-Diamond States

**Author:** Ric Fulop, Center for Bits and Atoms, Massachusetts Institute of Technology

Companion dynamical manuscript and reproducible symbolic/numerical verification suite inside the `ricfulop/fuzzy-null-cone-gravity` repository. [Paper I](../paper1/README.md) lives in `paper1/` (`../paper1/paper.tex`), while Paper II lives in this dedicated `paper2/` folder. Run all commands below from `paper2/`; each paper's source, PDFs, figures, reports, and verification scripts stay in its own folder.

## Read the draft and audit

- [Working draft v0.2 PDF](paper2_einstein_gravity.pdf), built from [paper2_einstein_gravity.tex](paper2_einstein_gravity.tex).
- [Paper I pre-submission interaction-and-symmetry audit](reports/paper1-action-audit.md).
- [Proposed common regulator for both papers](reports/common-regulator-bridge.md): exact product-mismatch identities and a prescription for evaluating one Galerkin action with complete Gaunt tensors. This is a follow-up proposal; the v0.2 manuscript and its archived calculations are unchanged. Check its algebra with `python3 verification/common_regulator_identities.py` (Python standard library only).
- [Active quantum research objective](GOAL.md): completion criteria for finite-rank one-loop amplitudes, dynamical corner states, and Euclidean--Lorentzian reconstruction from one theory. The first tracked prerequisite is the [lapse-constraint test on vacuum initial data](reports/quantum-bridge/lapse-constraint-on-vacuum-data.md). Tracked local CPU calculations share the recipe `bash verification/quantum_bridge/run.sh`, launched through the project's experiment runner; `verification/quantum_bridge/config.json` selects the calculation. Reproduce a historical test from its recorded source commit.

- [Verified constraint obstruction and its cutoff dependence](reports/quantum-bridge/constraint-result-and-next-calculation.md), [kinetic-completion result](reports/quantum-bridge/kinetic-completion-result.md), and [conditional corner/reconstruction derivations](reports/quantum-bridge/corners-and-reconstruction.md). These follow-up reports preserve the v0.2 manuscript and distinguish tested components from the three unfinished quantum deliverables.
- [Confirmed canonical mixed-bracket obstruction](reports/quantum-bridge/canonical-mixed-bracket.md) and [explicit uniform-corner quantum states, matching, and transfer kernels](reports/quantum-bridge/corner-quantization-and-transfer.md). The latter supplies a positive representation of the stated boundary sector, while proving that an area bound alone does not give a finite state count.
- [Continuum all-plus one-loop reference](reports/quantum-bridge/continuum-one-loop-reference.md): a dimensionally regulated quantum benchmark and its verification recipe. This calculation has no angular rank cutoff; the finite-$N$ Gaunt amplitude remains a separate unresolved task.
- [Reduction-first canonical Gaunt candidate](reports/quantum-bridge/reduced-canonical-gaunt.md): a declared regulator revision with reconstructed constraints and a canonical symplectic form, a ten-mode TT example, specific interaction coefficients, and a leading uniform-area/free-transfer construction. Its dependent geometry is not capped at the free angular rank; the complete interacting quantum goals remain open.
- [Cubic and quartic reduced interaction vertices](reports/quantum-bridge/reduced-vertices.md): the complete classical weak-field functionals, including the scalar and longitudinal constraint responses. Seven exact directions, held-out rotational invariants and 32 independent finite-metric evaluations passed their specified checks.
- [Sixth-order reduction and the finite-mode loop kernel](reports/quantum-bridge/reduced-sixth-order-and-loop-kernel.md): derived recursion identities for all fifth/sixth-order sectors and a precise four-point loop prescription retaining every vertex family. The higher gravitational coefficients and the actual loop integral remain unevaluated.
- [Independent sixth-order and loop audit](reports/quantum-bridge/sixth-order-independent-audit.md): checks the displayed recursion, propagator, normalization and combinatorial factors under their stated hypotheses. It restricts the static square-root check to four configuration legs at zero momentum and frequency.
- [Complete quartic tensor and finite-mode two-point calculation](reports/quantum-bridge/complete-quartic-and-two-point.md): all 27 quartic invariant coefficients passed their specified exact checks. After a preserved coarse-quadrature failure, the 64/128-point child passed its controls and evaluated twenty diagonal plus two mixed one-loop two-point cases, with maximum scaled analytic comparison error $8.49\times10^{-17}$. This is a fixed ten-mode coefficient, not four-graviton scattering or a physical large-N limit.
- [Local corner-area matching](reports/quantum-bridge/local-corner-matching.md): derives a nonzero local area bracket, shows why the naive matching constraints are not all first class, and constructs two explicit local reductions with different physical interpretations. These analytic results do not derive microscopic corner multiplicities or an interacting positive quantum theory.
- [Local existence and smooth ADM energy](reports/quantum-bridge/local-existence-and-adm-smoothness.md): a weighted elliptic and implicit-function argument establishes the local smooth reduced Hamiltonian for each fixed finite Schwartz TT basis, including equality of its volume energy to the ADM flux. Its neighborhood is not uniform in rank.
- [Constructive interacting quantum completion](reports/quantum-bridge/positive-local-quantum-completion.md): a chosen global extension of that local Hamiltonian admits direct Weyl quantization with a proved self-adjoint domain, lower bound, ground state and interacting Euclidean/Lorentzian transfer. It preserves the local formal loop prescription but does not determine the original microscopic theory or its large-rank limit.
- [Contracted four-point implementation](reports/quantum-bridge/contracted-four-point-implementation.md) and [verified higher-geometry/Green ingredients](reports/quantum-bridge/higher-geometry-green-validation.md): two geometry directions, 68 independent finite-metric evaluations and 27 first-Green checks passed, along with an independent pure-momentum sixth-order control. The full mixed higher Hamiltonian and complete gravitational four-point loop remain unevaluated.
- [Full fifth-order directional evaluation](reports/quantum-bridge/full-fifth-order-evaluation.md): two directions retaining every H5 momentum sector passed the declared checks after a preserved spatial-refinement failure. The maximum scaled spatial/refinement diagnostic is $7.58\times10^{-10}$ under the unchanged $10^{-6}$ gate. This does not supply the contracted loop contribution; the [analytic audit](reports/quantum-bridge/h5-pairing-audit.md) also derives an axial contraction reduction for that future work.
- [Sixth-order bilinear and tail reduction](reports/quantum-bridge/sixth-order-bilinear-and-tail-reduction.md): removes explicit second-response fields from the integrated Hamiltonian while retaining their scalar/York bilinears, and isolates the three-HEAT algebraic tail. The formulas have independent algebra checks; the complete mixed energy remains unevaluated.
- [Validated sixth-order source and integration ingredients](reports/quantum-bridge/sixth-order-backend-validation.md): nine Fourier controls, 72 synthetic tail checks, 680 source-component checks, exact momentum parity and two gravitational tails passed, including the required nonzero-monopole case. A preserved work-cap interruption was repaired without changing tolerances or the cap. This does not yet supply a mixed H6 value or a four-point result.

- [Complete sixth-order numerical result](reports/quantum-bridge/full-sixth-order-results.md): after a [preserved work-cap interruption](reports/quantum-bridge/full-sixth-order-pilot.md), all twenty mixed energy entries passed their numerical gates, giving H6=14796721.40459518 in the stated raw-profile units, with summed diagnostic 5.79e-5. Pure momentum passed its independent radial comparison. Independent full mixed-functional validation and the four-point loop remain pending.
- [Exact Gaussian-kernel consolidation](reports/quantum-bridge/gaussian-kernel-consolidation.md): the rational repair reduced the scalar/York kernels to 207/265 terms under unchanged quadrature, tolerance and evaluation limits. The [unreduced sixth-order reference design](reports/quantum-bridge/unreduced-sixth-order-reference.md) and [exact inverse fixtures](reports/quantum-bridge/spherical-green-control-design.md) specify the next independent spatial comparison; those numerical checks are pending.
- [Independent inverse-primitive protocol](reports/quantum-bridge/independent-inverse-primitive-protocol.md) and [initial result](reports/quantum-bridge/independent-inverse-initial-result.md): all 26 angular controls passed; two high-l scalar regular-third-derivative checks failed the initial resolution gate with roughly second-order refinement. The vector stage was not reached. A finer-grid child retains all tolerances. This remains a prerequisite for the independent mixed H6 reference.
- [Oblate quadratic operator and boundary preparation](reports/oblate-quadratic-preparation/README.md): formal Hamiltonian blocks, domain requirements, and bounded kinematic verification for the original Paper I Galerkin action.
- [Bounded oblate weak-form validation](reports/oblate-weak-curvature-validation/full-recheck-pass/README.md): all eight batches passed, including representative exact sampled Green checks and three fixed symbolic gauge controls; finite radial assembly and physical boundary choices remain open.

## Reproduce all calculations and figures

Install `numpy`, `scipy`, `sympy`, and `matplotlib`, then run from `paper2/`:

```bash
python3 verification/audit_paper1_action.py
python3 verification/bondi_sachs_4pt_scattering.py
python3 verification/fuzzy_scaling_resolution.py
python3 verification/projection_and_constraint_defects.py
python3 figs/generate_figures.py
```

Build the PDF with pdfLaTeX:

```bash
pdflatex paper2_einstein_gravity.tex
pdflatex paper2_einstein_gravity.tex
```

## Verification inventory

- `verification/audit_paper1_action.py` -> `audit_paper1_action_results.json`: Symbolic proof of the five structural discrepancies in Paper I's Appendix A Eq. (A.1) (commit `e541955495d0` / `49de5fa`) against Chandrasekaran Appendix B (`arXiv:2611.34924`).
- `verification/bondi_sachs_4pt_scattering.py` -> `bondi_sachs_4pt_results.json`: Exact nonlinear Bondi–Sachs identity check, fourth-order weak-field compilation (12/25/35 monomials), Schur-complement constraint reduction, independent four-graviton KLT matching across all helicity channels, factorization, $SO(3)$ foliation invariance, Ward projection, Weinberg soft limit, controlled ablations, and Moyal phase scans.
- `verification/fuzzy_scaling_resolution.py` -> `fuzzy_scaling_results.json`: Arbitrary-precision Racah Wigner-$6j$ vertex ratios up to $N = 512$ across fixed-$j$, mesoscopic $j = c\sqrt{N}$, fixed-$p_\perp$ saturated-cap, and decoupled sequences, plus comparison of Routes A, B, and C (band-limited commutative Gaunt product).
- `verification/projection_and_constraint_defects.py` -> `projection_defects_results.json`: Symbolic checks of Propositions 1–3, resolvent virtual-state self-energy on a 2,497-state Wigner-$6j$ occupation-budget Fock sector, and 1- and 2-body spin-2 Gel'fand–Naimark Lorentz boost compression defects.
- `figs/generate_figures.py` -> `figs/fig{1,2,3}_*.{pdf,svg}`: Publication vector figure generator using `figs/orx_figstyle.py`.
