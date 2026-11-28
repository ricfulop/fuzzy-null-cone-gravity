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
- [Cubic and quartic reduced interaction vertices](reports/quantum-bridge/reduced-vertices.md): the complete classical weak-field functionals, including the scalar and longitudinal constraint responses. The committed verification protocol evaluates both tensor parities and checks the geometry with a separate finite matrix exponential. Higher vertices, quantization and the finite-rank loop remain unfinished.

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
