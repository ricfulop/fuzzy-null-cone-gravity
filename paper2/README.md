# Paper II: Toward Interacting Einstein Gravity from Finite Causal-Diamond States

**Author:** Ric Fulop, Center for Bits and Atoms, Massachusetts Institute of Technology

Companion dynamical manuscript and reproducible symbolic/numerical verification suite inside the `ricfulop/fuzzy-null-cone-gravity` repository. [Paper I](../paper1/README.md) lives in `paper1/` (`../paper1/paper.tex`), while Paper II lives in this dedicated `paper2/` folder. Run all commands below from `paper2/`; each paper's source, PDFs, figures, reports, and verification scripts stay in its own folder.

## Read the draft and audit

- [Working draft v0.2 PDF](paper2_einstein_gravity.pdf), built from [paper2_einstein_gravity.tex](paper2_einstein_gravity.tex).
- [Paper I pre-submission interaction-and-symmetry audit](reports/paper1-action-audit.md).
- [Oblate quadratic operator and boundary preparation](reports/oblate-quadratic-preparation/README.md): formal Hamiltonian blocks, domain requirements, and bounded kinematic verification for the original Paper I Galerkin action.

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
