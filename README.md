# A Finite-State Model of Quantum Gravity in a Causal Diamond

**Author:** Ric Fulop ([ricfulop@mit.edu](mailto:ricfulop@mit.edu)), Center for Bits and Atoms, Massachusetts Institute of Technology

This repository contains the manuscript source and compiled PDF, verification scripts and archived data, Lean proof modules, and editable figure artwork. The reviewed revision adds explanatory figures, enforces the two-particle budget in the pair-creation diagnostic, and reports the B200 results with observable-specific sampling qualifications. [Revision details and figure specifications](reports/manuscript-review-2026-11-27.md).

---

## Summary of Results

1. **Classical Reduction and Axiomatic Bridge:** All 14 classical first-order null-cone identities (Sec. 3 and Appendices A–B of `arXiv:2611.34924`) hold identically in SymPy, and Appendix E.2's conditional Lorentzian-to-Euclidean equivalence theorem holds without gaps once Axioms 1–5 are met (checked in Z3 and Lean 4).
2. **Four Continuum Obstructions:** Direct continuum quantization encounters four independent obstructions: Sec. 5.2 non-entailment and circularity ($\mathbb{C}^2 \otimes \ell^2(\mathbb{N})$ counterexample), unrenormalized focusing-clock zero-point collapse ($P_{\mathscr{D}} \equiv 0$) and inverse-lapse divergence ($N^{-1} \to \infty$), Sec. 4.4 quadratic angular divergence ($\propto J^2$) at fixed $G_N$, and a two-loop quadratic sunset divergence ($\sim \Lambda_\parallel^2 \ln \Lambda_\parallel$) for any transverse-only regulator.
3. **Fuzzy-Sphere Null-Cone Regularization:**
   - **Fuzzy Monopole Bundles and Matrix Action $I_{\mathrm{FSG}}^{(N)}$:** Rectangular $SU(2)$-equivariant bundles $\mathcal{B}_{s_w}^{(N)} = \mathrm{Hom}(V_{(N-1-s_w)/2}, V_{(N-1+s_w)/2})$ of dimension $N^2 - s_w^2$, explicit $SU(2)$-equivariant embedding $\iota_{N \to M}$, closed Wigner-$6j$ parity-symmetrized cubic vertex $\mathcal{V}_3^{(N)}[F; Z_2, Z_3]$ (verified to $3.96 \times 10^{-14}$), $O(1/N^2)$ convergence to 4D GR, Regge–Wheeler/Zerilli isospectrality, sub-Planckian $SO(3,1)$ boost closure, and Ermakov–Pinney caustic regularization.
   - **Smoothly Smeared Massive-QEI and Exponential Clock:** $C^\infty_c$ bump smearing bounded uniformly across all 4D spin-2 harmonics via Fewster–Eveson massive quantum energy inequalities, monotone spectral-projector equivalence, and exact zero-inversion Balian–Brézin Fredholm determinants ($\lambda_{\mathrm{crit}}^{\mathrm{phys}} \in [17.0, 22.3] \gg 2$ inside $B(s)=s Q_*$).
   - **Area-Scaled Total Occupation Budget and Projected Dynamics:** $\sum_{k,j,m,p} n_{k,j,m,p}(k + t_j^2/k) \le B(s) = s Q_*$ with projected Hamiltonian $H_{\mathrm{phys}} = P_{B(s)} H_{\mathrm{FSG}} P_{B(s)}$ (zero sector leakage, $O(1/N^2)$ interacting amplitude and 1-loop self-energy convergence), yielding finite fixed-area sectors and finite projected perturbative sums while preserving the macroscopic $s^{-3/2}$ and $-\ln s/(s_{\min} s^2)$ pair-creation asymptotics via an $O(s^{-2})$ tail bound.
4. **Universal 4D One-Loop, Two-Loop, and Semiclassical Benchmarks:**
   - One-loop logarithmic corrections: Sen Schwarzschild mixed ($+\frac{106}{45}$, $0.09\%$) and singlet ($+\frac{77}{90}$, $0.26\%$), Benedetti–Casini flat-space graviton entanglement entropy ($-\frac{61}{45}$, $0.02\%$), and Anninos–Denef–Law–Sun de Sitter $S^4$ ($-\frac{571}{45}$, exact).
   - Two-loop Schwarzschild Goroff–Sagnotti shift $\frac{209}{1440}\frac{\ln(s/s_{\min})}{s}$, sub-Planckian pixel-size ceiling $s_{\min} \le \pi$ ($2 \le d_0 \le 22$), dynamical-rank $\mathrm{AdS}_4$ Hawking–Page transition, 2,048-state three-diamond corner fusion and generalized entropy additivity ($8.88 \times 10^{-16}$), and finite-dimensional realizations of Secs. 4.4, 5.4, and 6 (including Don Page's pure-graviton Hawking luminosity $\alpha_{\mathrm{grav}} = 3.8367 \times 10^{-6}$, consistent with Page’s printed precision).
5. **Matrix simulations on NVIDIA B200 GPUs**
   - The archived production set contains 28 ensembles with 4–8 independent chains. At rank 48, each chain has 331,392 coupled real variables. Four independent copies account for 1,325,568 variables in parallel.
   - Increasing the microscopic coupling from 0.5 to 2 reduces measured roughness by 40.45–43.47% across four ranks, at the same lattice and scalar self-interaction within each comparison. The analytic Gaussian reference gives a 47.69–48.10% reduction at coupling 2. Matched zero-coupling controls at the production stiffness remain pending in this revision.
   - The archived ordinary R-hat for roughness at coupling 2 is 1.000–1.0035. Local mean area has R-hat as high as 1.416. Its sampled distribution requires further mixing checks. Positive definiteness follows from the exponential parameterization.
   - The largest scaled-coupling propagator approaches its free value. A massless continuum limit requires fit errors and a controlled extrapolation in rank and volume.

---

## Repository Structure

```text
.
├── paper.tex              # APS Physical Review D (revtex4-2) manuscript (Overleaf-ready)
├── paper.pdf              # Compiled reviewed manuscript
├── figs/                  # Science-style vector PDF and editable SVG artwork
├── lean4/                 # 6 standalone Lean 4.36.0 modules (50 verified theorems, 0 sorryAx)
├── verification/          # Verification, plotting scripts, and archived numerical data
└── reports/               # Technical reports and revision notes
```

---

## Reproducing the Lean 4.36.0 Proofs (50 Theorems)

Each Lean 4 module is completely self-contained (zero external Mathlib dependencies; standard Lean core axioms only). With Lean 4.36.0 installed via `elan`:

```bash
for f in lean4/*.lean; do
  echo "=== Checking $f ==="
  lean "$f"
done
```

| Lean 4 Module | Theorems | Certified Logical or Algebraic Scope |
|---|---:|---|
| `lean4/lean4_verified_audit.lean` | 11 | Sec. 5.2 non-entailment, AdS sum iff, clock collapse triad, App. E.2 validity |
| `lean4/lean4_path2_breakthrough.lean` | 9 | Lichnerowicz $+\frac{106}{45}$ / $+\frac{77}{90}$, 4D isotropy, QEI clock, Type I Wedderburn |
| `lean4/lean4_path2_frontier.lean` | 8 | Fuzzy monopole $\mathcal{B}_{s_w}^{(N)}$, $SO(3,1)$ boost closure, parity symmetry, caustic bounce |
| `lean4/lean4_path2_uv_completion.lean` | 8 | 2-loop power count, Schwarzschild $W^3$, $\mathrm{AdS}_4$ rank sum, 3-diamond SSA & CEMZ |
| `lean4/lean4_original_problem_resolution.lean` | 7 | Sec. 4.4 denominator dominance, Sec. 5.4 $SU(2)$ factorization, Sec. 6 Eq. (6.13) cancellation |
| `lean4/lean4_v2_reviewer_fixes.lean` | 7 | Total occupation budget, $B(s)=s Q_*$ tail regions, $P_B H P_B$ zero leakage, massive QEI, $Z_{\mathscr{D},1}=45$ mod-8 |
| **Total** | **50** | **Zero `sorryAx`, zero compiler warnings** |

---

## Reproducing the Numerical, Symbolic, and GPU Monte Carlo Verification

Install standard Python scientific packages (`numpy`, `scipy`, `sympy`, `z3-solver`, `matplotlib`):

```bash
pip install numpy scipy sympy z3-solver matplotlib

# Matrix vertex, projected Hamiltonian, and B(s)=s*Q* tail-bound verification
python verification/v2_test1_matrix_vertex_and_projected_dynamics.py

# Fewster-Eveson massive QEI and Balian-Brezin exponential clock verification
python verification/v2_test2_fewster_eveson_and_exponential_clock.py

# Original continuum audit counterexamples and Z3 dependency checks
python verification/test1_section52_counterexamples.py
python verification/test5_formal_dependency_and_z3.py
```

To reproduce the 28-GPU NVIDIA B200 nonperturbative hybrid Monte Carlo production sweep on Modal (requires `modal` CLI):

```bash
modal run verification/modal_fsg_hmc.py --mode final
python verification/plot_modal_hmc.py
```

---

## Rebuilding the reviewed figures

These commands use completed archived runs and local CPU calculations.

```bash
python verification/v2_test1_matrix_vertex_and_projected_dynamics.py
python verification/plot_path2_frontier.py
python verification/plot_path2_breakthrough.py
python verification/plot_path2_uv_frontier.py
python verification/plot_original_problem_solved.py
python verification/plot_v2_unification_and_clock.py
python verification/plot_modal_hmc.py
python verification/plot_review_figures.py
```

The ten manuscript figures are exported at 18.4 cm width as PDF and SVG. Figure titles are in captions, panel letters are uppercase, and ordinary lettering is at least 7 pt at the artwork’s native size. The supplied LaTeX retains the existing REVTeX manuscript format.

## Compiling the Paper (Local TeX or Overleaf)

This repository is directly importable into **Overleaf** via *New Project → Import from GitHub*. Locally, with TeX Live (`revtex4-2`):

```bash
pdflatex paper.tex && pdflatex paper.tex
```

---

## License

MIT License — see [LICENSE](LICENSE).
