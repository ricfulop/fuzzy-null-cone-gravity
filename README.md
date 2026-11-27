# A Finite-State Model of Quantum Gravity in a Causal Diamond

**Author:** Ric Fulop ([ricfulop@mit.edu](mailto:ricfulop@mit.edu)), Center for Bits and Atoms, Massachusetts Institute of Technology

This repository contains the complete APS *Physical Review D* manuscript (`paper.tex` / `paper.pdf`), all **50 kernel-verified Lean 4.36.0 theorems** (zero warnings, zero `sorryAx`), all **29 self-contained Python symbolic, numerical, Z3, and GPU Monte Carlo scripts** with **32 exact JSON datasets**, **8 figure generators**, **6 publication vector PDF figures**, and **7 technical reports** auditing and regularizing Venkatesa Chandrasekaran's algebraic null-cone quantum gravity framework ([arXiv:2611.34924](https://arxiv.org/abs/2611.34924)).

---

## Summary of Results

1. **Classical Reduction and Axiomatic Bridge:** All 14 classical first-order null-cone identities (Sec. 3 and Appendices A–B of `arXiv:2611.34924`) hold identically in SymPy, and Appendix E.2's conditional Lorentzian-to-Euclidean equivalence theorem holds without gaps once Axioms 1–5 are met (checked in Z3 and Lean 4).
2. **Four Continuum Obstructions:** Direct continuum quantization encounters four independent obstructions: Sec. 5.2 non-entailment and circularity ($\mathbb{C}^2 \otimes \ell^2(\mathbb{N})$ counterexample), unrenormalized focusing-clock zero-point collapse ($P_{\mathscr{D}} \equiv 0$) and inverse-lapse divergence ($N^{-1} \to \infty$), Sec. 4.4 quadratic angular divergence ($\propto J^2$) at fixed $G_N$, and a two-loop quadratic sunset divergence ($\sim \Lambda_\parallel^2 \ln \Lambda_\parallel$) for any transverse-only regulator.
3. **Fuzzy-Sphere Null-Cone Regularization:**
   - **Fuzzy Monopole Bundles and Matrix Action $I_{\mathrm{FSG}}^{(N)}$:** Rectangular $SU(2)$-equivariant bundles $\mathcal{B}_{s_w}^{(N)} = \mathrm{Hom}(V_{(N-1-s_w)/2}, V_{(N-1+s_w)/2})$ of dimension $N^2 - s_w^2$, explicit $SU(2)$-equivariant embedding $\iota_{N \to M}$, closed Wigner-$6j$ parity-symmetrized cubic vertex $\mathcal{V}_3^{(N)}[F; Z_2, Z_3]$ (verified to $3.96 \times 10^{-14}$), $O(1/N^2)$ convergence to 4D GR, Regge–Wheeler/Zerilli isospectrality, sub-Planckian $SO(3,1)$ boost closure, and Ermakov–Pinney caustic regularization.
   - **Smoothly Smeared Massive-QEI and Exponential Clock:** $C^\infty_c$ bump smearing bounded uniformly across all 4D spin-2 harmonics via Fewster–Eveson massive quantum energy inequalities, monotone spectral-projector equivalence, and exact zero-inversion Balian–Brézin Fredholm determinants ($\lambda_{\mathrm{crit}}^{\mathrm{phys}} \in [17.0, 22.3] \gg 2$ inside $B(s)=s Q_*$).
   - **Area-Scaled Total Occupation Budget and Projected Dynamics:** $\sum_{k,j,m,p} n_{k,j,m,p}(k + t_j^2/k) \le B(s) = s Q_*$ with projected Hamiltonian $H_{\mathrm{phys}} = P_{B(s)} H_{\mathrm{FSG}} P_{B(s)}$ (zero sector leakage, $O(1/N^2)$ interacting amplitude and 1-loop self-energy convergence), yielding finite fixed-area sectors and all-loop finiteness while preserving the macroscopic $s^{-3/2}$ and $-\ln s/(s_{\min} s^2)$ pair-creation asymptotics via an $O(s^{-2})$ tail bound.
4. **Universal 4D One-Loop, Two-Loop, and Semiclassical Benchmarks:**
   - One-loop logarithmic corrections: Sen Schwarzschild mixed ($+\frac{106}{45}$, $0.09\%$) and singlet ($+\frac{77}{90}$, $0.26\%$), Benedetti–Casini flat-space graviton entanglement entropy ($-\frac{61}{45}$, $0.02\%$), and Anninos–Denef–Law–Sun de Sitter $S^4$ ($-\frac{571}{45}$, exact).
   - Two-loop Schwarzschild Goroff–Sagnotti shift $\frac{209}{1440}\frac{\ln(s/s_{\min})}{s}$, sub-Planckian pixel-size ceiling $s_{\min} \le \pi$ ($2 \le d_0 \le 22$), dynamical-rank $\mathrm{AdS}_4$ Hawking–Page transition, 2,048-state three-diamond corner fusion and generalized entropy additivity ($8.88 \times 10^{-16}$), and finite-dimensional realizations of Secs. 4.4, 5.4, and 6 (including Don Page's pure-graviton Hawking luminosity $\alpha_{\mathrm{grav}} = 3.8367 \times 10^{-6}$, within $0.017\%$).
5. **Nonperturbative Matrix-Gravity Hybrid Monte Carlo (NVIDIA B200):**
   - Exact $T^2 \times S^2_N$ Fourier-accelerated HMC (`verification/modal_fsg_hmc.py`) with exact collective monopole-shear scale updates across 28 converged dual IR/UV ensembles (all $\hat{R}(\phi_0) \le 1.005$, exact Schwinger–Dyson Ward identities satisfied to $0.05\%$ for $Z$ and $0.16\%$ for $\Phi$), plus 4 flat-measure collapse benchmarks.
   - **Macroscopic IR limit ($g_N = g_P/\sqrt{N}$, $N=8\dots 48$ at $\mu^2 \in \{0.5, 0.05\}$):** Up to $N=48$ ($s=2531$, $1.33 \times 10^6$ coupled degrees of freedom), all 46 spin-2 harmonics remain UV-finite and massless ($Z_{\mathrm{wf}}=1.0008$, $m_{\mathrm{grav}}^2 = +0.077$), the quantum area spectrum stays strictly positive and concentrates smoothly onto $\lambda=1$ ($\lambda \in [0.730, 1.390]$), and the intrinsic specific heat remains bounded ($C_V^{\mathrm{int}}/N_{\mathrm{dof}} \in [0.44, 0.51]$).
   - **Fixed-$g_N$ strong-coupling UV scan ($g_N \in \{0.5, 1.0, 1.5, 2.0\}$, $N=8\dots 32$):** Graviton backreaction suppresses traceless horizon anisotropy $Q_2/Q_{2,\mathrm{free}}$ by $-47.8\%$ at $g_N=2.0$ independently of $N$ ($0.5190, 0.5221, 0.5231, 0.5222$), with Gaussian Binder cumulant $U_4(Q_2)=0.02\pm 0.03$ and persistent $+2.5\%$ graviton wave-function dressing ($Z_{\mathrm{wf}}=1.0246$).
   - **Functional measure phase transition ($g=2.5$, $N=8\dots 32$):** Omitting the invariant DeWitt–Fujikawa spin-2 measure triggers a $180\times$ radiative horizon collapse ($\langle\widehat{\rho}^2\rangle \approx 0.006$, $C_V/N_{\mathrm{dof}} \propto N^2$ up to $366.3$), whereas including it stabilizes the smooth geometric phase ($\langle\widehat{\rho}^2\rangle \approx 1.0$, $C_V^{\mathrm{int}}/N_{\mathrm{dof}} \approx 0.46$).

---

## Repository Structure

```text
.
├── paper.tex              # APS Physical Review D (revtex4-2) manuscript (Overleaf-ready)
├── paper.pdf              # Compiled 17-page PDF (zero errors, zero overfull boxes)
├── figs/                  # Six publication vector PDF figures (Fig. 1 – Fig. 6)
├── lean4/                 # 6 standalone Lean 4.36.0 modules (50 verified theorems, 0 sorryAx)
├── verification/          # 29 verification & GPU HMC scripts, 8 plotters, 32 JSON datasets
└── reports/               # 7 technical reports documenting every stage of the study
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

## Compiling the Paper (Local TeX or Overleaf)

This repository is directly importable into **Overleaf** via *New Project → Import from GitHub*. Locally, with TeX Live (`revtex4-2`):

```bash
pdflatex paper.tex && pdflatex paper.tex
```

---

## License

MIT License — see [LICENSE](LICENSE).
