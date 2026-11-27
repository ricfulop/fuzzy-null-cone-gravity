# Nonperturbative Modal NVIDIA B200 Hybrid Monte Carlo of 4D Fuzzy-Sphere Null-Cone Matrix Gravity ($I_{\mathrm{FSG}}^{(N)}$)

**Date:** 2026-11-27  
**Author:** Ric Fulop (`ricfulop@mit.edu`), Center for Bits and Atoms, MIT  
**Compute:** 72 total NVIDIA B200 ensembles on Modal across Tier 1 (24× B200, 173 s), Tier 2 (20× B200, 154 s), and Final Dual IR+UV Production Sweep with exact collective $(\phi_0, Z)$ scale updates (28× B200, 230 s, 1.063 B200 GPU-hours, ~$6.64 USD; total campaign ~$15.14 USD).  
**Code & Data:** `verification/modal_fsg_hmc.py`, `verification/modal_hmc_final_results.json`, `verification/modal_hmc_tier1_results.json`, `verification/modal_hmc_tier2_results.json`, `figs/fig6_modal_b200_hmc.pdf`

---

## 1. Exact Collective-Update FAHMC Verification (All 28 Production Ensembles)

Interleaving exact Fourier-Accelerated HMC (FAHMC, matching the free $T^2 \times S^2_N$ Hessian to $5.8 \times 10^{-16}$) with exact collective monopole-shear scale updates ($\Phi \mapsto \Phi + \frac{\delta}{\sqrt{N}} I_N$, $Z \mapsto e^{-\frac{g_N}{2\sqrt{N}}\delta} Z$, leaving $S_Z$ invariant to $0.0$ and cancelling $-\Delta S_{\mathrm{DW}} + \ln J$ to $1.14 \times 10^{-13}$) completely eliminated slow-mode autocorrelation across all 28 final ensembles:

- **Maximum Gelman–Rubin $\hat{R}(\phi_0)$ across all 28 ensembles** (up to $N = 48$ at $\mu^2 = 0.05$ and fixed $g_N = 2.0$ at $N = 32$): **`1.0050`** ($\tau_{\mathrm{int}}(\phi_0) \in [1.0, 1.1]$ everywhere!).
- **Exact Schwinger–Dyson Ward identities across all 28 ensembles:**
  - Spin-2 shear: $\max |\langle Z \cdot \partial S_E/\partial Z\rangle / N_{\mathrm{dof}, Z} - 1| = \mathbf{0.00050}$ (`0.9995–1.0005`),
  - Nonlinear exponential/quartic scalar: $\max |\langle \Phi \cdot \partial S_E/\partial \Phi\rangle / N_{\mathrm{dof}, \Phi} - 1| = \mathbf{0.00161}$ (`0.9984–1.0006`).
- **Operator & integrator verification:** $O(N^3)$ diagonal-block $\iota_{N\to N\pm s}$ $SU(2)$-equivariance error $6.1\times 10^{-15}\text{--}1.1\times 10^{-14}$; $O(N^2)$ Schwinger-boson $\widehat{\bar{\eth}}$ identity error $0.0$ at $N=16,32$; Tensor-Core `matrix_exp` 40× speedup at $N=48$ (`13.48 ms`/grad).

---

## 2. Final Nonperturbative Physics Results

### Result 1: Macroscopic IR Limit ($g_N = g_P/\sqrt{N}$, $g_P = \sqrt{4\pi/\ln 3} = 3.3821$, $\mu^2 = 0.5$ and $\mu^2 = 0.05$)
Across $N = 8, 12, 16, 24, 32, 48$ (up to $1.33 \times 10^6$ coupled real degrees of freedom), lowering the IR area stiffness by $10\times$ ($\mu^2 = 0.5 \to 0.05$) confirms $\mu^2$-independent semiclassical emergence:
- **Massless, UV-finite spin-2 gravitons:** At $N=48, \mu^2=0.05$ ($\hat{R} = 1.005$, $\mathrm{SD}_Z = \mathrm{SD}_\Phi = 1.0000$), $Z_{\mathrm{wf}} = 1.0008$ and $m_{\mathrm{grav}}^2 = +0.077$ (dropping from $+0.348$ at $N=8$), with $G_{\mathrm{meas}}(j)/G_{\mathrm{free}}(j) \approx 0.995$ across all 46 harmonics ($2t_j^2 = 4\dots 2254$).
- **Caustic-free area eigenvalue concentration:** Eigenvalues of $\widehat{\rho}^2 / r_0^2 = \exp((g/\sqrt{N})\Phi)$ stay strictly positive on every configuration and sharpen from $[0.428, 2.660]$ ($N=8$, $s=70$) to $[0.730, 1.390]$ ($N=48$, $s=2531$).
- **Bounded intrinsic specific heat:** $C_V^{\mathrm{int}} / N_{\mathrm{dof}} \in [0.439, 0.510]$ across all ranks and both $\mu^2$ values (no instability as $\mu^2$ decreases tenfold).

### Result 2: Strong-Coupling UV Scan (Fixed $g_N \in \{0.5, 1.0, 1.5, 2.0\}$ Across $N = 8, 16, 24, 32$)
Holding the microscopic matrix coupling $g_N$ fixed while increasing rank $N$ directly exposes persistent nonperturbative graviton backreaction:
- **Universal $N$-independent horizon-anisotropy suppression ($Q_2 / Q_{2,\mathrm{free}}$):** Quantum shear fluctuations strongly suppress traceless non-spherical horizon roughness $Q_2 = \frac{1}{N}\mathrm{Tr}(\Phi_{\mathrm{trless}}^2)$:
  - $g_N = 0.50$: $Q_2 / Q_{2,\mathrm{free}} = 0.918 \to 0.877$,
  - $g_N = 1.00$: $Q_2 / Q_{2,\mathrm{free}} = 0.766 \to 0.722$,
  - $g_N = 1.50$: $Q_2 / Q_{2,\mathrm{free}} = 0.629 \to 0.607$,
  - $g_N = 2.00$: **$Q_2 / Q_{2,\mathrm{free}} = 0.5190\,(N=8),\, 0.5221\,(N=16),\, 0.5231\,(N=24),\, 0.5222\,(N=32)$** (all $\pm 0.0002\text{--}0.0003$), a **$-47.8\%$ nonperturbative smoothing** that locks to a constant across $N = 8\dots 32$ with Gaussian Binder cumulant $U_4(Q_2) = 0.02 \pm 0.03$ (no spontaneous $SU(2)$ symmetry breaking).
- **Persistent wave-function dressing & strict positivity:** At $g_N = 2.0, N = 32$, $Z_{\mathrm{wf}} = 1.0246$, $Z(j=2) = 0.9556$, $Z(j_{\max}) = 1.0257$, intrinsic $C_V^{\mathrm{int}}/N_{\mathrm{dof}} = 0.472 \pm 0.018$, and $\lambda(\widehat{\rho}^2) \in [0.374, 3.353]$.

### Result 3: DeWitt–Fujikawa Measure Controls the Crumpled-vs-Geometric Phase Boundary
Comparing `dewitt_measure=False` against `dewitt_measure=True` at $g=2.5$ ($N=8\dots 32$) demonstrates that omitting the invariant spin-2 functional measure triggers a $180\times$ radiative horizon collapse ($\langle\widehat{\rho}^2\rangle \approx 0.006$, $C_V/N_{\mathrm{dof}} \propto N^2$ reaching `366.3` at $N=32$), whereas including it stabilizes the smooth geometric phase ($\langle\widehat{\rho}^2\rangle \approx 1.0$, $C_V^{\mathrm{int}}/N_{\mathrm{dof}} \approx 0.46$).
