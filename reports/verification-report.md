# Technical Verification Report: `arXiv:2611.34924v1`
**Paper:** Venkatesa Chandrasekaran, *"An algebraic quantization of gravity from the null cone"* (Stanford Leinweber Institute for Theoretical Physics, Nov 25, 2026)  
**Audit Date:** November 27, 2026

---

## Executive Verdict

We executed all five independent verification tests (exact counterexamples, Williamson symplectic diagonalization, operator-algebra / Hamiltonian singularity checks, microstate & squashed-cone symbolic audits, and Z3 SMT proof-graph mechanization).

* **What passes:** The classical first-order null-cone phase space and local geometric identities (eqs. 3.26a, 4.21a–c, 5.67b, 5.76, D.52a), the internal moment consistency of the thermal clock distribution (Footnotes 23–24), the vanishing of $k^2$ on a marginally trapped surface in the Fursaev–Patrushev–Solodukhin squashed-cone formula (eq. 5.67b), and the conditional Lorentzian-to-Euclidean axiomatic equivalence in Appendix E.2 (assuming Axioms 1–5 hold).
* **What fails:** Every headline **nonperturbative** claim fails under exact verification:
  1. **Test 1 (FAIL — Circular + Counterexamples):** Both Section 5.2 finiteness proofs ($\dim \mathcal{H}_{\mathscr{D}} < \infty$ and AdS $Z_\beta < \infty$) misapply open-amplitude Cauchy–Schwarz (eq. 4.30) to non-trace-class closed traces and are refuted by explicit infinite-dimensional counterexamples (`C^2 (x) l^2(N)` and bounded-dilaton JT + massless scalar).
  2. **Test 2 (FAIL — Continuum Collapse $P_{\mathscr{D}} = 0$):** Exact Williamson symplectic diagonalization of the un-normal-ordered clock focusing operator $\widehat{\Xi}_y(v)$ (eq. 4.22) proves its bottom-of-spectrum zero-point energy diverges as $E_{0,\min} \sim 14.64\, N^2 J^2 / s$, collapsing $P_{\mathscr{D}} = 0$ identically past the Planck mode cutoff $N_\star \sim (s \ln(\ell/u_{\mathrm{cl}})/14.64)^{1/4}$.
  3. **Test 3 (FAIL — Unproved $U_{\mathcal{M}}$, Singular $H_c$, Toy 1D Check):** The same zero-point blowup drives the inverse lapse $1/N(\epsilon,y) = N(1,y)^{-1}\exp(\widehat{\Xi}_y(1)) \to +\infty$ inside the reduced Hamiltonian $H_c$ (eqs. 3.54–3.55); Kadison–Singer pure-state extensions leave area-changing transition amplitudes $\omega(U_{\mathcal{M}}) \in [-1, 1]$ completely undetermined; and Section 4.4's $\mathrm{Li}_2(e^{-2i/s})/(2s^2)$ calculation is a free 1D oscillator sum that becomes quadratically UV-divergent ($\propto J^2$) at finite $G_N$ once 4D angular harmonics are restored.
  4. **Test 4 (FAIL on Microscopic Claim; Nuanced Pass on MTS Algebra):** Perturbative graviton states in Appendix D.3 scale as $\log d(s) \approx 2.801\,(s Q_*)^{3/4} \ll s$ (an exponentially vanishing fraction of $e^s$), while the $+106/45\log s_b$ coefficient is imported from the covariant 1-loop heat kernel of Bastianelli–Bonezzi–Melis (`2206.13287`) with an $O(1)$-broad power-law ensemble ($\mathrm{Std}(\mathscr{A})/\langle\mathscr{A}\rangle = 26.9\%$). On a marginally trapped surface ($\Theta_{(\ell)} = 0$), the squashed-cone extrinsic curvature $k^2 = 0$ vanishes even with nonzero shear, but on a generic causal diamond corner (e.g., Minkowski sphere) the true 1-loop heat kernel and the Gauss–Bonnet proxy disagree by $-112\pi^2/5$, and the conical defect is off-shell where $R_{ab}^2$ is gauge-dependent.
  5. **Test 5 (MECHANIZED — Z3 SMT & Lean 4 Blueprint):** Z3 5.4.1 proves all six formal statements: non-entailment of $\dim\mathcal{H}_{\mathscr{D}} < \infty$ without `ClosedTraceConv`, redundancy of boosts/cycles once `ClosedTraceConv` is added, exact circularity of Section 5.2 $\leftrightarrow$ Appendix E.3, `UNSAT` of the continuum clock triad, and conditional validity of Appendix E.2.

![Verification summary figure](verification/verification_summary.svg)

---

## Detailed Ledger of Tests 1–5

### Test 1 — Section 5.2 Finiteness Proofs (`test1_section52_counterexamples.py`)

* **Test 1A ($\mathcal{H} = \mathbb{C}^2 \otimes \ell^2(\mathbb{N})$):** With $s \in \{0, 2\pi\}$, $C_b = \mathrm{diag}(1, e^{-\pi}) \otimes \hat{1}$ (satisfying eq. 5.30c: $e^{-2\pi}\hat{1} \le C_b^\dagger C_b \le \hat{1}$), $\eta = \pi$, and bounded observable $\|\widehat{\mathcal{O}}\| = 1$ having $\delta_\eta = 0.172533 > 0$:
  * Product exhaustion $P_M$: $Z_{P_M}[X] = +0.086428\,M \to +\infty$ ($M = 50 \implies +4.3214$).
  * Cycle exhaustion $P_{\mathrm{cyc}}(K, r=200)$: $\langle\Delta_\pi \widehat{\mathcal{O}}\rangle_P \le 0.00087 \ll \delta_\pi/2$ and $Z_P[X] \le -16.9944\,K \to -\infty$ ($K = 50 \implies -849.72$).
* **Test 1B (Continuous-spectrum JT + massless scalar analog):** Discretizing $s \in [0, 10]$ with $w(s) = e^{-s} f_b(s)$ and $\eta = 3.0$ gives $\delta_\eta = 0.262724$, $\Delta s = 0.1$, and cycle length $r = 488{,}620$ meeting eq. (5.32) ($\le \delta_\eta/2$) while $\dim \mathcal{H} = \infty$.
* **Test 1C (AdS boundary eqs. 5.44–5.45):** For $w_i = 1/(i+1)$ ($\sum w_i = \infty$) and a 2-state diamond ($p = 0.4$), balanced subsets give $S_{F_k} = 0$ identically, whereas adding $9{,}000$ unmatched indices yields $\pm 1.0814$. By unconditional convergence, $\lim_{F} S_F = 0 \iff \sum_i w_i < \infty$.
* **Test 1D (Open-state salvage fails):** Replacing the closed trace by a finite-norm vector $|\Phi\rangle \in \mathcal{H}\otimes\mathcal{H}^*$ (so eq. 4.29 holds) destroys eq. (5.30c) and eq. (D.25): with weight $1/2$ on $|1,1\rangle$ and $\frac{1}{2}2^{-n}$ on $|0,n\rangle$, $\mathrm{Var}_P(\hat{s}) \ge \frac{1}{2}s_*^2 = 19.7392$ for all $r \le 1000$.

### Test 2 — Continuum Clock Zero-Point Collapse (`test2_continuum_clock_uv.py`)

Using Williamson's theorem (exact symplectic diagonalization) on the paper's reduced symplectic form (eqs. 3.42b–c, A.14) and focusing operator (eq. 4.22) for one polarization $z^1$ on $v \in [0.1, 1]$:

| Radial cutoff $N$ | $\omega_{\min}$ | $\omega_{\max}$ | $s\,E_0^{(1\text{ harm.})}$ | $s\,E_0 / N^2$ |
|---:|---:|---:|---:|---:|
| 5 | 35.254 | 487.9 | 458.1 | 18.325 |
| 10 | 35.205 | 1,244.0 | 1,682.7 | 16.827 |
| 20 | 35.198 | 3,021.2 | 6,311.1 | 15.778 |
| 40 | 35.198 | 7,026.4 | 24,128.8 | 15.080 |
| 80 | 35.197 | 15,823.8 | 93,717.0 | 14.643 |

Summing over angular harmonics $j = 2 \dots J$ on a macroscopic diamond ($b = 10^5\ell_P$, $s = \pi \times 10^{10}$, clock cap $\ell(\mathscr{D}) = 10\,u_{\mathrm{cl}}$):

| Cutoff $N=J$ | Total modes | $E_{0,\min}$ (nats) | $u_{\min}/u_{\mathrm{cl}}$ | $P_{\mathscr{D}} \equiv 0$? | $\langle\Omega|\Psi_L\rangle$ upper bound | $\log_{10}\mathrm{Pr}_W$ (C.16d) |
|---:|---:|---:|---:|:---:|---:|---:|
| 10 | 1,170 | $6.27\times 10^{-6}$ | $1.000$ | No | $0.990$ | $0.0$ |
| 100 | $1.02\times 10^6$ | $4.75\times 10^{-2}$ | $1.049$ | No | $0.990$ | $0.0$ |
| 250 | $1.57\times 10^7$ | $1.835$ | $6.267$ | No | $0.941$ | $0.0$ |
| **500** | **$1.25\times 10^8$** | **$29.25$** | **$5.04\times 10^{12}$** | **Yes** | **$1.98\times 10^{-11}$** | **$-1.01\times 10^7$** |
| **1,000** | **$1.00\times 10^9$** | **$467.0$** | **$6.81\times 10^{202}$** | **Yes** | **$1.47\times 10^{-201}$** | **$-1.58\times 10^8$** |

* **Normal-ordering fails (Test 2D):** In a squeezed state $|r\rangle$ ($r = 1.0$), $\langle r|:\!(\hat{a}+\hat{a}^\dagger)^2\!:|r\rangle = e^{-2r}-1 = -0.8647$ (and $-4366.6$ across $K=100$ modes), destroying $\widehat{\Xi}_y \ge 0$ and the area upper bound eq. (5.21d).

### Test 3 — Interacting Algebra, $H_c$, and Section 4.4 Pair Creation (`test3_interacting_algebra_and_Hc.py`)

* **Test 3A (Section 4.4):** At $s = 1000$, $2s^2\mathrm{Re}(S) = 1.64179 \to \pi^2/6 = 1.64493$, and $\mathrm{Im}(S) = -7.2146\times 10^{-9}$ matches the asymptotic $-i(\log s + 1 - \log 2)/s^3$ to $7.0\times 10^{-8}$. Restoring the $2((J+1)^2 - 4)$ angular modes at fixed $s=100$ drives the amplitude to $1.62\times 10^4$ by $J = 10^4$ (quadratic UV divergence at finite $G_N$).
* **Test 3B ($H_c$ lapse singularity):** In eqs. (3.54)–(3.55), $1/N(\epsilon, y) = N(1, y)^{-1}\exp(\widehat{\Xi}_y(1))$ has operator lower bound $\ge 5.04\times 10^{12}$ at $N=J=500$ and $\ge 6.81\times 10^{202}$ at $N=J=1000$, violating Appendix C.2's bounded-generator assumption (eq. C.7).
* **Test 3C (Pure-state extension non-uniqueness):** On $\mathcal{A} = M_2(\mathbb{C}) \otimes M_2(\mathbb{C}) \supset \mathcal{A}_{\mathrm{in}} = M_2(\mathbb{C}) \otimes I_2$, pure extensions $\omega_\theta$ agree with $|\Omega_{\mathscr{A}}\rangle$ on $\mathcal{A}_{\mathrm{in}}$ to machine precision ($< 10^{-16}$) yet assign arbitrary values $\omega_\theta(U_{\mathcal{M}}) = \sin(2\theta) \in [-1, 1]$ to an area-changing unitary $U_{\mathcal{M}}$.

### Test 4 — Entropy Scaling, Ensemble Spread, and Squashed Cones (`test4_entropy_and_squashed_cone.py`)

* **Test 4A (Microscopic state count):** The exact Bose–Einstein count for eq. (D.45a) is $\log d(s) = \frac{4}{3}(\pi^4/5)^{1/4}(s Q_*)^{3/4} \approx 2.801\,(s Q_*)^{3/4}$, which is $0.89\%$ of $S_{\mathrm{BH}}$ at $s = 10^{10}$ and $0.0028\%$ at $s = 10^{20}$.
* **Test 4B (Ensemble spread):** With $\alpha = 257/90$, $q = 257/302$, Footnote 23 eq. (5.61) is exact, $\langle s/s_* \rangle = 257/347 \approx 0.7406$, fractional area std $= 26.86\%$, and fractional horizon-radius std $= 15.07\%$.
* **Test 4C (Squashed cones & interior cap):** Applying Fursaev–Patrushev–Solodukhin (`1306.4000`, eqs. 1.4–1.6) to Bastianelli–Bonezzi–Melis's (`2206.13287`, eq. 36) true 1-loop heat kernel $\mathrm{tr}[a_{2,\mathrm{tot}}] = \frac{53}{45} R_{abcd}^2 - \frac{361}{90} R_{ab}^2 + \frac{43}{36} R^2$ gives surface density $\frac{2\pi}{45}(212 R[q] + 63 \Theta_{(\ell)}\Theta_{(n)})$, vs the Gauss–Bonnet proxy $\frac{424\pi}{45} R[q]$ (difference $\frac{14\pi}{5}\Theta_{(\ell)}\Theta_{(n)}$):
  * On a marginally trapped surface ($\Theta_{(\ell)} = 0$), the difference is identically **$0$** even with nonzero shear.
  * On a round Minkowski sphere of radius $b$, the true BBM and Gauss–Bonnet integrals disagree by **$-112\pi^2/5$** ($2384\pi^2/45$ vs $3392\pi^2/45$).
  * The Schwarzschild interior Kretschmann integral diverges as $16\pi(\epsilon^{-3} - 1)$ as the cap $\Sigma_1$ approaches $r=0$.

### Test 5 — Formal SMT Mechanization & Kernel-Verified Lean 4.36.0 Proofs (`test5_formal_dependency_and_z3.py`, `lean4_verified_audit.lean`)

1. **Z3 5.4.1 (6/6 checks passed):**
   * `5A_sec52_without_trace_class`: **SAT** (countermodel without `ClosedTraceConv`).
   * `5A_sec52_with_trace_class_is_redundant`: **UNSAT** (once `ClosedTraceConv` is added, boosts/cycles are unnecessary).
   * `5B_circular_loop_5_2_and_E_3`: **SAT** (all-False model satisfies the circular chain `ClosedTraceConv -> DimHD_Finite -> Axiom5 -> ClosedTraceConv`).
   * `5C_clock_incompatibility_triad`: **UNSAT** (`ContinuumCCR` + `UnNormalOrderedPos` + `NonZeroClockSubspace` are mutually inconsistent).
   * `5D_appendix_E2_conditional_valid`: **UNSAT** of negation (Appendix E.2 is valid conditionally on Axioms 1–5).
   * `5E_real_arithmetic_trace_bound`: **UNSAT** of negation ($c > 0 \land c d \le T \implies d \le T/c$).
2. **Lean 4.36.0 (`~/.elan/bin/lean lean4_verified_audit.lean` — exit code 0, zero warnings, zero `sorryAx`):**
   * Certified all 11 core theorems directly in Lean 4's kernel (standard axioms only; `appendix_E2_valid_conditionally` uses zero axioms):
     * `uniform_lower_bound_finite_trace_implies_dim_bound` and `uniform_lower_bound_infinite_trace_unbounded`
     * `ads_544_unconditional_bounded_iff_weights_bounded` (exact `↔` for eqs. 5.44–5.45)
     * `cycle_partition_left_inv` & `cycle_partition_right_inv` (exact bijection partitioning `Bool × ℕ` into length-100 cycles)
     * `section52_counterexample_product_diverges_pos` ($Z_{P_M}[X] = 40M \to +\infty$) & `section52_counterexample_cycle_diverges_neg` ($401\,Z_{P_{\mathrm{cyc}}(K,100)}[X] \le -1{,}551{,}960\,K \to -\infty$)
     * `clock_continuum_incompatibility_triad`
     * `sec52_does_not_entail_finite_dim`, `circular_loop_admits_all_false_model`, and `appendix_E2_valid_conditionally`

### Test 6 — Exact Symbolic Audit of Appendices A & B (`test6a_appendix_A_audit.py`, `test6b_appendix_B_audit.py`, `test6b_appendix_B_extra.py`)

We constructed the full 4D metric (B.1), Christoffel symbols, Ricci tensor, scalar curvature, Green's functions, Lie derivatives, and Iyer–Wald/Komar charges from scratch in SymPy to verify whether the classical foundation has sign errors or missing boundary terms:

| Appendix / Equations | Statement Audited | Symbolic Residual | Result |
| :--- | :--- | :---: | :---: |
| **Sec. 3.3 / App. A:** `3.39 → 3.42a–d` | Smooth null-cone symplectic reduction ($\alpha$ coeff $1/4$, $-2\mathcal{K}(y)\delta J$ with $J$ coeff $1/16$) | `0` | **PASS** |
| **App. A.1:** `A.1` | Connection components $\beta_v = \partial_v\log N$, $\beta_A = w_A - \frac{1}{2}\mathring{\nabla}_A f$ | `0` | **PASS** |
| **App. A.2:** `A.23d–A.27` | Linearized Raychaudhuri ODE $L[h_1]=L[h_2]=0$ and Green's function jump $\partial_v G|_{v=w} = 1$ | `0` / `1` | **PASS** |
| **App. A.3:** `A.30b–A.33c` | $V$-coordinate constraint $\mathfrak{C}_v = (\partial_v V)^2 \mathcal{D}_{\varrho^2}(\varrho^2)$ and Green's function $G_{\varrho^2}$ | `0` / `1` | **PASS** |
| **App. A.3:** `A.34–A.36` | Second-class constraint operator $\Delta_{\mathfrak{C}}$ and inverse kernel ($+\delta(v-w)$ sign) | `0` | **PASS** |
| **App. A.3:** `A.37–A.41` | Coupled symplectic operator $\mathcal{D}_z U$ triangular reduction to $(y'', w, b'')$ | `0, 0, 0` | **PASS** |
| **App. B.1:** `B.3` | Null Raychaudhuri constraint $R_{vv}=0 \implies \partial_v f = -\sigma_{AB}\sigma^{AB}/\Theta_{(\ell)}$ | `0` | **PASS** |
| **App. B.1:** `B.4a` | Damour–Navier–Stokes angular constraint $\partial_v \Pi_A$ ($R_{ry}=0$) | `0` | **PASS** |
| **App. B.1:** `B.5a`, `B.5c` | Cross-focusing ($\frac{1}{2}q^{AB}R_{AB}=\Lambda$) and transverse connection $\partial_v\mu_c$ | `0`, `0` | **PASS** |
| **App. B.1:** `B.7` | Bulk Hamiltonian divergence identity reducing eq. (3.47) to pure surface charge (3.54) | `0` | **PASS** |
| **App. B.2:** `B.8d–f` | Displaced null boundary $s_a s^a$, null generator $m^a$, and exact inaffinity $\kappa_m$ | `0` | **PASS** |
| **App. B.2:** `B.10–B.11e` | Full 4D Einstein–Hilbert decomposition $\sqrt{-g}(R-2\Lambda) = L_0 + \Delta L_0 + \partial_a L^a$ | `0` | **PASS** |
| **App. B.2:** `B.12` | Boundary divergence $s_a L^a$ on displaced null surface $r = b(t, y)$ | `0` | **PASS** |
| **App. B.2:** `B.17` | Iyer–Wald / Komar Noether charge $\mathbf{Q}_\chi|_C$ on displaced cut $r = b(t, y)$ | `0` | **PASS** |

**Takeaway from Test 6:** The classical general-relativity foundation (Section 3 and Appendices A–B) is algebraically spotless. Every failure identified in this audit lies cleanly at the **nonperturbative quantization leap** (Sections 4–5 and Appendices C–E).

---

## Reproducible Verification Artifacts

* `verification/test1_section52_counterexamples.py` → `verification/test1_results.json`
* `verification/test2_continuum_clock_uv.py` → `verification/test2_results.json`
* `verification/test3_interacting_algebra_and_Hc.py` → `verification/test3_results.json`
* `verification/test4_entropy_and_squashed_cone.py` → `verification/test4_results.json`
* `verification/test5_formal_dependency_and_z3.py` → `verification/test5_z3_results.json`
* `verification/lean4_verified_audit.lean` (compiled cleanly with Lean 4.36.0, exit code 0, zero `sorryAx`)
* `verification/test6a_appendix_A_audit.py` → `verification/test6a_appendix_A_results.json`
* `verification/test6b_appendix_B_audit.py` → `verification/test6b_appendix_B_results.json`
* `verification/test6b_appendix_B_extra.py` → `verification/test6b_appendix_B_extra_results.json`
* `verification/plot_verification_summary.py` → `verification/verification_summary.svg`, `verification/verification_summary.pdf`
