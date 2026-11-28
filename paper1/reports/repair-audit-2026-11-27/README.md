# Action and clock repair calculations

**Date:** 27 November 2026. **Status:** implemented repair supplement. The manuscript, clock figures, and reproducibility scripts are revised together on the repair branch, based on citation-fix PR #5 (main commit `49de5fa`). The previous action is withdrawn. A new Einstein--Hilbert angular Galerkin bulk prescription is defined below; its boundary quantization, finite-rank constraint closure, and equivalence to the sampled matrix theory remain open.

## What worked, and what did not

| Question | Result |
|---|---|
| Does a minimal correction of Appendix A yield the Einstein equations? | No. Even its spherical, commuting sector fails the bulk variational test. |
| Can the spherical sector be repaired from Einstein-Hilbert? | Yes. The covariant spherical reduction gives Schwarzschild-de Sitter, including the mass-conservation constraint when variation precedes gauge fixing. |
| Can the axial and polar gravitational equations be obtained by eliminating constraints? | Yes, from the correct Einstein quadratic actions. The symbolic residuals against Regge-Wheeler and Zerilli are exactly zero for general angular momentum. |
| Can those quadratic actions be represented in the finite matrix bundles? | Yes, by an isometric harmonic map. Numerical bilinear checks at N=4,6,8 agree within 5.8e-12. This is a quadratic construction, not a derivation from the withdrawn nonlinear Appendix A action. |
| Can the clock lower bound be repaired without the disputed scalar-QEI attribution? | Yes. An explicit Fourier/trace-norm argument gives a mass- and radial-cutoff-independent lower bound for the specified oscillator model. |
| Does the true occupation budget ruin the favorable vacuum-clock estimate? | No. At s=5000 and B=2500 the recalculated second-moment bound gives vacuum acceptance at least 0.999998907949 for a threshold twice the classical clock. |
| Has the full nonlinear matrix theory been shown equivalent to Einstein gravity? | No. A definite nonlinear Einstein--Hilbert Galerkin bulk replacement is now supplied. Its finite-rank constraint algebra, null-boundary quantization, and relation to the simulated Euclidean model still need a common derivation. |

The constructive results are a nonlinear Galerkin bulk replacement, exact Einstein quadratic sectors in finite matrix coordinates, and a new conditional finiteness argument. The calculation does not establish a continuum quantum-gravity theory or a new derivation of the known black-hole potentials.

## 1. Spherical action: a decisive variational test

Use X=rho^2 for the area field, and a two-dimensional metric

$$
h_{ab}dx^a dx^b=-F\,dt^2-2N\,dt\,dr.
$$

Angular integration of the four-dimensional Einstein-Hilbert action gives, up to boundary terms,

$$
I_{\rm sph}=\frac{1}{4G}\int d^2x\,\sqrt{-h}
\left[X R[h]+\frac{(\nabla X)^2}{2X}+2-2\Lambda X\right].
$$

This expression must be varied covariantly before fixing the null gauge if all metric equations are to be retained. A direct Christoffel/Ricci calculation gives

$$
R[h]=-\frac1N\partial_r\left(\frac{F_r-2N_t}{N}\right).
$$

After one radial integration by parts the gauge-fixed bulk density is

$$
\mathcal L_{\rm sph}
=\frac{X_r(F_r-2N_t)}{N}
-\frac{X_tX_r}{X}
+\frac{F X_r^2}{2NX}
+2N-2\Lambda NX.
$$

The script derives all three Euler-Lagrange expressions and obtains exactly zero on

$$
X=r^2,\qquad N=1,\qquad F=1-\frac{2M}{r}-\frac{\Lambda r^2}{3}.
$$

At N=1 this density agrees exactly with the scalar part of **Chandrasekaran's Appendix B, Eqs. B.11a-b**. This check uses the original manuscript, not our Appendix A.

### Why changing a few factors in Appendix A is insufficient

Take the suggested correction rho^2 N R_q -> 2N, but leave the other spherical terms of Eq. A1 intact. With

$$
\vartheta=X_t/X-FX_r/(2X),
$$

its bulk density, modulo a boundary derivative, becomes

$$
\mathcal L_{\rm patch}=-2X_r\vartheta+2N-2\Lambda NX.
$$

On the proposed Schwarzschild-de Sitter background its variational residuals include

$$
\frac{\delta I_{\rm patch}}{\delta F}\propto4,\qquad
\frac{\delta I_{\rm patch}}{\delta N}\propto2-2\Lambda r^2.
$$

They do not vanish. Boundary terms cannot remove these interior residuals for compactly supported variations. The correct spherical action needs the radial-gradient and lapse terms displayed above.

There is a second useful check: the three gauge-fixed equations also admit M=M(t). The metric equation removed by premature null gauge fixing is recovered from the covariant action:

$$
\mathcal E_{ab}
=h_{ab}\Box X-\nabla_a\nabla_bX+
\frac{\partial_aX\partial_bX}{2X}
-\frac{h_{ab}(\nabla X)^2}{4X}
-h_{ab}(1-\Lambda X)=0.
$$

For the time-dependent-mass ansatz, the remaining component is

$$
\mathcal E_{tt}=-2\,\partial_t M.
$$

Thus the vacuum mass is constant. Merely showing the three gauge-fixed residuals vanish would have missed this constraint.

**Reproduction:** action_derivation.py; results in action_results.json.


## 1b. Nonlinear bulk reconstruction implemented in this revision

The failed matrix-product action cannot be repaired by asserting that its quadratic benchmarks imply a nonlinear completion. We instead specify a different regulator: the ordinary Einstein--Hilbert action pulled back to finitely many angular coefficient fields.

1. Expand the three components of a general two-dimensional metric $h_{ab}$ and the angular trace $X$ in scalar harmonics through $j=N-1$.
2. Expand two real angular shifts $A_a^A$ in spin-one harmonics and a real tracefree angular tensor $C_{AB}$ in spin-two harmonics. Their coefficients are stored in the existing rectangular bundles by the isometric map. Keep all ten real metric components before any gauge fixing.
3. Reconstruct

$$
ds_N^2=h_{ab}dx^a dx^b+(Xq_{AB}+C_{AB})(dy^A+A_a^A dx^a)(dy^B+A_b^B dx^b).
$$

Require $h$ to be Lorentzian and $Xq+C$ positive definite pointwise. Matrix positivity alone is insufficient for this domain restriction. Define

$$
I_{G,N}^{\mathrm{bulk}}=\frac1{16\pi G}\int d^4x\sqrt{-g_N}\,[R(g_N)-2\Lambda].
$$

All products, inverses, and determinants are pointwise operations on the reconstructed metric. Do not project their higher angular harmonics away before integration. This specifies nonlinear interactions at every order, although a numerical angular integral requires its own convergence checks. It is a finite angular field approximation, not a finite-dimensional many-body Hilbert space; the separate occupation restriction is still necessary.

For compactly supported variations, the coefficient equations are the Einstein residual contracted with the derivative of the reconstructed metric and integrated over the sphere. Every representable Einstein solution is stationary. Vanishing projected equations does not establish that the omitted harmonics vanish. The spherical reduction and missing mass equation above now follow from this bulk action. Around a spherical Einstein solution, the first variation vanishes and rotational invariance separates angular multiplets. Its Hessian is therefore the Einstein Hessian restricted to the retained complete metric multiplets, leading to the two derived master actions below.

**An independent nonlinear test:** the new `nonlinear_galerkin.py` directly builds Christoffel symbols and the Ricci tensor of $g_N=A(\theta)\bar g_{\mathrm{Schw}}$, with $A=1+\epsilon\sqrt5 P_2(\cos\theta)>0$. The metric coefficients fit every $N\ge3$. Without inserting a conformal-transformation formula, it derives

$$
R[g_N]=\frac1{r^2}\left[-\frac{3\Delta_\Omega A}{A^2}+\frac{3|\nabla_\Omega A|^2}{2A^3}\right].
$$

The integrated density, suppressing the overall factor $1/(4G)$, is

$$
\mathcal L_{\mathrm{ang}}=\frac32\left\langle\frac{|\nabla_\Omega A|^2}{A}\right\rangle
=9\epsilon^2-\frac{9\sqrt5}{7}\epsilon^3+O(\epsilon^4).
$$

The exact integral and its first variation agree with the direct Ricci calculation and Einstein trace at finite amplitudes $\epsilon=-0.2,0.1,0.3$. Matrix-to-symbol reconstruction agrees at $N=3,4,6,8$. A deliberately distinct fuzzy-product evaluation is recorded to demonstrate that nonlinear trace substitution is a different action. This is a selected nonlinear-family test, not a check of all constraints.

**Remaining scope:** finite angular truncation is not closed under all diffeomorphisms. A first-class constraint algebra at finite rank, a complete causal-diamond boundary and joint variational principle, boundary quantization, and a common measure with the Euclidean simulations have not been derived. The replacement defines and tests nonlinear bulk dynamics; it does not establish a nonlinear quantum-gravity completion. The older cubic vertex, projected interaction example, and 1,260-state gluing Hamiltonian remain separate specified models.

Reproduction: `nonlinear_galerkin.py`; results in `nonlinear_galerkin_results.json`. The complete updated derivation is in `sections/nonlinear_action.tex`.

## 2. Constraint elimination really produces the black-hole potentials

These are reproductions of established Einstein perturbation theory, with the algebra checked here. They replace the previous inference from a spectrum calculation that inserted the potentials.

### Axial sector

Let f=1-2M/r, L=(j-1)(j+2), and j>=2. Up to an overall positive harmonic normalization, the Einstein quadratic action in Regge-Wheeler gauge has density

$$
\mathcal L_-=\frac12\left[
\left(\dot h_1-h_0'+\frac{2h_0}{r}\right)^2+
\frac{L}{r^2f}h_0^2-\frac{Lf}{r^2}h_1^2
\right].
$$

Vary h_0 and h_1. A differential combination of the resulting equations gives the constraint

$$
\dot h_0=f(fh_1)'.
$$

Eliminating it and setting Psi_-=fh_1/r yields

$$
\left[-\partial_t^2+\partial_{r_*}^2-
f\left(\frac{j(j+1)}{r^2}-\frac{6M}{r^3}\right)\right]\Psi_-=0,
\qquad \partial_{r_*}=f\partial_r.
$$

The gauge-completed amplitudes

$$
k_0=h_0-\tfrac12\dot h_2,\qquad
k_1=h_1-\tfrac12h_2'+h_2/r
$$

are verified invariant under the standard odd-parity linear gauge transformations. Thus the quadratic result is not just a fit to a target potential.

### Polar sector

The input is the known Einstein quadratic density in alpha gauge, **Solomon (2023), Eq. (150)**. The script integrates out its auxiliary field, changes field variables, solves the remaining constraint, and integrates the mixed derivatives by parts. It does not start with the Zerilli potential.

With q=j(j+1), the resulting density in an intermediate normalization is

$$
\mathcal L_+=4q(q-2)\left[
\frac{\dot\psi^2}{f}-f(\psi')^2-\frac{V_Z}{f}\psi^2
\right],
$$

where the coefficient obtained by the elimination is

$$
V_Z=\frac{2f\left[\lambda^2(\lambda+1)r^3+
3\lambda^2Mr^2+9\lambda M^2r+9M^3\right]}
{r^3(\lambda r+3M)^2},\qquad \lambda=(j-1)(j+2)/2.
$$

The kinetic coefficients and the difference from the standard potential simplify to zero residuals for symbolic j (represented by q). This checks the constrained reduction for all j>=2, rather than only a numerical j=2 example.

**Limit:** These calculations start from Einstein-Hilbert perturbation theory. They do not prove that the withdrawn nonlinear matrix action has the same Hessian.

**Reproduction:** action_derivation.py and polar_derivation.py; action_results.json and polar_results.json.

## 3. A well-defined finite-matrix realization at quadratic order

For a spin-s rectangular matrix, use

$$
\langle X,Y\rangle_{s,N}
=\frac{\operatorname{Tr}(X^\dagger Y)}{\sqrt{N^2-s^2}}.
$$

The manuscript's harmonic map iota preserves this scalar normalized-trace pairing on the common harmonic subspace. A definite, type-correct scalar action on rectangular matrices is

$$
F\triangleright_s X
=\frac12\left[\iota_{N\to N+s}(F)X+
X\iota_{N\to N-s}(F)\right].
$$

For Hermitian F this map is self-adjoint in the bundle pairing and the identity acts as the identity. These properties, and the Hermiticity of S_N^(s)(X^dagger,X), are checked for s=1,2 and N=4,6,8.

Accordingly, bar-eth Z must use S_N^(1), not S_N^(2). For example, a possible definition of the shift contraction is

$$
w=-\frac12\left(\frac{N^{-1}\rho^2+\rho^2N^{-1}}2\right)
\triangleright_1 \partial_rV.
$$

This fixes dimensions and ordering for that expression. It is a specified choice, not a uniquely derived nonlinear gravitational prescription. Self-adjointness also does not imply that the scalar map preserves positivity.

After eliminating the Einstein constraints, the retained master harmonics can be mapped isometrically to the j>=2 subspace of the spin-2 bundle. Let L^2 act with eigenvalue j(j+1). Replacing the scalar angular eigenvalue in the two **derived** potentials by this finite self-adjoint operator gives

$$
I_N^{(2)}=\frac12\sum_{\sigma=\pm}\int dt\,dr_*
\left[\|\partial_t\Psi_\sigma\|_N^2-
\|\partial_{r_*}\Psi_\sigma\|_N^2-
\langle\Psi_\sigma,V_\sigma(L^2,r)\Psi_\sigma\rangle_N\right].
$$

The usual reality condition is imposed on each set of harmonic coefficients; arbitrary complex coefficients in the numerical checks simply test the complexified bilinear identity. This is the exact angular truncation of the two Einstein quadratic sectors. For every retained j, its equation is the corresponding master equation.

The radial inverse in V_Z is nonsingular for j>=2 and r>2M. Normalized-trace checks on independent random fields and derivatives agree with the harmonic coefficient action within 5.8e-12.

This provides a precise quadratic replacement. Extending it to a nonlinear action requires the shift/shear couplings, curvature of a general angular metric, constraints, and boundary terms to be specified together. Adding real traces to the old action is not enough.

**Reproduction:** matrix_realization.py; matrix_results.json.

## 4. An independent uniform lower bound for the quadratic clock

This argument applies to the oscillator matrices actually used in the manuscript's code. It requires no identification with the timelike scalar stress tensor of Fewster-Eveson.

Let f_g(u)=g(u)^2 for the stated smooth compact bump, and define

$$
c_n=\int_0^1 f_g(u)\cos(2\pi n u)\,du,\qquad
\omega_k=\sqrt{k^2+\mu^2}.
$$

For any finite radial truncation, write

$$
H_\mu^{\rm NO}=\frac12(q^TM_q q+p^TM_p p)
-\frac14\operatorname{Tr}(M_q+M_p),\qquad
K=(M_q-M_p)/2.
$$

The matrices M_q and M_p are positive Gram matrices because f_g>=0. The oscillator ground-energy formula and the Powers-Stormer trace inequality give

$$
\inf H_\mu^{\rm NO}
=\frac12\|\sqrt{M_q}\sqrt{M_p}\|_1
-\frac14\operatorname{Tr}(M_q+M_p)
\ \geq\ -\frac12\|K\|_1.
$$

For completeness, the trace step follows from

$$
\|\sqrt P-\sqrt Q\|_2^2\leq\|P-Q\|_1,\qquad
\|\sqrt P\sqrt Q\|_1\geq\operatorname{Tr}(\sqrt P\sqrt Q).
$$

The actual pair kernel can be written without a cancellation of large terms. Put D=omega_k omega_l. Then

$$
K_{kl}=
-\frac{\mu^2(k-l)^2}{\sqrt D(D+kl+\mu^2)}c_{|k-l|}
-\frac{D+kl-\mu^2}{\sqrt D}c_{k+l}.
$$

For l=k+d with d>=1,

$$
\frac{\mu^2d^2}{\sqrt D(D+kl+\mu^2)}
\leq\frac{\mu^2d^2}{2(k^2+\mu^2)^{3/2}},
\qquad
\sum_{k=1}^{\infty}\frac{\mu^2}{(k^2+\mu^2)^{3/2}}\leq1.
$$

The last inequality follows by comparison with the integral over x>=0, whose value is one; at mu=0 the difference term vanishes.

For the sum term, positivity and elementary algebra give

$$
0\leq \frac{D+kl-\mu^2}{\sqrt D}\leq k+l.
$$

One way to check the upper bound is to square it and use

$$
(D+kl-\mu^2)^2-(k+l)^2D
=D\,[2D-k^2-l^2-2\mu^2]-\mu^2(k+l)^2\leq0.
$$

Bounding the trace norm by the sum of absolute matrix entries therefore yields

$$
\|K\|_1\leq
\sum_{n=1}^{\infty}(2n^2-n)|c_n|=:2K_g,
\qquad
\boxed{H_\mu^{\rm NO}\geq-K_g\,\mathbf1}.
$$

This bound is uniform in both mu>=0 and the radial truncation. Since f_g is smooth with vanishing endpoint derivatives,

$$
|c_n|\leq \frac{C_4}{n^4},\qquad
C_4=\frac{\|f_g^{(4)}\|_{L^1(0,1)}}{(2\pi)^4},
$$

so K_g is finite. This proves a uniform lower bound, not merely finite-cutoff convergence. It is conservative and is not asserted to be optimal.

For the stated bump, numerical quadrature gives a weighted Fourier sum through n=80 of 1.1353307301. Adding the analytic tail expression

$$
C_4[2\zeta(2,81)-\zeta(3,81)]
$$

gives K_g approximately 0.61564. The unaugmented partial estimate is approximately 0.56767. The integrals and coefficients in these decimal evaluations use floating-point quadrature; the displayed analytic bound is the rigorous statement, and the decimals are not interval certificates.

Finite-cutoff trace-norm bounds are much tighter in the checked examples. At j=2, mu=sqrt(2)/(2pi), and 240 radial modes, the trace-norm lower bound is approximately -0.0554122. Numerical Williamson minima near -0.00878 involve poorly conditioned positive matrices and cancellation; they are diagnostic values, not certificates of the optimal infinite-mode floor.

**Reproduction:** clock_derivation.py; clock_results.json.

## 5. A conditional finite-area and finite-state consequence

Use exactly the quadratic model normalization

$$
\Xi_s=\frac1s\sum_{j=2}^{J(s)}\sum_{m,p}H_{\mu_j}^{\rm NO},
\qquad d_j=2(2j+1).
$$

If the angular-rank rule is J(s)<=floor(sqrt(s/s_min)), then

$$
D(J)=\sum_{j=2}^{J}d_j=2[(J+1)^2-4],\qquad
\frac{D(J)}s\leq\frac{8}{3s_{\min}}.
$$

The maximum of D(J)/(J^2 s_min) over integers J>=2 occurs at J=3, since
8/3-D(J)/J^2=2(J-3)^2/(3J^2). Therefore

$$
\Xi_s\geq-C_g,\qquad C_g=\frac{8K_g}{3s_{\min}}.
$$

Compression to the actual occupation-budget sector preserves this lower bound. Define on that sector

$$
\Xi_{s,B}=P_B\Xi_sP_B\big|_{\operatorname{ran}P_B},
\qquad \tau_{s,B}=u_{\rm cl}(s)e^{\Xi_{s,B}}.
$$

If the model's classical clock is chosen as u_cl(s)=c_0 sqrt(s), with c_0>0, its clock restriction tau<=ell implies

$$
s\leq(\ell/c_0)^2e^{2C_g}.
$$

**Consequently:** for a discrete area spectrum s_M=M s_min with finite corner degeneracy at each M, the clock allows only finitely many area sectors. Within each sector the imposed budget

$$
\sum_i n_i\left(k_i+\frac{t_{j_i}^2}{k_i}\right)\leq s_M Q_*
$$

allows finitely many modes and occupations. The resulting direct sum is finite-dimensional, and its closed finite-matrix traces are finite without first assuming convergence.

This strengthens the conditional model argument: a finite list of allowed areas need not be separately stipulated once this clock, area spectrum, and angular-rank rule are adopted. Their derivation from full nonlinear gravity is still an assumption, not a consequence of this theorem.

The conservative constant gives a much weaker area-bound multiplier than the manuscript's 1.0036. For s_min=ln(3), the all-area analytic expression evaluates to about 19.86. It establishes finiteness but does not justify the old precision claim.

There is also a uniform vacuum estimate. Since the normal-ordered vacuum mean is zero,

$$
\langle0|\Xi_{s,B}^2|0\rangle
\leq \frac{2D(J)K_g^2}{s^2}
\leq\frac{16K_g^2}{3s_{\min}s}.
$$

Chebyshev's spectral inequality thus gives an O(1/s) upper bound on rejection at any fixed positive threshold eta=ln(ell/u_cl). This conclusion follows from the new kernel bound and the explicit model normalization.

## 6. Calculating the actual occupation-projected clock

The chosen observable is exp(P_B Xi P_B) on ran(P_B), not P_B exp(Xi) P_B. The two generally differ.

### Exact small-budget exponential statistics

Every spin-2 quantum costs at least three, so four-particle states require B>=12. Below that threshold, the vacuum's invariant even sector contains only the vacuum and allowed same-channel two-particle states. The identical magnetic/polarization channels can be combined into permutation-symmetric combinations with sqrt(degeneracy) vacuum coupling. This is an exact symmetry reduction of the vacuum spectral measure.

For s=2B, Q_*=1/2, and the stated area-dependent angular cutoff:

| B | Mean exp(Xi_B) | Variance of exp(Xi_B) |
|---:|---:|---:|
| 4 | 1 exactly | 0 exactly |
| 6 | 1.0000629275 | 1.4563867e-4 |
| 8 | 1.0000575042 | 1.2064790e-4 |
| 9 | 1.0000566514 | 1.2062664e-4 |
| 10 | 1.0000505663 | 1.0583928e-4 |
| 11 | 1.0000657243 | 1.3626488e-4 |

The B=4 vacuum has no allowed pair partner, so the compressed clock leaves it exactly unchanged. These values are not a monotone scaling law: new states enter as B increases.

Independent full occupation-basis calculations with two channels at B=4,6,8 agree with the symmetry reduction to 2.3e-16.

### Large-budget second moments without a huge Hilbert-space matrix

Because Xi|0> contains only pairs, its projected norm can be evaluated directly:

$$
\langle0|\Xi_{s,B}^2|0\rangle
=\frac{1}{2s^2}\sum_j d_j\sum_{k,l}
|K_{kl}^{(j)}|^2\,
\mathbf1_{\omega_{k,j}+\omega_{l,j}\leq B}.
$$

Here the occupation cost omega_(k,j)=k+t_j^2/k differs from the oscillator frequency sqrt(k^2+mu_j^2) used in constructing K. The code keeps them separate.

| s | B | J_max | Variance of Xi_B | Vacuum acceptance lower bound at ell/u_cl=2 |
|---:|---:|---:|---:|---:|
| 20 | 10 | 4 | 9.8917818e-5 | 0.9997941155 |
| 80 | 40 | 8 | 4.1712313e-5 | 0.9999131813 |
| 240 | 120 | 14 | 1.0722458e-5 | 0.9999776826 |
| 5000 | 2500 | 67 | 5.2467912e-7 | 0.9999989079 |

The acceptance column is a **lower bound for the spectral-projector expectation in the vacuum**, not a measured sampling acceptance and not a guarantee for every state.

Smooth Fourier coefficients were evaluated by FFT quadrature; a doubled FFT grid changes the final s=5000 variance by about 1.1e-22. Low-frequency coefficients were independently compared with Gauss quadrature, and the small-budget variances agree with direct spectral calculations. This tests numerical quadrature, not gravitational modeling assumptions.

The favorable s=5000 second-moment result survives the proper many-particle restriction. The old large-area exponential means, variances, and Chernoff estimates have not been re-established by this calculation and must not be substituted for these second-moment results.

**Reproduction:** clock_large_budget.py; clock_large_budget_results.json.

## 7. Changes implemented in the manuscript

1. Replaced the incorrect Fewster-Eveson attribution with an explicit oscillator lower-bound proposition for the chosen quadratic clock.
2. Defined the compressed exponential and replaced unrestricted-occupation determinant estimates with exact small-budget and large-budget second-moment calculations.
3. Added the conditional finite-area and finite-state proposition, with its spectral, clock, angular-rank, and budget assumptions stated explicitly.
4. Withdrew the failed Appendix A expression and defined the new nonlinear bulk prescription. Added exact spherical and quadratic Einstein checks, a finite-matrix realization at quadratic order, and an anisotropic finite-amplitude check.
5. Distinguished this bulk reconstruction from the unresolved finite-rank constraint algebra, boundary quantization, and equivalence to the sampled model.
6. The existing 48% simulation result remains suppression of traceless matrix-field fluctuations. The new action/clock calculations do not change its meaning or establish physical area smoothing.

These changes are implemented in `paper.tex` and its shared `sections/` inputs. The old clock panels are replaced with the actual budget-projected statistics. No new Modal runs were launched and no simulation numbers were changed.

## Reproduction and provenance

Run all seven checks, in dependency order, from the repository root:

```bash
python verification/repair_audit_20261127/run_all.py
```

Dependencies: Python, NumPy, SciPy, SymPy. Each of the seven scripts in that directory can also run independently; the runner supplies their dependency order. They import the repository's existing harmonic-basis and oscillator-matrix builders. The proof-bound check includes 30 noncommuting positive-matrix tests, 21 mass/cutoff kernel checks, and the exact angular-rank identity. The newly derived arguments are explained above; successful assertions alone are not being offered as proofs.

Source hashes in source_hashes.json identify the exact manuscript and legacy scripts read during this work. The original audit used remote main `e541955495d0ba6e21a4a6505282436cac3e7e74`; the historical source hashes are preserved. This implementation starts from `49de5fa`, including the merged citation PR #5. The separate model-use/accounting draft PR #4 is not merged by this change. Current code hashes and run status are recorded in `reproducibility.json`.

Primary references:

- C. J. Fewster and S. P. Eveson, *Bounds on negative energy densities in flat spacetime*, [arXiv:gr-qc/9805024](https://arxiv.org/abs/gr-qc/9805024). Its two-dimensional Q1 factor exceeds one; the replacement proof above does not attribute our clock bound to this theorem.
- R. T. Powers and E. Stormer, *Free states of the canonical anticommutation relations*, Commun. Math. Phys. 16, 1-33 (1970), [DOI:10.1007/BF01645492](https://doi.org/10.1007/BF01645492). Square-root trace inequality.
- A. Aggarwal and M. Singh, *A generalization of Powers-Stormer inequality*, [arXiv:1606.03913](https://arxiv.org/abs/1606.03913). Explicit statement of the trace inequality used here.
- A. R. Solomon, *Off-Shell Duality Invariance of Schwarzschild Perturbation Theory*, Particles 6, 943-974 (2023), [DOI:10.3390/particles6040061](https://doi.org/10.3390/particles6040061). Starting Einstein quadratic actions; the polar elimination was independently recomputed here.
- K. Martel and E. Poisson, *Gravitational perturbations of the Schwarzschild spacetime*, Phys. Rev. D 71, 104003 (2005), [arXiv:gr-qc/0502028](https://arxiv.org/abs/gr-qc/0502028). Gauge-invariant master-field conventions.
- L. Lehner, R. C. Myers, E. Poisson, and R. D. Sorkin, *Gravitational action with null boundaries*, Phys. Rev. D 94, 084046 (2016), [arXiv:1609.00207](https://arxiv.org/abs/1609.00207). Continuum boundary and joint terms, whose pullback and quantization remain separate tasks.
- Chandrasekaran manuscript supplied with this project, Appendix B, Eqs. B.9-B.11. Used for the affine-gauge scalar-action cross-check.
