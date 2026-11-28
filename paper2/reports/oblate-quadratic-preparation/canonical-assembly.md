# Canonical quadratic operator on the original oblate background

## Status and fixed inputs

Preparation authorized by user message `msg_7c4562eb`. This is a derivation and implementation specification, not a newly executed gravitational calculation. The implementation is unexecuted. No spectrum, final constraint count, or nonlinear consistency result is claimed.

The physics input is the **original metric-coefficient Galerkin action** at Paper I commit `157ec0f04866ea36509b325ec65f0a68890f9b17`. It is not the separately considered Hamiltonian-first or reduction-before-cutoff candidate. The source is `paper1/sections/nonlinear_action.tex` at that commit. The canonical starting point is `paper1/reports/canonical-analysis-2026-11-27/README.md`, section 1; only its general Legendre formulas are used here, not its subsequently excluded rank-one witness. The oblate inputs are the committed reports `paper1/reports/oblate-constraint-rank-2026-11-27/README.md` and `paper1/reports/oblate-lapse-consistency-2026-11-27/README.md`.

Normalize the action by `(16 pi G)/(4 pi)`, use angular average `dOmega/(4 pi)`, and retain the radial canonical pairing with measure `dr`. Fix `Lambda=0`. Let

$$
x=\cos\theta,\quad d=1-x^2,\quad \rho=r^2+\epsilon,\quad
\Sigma=r^2+\epsilon x^2,
$$
$$
\bar\sigma=\operatorname{diag}(\Sigma/\rho,\Sigma/d,\rho d),\qquad
\bar a=-1,\quad\bar b=\bar V=\bar W=0,\quad\bar p=0.
$$

Use a compact shell bounded away from zero with positive `rho,Sigma`. Initially all radial integrations by parts are on compactly supported smooth fields. A nonzero boundary trace is not covered by that convention. Boundary realizations are treated separately in [the domain preparation](boundary-domain.md).

The actual coefficient fields are `q=(lambda,X,V,C)` and `z=(a,b,W)`, with

$$
\alpha^2=-a+b^2/\lambda,\qquad
\beta^r=b/\lambda,\qquad \beta^A=W^A-bV^A/\lambda.
$$

The coefficient cutoff applies to these fields, not independently to the derived ADM fields.

## 1. Spatial tangent and its kinetic Gram

Write a spatial perturbation as `h=T u`, where `u=delta q`. At the static background,

$$
h_{rr}=\delta\lambda,\qquad
h_{rA}=\bar\gamma_{AB}\delta V^B,\qquad
h_{AB}=\delta X\,q_{AB}+\delta C_{AB}.
$$

In particular, the vector columns of `T` depend on `r` and `epsilon`; their radial derivatives must be included when differentiating `h`. The nonlinear terms `gamma V^2` and `gamma V` in the full metric map do not add a bulk Hessian term at this stationary solution: the coefficient-map second derivative multiplies the first metric variation, which vanishes on compact tests. This statement does not remove possible boundary terms.

For tensor perturbations define

$$
G(h,k)=\int dr\left\langle\sqrt{\bar\sigma}
\left(h_{ij}k^{ij}-\operatorname{tr}h\operatorname{tr}k\right)
\right\rangle_\Omega,
\qquad S(r)=T^*G_rT.
$$

Here `G_r` denotes the angular form before the `dr` integration. `S` is indefinite but invertible in the stated spatial Legendre neighborhood; positivity is not assumed. In `(r,x,phi)` coordinates the volume density is `sqrt(det sigma)=Sigma`.

## 2. Full shift operator, including radial derivatives

Let `s=(b,W)` be a retained shift-coefficient perturbation and `U_s=Lie_v sigma`, with `v=(b/lambda,W)` at this background. Direct differentiation gives

$$
(U_s)_{rr}=2b'-b\,\partial_r\log\lambda+W^A\partial_A\lambda,
$$
$$
(U_s)_{rA}=\partial_Ab-b\,\partial_A\log\lambda+
\gamma_{AB}(W^B)',
$$
$$
(U_s)_{AB}=(b/\lambda)\partial_r\gamma_{AB}
+\mathcal L_W\gamma_{AB}.
$$

The terms `2 b'` and `gamma W'` are exactly retained tangent directions. Therefore

$$
U_s=U_0s+T J_1s',\qquad
J_1s'=(2b',0,W',0)
$$

in the schematic `(lambda,X,V,C)` ordering. Define

$$
C_0=T^*G_rU_0,\qquad J_0=S^{-1}C_0,\qquad
\mathcal J=J_0+J_1\partial_r,\qquad
\Pi=T S^{-1}T^*G_r.
$$

The normal response is `n_s=(1-Pi)U_0s`. Its quadratic form is algebraic in the radial shift profiles:

$$
A=-\tfrac12U_0^*G_r(1-\Pi)U_0
=\tfrac12(C_0^*S^{-1}C_0-U_0^*G_rU_0).
$$

No radial derivative was frozen in this cancellation. The derivative survives in `J` and in its adjoint. With the base coefficient pairing `dr`,

$$
\mathcal J^*p=J_0^*p-\partial_r(J_1^*p),\qquad
\langle p,\mathcal Js\rangle-\langle\mathcal J^*p,s\rangle
=[p^*J_1s]_{r_-}^{r_+}.
$$

The bracket identity is formal on compact tests until a common domain is chosen.

## 3. Actual canonical momentum and Legendre transform

The normalized quadratic kinetic Lagrangian is

$$
\mathcal L_{2,\mathrm{kin}}=\tfrac14G(T\dot u-U_s,T\dot u-U_s).
$$

Thus the momentum conjugate to the coefficient `u`, with canonical one-form `integral dr p^T delta u`, is

$$
p=\tfrac12S\dot u-\tfrac12T^*G_rU_s,
\qquad \dot u=2S^{-1}p+\mathcal Js.
$$

Substitution, including the full shift-only term, gives

$$
H_{2,\mathrm{kin}}=\langle p,S^{-1}p\rangle
+\langle p,\mathcal Js\rangle+\tfrac12\langle s,As\rangle.
$$

This fixes both the factor of two in the momentum block and the sign of `A`. Dropping the shift-only term would fail to recover the archived normal-shift Hessian. Paper I's symbol `K=2S^{-1}` is a **kinetic momentum block**. It is not the complete phase-space Hessian denoted `K` by the finite reducer.

## 4. Lapse-spatial and curvature blocks

Since the actual variable is `a=h_tt`, `delta alpha=-delta a/2` here. Define `V=-integral sqrt(sigma) R3` at unit lapse. With `mathcal V=V''` pulled back through `T`, its complete tensor form is derived in [the potential preparation](potential-hessian.md).

The lapse-spatial coupling is

$$
\langle f,\mathcal Lu\rangle
=\tfrac12\int dr\,d\Omega/(4\pi)\sqrt{\bar\sigma}\,
f\,\delta R[Tu],
$$
$$
\delta R[h]=\bar\nabla_i\bar\nabla_jh^{ij}
-\bar\Delta\operatorname{tr}h.
$$

On compact tests this is `G(Hess f,Tu)/2`, as in the lapse report. A weak form with only first derivatives is

$$
\langle f,\mathcal Lu\rangle
=-\tfrac12\int\sqrt{\bar\sigma}\,
\bar\nabla_if\left(\bar\nabla_jh^{ij}-\bar\nabla^i\operatorname{tr}h\right)
+\tfrac12\int_{\partial I\times S^2}\sqrt{\gamma_\partial}\,
f n_i\left(\bar\nabla_jh^{ij}-\bar\nabla^i\operatorname{tr}h\right),
$$

with the common angular-average normalization understood. The boundary flux is not generally zero for arbitrary boundary data.

At this static zero-momentum solution there is no quadratic `u-p` term, no `a-p` term, and no `a-a` or `a-s` temporal block. The `b`-dependence of the derived lapse multiplies the zero background curvature/kinetic density at this order. These simplifications are background-specific.

## 5. The complete formal block interface

With `w=(u,p)` and temporal perturbations `z=(delta a,s)`, the bulk quadratic form is

$$
H_2=\tfrac12\langle u,\mathcal Vu\rangle
+\tfrac12\langle p,2S^{-1}p\rangle
+\langle p,\mathcal Js\rangle
+\langle\delta a,\mathcal Lu\rangle
+\tfrac12\langle s,As\rangle.
$$

Consequently the reducer's blocks, as formal operators before radial discretization, are

$$
\boxed{K_{\rm full}=\begin{pmatrix}\mathcal V&0\\0&2S^{-1}\end{pmatrix},\quad
E=\begin{pmatrix}\mathcal L&0\\0&\mathcal J^*\end{pmatrix},\quad
D=\begin{pmatrix}0&0\\0&A\end{pmatrix}.}
$$

This is a concrete assembly prescription, not yet a finite matrix export. It derives from the retained Lagrangian's Legendre transform. In particular it does not replace that transform by independently truncating the continuum canonical ADM Hamiltonian.

For an algebraic complement of `ker A`, shift elimination changes only the momentum block to

$$
K_{pp,\mathrm{eff}}=2S^{-1}-\mathcal J_R A_R^{-1}\mathcal J_R^*.
$$

This can contain radial derivatives of `p`. Its definition requires the multiplier solution and adjoint domains to be compatible. The reduced lapse bracket is `L K_pp,eff L*`; its archived continuous realization alone does not prove those additional domain statements.

## 6. Harmonics and real canonical variables

Use the associated-Legendre amplitudes of the pinned `all_modes.py`: `Y_lm=P_l^m(x)`, with `exp(i m phi)` implicit, polar vector `(d Y', i m Y/d)`, axial vector `(i m Y,-Y')`, and its stated electric/magnetic tracefree tensors. No change to orthonormal spin-weighted coordinates is implicit.

For `m>0` choose real tensor harmonics `sqrt(2) Re(F)` and `sqrt(2) Im(F)`, where `F=T_lm(x) exp(i m phi)`. The same choice is made for scalar, vector, tensor, and temporal fields. If a complex sesquilinear block is `H=H_R+i H_I`, its real block in `(cos,sin)` order is

$$
\mathfrak R(H)=\begin{pmatrix}H_R&H_I\\-H_I&H_R\end{pmatrix}.
$$

This sign follows from the positive sine convention above. It differs from the commonly used realification for complex coefficient `q_cos+i q_sin`, because the field coefficients for this real basis are proportional to `q_cos-i q_sin`. Real momenta must be paired with these same real coordinates; simply duplicating eigenvalues is not a canonical export. The `m=0` fields are used once. Labels must encode species, degree, magnetic number, and real component.

## 7. Radial projection order is part of the calculation

Two distinct finite approximations must not be conflated:

1. **Canonical radial projection after the angular Legendre transform.** For a real radial basis matrix `F(r)` with `M=integral F^T F dr`, take `u=F Q` and `p=F M^{-1} P`. Then `integral p^T dot u dr=P^T dot Q`. Pull back every term of the already derived Hamiltonian to these coordinates and explicitly state the chosen temporal basis. The finite momentum block is the pullback of `2S(r)^{-1}`.
2. **Radial restriction of the Lagrangian before its Legendre transform.** First form the finite spatial mass matrix `W_M=integral F^T S F dr/2`, the corresponding shift coupling and shift norm, and then invert `W_M`. In general its inverse is not the matrix in item 1. Radial projection can introduce additional normal shift responses even when the angular response is tangent.

Neither choice is silently selected by this preparation. The first numerical protocol must identify its projection order and test the resulting radial approximation. Angular rank `N` and radial resolution are separate parameters.

An equal polynomial cutoff for every field need not contain `J_k`, the deformed gradient descriptors, or their radial derivatives. Therefore an unexpected small singular value in a radial matrix is not automatically the original field theory's rank defect. The domain/basis specification must provide inclusion maps, an enlarged compatible test space, or a controlled convergence interpretation.

## 8. Implementation and bounded validation status

The prepared `kinematic_blocks.py` supplies the complete general-r spatial tangent, derivative-free shift columns, exact derivative injection `J1`, angular density blocks, their post-integration canonical assembly, and the stated realification. It does not execute angular integrals, choose a boundary condition, or infer a numerical rank. Importing it launches no calculation. The implementation was unrun at the original preparation checkpoint; the subsequent [corrected kinematic check](kinematic-recheck/README.md) passed its bounded test surface.

The accompanying [interface manifest](interface.json) lists the required fields and explicitly records missing measured matrices. The previous exact-rational reducer must not be fed rationalized quadrature output; a numerical reduction method requires a separate precision and structural-null policy.

The original [first-check proposal](next-check.md) covered general-r tensor maps against direct Lie differentiation, density assembly against a full-matrix contraction at a deformed rational point and the pinned pairing in one exactly integrated round sector, the canonical Legendre identity and derivative-column cancellation in that sector, the r=1 specialization against the pinned angular source, and realification on an explicitly imaginary cross block. Its failed initial invocation and corrected passing rerun have separate immutable records.

The prepared `weak_forms.py` implements first-derivative potential and lapse densities, including the r-dependent basis connection and separate endpoint fluxes. It returns the four radial derivative slots before angular or radial integration. The later [eight-batch weak-form validation](../oblate-weak-curvature-validation/full-recheck-pass/README.md) passed representative `L=2`, `m=0,1,2` comparisons with independent curvature variations and Green identities at recorded rational points, three fixed symbolic gauge controls, and TT/ablation checks. This validates the stated local implementation checks without establishing an all-label Green theorem, angularly integrated operator, finite radial kernel preservation, common endpoint domain, or physical spectrum.

## 9. Input to the nonlinear preservation track

Paper I's parallel preparation requires an on-shell retained linearized seed, not merely a vector in a primary/secondary nullspace. The shared manifest therefore includes labelled coefficient functions or consistent time/radial jets, canonical-momentum-map evidence, and residuals for all retained linear equations and constraints. It also requires the **full continuum linear coefficient-Euler residual before angular projection**, including discarded temporal components. That residual is not recoverable from the retained `K,E,D` matrices alone; the tensor reconstruction maps must be preserved.

Paper I owns its `linear_metric`, `linearized_einstein`, `euler_pair_density`, and `candidate_source_density` implementations on its branch. At this background the spacetime perturbation is `h_tt=a1`, `h_tr=b1`, `h_tA=gamma0_AB W1^B`, with spatial components given in section 1. Its coordinate density uses `Sigma` with measure `dr dx/2` after azimuthal averaging, with no extra `r^2` factor. Pure continuum-gauge seeds are controls because their full linear Einstein residual vanishes; their finite coordinate-orbit leakage alone is not an on-shell nonlinear obstruction. No admissible seed or nonlinear compatibility value has been produced by this preparation.
