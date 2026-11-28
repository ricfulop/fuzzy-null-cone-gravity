# Second preservation excludes the rank-one vacuum witness

27 November 2026. Starting source: `dfdc3717ed98df97cc97e51ed509b87ad058d5f8`. The user approved this higher-order consistency calculation in message `msg_b2aec6b1`. It continues the [first preservation calculation](../dirac-stabilization-2026-11-27/README.md) for the [actual metric-coefficient action](../../sections/nonlinear_action.tex). No regulator change is made.

**Result.** For the specified flat, rank-one vacuum initial datum at $N=3,4$, the retained equations admit no sufficiently smooth continued solution. Although every secondary and its first time derivative can vanish initially, a retained angular-shift secondary necessarily obeys

$$
\boxed{\ddot{\mathcal C}_W[\xi]\big|_0
=-\frac{512}{1575}\kappa^3\int dr\,r^2w(r)\ne0,}
\qquad
\xi=w(r)D^AP_2\,\partial_A,
$$

for a nonnegative, nonzero smooth radial bump $w$ supported on a shell away from the origin and $\kappa\ne0$. All rates that solve first preservation give the same value. Second temporal rates cannot cancel it.

This is an obstruction to a **projected Galerkin trajectory through these fixed data**, stronger than the earlier obstruction to an exact continuum Einstein trajectory inside the ansatz. It does not exclude other initial data or classify the final constraint surface. In particular, the earlier lapse-block ranks 8/12 and full Poisson-matrix lower bounds 16/24 at this witness remain valid primary–secondary-surface results; this witness cannot promote them to ranks on a final dynamically admissible surface. The [full release goal](../../GOAL.md) remains open.

## 1. Setup and required smoothness

Use the same normalized action $\bar I=(16\pi G)/(4\pi)I$ and angular average $d\Omega/(4\pi)$ as the canonical checkpoints. The spatial and temporal coefficients are

$$
q=(\lambda,X,V,C),\qquad z=(a,b,W)=(h_{tt},h_{tr},A_t),
$$

with $\gamma=Xq_{S^2}+C$ and

$$
\sigma_{rr}=\lambda+\gamma_{AB}V^AV^B,\quad
\sigma_{rA}=\gamma_{AB}V^B,\quad \sigma_{AB}=\gamma_{AB},
$$

$$
\alpha=\sqrt{-a+b^2/\lambda},\qquad
\beta^r=b/\lambda,\qquad \beta^A=W^A-(b/\lambda)V^A.
$$

All spin-0, spin-1, and spin-2 coefficient multiplets through $L=N-1$ are retained. Products and inverses are pointwise after reconstruction. $P_L$ is the round spatial tangent projector, and $Q_L=1-P_L$. The proof is local in a radial shell and uses compactly supported variations, as in the bulk ADM action underlying the previous checkpoints.

The initial data are

$$
\sigma=dr^2+r^2d\Omega^2,\quad
a=-1,\quad b=W=V=C=0,\quad\lambda=1,\quad X=r^2,
$$

$$
K_{ij}=\kappa n_i n_j,\qquad n=(0,0,1),\qquad\kappa\ne0,
$$

in Cartesian coordinates. The momentum is fixed by the regular retained Legendre map from $\dot\sigma=2K$. All full continuum constraints, and hence every retained primary and secondary, vanish initially. The hypothesis fixes the temporal coefficients as well as the spatial coefficients and momenta.

Assume sufficient regularity to differentiate the retained secondary twice and the spatial Euler–Lagrange equation once. Smooth coefficient fields suffice; three time derivatives of the spatial coefficients, two of the temporal coefficients, and the required mixed radial derivatives are sufficient locally. The result does not analyze nonsmooth or distributional evolutions.

## 2. What first preservation already fixes

Let $T_I=\partial\sigma/\partial q^I$ and $B=T_{IJ}\dot q^I\dot q^J$. At the witness, the retained spatial equation leaves the following difference from the full Einstein spatial acceleration with the same initial temporal rates:

$$
R=\ddot\sigma_{\rm Gal}-\ddot\sigma_{\rm GR}
=Q_LB+Q_L(\mu_aK),\qquad\mu_a=\dot a.
$$

First preservation is solved by

$$
\mu_a=\frac{32\kappa}{9}P_2+u_{\rm low},\qquad
u_{\rm low}\in V_{L-2},
$$

with arbitrary retained $\mu_b,\mu_W$ and arbitrary radial profiles. Because $Q_L(u_{\rm low}K)=0$ and every shift-rate Lie derivative is tangent at the round point, **the entire tensor $R$ is the same throughout this affine family**. This is narrower than saying that $R$ is independent of every possible lapse rate before first preservation.

Set $x=\cos\theta$, $d=1-x^2$, and describe a spherical orthonormal tensor by

$$
(f,g,\widetilde v,\widetilde c)
=\left(U_{\hat r\hat r},\frac{U_{\hat\theta\hat\theta}+U_{\hat\phi\hat\phi}}2,
\frac{U_{\hat r\hat\theta}}{\sqrt d},
\frac{U_{\hat\theta\hat\theta}-U_{\hat\phi\hat\phi}}{2d}\right).
$$

For both $L=2,3$, the first-preserved residual is

$$
\frac{R}{\kappa^2}=
\left(
-\frac{8}{105}(35x^4-30x^2+3),
-\frac{8}{105}(35x^4-30x^2+3),
\frac{8}{21}x(7x^2-3),
\frac{8}{21}(7x^2-1)
\right).
$$

It is entirely discarded. Let $E_\sigma$ denote the unrestricted ADM spatial Euler derivative, with convention $\delta\bar I=\int E_\sigma[\delta\sigma]+\cdots$. At these data,

$$
E_\sigma[U]=-\frac12\int dr\,r^2\langle R,U\rangle_{\rm DW}.
$$

This follows from the kinetic Lagrangian $\langle\dot\sigma-\mathcal L_\beta\sigma,\alpha^{-1}(\dot\sigma-\mathcal L_\beta\sigma)\rangle_{\rm DW}/4$: changing the spatial acceleration by $R$ changes its Euler derivative by $-G_{\rm DW}(R,\cdot)/2$. Other initial values and first derivatives agree with the continuum comparison.

## 3. Differentiated angular-diffeomorphism identity

Use a fixed spatial angular descriptor

$$
\xi^x=3w(r)x(1-x^2),\qquad\xi^r=\xi^\phi=0.
$$

It is a retained polar $j=2$ vector for both cutoffs. Since $W=A_t$ enters the derived shift linearly with unit coefficient, its retained secondary is exactly the full momentum constraint tested against this fixed vector:

$$
\mathcal C_W[\xi]=\int dr\,d^2y\,\xi^A\mathcal H_A
$$

where $\mathcal H_A$ is the canonical momentum-constraint coordinate density for $\bar I$. In $(r,x=\cos\theta,\phi)$ coordinates the unit-sphere area density is one; the $r^2$ factor enters through the spatial volume density. There is no time-dependent testing weight from the coefficient-to-ADM map in this angular sector.

For clarity, derive the required identity from the unrestricted action, rather than assuming that this angular diffeomorphism is a symmetry of the truncated theory. Give the descriptor a time factor $\eta(t)$. Then

$$
\delta\sigma=\eta\mathcal L_\xi\sigma,\quad
\delta\alpha=\eta\mathcal L_\xi\alpha,\quad
\delta\beta=\dot\eta\,\xi+\eta[\xi,\beta].
$$

Write the unrestricted temporal Euler derivatives as $E_\alpha$ and $E_\beta$, so $E_\beta[\xi]=-\mathcal C_W[\xi]$. Integrating the $\dot\eta$ term by parts gives the off-shell identity

$$
\dot{\mathcal C}_W[\xi]
=-E_\sigma[\mathcal L_\xi\sigma]
-E_\alpha[\mathcal L_\xi\alpha]
-E_\beta[[\xi,\beta]].
$$

At the witness, $E_\alpha=E_\beta=0$ against **all** angular tests because the full continuum constraints vanish pointwise; also $\alpha=1$ and $\beta=0$. Consequently each of the last two terms and its first time derivative vanishes. This does not require the discarded continuum constraint derivatives to vanish. Their coefficients $\mathcal L_\xi\alpha$ and $[\xi,\beta]$ vanish initially.

At the round point $\mathcal L_\xi\sigma=T_Iv_\xi^I$ for a retained coefficient vector $v_\xi$. A putative continued retained solution must satisfy its spatial Euler–Lagrange equations and their derivatives. The actual coefficient chain rule is

$$
0=E_\sigma[T_I]+E_\alpha[\partial_I\alpha]
+E_\beta[\partial_I\beta].
$$

At fixed $a=-1,b=W=0$, the functions $\alpha=1$ and $\beta=0$ are independent of **every nearby spatial coefficient**. Thus $\partial_I\alpha=\partial_I\beta=0$ at the witness. Together with $E_\alpha=E_\beta=0$, this kills every temporal term in the differentiated chain rule, even when $\dot b$ or $\dot W$ is arbitrary. Holding the initial coefficient vector $v_\xi$ fixed therefore gives

$$
\dot E_\sigma[Tv_\xi]+E_\sigma[\dot T v_\xi]=0.
$$

Differentiating the off-shell identity now eliminates the unknown $\dot E_\sigma$, including dependence on the next time derivatives:

$$
\boxed{\ddot{\mathcal C}_W[\xi]
=-E_\sigma[U_\xi]
=\frac12\int dr\,r^2\langle R,U_\xi\rangle_{\rm DW},}
$$

$$
U_\xi=\mathcal L_\xi\dot\sigma-\dot T v_\xi.
$$

This argument differentiates equations that any smooth retained solution must obey. It does not presume a smooth extension of the pointwise temporal-Hessian kernel to a neighboring rank stratum.

## 4. Exact changing-tangent calculation

In $(r,x,\phi)$ coordinates the initial coefficient velocities are

$$
\dot\lambda=2\kappa x^2,\quad
\dot\gamma_{xx}=2\kappa r^2,\quad
\dot\gamma_{\phi\phi}=0,\quad
\dot V^x=\frac{2\kappa x(1-x^2)}r.
$$

The unrestricted angular transformation of these coefficient types is

$$
\delta_\xi\lambda=\xi^A\partial_A\lambda,\qquad
\delta_\xi\gamma=\mathcal L_\xi\gamma,\qquad
\delta_\xi V=\partial_r\xi+[\xi,V].
$$

At the witness $V=0$, the last expression contains $\partial_r\xi$, including $w'(r)$. These terms enter both $\mathcal L_\xi\dot\sigma$ and $\dot T v_\xi$ and cancel in their difference. The verifier checks the full tensor reconstruction before differentiating, with arbitrary symbolic $w(r)$; it does not set radial derivatives to zero.

The resulting coordinate tensor is

$$
(U_\xi)_{rr}=12\kappa w x^2d,\qquad
(U_\xi)_{rx}=0,\qquad
(U_\xi)_{xx}=12\kappa w r^2(1-3x^2),\qquad
(U_\xi)_{\phi\phi}=0.
$$

In the component convention of section 2,

$$
\frac{U_\xi}{\kappa w}=
\left(12x^2d,\ 6d(1-3x^2),\ 0,\ 6(1-3x^2)\right).
$$

An independent coefficient route differentiates $\delta_\xi q$ along the known initial $\dot q$ and applies $T_0$. It gives the same $U_\xi$; here $[\xi,\dot V]=0$ because the angular vector fields are proportional. A linear coefficient path is sufficient for these derivatives of the embedding and is not asserted to solve the evolution equations.

The round DeWitt pairing of components $(f,g,\widetilde v,\widetilde c)$ and $(F,G,\widetilde V,\widetilde C)$ is the angular average of

$$
-2fG-2gF-2gG+2d\widetilde v\widetilde V+2d^2\widetilde c\widetilde C.
$$

Direct polynomial integration with measure $dx/2$ yields

$$
\langle R,U_\xi\rangle_{\rm DW}
=-\frac{1024}{1575}\kappa^3w(r).
$$

Pairing with $Q_LU_\xi$ gives the same number because $R$ is discarded and the round projector is DeWitt self-adjoint. Substitution into the boxed Noether consequence proves the stated nonzero second derivative for both $N=3,4$.

## 5. Consequence and limits

A retained solution would have $\mathcal C_W[\xi](t)=0$ for every retained test and hence $\ddot{\mathcal C}_W[\xi](0)=0$. Choosing a nonnegative, nonzero compact radial bump contradicts this necessary condition. The obstruction does not depend on a particular choice among the 26/50 free first-step rates: all those rates leave $R$ unchanged, while $U_\xi$ depends only on the fixed initial spatial velocity. The differentiated identity also eliminates any prospective cancellation by second temporal rates.

The witness is therefore excluded by higher preservation. This is not a globally constructed tertiary constraint or a classification of strata: it is a necessary higher-order compatibility condition evaluated at specified data. Constraint surfaces may change rank near this symmetric point. The prior primary–secondary rank result is not false, but its witness cannot supply a rank result on a smooth dynamically admissible branch at these cutoffs.

The result does not exclude all projected solutions, exclude other temporal initial values or other spatial data, settle larger $N$, count physical degrees of freedom, or resolve boundary quantization and simulation matching. Nor does it demonstrate a contradiction of full general relativity: the original vacuum datum has continuum developments outside this finite coefficient ansatz. The manuscript and archived simulation claims are not promoted by this calculation.

## 6. Verification and next decision

Run the [local verifier](../../verification/higher_order_consistency_20261127/second_preservation.py) from `paper1/`:

```bash
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 \
  python verification/higher_order_consistency_20261127/second_preservation.py
```

The [exact results](checks.json), [execution log](verification.log), and [environment/source hashes](reproducibility.json) archive one successful verification run. It checks the reconstructed tensor transformation with arbitrary radial profiles, two derivations of its changing tangent, the initial ADM chain-rule zeros, the previously solved first-preservation equations, the all-$m$ homogeneous lapse kernel, and the exact rational obstruction. Independent angular quadrature at 24, 48, and 96 nodes agrees within $8.89\times10^{-16}$. The Noether argument is an analytic proof in this report, not a claim that numerical quadrature establishes an off-shell identity.

Two independent delegated reviews confirmed the coefficient, sign, radial cancellation, and rate independence. One additionally checked the coefficient-to-ADM chain rule explicitly. No simulation experiment, external compute, or GPU job was launched. The scientific verifier completed successfully on its first run.

The next proposed research step is to characterize a dynamically admissible branch of the actual retained system and evaluate the canonical brackets there. A new witness must survive stabilization; the excluded datum cannot be reused for a final degree-of-freedom claim. Whether to pursue that branch analysis or revise the formulation is a separate copilot decision. Neither direction is started by this checkpoint.
