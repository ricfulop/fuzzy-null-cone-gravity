# All-mode oblate constraint ranks; lapse continuation pending

27 November 2026. Starting source: 71a29e402e0941e2d85f9894c323d5440e0e98cb. The user approved the constraint-rank investigation in msg_ce3c51d4. This continues the [exact admissible oblate family](../admissible-branches-2026-11-27/README.md) of the [same Paper I action](../../sections/nonlinear_action.tex).

**Completed result.** On the exact static oblate solution, for sufficiently small nonzero deformation, the full temporal Hessian has rank **18 at $N=3$ and 26 at $N=4$**. Including the secondary–secondary brackets, the Poisson operator for the complete primary-plus-secondary set has local angular rank **36 and 52**, respectively. These are exact small-deformation results supported by all-magnetic-mode rational certificates, an explicit matching kernel, and independent finite-deformation numerical checks.

**Not completed.** These are not the ranks of a fully stabilized nonlinear Dirac system. A subsequent linearized lapse-preservation calculation reached a failed symbolic assertion in an exact candidate-gradient identity. Its partial outputs and failure log are preserved. Inspection, any correction, and a rerun await the user's copilot approval. No physical degree-of-freedom count or nonlinear constant-rank neighborhood is claimed.

## 1. Background and radial terms

Use $x=\cos\theta$, $d=1-x^2$, $\rho=r^2+\epsilon$, $\Sigma=r^2+\epsilon x^2$, and the exact static family

$$
\sigma=\frac{\Sigma}{\rho}dr^2+\frac{\Sigma}{d}dx^2+\rho d\,d\phi^2,\qquad
a=-1,\quad b=V=W=0,\quad p=0,\quad\Lambda=0.
$$

All actual spatial coefficients have $j\le2$, and the four-metric is Minkowski for every time. Work on a compact shell bounded away from $r=0$, with positive $\rho,\Sigma$ and compact radial smearings. Use the established normalized action $(16\pi G)/(4\pi)I$.

Let $T$ map coefficient velocities to spatial tensor variations. The retained spatial Gram matrix $T^*GT$ is invertible near round geometry, where $G$ is the spatially integrated DeWitt form. Write $\Pi$ for its tangent projector and $N=1-\Pi$. A retained shift coefficient vector $s=(b,W)$ induces

$$
v_s^r=b/\lambda,\qquad v_s^A=W^A,\qquad\lambda=\Sigma/\rho,
$$

and define its normal response

$$
n_s=N\mathcal L_{v_s}\sigma.
$$

At the static zero-momentum point the only nonzero temporal-Hessian block is

$$
A(s,t)=H''_{ss}[s,t]=-\frac12G(n_s,n_t).
$$

Lapse–lapse and lapse–shift blocks vanish there. This is a property of this point, not an assumption about arbitrary data.

The dependence of $n_s$ on radial profiles is algebraic. In the $rr$ component, the term containing $\partial_r b$ is $2\partial_r b$, represented exactly by a retained $\delta\lambda$. In the $rA$ components, the terms $\gamma_{AB}\partial_rW^B$ are represented exactly by $\delta V=\partial_rW$. Projecting these terms away is an exact tangent cancellation; the computation does not freeze radial derivatives. The remaining tangent coefficients used below in the secondary brackets still contain those derivatives.

## 2. Exact upper bound from tangent shifts

The angular metric can be written

$$
\gamma=\rho q_{S^2}-\epsilon\,dx\otimes dx.
$$

For an angular shift with $g=W(x)$,

$$
\mathcal L_W\lambda=\frac{2\epsilon}{\rho}xg,\qquad
\mathcal L_W\gamma=\rho\mathcal L_Wq_{S^2}
-\epsilon(dx\otimes dg+dg\otimes dx).
$$

Consequently every retained $W$ satisfying $W(x)\in V_{L-1}$ is exactly tangent at finite deformation. Here $V_j$ denotes scalar harmonics of degree at most $j$ and $L=N-1$.

The dimension of this angular subspace is $2L^2-2$. In a complex fixed-$m$ basis the polar vector's $x$ component is $d\,\partial_xY_{\ell m}$, with scalar degrees $\ell-1,\ell+1$, and the axial vector's $x$ component is $imY_{\ell m}$. Eliminating degrees $L+1$ and $L$ in $W(x)$ imposes two independent complex conditions in each nonnegative-$m$ block. The $m=0$ block is real; every $m>0$ block contributes twice its complex dimension. The [exact kernel certificate](exact-kernel.json) records these maps and nullvectors.

Independently take $v_s^r=f\in V_{L-2}$ and $W=0$, so the retained coefficient is $b=\lambda f$. Its $rr$ and angular variations stay in band. The mixed variation is represented by

$$
\delta V=\lambda\gamma^{-1}df
=\rho^{-1}\operatorname{grad}_{S^2}f
-\epsilon\rho^{-2}(\partial_\phi f)\partial_\phi.
$$

This is retained: the second term is a scalar harmonic of degree at most $L-2$ times a degree-one rotational Killing field. Arbitrary radial profiles and their derivatives are allowed. These radial directions have dimension $(L-1)^2$.

The two subspaces are independent and together give

$$
\dim\ker(s\mapsto n_s)\ \ge\
(2L^2-2)+(L-1)^2=3L^2-2L-1.
$$

There are $3(L+1)^2-2$ shift coefficients, hence both the normal-map rank and shift-Hessian rank are at most $8L+2$.

## 3. Exact lower bound and complete kernel

At first order in $\epsilon$, subtract the changing retained tangent before projecting with the round discarded projector. The response can only occupy degrees $L+1,L+2$, because the background deformation has degrees zero and two. The verifier computes the complete discarded coefficients in every magnetic block and forms their exact DeWitt Gram matrices.

| Cutoff | Complex ranks in $m=0,\ldots,L$ | Real total rank | Exact shift-kernel dimension |
|---|---|---:|---:|
| $L=2$ ($N=3$) | $4,4,3$ | 18 | 7 |
| $L=3$ ($N=4$) | $4,4,4,3$ | 26 | 20 |

Each lower bound is certified by a nonzero rational principal minor in [the exact leading matrices](exact-leading.json). All successive Sylvester minors of the selected blocks are positive. The leading normal Gram is therefore positive on its image; its rank equals that of the leading normal map. This explicitly handles the indefinite DeWitt form.

Analyticity of the finite angular integrals and inverse spatial Gram matrix promotes the nonzero leading minors to nonzero minors for every sufficiently small nonzero $\epsilon$. The exact tangent subspace above gives the matching upper bound. Thus

$$
\boxed{\operatorname{rank}A=\operatorname{rank}(s\mapsto n_s)=8L+2,}
$$

and both kernels are precisely the direct sum described in section 2. On a fixed compact shell the leading coefficients and analytic remainders are controlled uniformly after the appropriate powers of $r$ are restored.

All lapse directions are in the temporal-Hessian kernel at the static point. Its total nullity is consequently $4L^2$, namely 16 or 36. This is not a gauge-function count.

## 4. Including the secondary–secondary brackets

Let $C_a[f]$ denote the secondary obtained by varying the actual coefficient $a=h_{tt}$, and $C_s[s]$ a shift secondary. At the static point, the first variation of the lapse secondary is

$$
\delta C_a[f]=\frac12G(\nabla^2f,T\delta q),
$$

because the spatial metric is flat and compact integration by parts turns the linearized scalar curvature into the Hessian test. The factor $1/2$ follows from $\partial\alpha/\partial a=-1/2$. The momentum derivative of a shift secondary is the retained coefficient vector $J_s$ defined by

$$
TJ_s=\Pi\mathcal L_{v_s}\sigma.
$$

The lapse–lapse and shift–shift secondary brackets vanish at $p=0$, zero shift. The mixed bracket is

$$
B(f,s)=\{C_a[f],C_s[s]\}
=\frac12G(\nabla^2f,\Pi\mathcal L_{v_s}\sigma)
=-\frac12G(\nabla^2f,n_s).
$$

The last identity is the flat spatial divergence identity
$G(\nabla^2f,\mathcal L_v\sigma)=0$ after compact radial and angular integration by parts. It retains radial derivatives of the descriptors.

Since $\ker A=\ker n$, the mixed bracket annihilates every shift-kernel section, including arbitrary radial profiles. Choose an invertible complement of $A$. Its algebraic inverse removes the remaining mixed block by invertible block operations on the smeared constraint operator. Equivalently, the mixed operator factors through $A$; any radial derivatives remain in the factor multiplying it. There is no residual secondary bracket on the temporal-Hessian kernel.

It follows that the Poisson operator of the complete primary-plus-secondary set has local angular rank

$$
\boxed{\operatorname{rank}\mathbb P_{\rm primary+secondary}
=2(8L+2)=16L+4=36,\ 52.}
$$

This is a rank statement for the radial field theory under the specified compact-smearing realization, not a finite total rank of an infinite-dimensional operator. It does not include constraints that may arise at subsequent stabilization stages.

## 5. What has been derived for the next stage

The stationary solution has zero canonical vector field, so evaluating first preservation only at that point gives $A\mu_s=0$ and leaves all lapse rates free. Those zero values do not establish closure nearby.

A linearized continuation analysis has been set up. Choose a complement $R$ to the shift kernel. Let

$$
K=2(T^*GT)^{-1},\qquad
K_{\rm eff}=K-J_RA_R^{-1}J_R^*,\qquad
L_fh=\frac12G(\nabla^2f,Th).
$$

Eliminating the nondegenerate shift pairs gives candidate next lapse conditions

$$
\mathcal T_f(p)=L_fK_{\rm eff}p=0.
$$

For a kernel shift $k$, exact tangency and flatness give $L_fJ_k=0$ and $V''J_k=0$, where $V''$ is the Hessian of the unit-lapse spatial-curvature potential. Thus the null shift secondaries have identically vanishing preservation in this linearized system. This does not make them nonlinear first-class generators.

Let $P_{\mathcal N}$ project, in the DeWitt form, onto the nondegenerate image of $s\mapsto n_s$, and put

$$
Z_f=(1-P_{\mathcal N})N\nabla^2f.
$$

The candidate lapse–tertiary bracket reduces to

$$
\{L_fh,\mathcal T_g\}_{\rm lin}
=-\frac12G(Z_f,Z_g).
$$

This uses the integrated flat identity $G(\nabla^2f,\nabla^2g)=0$ for compact tests. It is an analytic reduction, not a numerical replacement of radial derivatives.

The partial [lapse outputs](lapse-leading.json) contain exact leading coefficients for separate radial-profile and first-derivative slots. They suggest eight active scalar channels at $N=3$ and twelve at $N=4$. The proposed complementary exact gradient families still require a fully passing verification record; the extended checker stopped before completing that certification. Accordingly no completed lapse-kernel theorem, stabilized nonlinear rank, or degree-of-freedom count is asserted by this checkpoint.

## 6. Passing checks and the failed follow-up

The [all-mode verifier](../../verification/oblate_constraint_rank_20261127/all_modes.py) passed. Run from the paper1 folder:

    uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 python verification/oblate_constraint_rank_20261127/all_modes.py --verify

Its [log](verification.log), [results](checks.json), and [source/environment manifest](reproducibility.json) include exact leading minors, exact kernel dimensions, and complete spatial and shift matrices in every magnetic sector. At each of $\epsilon/r^2=0.01,0.05,0.1$, 64- and 128-node quadratures agree. Numerical shift ranks are unchanged across absolute thresholds from $10^{-8}$ to $10^{-14}$. The complete spatial signatures are $(35+,9-)$ at $N=3$ and $(70+,16-)$ at $N=4$. These numerical points support the calculation; the exact rank theorem uses the small-deformation proof.

The later [lapse checker](../../verification/oblate_constraint_rank_20261127/lapse_prolongation.py) exited with an assertion failure while comparing a proposed $m=1$ gradient vector with its retained polar/axial representation. Its [failure log](lapse-failed-check.log) is preserved separately from the passing rank log. All preceding leading-Gram calculations printed successfully, but partial output is not a substitute for a completed verification.

An independent analytic derivation supports the proposed gradient identity. Algebraic simplification or a basis-convention mismatch is the likely cause of the assertion failure; the residual has not yet been inspected, so this diagnosis is provisional. No correction or further rerun is authorized by this checkpoint. The user's copilot approval is required before proceeding after the failed run.

The failure concerns the newly added lapse-gradient identity and does not alter the separately passing all-mode rank certificate. Two independent reviewers support the shift-kernel argument and its rank. The mixed-secondary reduction and the analytic formulation of further preservation were checked independently as well. No external/GPU compute or simulation experiment was launched.

The next action proposed to the user is limited to inspecting the failed symbolic residual, correcting the identity or its implementation if warranted, and rerunning the local lapse checks. The full nonlinear release goal remains open.
