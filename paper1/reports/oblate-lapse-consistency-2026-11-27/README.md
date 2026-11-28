# Oblate lapse preservation: corrected verification and linearized result

27 November 2026. Parent source: 6e31a5e3fa1e0f0b9f7b54fb261df34de9b7742d. User approval msg_040d0d1b covers inspection of the failed identity, warranted correction, and the local rerun. This continues the [all-mode oblate constraint analysis](../oblate-constraint-rank-2026-11-27/README.md) of the [actual Paper I action](../../sections/nonlinear_action.tex).

**Completed checkpoint.** The failure was a symbolic comparison defect, not a failure of the proposed lapse descriptor. The corrected checker passed every magnetic block, the exact retained-gradient identities, and the leading radial forms. At $N=3$ there is one exact redundant lapse-profile family and eight complementary real angular channels; at $N=4$ there are four such families and twelve complementary channels. The complementary lapse–tertiary bracket is nondegenerate for sufficiently small nonzero oblateness in the continuous radial form extension specified below.

**Scope.** These are linearized preservation results about an exact stationary solution. They do not classify the surrounding nonlinear constraint surface, establish a physical degree-of-freedom count, or settle the boundary quantum theory and simulation connection. The full nonlinear release goal remains open.

## 1. Diagnosis, correction, and provenance

The original assertion compared each component of $W_1-c_1P_1-c_2A_2$, after applying SymPy factor, structurally with the integer zero. The [saved diagnostic](residual-inspection.json) records:

| Component | Internal factored expression | Structural equality with zero | Symbolic is-zero property |
|---|---|---|---|
| $x$ | Integer(0) | true | true |
| $\phi$ | Mul(Integer(0), I) | false | true |

Thus the second component was the unevaluated expression $0\,i$. An independent exact basis solve reproduces the original coefficients

$$
c_1=-\frac{r^2+5\epsilon/6}{r^2+\epsilon},
\qquad
c_2=-\frac{i\epsilon}{18(r^2+\epsilon)}.
$$

The correction divides the first residual component by $\sqrt{1-x^2}$, multiplies the second by that factor, and cancels the resulting rational expressions exactly. Both are identically zero on $-1<x<1$. Equality of the smooth scalar-generated vector field and retained vector harmonics extends across the coordinate poles. No numerical tolerance, altered descriptor, or altered coefficient is used.

The [original failure record](../oblate-constraint-rank-2026-11-27/lapse-failure.json), [failed log](../oblate-constraint-rank-2026-11-27/lapse-failed-check.log), and partial results are unchanged. The original verifier is pinned by the parent commit. The corrected [source](../../verification/oblate_constraint_rank_20261127/lapse_prolongation.py) writes new outputs to this separate report directory.

The approved rerun exited successfully. Its [log](verification.log), [summary](lapse-summary.json), [exact coefficients](lapse-leading.json), and [environment/source hashes](lapse-reproducibility.json) are archived here. The checker took 1.717 seconds internally; the local command took about 2.31 seconds. These are observed timings, not a performance guarantee. No GPU, remote compute, or tracked simulation experiment was launched.

Reproduce from the paper1 folder:

    uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 python verification/oblate_constraint_rank_20261127/lapse_prolongation.py

The checker consumes the existing exact leading shift certificate; it does not regenerate or overwrite that archive.

## 2. Linearized reduction at the exact solution

Use $L=N-1$, $x=\cos\theta$, $d=1-x^2$, $\rho=r^2+\epsilon$, and $\Sigma=r^2+\epsilon x^2$. The static spatial metric is

$$
\sigma=\frac{\Sigma}{\rho}dr^2+\frac{\Sigma}{d}dx^2+\rho d\,d\phi^2,
\qquad a=-1,\quad b=V=W=0,\quad p=0.
$$

It is flat and belongs to the actual retained ansatz. Work on a fixed compact shell $I=[r_-,r_+]$ with $r_->0$ and positive $\rho,\Sigma$. Initially all integrations by parts use smooth compactly supported radial tests. Use the established normalized action and its integrated DeWitt bilinear form $G$.

Let $T$ map retained spatial coefficient perturbations to tensor perturbations, $\Pi$ be its $G$-orthogonal projector, and $N_\perp=1-\Pi$. For a retained shift coefficient test $s$ define

$$
TJ_s=\Pi\mathcal L_{v_s}\sigma,\qquad
n_s=N_\perp\mathcal L_{v_s}\sigma,\qquad
A(s,t)=-\frac12G(n_s,n_t).
$$

The previous checkpoint proves ranks $R=18,26$ and equality of $\ker A$ with $\ker n$ at small nonzero deformation. Choose a smooth complement $R$ to that kernel. Elimination of these nondegenerate shift primary–secondary pairs gives

$$
K=2(T^*GT)^{-1},\qquad
K_{\rm eff}=K-J_RA_R^{-1}J_R^*.
$$

Here adjoints include radial integration by parts. Although $n_s$ and $A$ contain no radial derivatives of shift profiles, $J_s$ does contain them.

For a retained test $f$ of the actual coefficient $a=h_{tt}$, the linearized lapse secondary and its preservation condition are

$$
L_fh=\frac12G(\nabla^2f,Th),\qquad
\mathcal T_f(p)=L_fK_{\rm eff}p.
$$

The factor $1/2$ follows from $\alpha_{,a}=-1/2$ at this solution. A shift-kernel section $k$ obeys $TJ_k=\mathcal L_{v_k}\sigma$ exactly. Flatness and spatial diffeomorphism invariance give $L_fJ_k=0$ and $V''J_k=0$, with $V''$ the unit-lapse curvature-potential Hessian. Consequently null shift secondaries have identically vanishing preservation in this linearized system. This is not a nonlinear first-class assertion.

Let $P_{\mathcal N}$ be the $G$-orthogonal projector onto the nondegenerate normal shift image. Set

$$
Z_f=(1-P_{\mathcal N})N_\perp\nabla^2f.
$$

The exact reduced bracket is

$$
\boxed{\{L_fh,\mathcal T_g\}_{\rm lin}
=L_fK_{\rm eff}L_g^*
=-\frac12G(Z_f,Z_g).}
$$

To obtain this identity, the uneliminated kinetic term supplies $\tfrac12G(\Pi\nabla^2f,\Pi\nabla^2g)$, while shift elimination supplies $\tfrac12G(P_{\mathcal N}N_\perp\nabla^2f,P_{\mathcal N}N_\perp\nabla^2g)$. The compact-support flat identity $G(\nabla^2f,\nabla^2g)=0$ then gives the displayed expression. The indefinite trace term and all radial derivatives are retained.

## 3. Exact redundant families

For any retained scalar $f$ whose spatial gradient is supplied by retained temporal shift coefficients, use

$$
b=\partial_rf,\qquad W=\gamma^{-1}d_\Omega f.
$$

Then $v_s=\operatorname{grad}_\sigma f$, so $\mathcal L_{v_s}\sigma=2\nabla^2f$. It follows that $Z_f=0$ as a tensor. More precisely,

$$
\mathcal T_f=\frac12J_{s_{\ker}}^*p
$$

after splitting $s$ into the chosen shift complement and kernel. Thus this tertiary is redundant modulo the remaining shift secondaries. The lapse secondary $L_fh=0$ must still be imposed.

At both cutoffs, every angle-independent $f=F_0(r)$ is such a descriptor. At $L=3$ there are three further real families, obtained from

$$
f_z=F_z(r)\left(r^2x+\frac{\epsilon x^3}{3}\right),
$$

$$
f_\pm=F_\pm(r)\sqrt d
\left(r^2+\frac{2\epsilon}{3}+\frac{\epsilon x^2}{3}\right)e^{\pm i\phi}.
$$

Real and imaginary parts of the $m=1$ template supply two real families. Each radial coefficient is arbitrary smooth and compactly supported. Differentiating the radial coefficient does not change angular degree.

For $f_z$, $\partial_x(f_z/F_z)=\Sigma$ and $W=F_zd\,\partial_x$. With $Y_\pm=\sqrt d\,e^{\pm i\phi}$, the other fields satisfy

$$
W=F_\pm\left(\operatorname{grad}_{S^2}Y_\pm
\mp\frac{i\epsilon}{3\rho}Y_\pm\partial_\phi\right).
$$

The second term is a scalar degree-one harmonic times a rotational Killing vector and has degree at most two. The corrected checker also verifies its exact polar-degree-one/axial-degree-two decomposition in the stated associated-Legendre convention. The descriptor itself has degree at most three. Hence all coefficients are retained at $N=4$.

These give one real radial-profile function at $N=3$ and four at $N=4$. They are not counts of nonlinear gauge functions.

## 4. Complementary radial brackets

The verifier differentiates the covariant Hessian and subtracts the changing spatial tangent before round projection. The $F''$ contribution is $Y_{\ell m}\,dr^2$ and is exactly tangent at every deformation. The quotient field $Z_f$ therefore depends only on $F,F'$.

In each magnetic sector the checker projects the leading Hessian defect away from the independent leading normal shift columns, using their exact indefinite DeWitt Gram inverse. The inherited shift certificate makes this inverse well defined. It verifies the resulting orthogonality, every cross block, and the radial rank-one forms.

For $f=F(r)Y_{\ell m}$, the leading form is

$$
G(Z_f,Z_g)=c_{\ell m}\epsilon^2
\int_I r^{-6}(aF+brF')^*(aG+brG')\,dr+O(\epsilon^3).
$$

The bracket has the additional factor $-1/2$. In the unnormalized associated-Legendre convention the exact coefficients are:

| $N$ | $\ell$ | $(a,b)$ | $c_{\ell m}$ for $m=0,\ldots,\ell$ |
|---|---:|---|---|
| 3 | 1 | $(1,-1)$ | $-16/525,\ -64/1575$ |
| 3 | 2 | $(0,1)$ | $-4/245,\ 256/40425,\ 64/1323$ |
| 4 | 2 | $(0,1)$ | $-4/245,\ -4/49,\ -8/49$ |
| 4 | 3 | $(1,1)$ | $-20/1617,\ 43200/2523367,\ 81600/220759,\ 332275200/280832783$ |

These are complex Hermitian sector coefficients; each $m>0$ sector supplies two real angular channels. Every displayed coefficient is nonzero. The two top multiplets therefore supply $4L=8,12$ real radial-profile channels. The summary JSON field named leading_active_scalar_multiplets_real counts these channels, not irreducible multiplets.

The mixed signs are permitted by the DeWitt form. They do not imply a zero direction within any individual radial block, and they are not a positivity or stability theorem for the physical Hamiltonian.

## 5. Exact small-deformation statement and radial domain

The previous section alone is a leading-order certificate. The following analytic argument promotes it to a statement about the exact reduced form for sufficiently small nonzero $\epsilon$ on the fixed shell.

At $L=2$ take the complement of the monopole. At $L=3$, the degree-one projection of the three deformed dipole templates is an invertible smooth radial matrix near $\epsilon=0$: its limit is $r^2$ times an invertible constant matrix. Thus a complement to the four exact families is the sum of the two undeformed top multiplets.

The background and the two projections are axisymmetric and invariant under equatorial reflection. At fixed $m$ these two top degrees have opposite reflection parity. Each complementary scalar channel is therefore separate from the other, also at finite deformation. Invariance under $\phi\mapsto-\phi$ also makes the cosine and sine copies identical and uncoupled.

For one channel, divide the reduced form by $\epsilon^2$ and multiply by the sign of $c_{\ell m}$. The leading term controls the $H^1_0(I)$ norm. Indeed $aF+brF'$ is proportional to $r(F'+\eta F/r)$, with $\eta=-1,0,1$. Multiplication by $r^\eta$ reduces the estimate to the Dirichlet Poincare inequality; all radial weights and their inverses are bounded on this shell.

For the remainder, divide the independent exact normal shift columns by $\epsilon/r^2$ before forming their Gram inverse. Their limiting Gram is nonsingular by the earlier exact minors, so this inverse and the quotient projection depend smoothly and uniformly on $\epsilon$ on the compact shell. Since $Z_f$ depends only on $F,F'$, its finite angular coefficient matrices give

$$
\left|\epsilon^{-2}G(Z_f,Z_f)-G_0^{\rm red}(F,F)\right|
\le C|\epsilon|\,\|F\|_{H^1(I)}^2.
$$

For sufficiently small nonzero $\epsilon$, each signed complementary block remains coercive. It is consequently an invertible weak operator $H^1_0(I)\to H^{-1}(I)$. Taken together, the complementary blocks are nondegenerate. The kernel of the reduced operator on the retained scalar space with $H^1_0(I)$ radial profiles consists exactly of the section 3 templates with profiles in $H^1_0(I)$, including the limits of compactly supported smooth profiles.

This is a statement about the continuous extension of the reduced form initially derived on smooth compactly supported tests. It does not apply the unrestricted Hessian integration-by-parts identity to arbitrary boundary traces. In particular, defining this weak form does not automatically supply a common domain for the original canonical generators, their adjoints, and physical evolution. Unrestricted radial profiles, boundary charges, homogeneous solutions, and solvability with physical boundary data require separate treatment. No explicit numerical bound on the allowed deformation is claimed.

## 6. Constraint interpretation and next question

After eliminating the 18/26 nondegenerate shift pairs, the complementary linearized lapse chains have the structure

$$
\pi_a\ \longrightarrow\ Lh\
\longrightarrow\ \mathcal T(p)\
\longrightarrow\ \mathcal Q(h,\delta a).
$$

The derivative of $\mathcal Q$ with respect to the complementary lapse perturbation is, up to the constraint-sign convention, the nondegenerate bracket in section 5. Thus it determines that lapse perturbation in a radial realization admitting the required inhomogeneous solutions; subsequent preservation determines its rate. The exact low lapse tertiaries are combinations of null shift secondaries.

This establishes additional independent linearized lapse preservation conditions: the complementary bracket cannot arise from combinations of null shift secondaries because every $L_f$ commutes with those secondaries. It does not turn the remaining background-specific null directions into nonlinear gauge symmetries. A full degree-of-freedom count also needs a common boundary realization and independence of the complete stabilized set.

The approved correction and rerun are complete. The next proposed research step is to test the nonlinear preservation of these remaining null directions near the oblate solution, on the actual constraint surface with radial dependence retained. That is a new stage, not a conclusion of this verifier. It requires the user's next copilot decision.
