# Sixth-order energy as nested bilinears and first-response integrals

This is an analytic implementation reduction of the same local, reduction-first
Hamiltonian, not a measured sixth-order coefficient. It uses the conventions
and hypotheses of [the sixth-order recursion](reduced-sixth-order-and-loop-kernel.md)
and [the contracted four-point plan](contracted-four-point-implementation.md).
The fifth-order validation is a separate tracked calculation. No finite-rank
four-graviton amplitude or physical cutoff limit is established here.

## 1. Remove the explicit second York field from the energy

All fields below are physical spatial fields, including their Gaussian
envelopes. Let

$$
W=\mathcal A_0^{-1}b_1,\quad T=L_0W,\quad T_1=L_1W,
\quad r=b_2-\mathcal A_1W,\quad W_2=\mathcal A_0^{-1}r,
$$
$$
u=\mathcal C(s_2-R_2),\quad v=\mathcal C(R_2+7s_2),
\quad w=7u+v,\quad \mathcal C=(-\Delta)^{-1}/8.
$$
The canonical Gram correction is

$$
g_\alpha=-\frac1{12}\int\operatorname{tr}([Q,e_\alpha][Q,S])
+2\langle b_1(Q,e_\alpha),W\rangle,
\qquad E_c=\sum_\alpha e_\alpha(\mathsf G_0^{-1}g)^\alpha.
$$

Writing $[Q,S]=QS-SQ$, the pointwise expansion gives

$$
U_1=-\tfrac12[Q,S]-T,
\qquad
U_2=\tfrac16[Q,[Q,S]]+QT-T_1-L_0W_2-E_c.
$$

The trace of a product of an antisymmetric and a symmetric matrix vanishes.
Also $\operatorname{tr}(S[Q,[Q,S]])=-\operatorname{tr}[Q,S]^2$.
Consequently the complete fourth-order momentum scalar is

$$
s_4=s_4^{(0)}-2S:L_0W_2,
$$
$$
s_4^{(0)}=T:T+2\operatorname{tr}(SQT)-2S:T_1
-\frac1{12}\operatorname{tr}[Q,S]^2-2S:E_c.                 \tag{1}
$$

The weight of $s_4$ in the adjoint-reduced energy is $-w$. Since $S$ is
flat transverse and traceless, define $t_j=S_{ij}\partial_iw$ and integrate
by parts to obtain

$$
2\int wS:L_0W_2=-4\langle t,W_2\rangle.
$$

There is no surviving boundary term at infinity: $S$ and its derivatives
are polynomial-Gaussian, while $w$ and $W_2$ have decaying Green tails.
Combining this term with the fourth-order longitudinal norm gives

$$
-2\langle r,\mathcal A_0^{-1}r\rangle
-4\langle t,\mathcal A_0^{-1}r\rangle
=-2\langle r+2t,\mathcal A_0^{-1}r\rangle.                \tag{2}
$$

Thus production needs a vector cross bilinear, rather than the values and
derivatives of $W_2$ throughout space. Completing a square would also be
valid, but introduces a cancelling $t$ self-pairing that is unnecessary.

## 2. Localized sources for the two remaining inverse pairings

For the stated flat TT $Q$, symmetry and transversality reduce the first
York variation to

$$
(\mathcal A_1W)_i
=\partial_j(T_1)_{ij}-Q_{jk}\partial_kT_{ij}
-\tfrac12(\partial_iQ_{jk})T_{jk}.                         \tag{3}
$$

The contracted first connection vanishes, and two of the remaining three
Christoffel terms cancel against each other. Equation (3) retains the
contravariant convention for $W$ and covector convention for the output.

Similarly $D_1f=-Q_{ij}\partial_i\partial_jf$. The scalar test and source
are therefore

$$
F=-2S:T-R_3-8Q:\nabla^2u,
\qquad K=R_3-14S:T-8Q:\nabla^2v.                           \tag{4}
$$

Every term in $F,K,r,t$ is a polynomial-Gaussian or a
polynomial-Gaussian times one first-response HEAT field. Their Gaussian
multiplier has positive exponent. No projection back onto the free $j=2$
space is made. Products must retain their full auxiliary angular support.

The resulting energy is

$$
\begin{aligned}
H_6={}&-\mathcal M_4+g^{\mathsf T}\mathsf G_0^{-1}g
-\int\mathcal L_6-4\langle b_3,W\rangle
+2\langle W,\mathcal A_2W\rangle\\
&+\int R_4(v-u)-\int w s_4^{(0)}
-2\langle r+2t,\mathcal A_0^{-1}r\rangle-\langle K,\mathcal CF\rangle\\
&+4\int(\nabla v)^{\mathsf T}Q^2\nabla u
+\int v(R_2+7s_2)u+28\int s_2u^2.                         \tag{5}
\end{aligned}
$$

Here $\int R_6=\int\mathcal L_6$ uses the same decaying surface-term
identity as the earlier local-curvature calculation, and
$\mathcal M_4=\frac1{360}\int\operatorname{tr}([Q,[Q,S]]^2)$.
For the second York variation, put $T_2=L_2W$. Its weak form is

$$
2\langle W,\mathcal A_2W\rangle
=-\int\big[T_1:T_1+2T:T_2-4\operatorname{tr}(QTT_1)
+\operatorname{tr}(Q^2T^2)+\operatorname{tr}(QTQT)\big].     \tag{6}
$$

Equation (6) follows by expanding
$-\operatorname{tr}(e^{-Q}(L_{\bar\gamma}W)e^{-Q}(L_{\bar\gamma}W))$
through order $Q^2$, holding $W$ fixed. It avoids second derivatives of
$W$ in this integrated family.

At $Q=0$, (5) reduces to $77\int s_2(\mathcal Cs_2)^2$, with the
contributions $49+28$ retained separately. This is an algebraic control,
not a new numerical evaluation.

## 3. The single local term with an algebraic tail

In (5), every local first-response integral except

$$
-\int w\,T:T                                                     \tag{7}
$$

has a positive Gaussian multiplier and at most two HEAT factors. Term
(7) has three first HEAT factors and no Gaussian multiplier. Generically
$w=O(r^{-1})$ and $T=O(r^{-2})$, so its radial density is $O(r^{-3})$
and the omitted integral is $O(R^{-2})$. A Gaussian tail bound for the
whole of (7) would be incorrect.

For the present width-one basis all these first HEAT factors have rate one.
An exact angular collapse produces radial terms of the form

$$
4\pi c_{D m n p}\,r^{D+2}I_m(r^2)I_n(r^2)I_p(r^2),
\qquad I_m(z)=\int_0^1t^{2m}e^{-zt^2}\,dt .
$$

For $r>0$ the recurrence gives the exact identity

$$
I_m(r^2)=c_m\sqrt\pi\,\frac{\operatorname{erf}r}{r^{2m+1}}
-e^{-r^2}\sum_{j=1}^{m}d_{mj}r^{-2j},
$$
$$
c_m=\frac{(2m-1)!!}{2^{m+1}},\qquad
d_{mj}=\frac{(2m-1)!!}{2^j[2(m-j)+1]!!}.                    \tag{8}
$$

The sum is empty for $m=0$ and $(-1)!!=1$. Merge coefficients exactly
before evaluating the resulting exterior radial density:

$$
E(r)=4\pi\sum_{k=0}^{3}\pi^{k/2}(\operatorname{erf}r)^k
e^{-(3-k)r^2}\sum_\nu A_{k\nu}r^\nu.                       \tag{9}
$$

The $A_{k\nu}$ are rational. Any nonintegrable coefficient in the
$k=3$ part must cancel exactly; the physical falloff above requires
$A_{3\nu}=0$ for $\nu>-3$. This is an acceptance condition on the
assembled integrand, not permission to delete terms by hand.

Its algebraic exterior integral can then be included explicitly:

$$
\int_R^\infty E_{\rm alg}(r)\,dr
=4\pi^{5/2}\sum_{\nu\le-3}
\frac{A_{3\nu}R^{\nu+1}}{-\nu-1}.                         \tag{10}
$$

The remaining exterior error is exponentially bounded. Define
$J_\nu(a,R)=\int_R^\infty r^\nu e^{-ar^2}\,dr$.
Using $1-\operatorname{erf}(r)^3\le3\operatorname{erfc}(r)
\le3e^{-r^2}/(\sqrt\pi r)$ gives

$$
\begin{aligned}
\left|\int_R^\infty(E-E_{\rm alg})\,dr\right|
\le{}&12\pi^2\sum_\nu|A_{3\nu}|J_{\nu-1}(1,R)\\
&+4\pi\sum_{k=0}^{2}\pi^{k/2}
\sum_\nu|A_{k\nu}|J_\nu(3-k,R).                           \tag{11}
\end{aligned}
$$

For any integer $\nu$ and positive denominator, a convenient bound is

$$
J_\nu(a,R)\le
\frac{R^{\nu-1}e^{-aR^2}}
{2a-\max(\nu-1,0)/R^2}.                                   \tag{12}
$$

It follows by differentiating $r^{\nu-1}e^{-ar^2}$; for $\nu>1$ bound
$J_{\nu-2}\le R^{-2}J_\nu$, while for $\nu\le1$ drop the nonpositive
remainder. At the origin, use the original regular $I_m$ representation;
the singular pieces of (8) should never be evaluated separately there.
Floating evaluation of (10)--(12) is not automatically an interval bound.
Coefficient cancellation and numerical roundoff still require independent
precision checks before attaching an accuracy claim.

A separate leading-tail check uses the exact source moments
$M=\int b_1$ and $M_h=\int(14s_2-6R_2)$. The first-response asymptotics give

$$
w\sim\frac{M_h}{32\pi r},\qquad
\int_{S^2}T:T\,d\Omega\sim\frac{15|M|^2}{32\pi r^4}.
$$

Thus the exterior integral of (7) has the leading term

$$
-\frac{15M_h|M|^2}{2048\pi^2R^2},
\qquad
A_{3,-3}=-\frac{15M_h|M|^2}{4096\pi^{9/2}}.                \tag{13}
$$

In these conventions
$M_h=14\int S:S+\tfrac32\int|\nabla Q|^2$ is positive for nonzero data.
A required tail-control direction must have an exactly verified $M\ne0$.
The historical electric-axis direction has $M=0$; the historical
`noncommuting_labels` direction also has a vanishing monopole by its
electric/magnetic label pairing. Those two directions alone cannot
validate the generic nonzero coefficient in (13).

## 4. Remaining implementation and validation

The rearrangement and equations (5)--(12) have independent read-only
algebra checks. Equation (13) supplies a separate hand-derived leading-tail
control. No new scientific run has tested this representation. The existing
Gaussian pairing backend does not support the three-HEAT, zero-exponent
case, and its guard must not be bypassed.

Required checks include scalar and York adjointness, exact Gaussian
normalizations, the York negative-form spectral bounds, source localization,
all canonical Gram terms, held-out mixed directions, an explicit unreduced
$W_2/u_4$ comparison, momentum reversal and homogeneous scaling. The
algebraic exterior formula requires its own moment and tail controls.
Directional-jet contractions and all five loop families remain subsequent
work. These steps do not supply a physical large-$N$ limit, microscopic
corner-state dynamics, or equivalence to the archived Paper I regulator.
