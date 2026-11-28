# Independent mixed-H6 bilinear and exterior-tail audit

This note contains hand derivations and a read-only audit of
[the bilinear/tail reduction](sixth-order-bilinear-and-tail-reduction.md),
using the definitions in
[the sixth-order recursion](reduced-sixth-order-and-loop-kernel.md).
No numerical field evaluation, scientific run, or executable
implementation validation was performed to obtain these results.
All constants below use the raw width-one TT fields and $16\pi G=1$.
They do not include the factors converting normalized canonical coordinates
to those fields.

The proposed reduction is algebraically complete: it leaves two inverse
cross pairings and one local three-HEAT integral without a Gaussian
multiplier. The inverse cross pairings can be evaluated without reconstructing
their inverse fields. The three-HEAT integral has a genuine algebraic tail.
The formulas remain conditional on the local reduced construction and its
stated decaying inverses. They establish neither a measured mixed-H6
coefficient nor a gravitational four-point amplitude or a uniform-rank limit.

## 1. Pointwise fourth-order momentum scalar

Use physical spatial tensors, including their Gaussian envelopes, and write
$$
W=\mathcal A_0^{-1}b_1,\qquad T=L_0W,\qquad T_1=L_1W,
\qquad r=b_2-\mathcal A_1W,\qquad W_2=\mathcal A_0^{-1}r.
$$
The vector components of $W,W_2$ are contravariant, whereas $b_1,r$ are
covectors. All integrated pairings use the flat coordinate volume.
Set
$$
g=\mathsf G_2a,\qquad a=\mathsf G_0^{-1}p,\qquad
E_c=\sum_\alpha e_\alpha(\mathsf G_0^{-1}g)^\alpha.
$$
The raw tangent terms are
$$
E_1=\tfrac12(QS+SQ),\qquad
E_2=\tfrac16(Q^2S+QSQ+SQ^2).
$$
Substitution into the original momentum expansion gives
$$
U_1=-\tfrac12[Q,S]-T,
$$
$$
U_2=\tfrac16[Q,[Q,S]]+QT-T_1-L_0W_2-E_c.
$$
For symmetric $Q,S,T$, the commutator $[Q,S]$ is antisymmetric, so
$\operatorname{tr}([Q,S]T)=0$. Cyclicity also gives
$$
\operatorname{tr}(S[Q,[Q,S]])=-\operatorname{tr}[Q,S]^2.
$$
Therefore $\operatorname{tr}(U_1^2)+2\operatorname{tr}(SU_2)$ is
$$
\boxed{
s_4=T:T+2\operatorname{tr}(SQT)-2S:T_1
-\frac1{12}\operatorname{tr}[Q,S]^2
-2S:L_0W_2-2S:E_c.}
\tag{A1}
$$
In particular, the raw coefficient is $1/4-1/3=-1/12$.
It would be wrong to replace $\operatorname{tr}(U_1^2)$ with the
Frobenius square of $U_1$, which is generally nonsymmetric.

The Gram contraction can be obtained from the first response alone:
$$
g_\alpha=-\frac1{12}\int
\operatorname{tr}([Q,e_\alpha][Q,S])
+2\langle b_1(Q,e_\alpha),W(Q,S)\rangle.
\tag{A2}
$$
Both $g^{\mathsf T}\mathsf G_0^{-1}g$ in the integrated sixth-order
momentum energy and $-2S:E_c$ in (A1) are required. Keeping one does not
replace the other.

## 2. Removal of the second York field

Put
$$
u=\mathcal C(s_2-R_2),\quad v=\mathcal C(R_2+7s_2),\quad
w=7u+v,\quad t_j=S_{ij}\partial_iw,\quad
\mathcal C=\tfrac18(-\Delta)^{-1}.
$$
The coefficient multiplying $s_4$ in the adjoint energy is $-w$.
Flat tracelessness and transversality of $S$ imply
$$
\int w S:L_0W_2
=2\int wS_{ij}\partial_iW_2^j
=-2\int S_{ij}(\partial_iw)W_2^j.
$$
The boundary term vanishes because $S$ is polynomial-Gaussian.
Consequently
$$
2\int w S:L_0W_2=-4\langle t,W_2\rangle.
$$
Combining this with the fourth-order York norm proves
$$
\boxed{
-2\langle r,\mathcal A_0^{-1}r\rangle
-4\langle t,\mathcal A_0^{-1}r\rangle
=-2\langle r+2t,\mathcal A_0^{-1}r\rangle.}
\tag{A3}
$$
Only one vector inverse cross pairing remains. Neither values nor spatial
derivatives of $W_2$ are needed by this energy formula. The cross pairing
does not have a definite sign merely because the York self form is negative.

For completeness, the compact first York variation is also correct:
$$
(\mathcal A_1W)_i
=\partial_j(T_1)_{ij}-Q_{jk}\partial_kT_{ij}
-\tfrac12(\partial_iQ_{jk})T_{jk}.
\tag{A4}
$$
The contracted connection $\Gamma_1{}^\ell{}_{jj}$ vanishes for flat-TT
$Q$. In the other connection term, the first and third derivatives of
$Q$ cancel after exchanging the two indices contracted with symmetric
$T$. The surviving derivative is the last term in (A4).

## 3. Scalar adjoint reduction and the complete energy

Since $D_1f=-\partial_i(Q_{ij}\partial_jf)$ and $\partial_iQ_{ij}=0$,
$$
F=s_3-R_3+8D_1u=-2S:T-R_3-8Q:\nabla^2u,
$$
$$
K=R_3+7s_3+8D_1v=R_3-14S:T-8Q:\nabla^2v.
\tag{A5}
$$
Eliminating $u_4$ by scalar and $D_1$ adjointness leaves
$-\langle K,\mathcal CF\rangle$. Further,
$$
-8\langle v,D_2u\rangle
=4\int(\nabla v)^{\mathsf T}Q^2\nabla u.
$$
Let $s_4^{(0)}$ denote the right-hand side of (A1) with the
$-2S:L_0W_2$ term omitted. The complete result is
$$
\begin{aligned}
H_6={}&-\mathcal M_4+g^{\mathsf T}\mathsf G_0^{-1}g
-\int\mathcal L_6-4\langle b_3,W\rangle
+2\langle W,\mathcal A_2W\rangle\\
&+\int R_4(v-u)-\int w s_4^{(0)}
-2\langle r+2t,\mathcal A_0^{-1}r\rangle-\langle K,\mathcal CF\rangle\\
&+4\int(\nabla v)^{\mathsf T}Q^2\nabla u
+\int v(R_2+7s_2)u+28\int s_2u^2,
\end{aligned}
\tag{A6}
$$
where
$$
\mathcal M_4=\frac1{360}\int\operatorname{tr}([Q,[Q,S]]^2).
$$
Thus equation (5) of the owning reduction is confirmed, including its
canonical correction, both inverse cross pairings, and all scalar signs.
It retains $q^6,q^4p^2,q^2p^4,p^6$.

The weak $\mathcal A_2$ identity, equation (6) there, is also correct.
Holding the contravariant $W$ fixed, set
$$
J_1=T_1-QT,\qquad J_2=T_2-QT_1+\tfrac12Q^2T,\qquad T_2=L_2W.
$$
The order-two term of the exact negative norm gives
$$
2\langle W,\mathcal A_2W\rangle
=-\int\operatorname{tr}(J_1^2+2TJ_2)
$$
$$
=-\int\left[
T_1:T_1+2T:T_2-4\operatorname{tr}(QTT_1)
+\operatorname{tr}(Q^2T^2)+\operatorname{tr}(QTQT)
\right].
\tag{A7}
$$
The last two traces are distinct. Matrix order must be preserved.

Two reductions check the scalar bookkeeping. At $Q=0$, all York and
curvature terms disappear, $v=7u$, and (A6) gives
$(49+28)\int s_2(\mathcal Cs_2)^2$. At $S=0$, $u=-v$, $F=-K$,
and it becomes the stated configuration Schur-complement expression,
including $-\int R_2v^2$ and
$-4\int(\nabla v)^{\mathsf T}Q^2\nabla v$.

## 4. Smallest production evaluator

The following separation requires no production reconstruction of $W_2$,
$u_3$, or $u_4$.

| Contribution | Required representation |
| --- | --- |
| $\mathcal M_4,\int\mathcal L_6$ | Local polynomial-Gaussian moments |
| $g$ and its finite Gram contraction | First-field pairings and the retained-basis matrix |
| $b_3W$, weak $\mathcal A_2$, $R_4(v-u)$, $-w(s_4^{(0)}-T:T)$ | Positive Gaussian multiplier; at most two first HEAT fields |
| $4\nabla v\,Q^2\nabla u$, $v(R_2+7s_2)u$, $28s_2u^2$ | Positive Gaussian multiplier; two first HEAT fields |
| $-\langle K,\mathcal CF\rangle$ | One scalar inverse cross bilinear |
| $-2\langle r+2t,\mathcal A_0^{-1}r\rangle$ | One York inverse cross bilinear |
| $-\int w\,T:T$ | Three first HEAT fields; no Gaussian multiplier |

The inverse endpoints have a particularly simple structure. Each of
$F,K,r,t$ is a finite sum of a polynomial times $\exp(-3r^2/2)$ and
terms $P(x)\exp(-r^2/2)H(x)$, where $H$ is a first HEAT field of rate one.
Its defining representation
$$
H(x)=\sum_{p,m}c_{pm}x^p\int_0^1t^{2m}e^{-t^2r^2}\,dt
$$
turns each endpoint into a Gaussian mixture with rates
$$
\alpha(t)=\tfrac12+t^2\in[\tfrac12,\tfrac32].
$$
The inverse can be contracted analytically between each pair of
polynomial-Gaussians, before performing at most two heat-parameter
integrations. The positive lower bound on these rates makes the parameter
integrands smooth; no singular heat endpoint or extra inverse parameter
is required. This procedure is an all-space inverse pairing, not an
approximation imposing boundary conditions on an intermediate potential.

With $\widehat f(k)=\int e^{-ik\cdot x}f(x)\,dx$, the bilinear conventions are
$$
B_C(f,g)=\frac1{8(2\pi)^3}\int
\frac{\widehat f(-k)\widehat g(k)}{k^2}\,d^3k,
$$
$$
B_Y(b,c)=\frac1{(2\pi)^3}\int
\widehat b_i(-k)\left(-\frac{\delta_{ij}}{k^2}
+\frac{k_ik_j}{4k^4}\right)\widehat c_j(k)\,d^3k.
\tag{A8}
$$
For real sources these use the sum of real-real and imaginary-imaginary
Fourier products. For future complex polarizations, the continuation is
bilinear: do not conjugate the source coefficients.

For $P(x)e^{-\alpha r^2}$ the Fourier polynomial is obtained from
$$
\widehat{x^p e^{-\alpha r^2}}
=i^{|p|}\partial_k^p\left[
(\pi/\alpha)^{3/2}e^{-k^2/(4\alpha)}\right].
$$
For two rates $\alpha,\beta$, put
$\gamma=(\alpha+\beta)/(4\alpha\beta)$. If every $p_i$ is even and
$|p|=2n$, the needed moment is
$$
\int\frac{k^p e^{-\gamma k^2}}{(k^2)^s}\,d^3k
=2\pi\,\frac{\prod_i(p_i-1)!!}{(2n+1)!!}
\Gamma(n+\tfrac32-s)\,
\gamma^{-(n+3/2-s)}.
\tag{A9}
$$
Odd moments vanish. For $s=1$, $n\ge0$ is allowed; for $s=2$, $n\ge1$
is necessary. The two divergence factors supply that minimum degree in
the York longitudinal term. A constant numerator over $k^4$ is not an
integrable moment and must not be silently discarded.

The Fourier Gaussian prefactor is $1/[8(\alpha\beta)^{3/2}]$; the scalar
$\mathcal C$ contributes a further factor $1/8$.
After cancellation of half-integer powers, every polynomial kernel is
$\pi^{3/2}/\sqrt{\alpha+\beta}$ times a rational Laurent polynomial in
$\alpha,\beta,(\alpha+\beta)^{-1}$. An independent Decimal implementation
therefore needs polynomial operations and positive square roots, rather
than the HEAT moment evaluator.

Useful normalization targets, for $g_\alpha=e^{-\alpha r^2}$, are
$$
\langle g_\alpha,\mathcal Cg_\beta\rangle
=\frac{\pi^{3/2}}{16\alpha\beta\sqrt{\alpha+\beta}},
$$
$$
\langle c g_\alpha,\mathcal A_0^{-1}(d g_\beta)\rangle
=-\frac{11\pi^{3/2}}{24\alpha\beta\sqrt{\alpha+\beta}}\,c\cdot d.
\tag{A10}
$$
The second follows from the angular multiplier $-1+1/12=-11/12$.
Derivative controls are
$$
\sum_i\langle\partial_if,\mathcal C\partial_ig\rangle
=\tfrac18\int fg,\qquad
\langle\nabla f,\mathcal A_0^{-1}\nabla g\rangle
=-\tfrac34\int fg.
$$
For a vector self pairing, with
$E=\sum_i\langle b_i,(-\Delta)^{-1}b_i\rangle$,
$$
-E\le \langle b,\mathcal A_0^{-1}b\rangle\le-\tfrac34E.
$$
These bounds do not apply as sign restrictions on arbitrary cross pairings.

## 5. Audit of the exact exterior representation

Only $-\int wT:T$ lacks a Gaussian multiplier. The following confirms
equations (8)--(12) of the owning reduction; it does not validate a numerical
implementation of them.

For $r>0$, repeated use of
$2r^2I_m=(2m-1)I_{m-1}-e^{-r^2}$ gives
$$
I_m(r^2)=c_m\sqrt\pi\,\frac{\operatorname{erf}r}{r^{2m+1}}
-e^{-r^2}\sum_{j=1}^m d_{mj}r^{-2j},
$$
$$
c_m=\frac{(2m-1)!!}{2^{m+1}},\qquad
d_{mj}=\frac{(2m-1)!!}{2^j[2(m-j)+1]!!}.
\tag{A11}
$$
For example,
$$
I_2(r^2)=\frac{3\sqrt\pi\,\operatorname{erf}r}{8r^5}
-e^{-r^2}\left(\frac1{2r^2}+\frac3{4r^4}\right).
$$
Thus both the recurrence indexing and its constants are correct.

The angular collapse uses sphere averages, so the signed radial density is
$$
E(r)=4\pi\sum_{k=0}^3\pi^{k/2}(\operatorname{erf}r)^k
e^{-(3-k)r^2}\sum_\nu A_{k\nu}r^\nu.
\tag{A12}
$$
Its coefficients are rational at the present common rate one. Merge them
exactly. Since $w=O(r^{-1})$ and $T=O(r^{-2})$, the complete radial density
is $O(r^{-3})$. Consequently $A_{3\nu}=0$ for $\nu>-3$ is an appropriate
exact acceptance condition. The genuine $\nu=-3$ coefficient must remain.
No failed coefficient may be deleted to manufacture the expected falloff.

Replacing the three error functions by one gives the algebraic exterior:
$$
\int_R^\infty E_{\rm alg}(r)\,dr
=4\pi^{5/2}\sum_{\nu\le-3}
\frac{A_{3\nu}R^{\nu+1}}{-\nu-1}.
\tag{A13}
$$
For $J_\nu(a,R)=\int_R^\infty r^\nu e^{-ar^2}\,dr$,
$$
\left|\int_R^\infty(E-E_{\rm alg})\,dr\right|
\le12\pi^2\sum_\nu|A_{3\nu}|J_{\nu-1}(1,R)
+4\pi\sum_{k=0}^2\pi^{k/2}
\sum_\nu|A_{k\nu}|J_\nu(3-k,R).
\tag{A14}
$$
The first prefactor is $4\pi\pi^{3/2}\cdot3/\sqrt\pi=12\pi^2$,
using $1-\operatorname{erf}(r)^3\le3\operatorname{erfc}(r)$.
The other terms use $|\operatorname{erf}r|^k\le1$.

Integration by parts gives
$$
2aJ_\nu=R^{\nu-1}e^{-aR^2}+(\nu-1)J_{\nu-2}.
$$
For $\nu>1$, $J_{\nu-2}\le R^{-2}J_\nu$; for $\nu\le1$, the
second term is nonpositive. Therefore
$$
J_\nu(a,R)\le
\frac{R^{\nu-1}e^{-aR^2}}
{2a-\max(\nu-1,0)/R^2}
\tag{A15}
$$
is correct for $R>0$ and a positive denominator.

Use the original regular $I_m$ representation near the origin. For the
exterior, (A13) plus the bound (A14) avoids replacing an algebraic tail
with a Gaussian estimate. Floating evaluation of the rational coefficients,
signed algebraic sum, and logarithmic bounds introduces additional error;
the analytic inequality alone does not certify those floating operations.

## 6. Independent leading-tail normalization

Let
$$
M_i=\int(b_1)_i,\qquad
M_w=\int(14s_2-6R_2).
$$
The compact flat-TT expression for $b_1$ gives
$$
M_i=-\tfrac12\int S_{jk}\partial_iQ_{jk}.
\tag{A16}
$$
The other two derivative terms integrate to zero by transversality.
The scalar and York kernels then imply, with $n=x/r$,
$$
w=\frac{M_w}{32\pi r}+O(r^{-2}),
$$
$$
W_i=-\frac{7M_i+n_i(n\cdot M)}{32\pi r}+O(r^{-2}),
$$
$$
T_{ij}=\frac3{16\pi r^2}
\left[n_iM_j+n_jM_i+(n_in_j-\delta_{ij})(n\cdot M)\right]
+O(r^{-3}).
\tag{A17}
$$
The square of the tensor in brackets is $2|M|^2+4(n\cdot M)^2$.
Hence
$$
\int_{S^2}T:T\,d\Omega
=\frac{15|M|^2}{32\pi r^4}+O(r^{-5}),
$$
and the leading exterior contribution is
$$
\boxed{
-\int_{r>R}wT:T\,d^3x
=-\frac{15M_w|M|^2}{2048\pi^2R^2}+O(R^{-3}).}
\tag{A18}
$$
The exact finite expansion above determines subsequent powers.
In the normalization of (A12),
$$
\boxed{
A_{3,-3}=-\frac{15M_w|M|^2}{4096\pi^{9/2}}.}
\tag{A19}
$$
This is an independent leading-coefficient test for the signed term
$-\int wT:T$, not a prescription to insert its expected answer.

Since $\int R_2=-\tfrac14\int|\nabla Q|^2$,
$$
M_w=14\int S:S+\tfrac32\int|\nabla Q|^2>0
$$
for nonzero data. Thus the generic negative $R^{-2}$ tail cannot cancel
against the Gaussian-localized families. Rearrangement can integrate it
analytically, but cannot justify dropping it.

## 7. Guaranteed nonzero-monopole direction and its sign

The electric and magnetic profiles here are exactly those defined in
[the raw basis implementation](../../verification/quantum_bridge/reduced_gaunt_exact.py).
Let $E(A)$ be the electric TT profile generated from a real STF label $A$
and $f(x)=e^{-r^2/2}$. Define
$$
B(B)_{ij}=\tfrac12\left(
\epsilon_{iak}\partial_aE(B)_{kj}
+\epsilon_{jak}\partial_aE(B)_{ik}\right),
\qquad \epsilon_{123}=+1.
\tag{A20}
$$
The argument of $B(\cdot)$ is a label matrix; the function denotes the
magnetic profile. This fixes the symmetric-curl sign.
Define the label vector
$$
V_i(A,B)=\epsilon_{ijk}A_{j\ell}B_{\ell k}.
\tag{A21}
$$
It satisfies $V(B,A)=-V(A,B)$ and $V(A,A)=0$.

Here is the normalization derivation without evaluating the spatial
polynomial dictionaries. For $n=k/|k|$, let $P=I-nn^{\mathsf T}$ and
$$
\Pi_n A=PAP-\tfrac12P\,\operatorname{tr}(PA).
$$
The Fourier electric field is
$$
\widehat{E(A)}(k)=|k|^4\Pi_nA\,(2\pi)^{3/2}e^{-|k|^2/2}.
$$
For $N(n)_{ik}=\epsilon_{iak}n_a$, the symmetric-curl symbol on a TT
tensor is
$$
\mathscr C_nX=\tfrac12(NX-XN)=NX.
$$
Thus
$$
\widehat{B(B)}(k)=i|k|^5\mathscr C_n(\Pi_nB)
(2\pi)^{3/2}e^{-|k|^2/2}.
$$
The transverse tensor contraction reduces to
$$
D(n)=(\Pi_nA):\mathscr C_n(\Pi_nB)
=-n\cdot V(A,B)+n\cdot(An\times Bn).
\tag{A22}
$$
Using the second and fourth sphere moments,
$$
\int_{S^2}n_iD(n)\,d\Omega
=-\frac{4\pi}{3}V_i+\frac{4\pi}{15}V_i
=-\frac{16\pi}{15}V_i.
$$
Parseval cancels the $(2\pi)^3$ from the squared seed transform.
The remaining radial moment is
$$
\int_0^\infty k^{12}e^{-k^2}\,dk
=\tfrac12\Gamma(13/2)=\frac{10395\sqrt\pi}{128}.
$$
It follows that
$$
\int B(B):\partial_iE(A)
=-\frac{693}{8}\pi^{3/2}V_i(A,B),
$$
and (A16) gives the raw monopole formula
$$
\boxed{
Q=E(A),\ S=B(B)
\quad\Longrightarrow\quad
M=\frac{693}{16}\pi^{3/2}V(A,B).}
\tag{A23}
$$

Now choose
$$
A=\operatorname{diag}(-\tfrac12,-\tfrac12,1),\qquad
X=\begin{pmatrix}0&0&1\\0&0&0\\1&0&0\end{pmatrix},
\qquad Q=E(A),\quad S=B(X).
$$
Then $V(A,X)=(0,3/2,0)$ and
$$
\boxed{M=(0,\tfrac{2079}{32}\pi^{3/2},0).}
\tag{A24}
$$
In particular its sign is positive along $y$. As a direct sign check,
equation (A22) is
$D(n)=-(3/2)n_y(1-n_z^2)$ for these two labels.

The raw electric Gram normalization is
$$
\int E(A):E(B)=\frac{189}{8}\pi^{3/2}\operatorname{tr}(AB).
$$
Its angular factor is $(8\pi/5)\operatorname{tr}(AB)$ and its radial
factor is $945\sqrt\pi/64$, providing a check on the same Fourier
normalization. The squared-gradient and magnetic norms are $11/2$ times
the corresponding electric norm. Therefore this direction has
$$
\int S:S=\frac{2079}{8}\pi^{3/2},\qquad
\int|\nabla Q|^2=\frac{6237}{32}\pi^{3/2},
$$
$$
\boxed{M_w=\frac{251559}{64}\pi^{3/2}.}
\tag{A25}
$$
The independent rational coefficient target for the signed tail is thus
$$
\boxed{
A_{3,-3}
=-\frac{15}{4096}\,
\frac{251559}{64}\left(\frac{2079}{32}\right)^2.}
\tag{A26}
$$
All of (A24)--(A26) are analytic targets, not measured results.
For amplitudes $Q=qE(A),S=sB(X)$, $M$ acquires the factor $qs$ and
$$
M_w=\frac{18711q^2+232848s^2}{64}\pi^{3/2}.
$$
Canonical-coordinate factors must be applied separately when that
normalization is desired.

### Historical directions do not test this leading tail

For general raw labels
$Q=E(A_q)+B(B_q)$ and $S=E(A_s)+B(B_s)$, the same-parity
overlaps in (A16) vanish by spatial parity, and the cross terms give
$$
M=\frac{693}{16}\pi^{3/2}
\big[V(A_q,B_s)+V(B_q,A_s)\big].
\tag{A27}
$$
The electric-axis direction has $M=0$. The historical
noncommuting-label direction has
$(A_q,B_q,A_s,B_s)=(A,X,X,A)$, so it also has
$M\propto V(A,A)+V(X,X)=0$. Noncommuting labels alone do not guarantee
a nonzero York monopole. Those directions can check other mixed terms,
but cannot establish the generic $R^{-2}$ tail normalization.
The proposed $(E(A),B(X))$ direction supplies the missing test.

## 8. Independent controls and remaining limits

Before accepting a mixed-H6 value:

1. Check the unequal-rate scalar/York kernels, odd Fourier sectors,
   gradient identities, symmetry, and the York self-form bounds in
   (A8)--(A10). Refine each heat parameter separately and compare selected
   blocks with independent precision; a shared moment routine is not an
   independent arithmetic check.
2. Check (A1) against the original $U_1,U_2$ expression, retain both
   canonical Gram contributions, and verify
   $\int s_3=0$ and $\int s_4=-a^{\mathsf T}\mathsf G_2a$.
3. On at least one genuinely mixed direction, compare the rearranged
   energy with explicit $W_2,u_4$ reconstruction. A radial-harmonic
   implementation is useful for this check, even if production uses
   Gaussian-mixture bilinears. Do not impose a decaying harmonic boundary
   condition on an auxiliary biharmonic potential; only the physical
   inverse fields have the prescribed decay.
4. Verify the exact cancellations in (A12), the nonzero target (A26),
   and matching of the regular and exterior formulas away from the origin.
   Refine the core quadrature and splitting radius independently, keeping
   the signed algebraic exterior and Gaussian remainder bound separate.
5. Recover the $77$ identity and configuration Schur complement; check
   momentum reversal, joint sixth-degree scaling, and rotations on
   held-out mixed directions.

Three small analytic targets distinguish tail retention from exact
cancellation. A single angular bin $(D,m,n,p)=(0,0,0,1)$ with coefficient
one has $A_{3,-3}=1/16$ and algebraic exterior
$\pi^{5/2}/(8R^2)$. The pair of bins
$(2,0,0,1)$ with coefficient one and $(0,0,0,0)$ with coefficient $-1/2$
instead uses $r^2I_1=(I_0-e^{-r^2})/2$ to give the exact radial density
$-(\pi^2/2)e^{-r^2}\operatorname{erf}(r)^2$.
Only $A_{2,0}=-1/8$ survives; all algebraic coefficients cancel.
The uncancelled single bin $(0,0,0,0)$ has power $\nu=-1$ and must fail
the H6 falloff gate.

For the declared $j=2$ free profiles, cubic scalar inverse endpoints have
total $j\le6$, vector endpoints have total $j\le6$ and Cartesian orbital
support at most $\ell=7$, and local sixth-order scalar densities have
$j\le12$. These are auxiliary support bounds, not authorization to
project products to the free $j=2$ space or a claim about a physical
large-$N$ limit. A finite number of checked directions is not a reconstructed
sixth-order tensor. The subsequent directional jets, loop diagram sum,
connected observable, and physical cutoff limits remain separate work.
