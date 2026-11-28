# Sixth-order constraint reduction and the finite-mode loop calculation

## Status

The [cubic and quartic calculation](reduced-vertices.md) has passed its
specified checks. This note derives the higher-order reduction identities
needed to continue to a complete one-loop four-point function of the same
candidate. The general fifth/sixth-order Gaunt coefficients and the
gravitational loop integral have **not yet been evaluated**.

All statements retain the local smoothness, invertible Gram matrix,
decaying York inverse and positive conformal solution hypotheses of the
[reduction-first construction](reduced-canonical-gaunt.md). They do not
establish equivalence to the historical finite-rank actions.

## 1. Only three scalar responses are needed through sixth order

Scale both canonical fields by a formal parameter, $(q,p)\mapsto
(\lambda q,\lambda p)$, and use homogeneous orders in $\lambda$.
Write

$$
\bar R=\sum_{n\geq2}R_n,\qquad
\sigma^2=\sum_{n\geq2}s_n,\qquad
\psi=1+\sum_{n\geq2}u_n,\qquad
\bar\Delta=\Delta+\sum_{n\geq1}D_n.
$$

The pointwise momentum scalar is denoted $s_n$ to distinguish it from
the extrinsic-curvature tensor. In particular $s_2=S:S$, where
$S=\sum_\alpha(\mathsf G_0^{-1}p)^\alpha e_\alpha$.
Since the conformal metric is unimodular,

$$
D_n f=\partial_i\left[\frac{(-Q)^n{}^{ij}}{n!}\partial_j f\right],
\qquad \mathcal C=\frac18(-\Delta)^{-1},\qquad B_n=R_n+7s_n.
$$

Expanding the scalar constraint gives

$$
\begin{aligned}
u_2&=\mathcal C(s_2-R_2),\\
u_3&=\mathcal C(s_3-R_3)+8\mathcal C D_1u_2,\\
u_4&=\mathcal C(s_4-R_4)
+8\mathcal C(D_1u_3+D_2u_2)-\mathcal C(B_2u_2).
\end{aligned}
$$

These signs follow directly by rewriting the equation as

$$
u=\mathcal C(s-\bar R)+8\mathcal C(\bar\Delta-\Delta)u
-\mathcal C[(\bar R+7s)u]+28\mathcal C(su^2)+\cdots.
$$

Let $I_n=\int s_n$. The exact energy
$H=\int[s\psi^{-7}-\bar R\psi]$ then yields

$$
\boxed{H_5=I_5-\int R_5-\langle B_3,u_2\rangle-\langle B_2,u_3\rangle,}
$$

$$
\boxed{
H_6=I_6-\int R_6-\langle B_4,u_2\rangle
-\langle B_3,u_3\rangle-\langle B_2,u_4\rangle
+28\int s_2u_2^2.}
$$

Neither $u_5$ nor $u_6$ enters: every conformal response in the energy
already multiplies a scalar beginning at order two. It is nevertheless
essential to retain the pointwise $s_3,s_4,R_2,R_3,R_4$; integrated
surface-derivative identities cannot replace these sources.

## 2. Integrated momentum terms from the inverse Gram matrix

Put $a=\mathsf G_0^{-1}p$ and expand
$\mathsf G=\mathsf G_0+\mathsf G_2+\mathsf G_3+\mathsf G_4+\cdots$.
The previously derived $\mathsf G_1=0$ gives

$$
\begin{aligned}
I_2&=p^{\mathsf T}\mathsf G_0^{-1}p,& I_3&=0,\\
I_4&=-a^{\mathsf T}\mathsf G_2a,&
I_5&=-a^{\mathsf T}\mathsf G_3a,\\
I_6&=-a^{\mathsf T}\mathsf G_4a+
a^{\mathsf T}\mathsf G_2\mathsf G_0^{-1}\mathsf G_2a.
\end{aligned}
$$

Thus the energy needs no pointwise construction of $s_5$ or $s_6$.
The last term in $I_6$ is a finite retained-basis matrix contraction;
discarding it would change the canonical momentum map.

## 3. The York response can also be reduced to two inverse solves

For the fixed direction $S$, define

$$
E_S=D\exp_Q[S]
=\sum_{n\geq0}\frac1{(n+1)!}\sum_{a=0}^nQ^aSQ^{n-a}.
$$

Let $\mathcal D_{\bar\gamma}T$ denote $\bar D^jT_{ij}$, with a covector
output, and let $L_{\bar\gamma}$ act on a contravariant vector.
Because $\det\bar\gamma=1$, its covariant tensor components can be
written without a volume-derivative term:

$$
(L_{\bar\gamma}W)_{ij}
=W^k\partial_k\bar\gamma_{ij}
+\bar\gamma_{ik}\partial_jW^k
+\bar\gamma_{jk}\partial_iW^k
-\frac23\bar\gamma_{ij}\partial_kW^k.
$$

Expand

$$
\mathcal A=\mathcal D_{\bar\gamma}L_{\bar\gamma}
=\mathcal A_0+\mathcal A_1+\mathcal A_2+\cdots,\qquad
b=\mathcal D_{\bar\gamma}E_S=b_1+b_2+b_3+\cdots.
$$

Here $\mathcal A_0=\Delta_L$. Flat transversality sets $b_0=0$.
The first two responses are

$$
W_1=\mathcal A_0^{-1}b_1,\qquad
r_2=b_2-\mathcal A_1W_1,\qquad
W_2=\mathcal A_0^{-1}r_2.
$$

All vector/covector pairings below use the coordinate volume $d^3x$.
Integration by parts gives

$$
\langle V,\mathcal A W\rangle
=-\frac12\langle L_{\bar\gamma}V,L_{\bar\gamma}W\rangle_{\bar\gamma}.
$$

Consequently $\mathcal A$ is symmetric in this pairing on the stated
decaying domain, and the subtracted longitudinal norm is

$$
\mathcal N=\|L_{\bar\gamma}\mathcal A^{-1}b\|_{\bar\gamma}^2
=-2\langle b,\mathcal A^{-1}b\rangle.
$$

Expansion of this scalar inverse, using symmetry, gives its directional
coefficients:

$$
\begin{aligned}
\mathcal N_2&=-2\langle b_1,W_1\rangle,\\
\mathcal N_3&=-4\langle b_2,W_1\rangle
+2\langle W_1,\mathcal A_1W_1\rangle,\\
\mathcal N_4&=-2\langle r_2,\mathcal A_0^{-1}r_2\rangle
-4\langle b_3,W_1\rangle
+2\langle W_1,\mathcal A_2W_1\rangle.
\end{aligned}
$$

The third vector response cancels from this expression. It is therefore
unnecessary to solve for $W_3$ when computing $\mathsf G_4$.
This saves a dependent-field inverse without dropping an interaction.

For the raw, unprojected metric norm $\mathcal M$, the matrix-exponential
identity gives

$$
\begin{aligned}
\mathcal M_2[S,S]&=-\frac1{12}\int\operatorname{tr}[Q,S]^2,\\
\mathcal M_3[S,S]&=0,\\
\mathcal M_4[S,S]&=\frac1{360}\int
\operatorname{tr}\big([Q,[Q,S]]^2\big).
\end{aligned}
$$

For example, diagonalizing the symmetric matrix $Q$ reduces the norm
factor for each component of $S$ to
$[\sinh(x/2)/(x/2)]^2=1+x^2/12+x^4/360+\cdots$.
Then

$$
\mathsf G_2=\mathcal M_2-\mathcal N_2,\qquad
\mathsf G_3=-\mathcal N_3,\qquad
\mathsf G_4=\mathcal M_4-\mathcal N_4,
$$

with matrix elements obtained by polarization of the directional
forms. Odd raw-norm coefficients vanish, but the York subtraction does
not justify setting the full $\mathsf G_3$ to zero.

## 4. The required pointwise momentum scalars

Let $E_{S,n}$ be the order-$Q^nS$ tangent term and $L_n$ the
order-$Q^n$ part of the longitudinal operator. Define

$$
Z_1=E_{S,1}-L_0W_1,
$$

$$
Z_2=E_{S,2}-L_0W_2-L_1W_1
-\sum_\alpha e_\alpha
(\mathsf G_0^{-1}\mathsf G_2a)^\alpha.
$$

The covariant conformal momentum tensor is
$\sigma_{ij}=S_{ij}+(Z_1)_{ij}+(Z_2)_{ij}+\cdots$.
With

$$
U_1=Z_1-QS,\qquad
U_2=Z_2-QZ_1+\frac12Q^2S,
$$

the pointwise sources are

$$
s_3=2\operatorname{tr}(SU_1)
=-2S:L_0W_1,\qquad
s_4=\operatorname{tr}(U_1^2)+2\operatorname{tr}(SU_2).
$$

The required consistency identities are $\int s_3=0$ and
$\int s_4=-a^{\mathsf T}\mathsf G_2a$. The first does not permit
discarding $s_3$ from $u_3$: it generally vanishes only after integration.

Together, Sections 1--4 give a finite evaluation prescription for the
full $H_5,H_6$, not only a momentum or configuration direction.
It still requires actual spatial integrations and retained-basis
contractions before numerical vertices are available.

## 5. Two analytic reductions that check the bookkeeping

### Pure momentum

At $q=0$, all curvature and higher Gram terms vanish:

$$
u_2=\mathcal C s_2,\quad u_3=0,\quad
u_4=-7\mathcal C(s_2u_2),\quad H_5=0.
$$

Green-function symmetry in the general expression then yields

$$
H_6=77\int s_2(\mathcal C s_2)^2,
$$

which recovers the coefficient independently evaluated in the preceding
constraint run. This is a normalization control for the general formula,
not a newly computed collection of mixed sixth-order vertices.

### Time-symmetric configuration data

At $p=0$, define the conformal Laplacian
$\mathcal L_Q=-8\bar\Delta+\bar R$ with a decaying inverse on sources.
Writing $\psi_0=1+w$ gives
$\mathcal L_Qw=-\bar R$. The exact configuration potential is therefore

$$
\boxed{V(q)=H(q,0)
=-\int\bar R+\langle\bar R,\mathcal L_Q^{-1}\bar R\rangle.}
$$

Put $v=\mathcal C R_2$. Expanding this single inverse produces

$$
V_5=-\int R_5+2\langle R_3,v\rangle
+8\langle v,D_1v\rangle,
$$

$$
\begin{aligned}
V_6={}&-\int R_6+2\langle R_4,v\rangle\\
&+\langle R_3+8D_1v,\mathcal C(R_3+8D_1v)\rangle
+8\langle v,D_2v\rangle-\int R_2v^2.
\end{aligned}
$$

Here $\langle v,D_1v\rangle=\int Q^{ij}\partial_i v\partial_jv$,
whereas $\langle v,D_2v\rangle=-\frac12\int
(Q^2)^{ij}\partial_i v\partial_jv$.
These independent Schur-complement expressions are checks on the
general scalar expansion.

## 6. Exact angular support at any fixed perturbative order

Suppose all free tensors have total angular momentum at most
$L=N-1$. Tensor products obey the Clebsch--Gordan upper bound; flat
derivatives and the decaying flat scalar/vector Green operators commute
with rotations. Induction on the above recursion therefore bounds the
total angular momentum of an order-$n$ dependent tensor by $nL$.
This statement concerns total angular momentum in the appropriate
tensor representation, not the orbital degree of a single Cartesian
component.

In particular, the required scalar responses $u_2,u_3,u_4$ have
support at most $2L,3L,4L$, and the vector responses $W_1,W_2$ at
most $2L,3L$. Complete angular contractions through sixth order are
finite. For the selected $N=3$ block, the scalar Green solves need
total angular momentum through eight and the vector solves through six.
The radial functions are not finitely many Gaussian polynomials after
inverse solves; radial accuracy must still be controlled separately.

This perturbative support bound does not cap the exact nonlinear
geometry at finite amplitude. It also does not justify reinserting an
intermediate projection at the free cutoff $L$.

## 7. A precise finite-mode one-loop kernel

Choose, explicitly, a midpoint/Weyl perturbative quantization in the
canonical variables, with canonical phase-space measure. Additional
order-$\hbar$ counterterms are separate specified inputs; they must
not be silently fitted to the desired amplitude.
This prescription is a quantum choice, not a consequence of the
classical symplectic embedding alone.

Normalize the free canonical modes as
$H_2=\frac12\sum_i(P_i^2+\omega_i^2Q_i^2)$ with $\omega_i>0$.
The Euclidean phase action has kinetic term $-iP\dot Q$.
With internal Fourier convention $e^{i\nu\tau}$, each free inverse
Hessian block is

$$
\mathcal G_i(\nu)=\frac1{\nu^2+\omega_i^2}
\begin{pmatrix}1&\nu\\-\nu&\omega_i^2\end{pmatrix}
$$

in the order $(Q_i,P_i)$. This phase-space matrix is not being
claimed to be a positive coordinate-space sampling covariance.

Give the four external legs distinct formal parameters and backgrounds
$z(\tau)=\sum_{r=1}^4\eta_r z_r e^{-i\Omega_r\tau}$ with
$\sum_r\Omega_r=0$. For each nonempty subset $S$ of the legs, let

$$
\mathcal B_S
=\left[\prod_{r\in S}\partial_{\eta_r}\,
\nabla_z^2H\!\left(\sum_r\eta_rz_r\right)\right]_{\eta=0}.
$$

This block uses $H_{|S|+2}$. If $\operatorname{OP}_k$ is the set of
ordered partitions of the four distinct legs into $k$ nonempty blocks,
the one-loop 1PI four-point coefficient, after stripping the overall
frequency-conservation delta function, is

$$
\boxed{
\Gamma^{(1)}_{1234}
=\frac{\hbar}{2}\sum_{k=1}^{4}\frac{(-1)^{k+1}}{k}
\sum_{(S_1,\ldots,S_k)\in\operatorname{OP}_k}
\int\frac{d\nu}{2\pi}\,
\operatorname{tr}\prod_{a=1}^{k}
\left[\mathcal G(\nu_a)\mathcal B_{S_a}\right],
}
$$

where $\nu_a=\nu+\sum_{b<a}\sum_{r\in S_b}\Omega_r$ in this
Fourier convention. Distinct external labels keep all differentiation
and cyclic multiplicities explicit.

The ordered-partition counts are $1,14,36,24$ for $k=1,2,3,4$.
They retain the sixth-order tadpole, fifth/cubic pair, quartic pair,
quartic/cubic/cubic terms and cubic box. No continuum scaleless-integral
argument is used to discard a finite-mode term.

The general Gaussian determinant framework and its spectral assumptions
are discussed by
[Vassilevich, PDF pp. 5, 7–8](https://www.alphaxiv.org/pdf/hep-th/0306138?page=7).
The displayed phase-space blocks, external-mode choice and finite-rank
Hamiltonian must be supplied and validated for this candidate; that
reference does not provide them.

The formula above is an **evaluation prescription**, not an evaluated
gravitational amplitude. In addition, a 1PI vertex alone is not the
one-loop scattering observable. The corrected two- and three-point
functions, exchange contributions, external-leg normalization,
analytic continuation, wave-packet/asymptotic limit and renormalization
must be assembled consistently.

## 8. Static-background check using the same constrained geometry

There is a useful independent check of the loop implementation. At fixed
$q$, expand in small momenta about the time-symmetric solution:

$$
H(q,p)=V(q)+\frac12p^{\mathsf T}\mathsf M(q)p+O(p^4).
$$

Let

$$
\mathsf J_{\alpha\beta}
=\int \psi_0^{-6}\,T_{\alpha\,ij}T_\beta^{ij},\qquad
\boxed{\mathsf M(q)=2\mathsf G^{-1}\mathsf J\mathsf G^{-1}.}
$$

To derive this, write $\psi=\psi_0+t^2\chi+O(t^4)$ for
$p\mapsto tp$.
The linearized scalar equation is
$(-8\bar\Delta+\bar R)\chi=\sigma_1^2\psi_0^{-7}$.
Multiplying by the homogeneous solution $\psi_0$, subtracting the
adjoint equation and integrating gives
$-8\oint\partial_r\chi=\int\sigma_1^2\psi_0^{-6}$.
This is exactly the quadratic momentum contribution to the ADM energy.

On a neighborhood where $\mathsf M$ and
$\mathsf K=\nabla_q^2V$ are positive, the static Gaussian loop is

$$
V^{(1)}(q)
=\frac{\hbar}{2}
\operatorname{tr}\sqrt{\mathsf M^{1/2}\mathsf K\mathsf M^{1/2}},
$$

up to the chosen background-independent vacuum subtraction and explicit
counterterms. The determinant of the phase Hessian is
$\det(\nu^2I+\mathsf M\mathsf K)$; this shows why an extra
field-dependent determinant must not be added twice when momenta are
integrated out.

Four configuration derivatives of this expression check only the
$qqqq$ 1PI component at $p=0$ and zero external frequency. They do not
check general phase-space components: for example the $qp^4$, $q^2p^4$
and $p^6$ vertices are invisible to this static check. They do not provide
the nonzero-frequency scattering amplitude or prove interacting
reflection positivity. The full
Hamiltonian's domain, ordering, positivity and limits remain part of the
original goal.

## 9. Next implementation and acceptance checks

The next implementation must supply the complete dependent-field
contractions in Sections 1--4, with all retained physical intermediate
modes. Required checks include recovery of the verified $H_3,H_4$,
the pure-momentum sixth coefficient and the independent static
Schur-complement identities; separate angular support and radial
convergence checks; and comparisons against finite-amplitude constrained
data with both configuration and momentum excited.

The loop integrator should first reproduce polynomial-oscillator controls
with all five vertex families present, including canonical coordinate
changes. Those controls validate bookkeeping only. The gravitational
result requires insertion of the candidate's full vertices and comparison
with its own tree/propagator calculation before attempting the Einstein
benchmark and the physical large-$N$ limit.

The complete local corner states, microscopic state counting,
interacting Euclidean--Lorentzian reconstruction and the observable map
to the historical simulations remain unresolved. No step here converts
the quadratic far-sphere area construction into those missing results.
