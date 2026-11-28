# Independent H5 pairing and axial-contraction audit

This is an analytic audit of the proposed first-Green-field H5 assembly.
It uses the definitions in [the higher-order reduction](reduced-sixth-order-and-loop-kernel.md),
the cubic invariant in [the earlier vertex derivation](reduced-vertices.md),
and the explicit STF basis in
[reduced_gaunt_exact.py](../../verification/quantum_bridge/reduced_gaunt_exact.py).
No numerical calculation, experiment, or new coefficient evaluation was
performed. This note does not validate an implementation or complete the
gravitational four-point loop.

## 1. The proposed H5 formula is correct

All contractions in this section use flat Cartesian indices unless a
conformal metric is written explicitly. Put

$$
W=\mathcal A_0^{-1}b_1,\quad T_0=L_0W,\quad T_1=L_1W,\quad
u=\mathcal C(s_2-R_2),\quad v=\mathcal C(R_2+7s_2),\quad s_2=S:S.
$$

Here W is a contravariant vector when inserted into the curved
longitudinal operator. In particular,

$$
(T_1)_{ij}=W^k\partial_kQ_{ij}
+Q_{ik}\partial_jW^k+Q_{jk}\partial_iW^k
-\frac23Q_{ij}\partial_kW^k.                         \tag{1}
$$

The exact proposed assembly is

$$
\begin{aligned}
H_5={}&-\int\mathcal L_5-4\int b_2\cdot W
-2\int T_0:T_1+2\int\operatorname{tr}(QT_0^2)\\
&+\int R_3(v-u)+2\int(S:T_0)(7u+v)
-8\int(\nabla v)^TQ\nabla u .                       \tag{2}
\end{aligned}
$$

Every sign and factor in (2) follows from the existing reduced identities
as follows.

### York contribution

Since the odd raw norm vanishes, $\mathsf G_3=-\mathcal N_3$, and
the inverse-Gram momentum energy gives

$$
I_5=\mathcal N_3[S,S]
=-4\langle b_2,W\rangle+2\langle W,\mathcal A_1W\rangle.       \tag{3}
$$

Hold a test vector W fixed while differentiating the operator identity

$$
\langle W,\mathcal A W\rangle
=-\frac12\int\operatorname{tr}
\big[(e^{-Q}L_{e^Q}W)^2\big].                                \tag{4}
$$

Only after that variation insert the particular response
$W=\mathcal A_0^{-1}b_1$. With $Q\mapsto\lambda Q$,

$$
e^{-\lambda Q}L_{e^{\lambda Q}}W
=T_0+\lambda(T_1-QT_0)+O(\lambda^2).
$$

Unimodularity removes any volume variation. Taking the linear term in
(4) proves

$$
2\langle W,\mathcal A_1W\rangle
=-2\int T_0:T_1+2\int\operatorname{tr}(QT_0^2).                \tag{5}
$$

There is no W2 term in this identity. Differentiating the dependence of
W itself during this step would compute a different expression.

### Scalar contribution

The preceding reduction gives

$$
H_5=I_5-\int R_5-\langle R_3+7s_3,u\rangle
-\langle v,s_3-R_3\rangle-8\langle v,D_1u\rangle.             \tag{6}
$$

Use the pointwise identity $s_3=-2S:T_0$ and
$D_1u=-\partial_i(Q^{ij}\partial_j u)$. Integration by parts yields

$$
\langle v,D_1u\rangle=\int(\nabla v)^TQ\nabla u .             \tag{7}
$$

Combining (3), (5), (6), (7), and $\int R_5=\int\mathcal L_5$
gives (2). In particular, the minus sign in the final gradient pairing
is correct. R3 must remain the pointwise curvature; replacing it with
the integrated-density representative L3 would change (2).

## 2. What the assembly requires and useful independent checks

Only the first inverse fields W, u, v occur. H5 needs neither W2 nor
u3, nor the G2a correction, s4, or b3 required at sixth order.
The local ingredients are L5, pointwise R2/R3, b1/b2, Q and S.
L5 needs the connection through order four; b2 needs the tangent and
connection through order two.

Each nonlocal term has at most two first-Green factors and a
polynomial-Gaussian multiplier. In particular T1 in (1) contains Q or
a derivative of Q in every term. Thus the all-space integrands in
(2) have Gaussian suppression despite the individual Green fields'
algebraic tails. This statement is specific to this H5 representation;
it does not remove the previously identified mixed H6 algebraic tails.

The allowed sectors remain q5, q3p2 and qp4. Under $S\mapsto-S$,
W, T0, T1 and b2 all change sign, while u and v do not. Every term in
(2) is consequently even in canonical momentum. At q=0 all terms vanish.
At p=0, $W=0$, $u=-v$, so

$$
H_5(q,0)=-\int\mathcal L_5+2\int R_3v
+8\int(\nabla v)^TQ\nabla v,
$$

which agrees with the independent configuration Schur-complement formula.

Useful implementation checks are the scalar adjoint identity (7) and
the independently evaluated operator pairing (5). A local auxiliary
check on (1) is

$$
\operatorname{tr}T_1=\operatorname{tr}(QT_0),
$$

obtained by varying conformal tracelessness. It catches a lost metric
term or confusion between a fixed covector and a fixed contravariant
vector in L1.

## 3. Axial cubic Hessian: explicit diagonal form

Distinguish the spatial tensor Q from its STF electric/magnetic label
matrices E and B. Write the known cubic polynomial as

$$
H_3=\kappa\,[c_E\operatorname{tr}E^3+c_B\operatorname{tr}(EB^2)],
\quad \kappa=(2\pi/3)^{3/2},\quad
c_E=446/3,\quad c_B=25004/9.                                  \tag{8}
$$

Consider an external direction whose configuration label is
$E_{\rm ext}=\alpha A$, $B_{\rm ext}=0$, with
$A=\operatorname{diag}(-1/2,-1/2,1)$. Any external momentum component
does not contribute to this cubic derivative because H3 contains no p.
For arbitrary symmetric STF internal label directions U,V,

$$
\begin{aligned}
C_{EE}(U,V)&=6\kappa c_E\alpha\operatorname{tr}(AUV),\\
C_{BB}(U,V)&=2\kappa c_B\alpha\operatorname{tr}(AUV),\\
C_{EB}(U,V)&=0.                                               \tag{9}
\end{aligned}
$$

Transpose and cyclicity make $\operatorname{tr}(AUV)$ symmetric in
U,V. In the repository's ordered raw STF basis

$$
(F_0,F_1,F_2,F_3,F_4)=(D,A,X_{xy},X_{xz},X_{yz}),
\quad D=\operatorname{diag}(1,-1,0),
$$

where each X has two unit off-diagonal entries, direct matrix
multiplication gives

$$
\operatorname{tr}(AF_iF_j)=\delta_{ij}\tau_i,
\qquad \tau=(-1,3/4,-1,1/2,1/2).                              \tag{10}
$$

The norms are $n_i=\operatorname{tr}F_i^2=(2,3/2,2,2,2)$.
Let $f_{E,i}=\sqrt{2/g_{E,i}}$ and
$f_{B,i}=\sqrt{2/g_{B,i}}$ be the configuration factors from
[complete_quartic.py](../../verification/quantum_bridge/complete_quartic.py).
For normalized canonical internal legs the diagonal entries are

$$
C_{E,i}=6\kappa c_E\alpha f_{E,i}^2\tau_i,\qquad
C_{B,i}=2\kappa c_B\alpha f_{B,i}^2\tau_i.                     \tag{11}
$$

If the external leg is the unit normalized electric-axis coordinate,
then $\alpha=f_{E,1}$; do not omit this third normalization factor.
An external scalar multiple multiplies alpha by the same scalar.

The equal-norm pairs (D,Xxy) and (Xxz,Xyz) have identical entries after
canonical normalization. Their normalized angular factors are

$$
\tau_i/n_i=(-1/2,1/2,-1/2,1/4,1/4).                           \tag{12}
$$

As an independent consequence, each parity block has zero trace:
$2(-1/2)+(1/2)+2(1/4)=0$. Since the free frequency is constant within
each parity multiplet, its cubic one-point loop contraction vanishes
separately in each parity block.

### Symmetry explanation and limits

Axial rotations split each real spin-two multiplet into m=0, |m|=1,
and |m|=2 subspaces. A symmetric invariant bilinear form is a scalar on
each two-dimensional angular subspace, explaining the diagonal Hessian
and pair degeneracies. Spatial parity forbids electric/magnetic mixing
for this electric external direction. The chosen axis must be aligned
with the displayed basis; a rotated axis requires the corresponding
known spin-two basis rotation.

**Axiality alone does not imply the same-parity statement.** For a pure
magnetic axial external direction $B_{\rm ext}=\beta A$, the nonzero
cubic block instead is

$$
C_{EB}(U,V)=2\kappa c_B\beta\operatorname{tr}(AUV),              \tag{13}
$$

with the appropriate electric and magnetic canonical factors. It
connects opposite parities. Additional radial copies can also mix within
the same angular representation. The reduction here uses the specified
one-profile-per-parity model.

## 4. Analytical two-dimensional loop blocks replace eigendecomposition

For the electric axial external direction, (11) makes the cubic C
diagonal in the ten configuration modes. The integrated matrix
$K_r=\int G(\nu)C^{(r)}G(\nu+\Omega_r)d\nu/(2\pi)$ consequently
splits into ten independent two-dimensional phase-space blocks:

$$
(K_r)_i=\frac{C^{(r)}_i}{\omega_i(4\omega_i^2+\Omega_r^2)}
\begin{pmatrix}1&\Omega_r/2\\\Omega_r/2&-\omega_i^2\end{pmatrix}.
                                                                    \tag{14}
$$

The off-diagonal entries generally survive at nonzero frequency. To
contract (14) against a symmetric H5 Hessian without numerical
eigendecomposition, set

$$
\begin{aligned}
\theta_i&=\tfrac12\arctan\!\frac{\Omega_r}{1+\omega_i^2},\\
v_{i,+}&=\cos\theta_i\,e_{Q_i}+\sin\theta_i\,e_{P_i},\\
v_{i,-}&=-\sin\theta_i\,e_{Q_i}+\cos\theta_i\,e_{P_i},\\
\mu_{i,\pm}&=\frac{1-\omega_i^2
\pm\sqrt{(1+\omega_i^2)^2+\Omega_r^2}}2.
\end{aligned}                                                     \tag{15}
$$

Then (14) is its scalar prefactor times
$\mu_{i,+}v_{i,+}v_{i,+}^T+\mu_{i,-}v_{i,-}v_{i,-}^T$.
All variables here use the fixed normalized canonical units of the
loop prescription. No positivity of K is assumed; the two eigenvalues
of the displayed matrix have opposite signs.

For a general external triple this still needs at most twenty H5
directional contractions per r. Its benefit is a known sparse basis,
not a reduction in that worst-case count. At zero frequency the two
directions are simply Q and P.

## 5. Additional reduction when every external leg shares the electric axis

Suppose all four external directions lie in
$\operatorname{span}\{e_{Q_{E,1}},e_{P_{E,1}}\}$. This includes the
frequency-dependent free phase directions induced by four sources for
the electric-axis configuration mode. Every direction is fixed by the
full stabilizer of A: axial rotations and a half-turn about an
equatorial axis. The latter makes the stabilizer an embedded O(2),
not merely SO(2). Spatial parity also fixes every external direction.

Therefore the H5 Hessian after the other three external differentiations,
and the H6 Hessian after all four, have identical phase blocks on the two
members of each |m| pair and do not mix parity sectors. The equatorial
half-turn rules out an angular antisymmetric commutant in mixed Q/P
blocks; SO(2) alone would not suffice for that stronger statement.

One can retain only the representative labels D, A, Xxz in each parity,
with multiplicities respectively 2,1,2. For the fifth/cubic contraction
this gives six representative modes times two directions (15), hence
at most **twelve H5 contractions per external singleton**, or 48 over
the four singletons. The sixth-order tadpole similarly needs twelve
representative QQ/PP contractions. Thus this common-axis choice can
reduce the nominal 100 higher-vertex contractions to at most 60 while
retaining all ten physical internal modes through exact multiplicities.

The prefactor C in (14) must still be evaluated for its particular
representative and external leg; the m=0, |m|=1 and |m|=2 coefficients
are not equal to one another. If any other external leg breaks the
common electric-axis stabilizer, this representative replacement is
not justified, even when one individual cubic insertion remains diagonal.

For three purely configurational external H5 legs at zero momentum,
momentum-reversal parity additionally forces the internal QP Hessian
to vanish. This simplification does not apply to generic physical
configuration-source phase directions at nonzero frequencies, which
contain external P components. Keeping (14)--(15) avoids that omission.

These are exact algebraic reductions for the declared finite basis.
They are not evaluated H5/H6 coefficients, a completed loop amplitude,
or a continuum/scattering-limit argument.
