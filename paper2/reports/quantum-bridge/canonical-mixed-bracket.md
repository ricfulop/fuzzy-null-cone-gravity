# Mixed constraint bracket of a Hamiltonian-first Gaunt candidate

## Question and status

The [kinetic-completion experiment](kinetic-completion-result.md) repaired
one lapse derivative by restoring discarded kinetic channels. It did not
show that a finite canonical restriction preserves the gravitational
constraint algebra. This report defines that candidate separately and
derives a mixed-bracket test before a tracked independent verification.

**Status at this source snapshot:** analytic derivation and verification
code prepared; the new experiment has not yet been interpreted. Numerical
values below are exact analytical predictions, not reported run results.

The prediction is that a radial spatial generator and a lapse generator
have a nonzero canonical bracket on vacuum constraint data when their
harmonic number reaches the angular cutoff. If confirmed, this would rule
out simply importing the continuum first-class algebra for this candidate.
It would not complete the Dirac algorithm or rule out other finite-rank
theories.

## 1. Declare the changed theory

Use the continuum ADM constraints, with $16\pi G=1$ and $\Lambda=0$,

$$
 C[\alpha]=\int d^3x\,\alpha\left[
 \frac{\pi^{ij}\pi_{ij}-\tfrac12\pi^2}{\sqrt\gamma}
 -\sqrt\gamma\,{}^{(3)}R\right],
 \qquad
 D[\beta]=\int d^3x\,\pi^{ij}{\cal L}_{\beta}\gamma_{ij}.
$$

The second expression equals the usual momentum constraint after an
integration by parts with the stated boundary conditions. The continuum
form and Hamilton equations are the starting point, not an assumption
about their projected algebra. See
[Arnowitt, Deser and Misner, PDF pp. 7--8, Eqs. (3.14)--(3.15)](https://www.alphaxiv.org/pdf/gr-qc/0405109v1?page=7).

Fix a flat background metric $\bar\gamma$ on a spherical annulus. In a
fixed background frame, retain scalar, vector and trace-free tensor
harmonics through $j=L=N-1$ in both canonical tensor variables:

$$
 \gamma=\bar\gamma+q_L,\qquad
 \pi/\sqrt{\bar\gamma}=p_L,\qquad P_Lq_L=q_L,\quad P_Lp_L=p_L.
$$

The $rr$ component and angular trace have spin zero; the radial-angular
components have spin one; the angular trace-free components have spin two.
The projector is orthogonal in the background tensor $L^2$ pairing.
Radial dependence remains continuous.

The candidate action is the pullback of
$\int(\pi^{ij}\dot\gamma_{ij}-C[\alpha]-D[\beta])\,dt$ to these canonical
variables, with independent retained lapse and shift smearings. Products,
inverse metrics, densities and constraints are evaluated pointwise and
integrated with complete Gaunt contractions. No binary intermediate
projection is inserted into the nonlinear constraint densities.

This is a Hamiltonian-first revision. It is not the Legendre transform
of Paper I's metric-coefficient Galerkin Lagrangian, and it is not silently
adopted as the common theory. In particular, independently retaining
canonical momenta and lapse/shift coefficients changes the regulator.

Use the normalized angular integral
$\langle v\rangle=(4\pi)^{-1}\int d\Omega\,v$ in both the symplectic
form and the generators. All coefficients below use this convention.
The canonical bracket is
$\{F,G\}=\langle F_q,G_p\rangle-\langle F_p,G_q\rangle$,
including the radial integral $r^2dr$. A positive common normalization
does not change whether the bracket vanishes.

## 2. Initial data on every retained constraint

Take

$$
 \gamma=\bar\gamma,\qquad
 K_{ij}=k\,\delta_{iz}\delta_{jz},\qquad
 p^{ij}=K^{ij}-\delta^{ij}K
       =k\,\mathrm{diag}(-1,-1,0)
$$

in Cartesian components. Then
$K_{ij}K^{ij}-K^2=0$, ${}^{(3)}R=0$, and
$\partial_jp^{ij}=0$ pointwise. Thus every retained Hamiltonian and
momentum constraint vanishes, not just the two particular smearings used
below. Set the primary momenta conjugate to lapse and shift to zero as
well. The data are in the retained space for every $L\geq2$.

At these data,

$$
 C_p=2\alpha K,\qquad
 C_q^{\mathrm{kin}}=0,\qquad
 C_q^{\mathrm{pot}}=-\operatorname{Hess}\alpha
                    +\bar\gamma\,\Delta\alpha .
$$

The kinetic metric gradient vanishes because
$2p^2-(\operatorname{tr}p)p=0$.
For a fixed spatial smearing,

$$
 D_p={\cal L}_\beta\bar\gamma,\qquad D_q=-{\cal L}_\beta p
$$

where $p$ is understood as the density relative to the background volume.
In a Cartesian frame at this constant momentum,
$D_q=Mp+pM^T-(\operatorname{tr}M)p$, with $M^i{}_j=\partial_j\beta^i$.

Set $x=\cos\theta$, $u=P_\ell(x)$, and choose

$$
 \alpha=f(r)u,\qquad \beta=f(r)u\,e_r,\qquad
 f(r)=(r-1)^2(2-r)^2,\quad 1\leq r\leq2 .
$$

These are test functions for the generators, not a positive background
lapse. Both $f$ and $f'$ vanish at the radial boundaries. All boundary
terms used here therefore vanish. The angular fields are smooth at the
poles. The normalized radial integral is

$$
 R_f=\int_1^2 r f^2\,dr=\frac1{420},\qquad
 \int_1^2r^2 f f'\,dr=-R_f.
$$

More generally the same calculation works for any nonzero smooth
interior-supported $f$, replacing $1/420$ by its positive $R_f$.

## 3. Exact angular calculation

Represent an axisymmetric symmetric tensor by four polynomial components
$(a,b,c,d)$:

$$
 T_{rr}=a,\quad T_{\theta\theta}+T_{\phi\phi}=b,\quad
 T_{r\theta}=\sqrt{1-x^2}\,c,\quad
 (T_{\theta\theta}-T_{\phi\phi})/2=(1-x^2)d.
$$

The other two components vanish on the specified data. They remain in the
candidate theory; axial symmetry makes their first functional derivatives
vanish in this test. The positive tensor pairing is

$$
 \langle T:S\rangle=
 \langle a_Ta_S+\tfrac12b_Tb_S+
 2(1-x^2)c_Tc_S+2(1-x^2)^2d_Td_S\rangle.
$$

This is the Frobenius pairing, not the DeWitt kinetic pairing.
Write $W=1-x^2$ and set $k=1$ within the angular formulas. Then

$$
 K=(x^2,W,-x,\tfrac12),\quad T=uK,
$$
$$
 D_r=(-Wu,(1+x^2)u,0,-u/2),\qquad
 D_t=(2W(u+xu'),0,xu+x^2u',0).
$$

The spatial generator's metric gradient is
$D_q=k[f'D_r+(f/r)D_t]$.
Both $D_r:T$ and $D_t:T$ vanish pointwise before projection.

The potential gradient and $D_p$ have only harmonic number $\ell$.
Their projected contraction therefore equals their continuum contraction
for $\ell\leq L$; its radial integral vanishes. Consequently the complete
projected bracket at this point is

$$
 \{D_L,C_L\}
 =-2k^2\int_1^2r^2 f\left[
 f'\langle Q_LD_r:Q_LT\rangle+
 \frac f r\langle Q_LD_t:Q_LT\rangle\right]dr.
$$

Define $C_r=\langle Q_LD_r:Q_LT\rangle$ and
$C_t=\langle Q_LD_t:Q_LT\rangle$. Radial integration gives

$$
 \frac{\{D_L,C_L\}}{k^2R_f}=2(C_r-C_t).
$$

For $\ell\leq L-2$, all these tensors are retained, so this bracket is
zero. For $\ell=L-1,L$, only their highest common harmonic
$j=\ell+2$ contributes to the discarded pairing. Put

$$
 c_\ell=\frac{(\ell+1)(\ell+2)}{(2\ell+1)(2\ell+3)}.
$$

The highest component of $T$ is
$$
 \left(c_\ell P_j,-c_\ell P_j,
 -\frac{c_\ell}{\ell+2}P_j',
 \frac{c_\ell}{2(\ell+1)(\ell+2)}P_j''\right).
$$

Using
$\langle(1-x^2)^s(P_j^{(s)})^2\rangle
 =(j+s)!/[(2j+1)(j-s)!]$ gives

$$
 C_r=-\frac{c_\ell^2}{(\ell+1)(\ell+2)},\qquad
 C_t=-\frac{2c_\ell^2(\ell+1)}{\ell+2}.
$$

The predicted mixed bracket is therefore

$$
 \boxed{\frac{\{D_L,C_L\}}{k^2R_f}
 =\frac{2(\ell+1)(\ell+2)\,[2(\ell+1)^2-1]}
 {(2\ell+1)^2(2\ell+3)^2}>0,\qquad \ell=L-1,L .}
$$

For $L=2$, the normalized coefficients are $28/75$ at $\ell=1$
and $408/1225$ at $\ell=2$. With the specified radial profile, the
brackets at $k=1$ are $1/1125$ and $34/42875$, respectively.
These are predictions to be checked by the independent flow calculation.

For this sequence of unnormalized Legendre smearings the coefficient tends
to $1/4$ as $\ell$ grows with the cutoff. Unit-$L^2$ normalization of
both smearings multiplies it by $2\ell+1$. Their spatial derivatives also
grow; neither sequence is a fixed-physical-momentum scattering limit.
For any fixed $\ell$, this witness vanishes once $L\geq\ell+2$.

## 4. Independent full Hamiltonian-flow check

The verification does not infer the full bracket from a kinetic norm.
It separately computes the canonical flow
$\dot\gamma=P_L C_p=2fP_L(uK)$ and
$\dot p=-P_L C_q$, and differentiates the momentum generator:

$$
 \dot D[\beta]=\int r^2dr\,
 \langle p:{\cal L}_\beta\dot\gamma+
 \dot p:{\cal L}_\beta\bar\gamma\rangle .
$$

In the background orthonormal spherical frame,

$$
 M=\begin{pmatrix}
 f'u & -(f/r)\sqrt W\,u' & 0\\
 0 & (f/r)u & 0\\
 0 & 0 & (f/r)u
 \end{pmatrix},\qquad
 {\cal L}_\beta\dot\gamma
 =fu\,\partial_r\dot\gamma+M^T\dot\gamma+\dot\gamma M.
$$

The script evaluates these three-by-three matrices directly at quadrature
points. It also evaluates the projected-gradient definition of the
Poisson bracket and the unprojected continuum control.

The curvature-potential term is included. In particular,
$$
 \int d\Omega/(4\pi)\,
 ({\cal L}_\beta\bar\gamma):C_q^{\mathrm{pot}}
 =\frac4{2\ell+1}
 \left[\frac{f'^2+ff''}{r}
 +(1-\ell(\ell+1))\frac{ff'}{r^2}\right].
$$
Its integral with $r^2dr$ vanishes by the boundary conditions, although
the local integrand is not zero. The numerical check records both facts
and cannot pass by simply omitting that term.

The planned run uses exact rational projections for all 88 pairs
$L=2,\ldots,12$, $\ell=0,\ldots,L$, plus ten independent matrix-flow
cases at 60-digit precision. The two quadratures use 24/32 angular nodes
and 6/8 radial nodes. An additional control enlarges the *canonical*
space to $J\geq\ell+2$, which removes this particular discarded tensor.
That control changes the phase space; it is not just better quadrature
or auxiliary padding of the same finite theory.

## 5. Scope of a possible confirmed obstruction

A nonzero bracket at a point where every proposed primary and secondary
constraint vanishes prevents this set from being first class there.
It does not establish whether a later Dirac algorithm eliminates that
point, produces further constraints, fixes some multipliers, or leads to
a usable second-class theory. Those possibilities require their own
analysis and measure.

This test concerns a declared alternative prescription. It neither
duplicates a full Dirac analysis of Paper I's original Lagrangian nor
equates the two regulators. A repair that changes only the kinetic
matrix or makes the lapse linear has not thereby repaired the constraints,
ghosts and measure.

The full one-loop amplitude and its controlled limit remain the target.
The test determines whether continuum gauge reduction may be used for
this candidate without additional terms. The dynamical corner-state
count and Euclidean--Lorentzian reconstruction still require the physical
state space identified in the [separate analysis](corners-and-reconstruction.md).
