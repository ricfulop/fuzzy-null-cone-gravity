# A lapse-constraint test on vacuum initial data

**Status:** analytical derivation and predictions; an independent tracked
calculation is prepared. Numerical validation is not claimed by this source
snapshot. This tests a prerequisite for the requested loop amplitude, not the
amplitude itself. All three quantum-gravity objectives remain open.

## Why this test is needed

The common-action proposal adopts Paper I's full metric-coefficient Galerkin
pullback of the Einstein action. A complete Gaunt evaluation removes premature
product truncation, but the independent velocity space is still angularly
truncated. Legendre transformation in that space is not the same as simply
restricting the continuum Hamiltonian. Paper I established this distinction
for a positive kinetic block off the full constraint surface.

Here the DeWitt trace term, scalar--tensor cancellations, the actual metric
coefficient map, and initial data satisfying *both* continuum constraints are
included. The calculation uses an ADM representative of the Einstein action
with its usual temporal total derivative removed. All lapse test functions
can be supported away from radial boundaries; a gravitational corner
variational principle is not assumed or derived by this bulk calculation.

The continuum starting point is the canonical Einstein action. ADM state:
“Their variation yields the four constraint equations” in their discussion
of lapse and shift multipliers. The first-order action and reduction are on
[pp. 8 and 10](https://www.alphaxiv.org/pdf/gr-qc/0405109v1?page=10).
Our finite-rank result below is a new inference/calculation for this specified
regulator; the cited continuum formulation does not establish it.

## 1. Constraint-satisfying data inside the metric ansatz

Set the cosmological constant to zero. Work locally on a spatial annulus in
flat Cartesian space with

$$
\gamma_{ij}=\delta_{ij},\qquad K_{ij}=k\,\delta_{iz}\delta_{jz},
$$

where k is constant and nonzero. Both continuum vacuum constraints vanish:

$$
{}^{(3)}R+K^2-K_{ij}K^{ij}=0,\qquad
D_j(K^{ij}-\gamma^{ij}K)=0.
$$

The first equality is k^2-k^2=0 and the second is the divergence of the
constant Cartesian tensor k diag(-1,-1,0). No off-shell scalar proxy is used.
These are local data, not an assertion of asymptotic flatness with constant
extrinsic curvature at infinity.

Write x=cos(theta). In an orthonormal spherical frame the nonzero components
of K/k are

$$
K_{rr}/k=x^2,\quad K_{r\theta}/k=-x\sqrt{1-x^2},\quad
K_{\theta\theta}/k=1-x^2,\quad K_{\phi\phi}=0.
$$

Their scalar components have degrees 0 and 2; their spin-one and spin-two
components have degree 2. Thus the complete velocity and momentum data fit
the retained bands for L=N-1 >=2. At this point the common 2+2 metric map has
h_rr=1, X=r^2, C=0 and A_r=0. Its tangent map onto these six spatial metric
components is invertible within the respective retained bands. It has no
lapse dependence, so the fixed canonical momenta used below are legitimate
momenta of that coefficient map at this point.

Because both the lapse test and the momenta are axisymmetric, multiplication
by the lapse preserves m. The m!=0 kinetic blocks have zero momenta and
decouple from this Hamiltonian evaluation. The unused azimuthal tensor
components also decouple at this spherical spatial background. They have
not been discarded from the definition of the gravitational action.

Vary a *retained metric coefficient*, not an independently bandlimited lapse:

$$
h_{tt}=-1-2\epsilon f(r)P_\ell(x),\qquad h_{tr}=0,
\qquad \alpha=\sqrt{1+2\epsilon f(r)P_\ell(x)}.
$$

Here 0<=ell<=L, f is smooth and compactly supported inside the radial
annulus, and |2 epsilon f|<1 keeps the metric Lorentzian. All other metric
coefficients are held fixed. The derived lapse may have infinitely many
angular harmonics, exactly as required by the common action.

## 2. Perform the finite-dimensional Legendre transform

Suppress the positive common density r^2/(4G), and set f=k=1 for the local
angular calculation. Let P act with the corresponding scalar, vector and
tracefree tensor harmonic bands. At this spherical spatial metric P commutes
with the algebraic DeWitt operator D, whose bilinear form is

$$
\langle v,Dw\rangle
=\left\langle v_{ij}w^{ij}-\operatorname{tr}v\operatorname{tr}w\right\rangle,
\qquad\langle a\rangle=\frac12\int_{-1}^1a(x)\,dx.
$$

With zero shift the retained kinetic action and momenta are

$$
\mathcal L_{N,\mathrm{kin}}=\frac14\langle v,A_\epsilon v\rangle,
\qquad A_\epsilon=P\alpha^{-1}DP,
\qquad p=\frac12 A_\epsilon v.
$$

D is invertible on spatial symmetric tensors, including its negative trace
direction. In a neighborhood of epsilon=0, A_epsilon is invertible on the
retained velocity space. At epsilon=0, v=2K and p=DK. Holding these p fixed,

$$
\mathcal H_N(\epsilon)=\langle p,A_\epsilon^{-1}p\rangle.
$$

The curvature potential vanishes at the flat spatial metric. The shift
equations are the projected continuum momentum constraints at epsilon=0
and vanish on the chosen data. The lapse equations also vanish there because
the continuum Hamiltonian constraint vanishes pointwise. Primary momenta of
lapse/shift coefficients vanish. This is therefore a point on the
**primary--secondary constraint surface** of the actual pullback.

For u=P_ell, alpha^-1=1-epsilon u+(3/2)epsilon^2 u^2+..., and inversion of the
retained kinetic matrix gives

$$
\mathcal H_N(\epsilon)
=\langle K,DK\rangle+\epsilon\langle K,DuK\rangle
+\epsilon^2\left[
\langle PuK,DPuK\rangle-\frac32\langle uK,DuK\rangle
\right]+O(\epsilon^3).
$$

The first two terms and the pointwise full quadratic contraction vanish for
these rank-one data. With Q=1-P,

$$
\boxed{\mathcal H_N''(0)
=-2\langle Q(uK),D Q(uK)\rangle.}
$$

Since tr(K)=1 and u is retained, tr(Q(uK))=Q(u)=0. Consequently the right
side is minus twice a nonnegative ordinary tensor norm. This positivity
follows here from the *full constraint-satisfying data*, rather than from
discarding the negative trace direction of gravity.

## 3. An explicit N=3 prediction

At L=2, u=P_1=x, the discarded components are

$$
Q(uK_{rr})=\frac25P_3,\qquad
Q[u(K_{\theta\theta}+K_{\phi\phi})]=-\frac25P_3,
$$

$$
Q(uK_{r\theta})=-\frac{2}{15}\sqrt{1-x^2}\,P_3',\qquad
Q\!\left[u\frac{K_{\theta\theta}-K_{\phi\phi}}2\right]
=\frac1{30}(1-x^2)P_3''.
$$

Their squared tensor norm is 2/15. Thus

$$
\mathcal H_{N=3}''(0)=-\frac4{15}k^2 f(r)^2,
$$

before restoring the positive radial factor r^2/(4G). For u=P_2 the
corresponding predicted derivative is -24/175. A linear change from epsilon
to the h_tt harmonic coefficient only rescales these nonzero values.

Equivalently, a Poisson bracket between a lapse primary momentum and its
secondary equation is nonzero at this point. The continuum primary and
secondary first-class pattern is **not automatically inherited**. This
does not yet classify all final constraints: stabilization could determine
multipliers or impose additional constraints, and the full degree count and
positivity remain to be calculated.

## 4. Independent checks and cutoff dependence

The calculation in `paper2/verification/quantum_bridge/run.py` will:

1. Verify the vacuum constraint polynomials and the normalized spin-harmonic
   orthogonality with exact rational arithmetic.
2. Compute the entire axisymmetric lapse Hessian for L=2,...,12. The analytic
   prediction is a diagonal block with two negative entries at ell=L-1,L;
   lower fixed harmonics have zero second derivative on these data.
3. Independently integrate the full, finite-amplitude factor
   (1+2 epsilon P_ell)^(-1/2), invert its finite matrices, and differentiate
   the resulting Hamiltonian. Use 60-digit decimal arithmetic, independent
   96- and 128-node Gauss--Legendre quadratures, three decreasing difference
   steps, and Richardson extrapolation.

The numerical Hamiltonian is assembled from all nonzero tensor blocks. Write
a=K_rr, b=K_thetatheta+K_phiphi, c=K_rtheta, d=(K_thetatheta-K_phiphi)/2.
If M_s is the Gram matrix weighted by alpha^-1 in the spin-s basis and
a_*,b_*,c_*,d_* are the unweighted basis moments of these components, then

$$
\mathcal H_N
=-2a_*^TM_0^{-1}b_*-\tfrac12b_*^TM_0^{-1}b_*
+2c_*^TM_1^{-1}c_*+2d_*^TM_2^{-1}d_*.
$$

This expression is evaluated without substituting the leakage formula.
The exact result, finite-difference convergence, and quadrature convergence
must agree before the computational result is accepted.

The fixed ell=1 test loses its second-order defect once L>=3. A test whose
lapse harmonic moves to ell=L is a different sequence. Unit-L2-normalized
results will be reported separately from unnormalized Legendre amplitudes.
Neither sequence by itself proves or disproves a physical large-N scattering
limit. Internal modes near the cutoff are relevant to a loop calculation,
so agreement of a few fixed low harmonics is insufficient to import the
continuum gauge measure.

## 5. Consequence for the research program

A confirmed nonzero bracket requires completing the finite-rank Dirac
algorithm before interpreting its degrees of freedom and measure. An
alternative is a declared change of construction: reduce or integrate out
the appropriate constrained variables before imposing the angular and state
cutoffs, and derive the induced interactions and determinants. Simply
inserting continuum ghosts or dropping the additional lapse dependence does
not establish equality to the common metric-coefficient action.

The corner and Euclidean constructions must use the same resulting physical
theory. A separately chosen finite Hermitian Hamiltonian provides a useful
conditional transfer-operator example but cannot fill either missing link.
This calculation therefore identifies a specific decision needed before the
three main deliverables can be claimed, rather than a universal no-go for
discrete gravity or an assertion that the original goal is now complete.
