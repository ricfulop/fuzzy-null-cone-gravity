# Bounded independent reference for the mixed sixth-order coefficient

This is an analytic implementation design, **not numerically tested**. It
authorizes no calculation or additional CPU job. The [production pilot](full-sixth-order-pilot.md) evaluates
the [contracted H6 formula](sixth-order-bilinear-and-tail-reduction.md); its
independent mixed-functional validation remains pending. The recommendation
below checks one mixed coefficient through explicit spatial dependent fields,
using the [unreduced recursion](reduced-sixth-order-and-loop-kernel.md),
sections 1–4. No source code or frozen calculation is changed by this note.
A separate read-only algebra audit checked the radial Green formulas, Kelvin
normalization, exterior-source effects and finite-ball flux identities, with
no remaining concrete mathematical objection. Numerical validation is pending.

**Target and independence.** Keep the pilot's raw tensors
$Q=E(A)$, $S=B(XZ)$, with $A=\operatorname{diag}(-1/2,-1/2,1)$ and
$XZ_{13}=XZ_{31}=1$ (one-based spatial indices). Both include their physical
$e^{-r^2/2}$ envelope. Retain all ten raw basis vectors in the canonical
correction. These are raw packet amplitudes, not normalized oscillator $Q_i,P_i$.
All reported energies include the full spatial measure and all powers of $\pi$.

Use a real spherical-harmonic radial Green implementation, with all $m$, to
recompute $u_2,W_1,W_2,u_3,u_4$. The inverse path must not call
`GaussianHeatSource`, the Fourier bilinear compiler, `PairingAccumulator`,
`collapse_triple_heat`, or the production triple integrator. It must not consume
their compiled $F,K,r,t$ sources or their tail coefficients. The same raw modes
and separately checked local geometry may be shared. This gives an independent
evaluation of the unreduced functional, conditional on those anchored local
ingredients; it is not a finite-amplitude solution of the full nonlinear theory.

**Minimal code.** Two new modules suffice in a later tracked implementation:

- `spherical_green_reference.py`: full-$m$ harmonic projection/reconstruction,
  regular radial scalar solves, spatial jets, a York factorization, and explicit
  exterior Laurent sources. Extend the existing second-order radial strategy in
  [reduced_gaunt_completion.py](../../verification/quantum_bridge/reduced_gaunt_completion.py)
  in the new module; preserve the archived axisymmetric implementation.
- `sixth_unreduced_reference.py`: build the following physical sources, perform
  spatial quadrature, and report separate unreduced energy and boundary ledgers.
  The tracked driver owns levels, caps, tolerances, checkpoints and authorization.

The existing axisymmetric Legendre routine cannot be used unchanged: $B(XZ)$
breaks axial symmetry. Do not remove magnetic-number copies on symmetry grounds
before checking their projections. No dense H5/H6 tensor or finite-amplitude
nonlinear York solve is required for this one reference direction.

**Unreduced evaluation sequence.** Write $\mathcal P=(-\Delta)^{-1}$,
$\mathcal C=\mathcal P/8$, $\mathcal Y=\Delta_L^{-1}$ and
$D_nf=\partial_i[(g^{-1})_n^{ij}\partial_jf]$, where
$(g^{-1})_1=-Q$, $(g^{-1})_2=Q^2/2$. Evaluate these differential operators
from physical metric jets, including derivatives of the inverse coefficients.

1. Solve $u=u_2=\mathcal C(s_2-R_2)$ and $W_1=\mathcal Yb_1$ independently.
   Form $T=L_0W_1$, $T_1=L_1W_1$, $T_2=L_2W_1$ and
   $\mathcal A_1W_1,\mathcal A_2W_1$ with the full metric/inverse/connection
   jet formulas. Here $T_2$ does **not** mean $L_0W_2$.
2. Form $r=b_2-\mathcal A_1W_1$ directly, then solve $W_2=\mathcal Yr$.
   Recompute the ten components
   $$g_\alpha=-\frac1{12}\int\operatorname{tr}([Q,e_\alpha][Q,S])
                +2\int b_1(Q,e_\alpha)\cdot W_1.$$
   These need ten spatial pairings, not ten further York inverses. Use the
   physical Gram matrix to form $E_c=\sum e_\alpha(\mathsf G_0^{-1}g)^\alpha$.
   Compare $g$ componentwise with the archived-quartic targets; do not import
   the production $g$ into the reference or add a separate $c$ term below.
3. Use raw tangents and matrix products rather than the simplified $s_4$ formula:
   $$E_1=(QS+SQ)/2,\qquad E_2=(Q^2S+QSQ+SQ^2)/6,$$
   $$Z_1=E_1-T,\qquad Z_2=E_2-L_0W_2-T_1-E_c,$$
   $$U_1=Z_1-QS,\qquad U_2=Z_2-QZ_1+Q^2S/2,$$
   $$s_3=2\operatorname{tr}(SU_1),\qquad
     s_4=\operatorname{tr}(U_1^2)+2\operatorname{tr}(SU_2).$$
   $\operatorname{tr}(U_1^2)$ is not the Frobenius norm of $U_1$.
4. Set $B_j=R_j+7s_j$ and solve
   $$F_3=s_3-R_3+8D_1u,\qquad u_3=\mathcal CF_3,$$
   $$F_4=s_4-R_4+8D_1u_3+8D_2u-B_2u,\qquad u_4=\mathcal CF_4.$$
   $F_3,r$ are localized, but $F_4$ has a bare algebraic $T:T$ tail.
5. Integrate the original strong-York and scalar expression:
   $$I_6=-\mathcal M_4-2\int r\cdot W_2-4\int b_3\cdot W_1
             +2\int W_1\cdot\mathcal A_2W_1
             +g^{\mathsf T}\mathsf G_0^{-1}g,$$
   $$H_6^{\rm ref}=I_6-\int\mathcal L_6
             -\int B_4u-\int B_3u_3-\int B_2u_4
             +28\int s_2u^2,$$
   where $\mathcal M_4=\int\operatorname{tr}([Q,[Q,S]]^2)/360$.
   No adjoint elimination of $u_4$ or $W_2$ is used. The $q^6,q^4p^2,q^2p^4,p^6$
   sectors and the canonical correction remain present. Pointwise $R_2,R_3,R_4$
   cannot be replaced by integrated Gamma–Gamma densities.

**Scalar radial Green data.** Use real harmonics with
$\int Y_{\ell m}Y_{\ell' m'}d\Omega=\delta_{\ell\ell'}\delta_{mm'}$.
For one source coefficient $f_{\ell m}$,
$$
(\mathcal Cf)_{\ell m}(r)=\frac1{8(2\ell+1)}
 \left[r^{-\ell-1}\int_0^r s^{\ell+2}f_{\ell m}(s)ds
       +r^\ell\int_r^\infty s^{1-\ell}f_{\ell m}(s)ds\right]. \tag{1}
$$
Regularity at the origin is $r^\ell$ times a smooth even function. Preserve that
form in origin interpolation and derivative limits; do not divide noisy angular
projections by $r^\ell$. The finite-radius boundary datum is
$$
R u'_{\ell m}(R)+(\ell+1)u_{\ell m}(R)
 =\frac{R^\ell}{8}\int_R^\infty s^{1-\ell}f_{\ell m}(s)ds. \tag{2}
$$
A homogeneous Robin boundary is justified only after bounding this exterior
source datum. It is wrong for $F_4\sim r^{-4}$.

Retain every nonzero exterior power, including subleading multipoles. If
$f_{\ell m}(s)=\sum_p d_{\ell m,p}s^{-p}+\epsilon_{\ell m}(s)$ for $s\ge R$,
with $p+\ell>2$, the source outside $R$ contributes inside the ball
$$
\delta u_{\ell m}(r)=\sum_p
 \frac{d_{\ell m,p}r^\ell}
 {8(2\ell+1)(p+\ell-2)R^{p+\ell-2}},\quad r<R. \tag{3}
$$
The omitted remainder is bounded by (1) with $|\epsilon_{\ell m}|$.
Its effect on $-\int B_2u_4$ must be weighted by $|B_2|$ and included in the
physical energy error budget; a small far-field residual alone is insufficient.

For $r\ge R$, put $A_{\ell m}(R)=\int_0^R s^{\ell+2}f_{\ell m}(s)ds$.
The exact algebraic-source continuation is
$$
u_{\ell m}(r)=\frac1{8(2\ell+1)}\left\{
r^{-\ell-1}\left[A_{\ell m}(R)+\sum_p d_{\ell m,p}J_{\ell p}(r,R)\right]
 +\sum_p\frac{d_{\ell m,p}r^{2-p}}{p+\ell-2}\right\}, \tag{4}
$$
$$
J_{\ell p}(r,R)=
\begin{cases}
(r^{\ell+3-p}-R^{\ell+3-p})/(\ell+3-p),&p\ne\ell+3,\\
\log(r/R),&p=\ell+3.
\end{cases}
$$
Add the bounded Gaussian-remainder continuation. Do not replace $A(R)$ by a
convergent all-space moment without checking convergence: $\ell=2,p=4$
already makes that moment diverge. Nonresonant particular terms have coefficient
$d/[8(\ell(\ell+1)-(p-2)(p-3))]$ multiplying $r^{2-p}$.
Retain logarithms when $p=\ell+3$; absence of the leading resonance is not a
license to omit possible subleading resonances.

**York solve and its boundary conditions.** A componentwise scalar construction is
$$
V=\mathcal P f,\qquad h=\nabla\cdot V,\qquad
\chi=\mathcal P h,\qquad \mathcal Yf=-V-\tfrac14\nabla\chi. \tag{5}
$$
Use (1)–(4) without the factor $1/8$ for $\mathcal P$. This has Fourier symbol
$-\delta_{ij}/k^2+k_ik_j/(4k^4)$ and coordinate kernel
$$K_{ij}(x)=-\frac{7\delta_{ij}+n_in_j}{32\pi r}. \tag{6}$$
Only the York vector must decay. For source moment $M=\int f$,
$\chi\sim-(M\cdot n)/(8\pi)$ is a bounded, nondecaying dipole. Imposing
$\chi\to0$, or harmonic exterior data for its algebraic source $h$, changes
the York inverse. Check that the nonintegrable **algebraic $r^{-2},\ell=0$
coefficient** of $h$ vanishes; do not silently remove a numerical violation.
This does not set the full $h_{00}(r)$ or its integral to zero. For example,
$f=\nabla g$ gives $V=\nabla\mathcal Pg$, $h=-g$, $\chi=-\mathcal Pg$ and
$\mathcal Yf=-3\nabla\mathcal Pg/4$. A gradient-Gaussian control must preserve
this localized monopole. The additive constant in $\chi$ does not affect $W$.

Obtain first radial derivatives by the declared convergent grid stencil and
second/third derivatives from the radial ODE and differentiated source, with
regular-origin limits. Precompute Cartesian derivatives of real solid
harmonics; this avoids polar coordinate singularities in $W,dW,ddW$.
The existing `york_series` can apply the full two-connection operator to these
independent jets; finite-metric samples below anchor its variations separately.

**Complete exterior data, not just a leading subtraction.** Derive the algebraic
$W_1$ expansion from the localized $b_1$ using source moments and derivatives
of (6),
$$W_{1i}^{\rm alg}(x)=\sum_\alpha
 \frac{(-1)^{|\alpha|}}{\alpha!}
 \left(\int y^\alpha b_{1j}(y)d^3y\right)\partial^\alpha K_{ij}(x). \tag{7}$$
Collect moments into STF tensors and their traces. The Kelvin kernel is
biharmonic, not componentwise harmonic: one trace level survives because
$\Delta^2K=0$ away from the origin. Even an isotropic Gaussian vector source
has an $r^{-3}$ York contribution beyond its $r^{-1}$ monopole. Keeping only
scalar STF multipoles of $W$ would miss it. The finite angular support of the
polynomial-Gaussian $b_1$ and the one surviving trace level give a finite
algebraic expansion; bound its Gaussian remainder using Gaussian radial moments.
Do not use production HEAT-bin coefficients to determine (7).

Differentiate the entire expansion to obtain $T^{\rm alg}$, square it with
full tensor contractions, and project every resulting power into real scalar
harmonics. This supplies all $d_{\ell m,p}$ for $F_4$ in (3)–(4). The other
$F_4$ terms have positive Gaussian factors even though $u_3,W_2$ have tails.
For $r=b_2-\mathcal A_1W_1$, source moments are obtained by independently
refined radial quadrature with Gaussian remainder bounds; apply the same York
exterior construction to $W_2$. Source-moment uncertainty must propagate into
its exterior boundary data and final pairings.

For an explicit leading check let $M=\int b_1$, $n=x/r$ and
$\mathcal A=15|M|^2/(32\pi)$. Direct differentiation of (6) gives
$$T:T=\frac{9\,[|M|^2+2(n\cdot M)^2]}{128\pi^2r^4}+O(r^{-5}),
\qquad \int_{r>R}s_4=\frac{\mathcal A}{R}+O(R^{-2}). \tag{8}$$
The $r^{-4}$ source has only $\ell=0,2$: its possible $\ell=1$ resonance is
absent by parity. The corresponding particular term in $u_4$ is
$$\frac{9[(n\cdot M)^2-2|M|^2]}{2048\pi^2r^2}, \tag{9}$$
in addition to homogeneous multipoles. Formula (3) includes an interior
constant $15|M|^2/(2048\pi^2R^2)$ from the omitted leading monopole source.
Truncating that source biases $-\langle B_2,u_4\rangle$ by
$+(\int B_2)15|M|^2/(2048\pi^2R^2)$ at leading order.
These are analytic checks, not fitted tail coefficients. Parity may cancel
some subleading spherical integrals; derive those cancellations before using them.

Write $u\sim\mu_u/r$ and $v=\mathcal CB_2\sim\mu_v/r$. The original energy's direct
exterior contribution is $-7\mu_u\mathcal A/(2R^2)$; the production triple term
has $-(7\mu_u+\mu_v)\mathcal A/(2R^2)$. The additional $-\mu_v$ contribution is transmitted
through the exterior source of $u_4$. Appending the production triple tail to
the unreduced density would double-count this effect. Both designs also need
their full subleading exterior powers and bounded Gaussian remainders.

**Source anchors and finite-ball identities.** Before interpreting a mixed
comparison, gate the following independently of its final total:

- Exact TT identities for the raw basis; envelope rates; $E_1,E_2$ matrix
  coefficients. At held-out spatial points compare $R_2$–$R_4$, $\mathcal L_6$,
  $b_1$–$b_3$ with the finite-exponential metric evaluator for this same mixed
  $Q,S$. Extract $\mathcal A_1W,\mathcal A_2W$ from finite-metric York samples
  holding independent $W$ jets fixed. Reuse the previously declared Decimal
  precision/refinement method, not the expected Taylor series as a reference.
- Check inverse units directly: $\mathcal C e^{-r^2}=\sqrt\pi\operatorname{erf}(r)/(32r)$,
  with origin limit $1/16$, and $[\mathcal Y(c e^{-r^2})](0)=-11c/24$.
  Also require the gradient-Gaussian identity after (6), including its sign
  and nonzero localized monopole. These are operator normalization controls.
- Verify $-8\Delta u= s_2-R_2$, $\Delta_LW_1=b_1$,
  $\Delta_LW_2=r$, $-8\Delta u_3=F_3$, $-8\Delta u_4=F_4$ with a different
  derivative stencil and held-out radii. Check angular projection/reconstruction
  residuals and regular-origin limits. Residuals alone do not certify boundary data.
- Recover gravitational H3/H4 from this spatial route and the raw archived
  quartic polynomial. Obtain, then test, $\int F_3=H_3$ and $\int F_4=H_4$;
  consequently the scalar monopoles are $H_3/(32\pi r)$ and $H_4/(32\pi r)$.
  **Do not impose those archive values as boundary conditions and call the
  resulting agreement a check.** Include the algebraic source tail in $\int F_4$.
- Check $\int s_3=0$, $\int s_4=-a^{\mathsf T}g$ and the ten canonical
  component targets. The $s_4$ integral needs its $O(R^{-1})$ tail (8).
  Here $a=\mathsf G_0^{-1}p$ denotes raw S coefficients.
- Recover the pure-momentum $77\int s_2(\mathcal Cs_2)^2$ reference with this
  backend. Use the archived radial-Legendre value as the independent target,
  not the archive's HEAT-derived `pure_p_H6` convenience field.

For $T=L_0W_1$, the finite-ball normalization is
$$\int_{B_R}T:T=-2\int_{B_R}b_1\cdot W_1
                  +2\oint_{S_R}W_1\cdot(Tn). \tag{10}$$
The last term starts at $-\mathcal A/R$, and cannot be discarded in a norm
check because the source pairing itself has Gaussian integration tails.
Apply the analogous identity to $W_2$ with its own boundary flux.

If $\mathcal D_{\rm weak}$ denotes the five physical weak-A2 densities, then
$$\int_{B_R}(2W_1\cdot\mathcal A_2W_1-\mathcal D_{\rm weak})
=2\oint_{S_R} n_a W_1^i
 [\delta^{ja}(T_2)_{ij}-Q^{ja}(T_1)_{ij}+\tfrac12(Q^2)^{ja}T_{ij}]. \tag{11}$$
This tests the weak replacement by a strong spatial evaluation. Do not gate
pointwise equality. The flux has Gaussian metric factors at this order.
An extra independent scalar solve for $v=\mathcal CB_2$ permits the check
$$\int_{B_R}B_2u_4-\int_{B_R}vF_4
=8\oint_{S_R}(v\partial_nu_4-u_4\partial_nv). \tag{12}$$
Its leading value is $15(\int B_2)|M|^2/(2048\pi^2R^2)$, matching the missing
interior-source effect above. Evaluate both volume and surface sides; do not
infer one from the other or subtract a fitted mismatch.

**Bounded discretization and acceptance.** For the stated $j=2$ free modes the
scalar responses need $\ell\le4,6,8$; vector total angular momentum is at most
four for $W_1$ and six for $W_2$. Cartesian vector components can require
orbital $\ell\le5,7$ respectively. This is not a free-$j=2$ intermediate cutoff.
Retain all $m$ and a guard band; test coefficients beyond the predicted support.
The final scalar sixth-order density has total support at most twelve.

A concrete initial finite-difference schedule has three separate slices:

- Radial: $R=24$, intervals $128,256,512$, fixed fine angular grid $(20,48)$.
- Angular: $R=24$, fixed $512$ intervals, grids $(12,32),(16,40),(20,48)$.
- Exterior: $R=18,24,30$, fixed spacing $3/64$, hence $384,512,640$ intervals,
  and fixed fine angular grid $(20,48)$.

Include these cross levels in the work preflight, counting the common finest
$R=24$ level only once. Use scalar guard modes through $\ell=10$ and Cartesian
vector guard modes through $\ell=9$ when checking the predicted support.
These are proposed bounded levels, not a claim of sufficient accuracy.
Reuse overlapping radial geometry and harmonic tables. Upgrade neither order nor
work caps silently if the declared levels fail. The origin uses regular solid
harmonic limits; the exterior uses (2)–(7), not an empirical fitted Robin value.

Separate radial, angular and exterior-radius changes. For a second-order radial
solver, record three-level behavior before trusting Richardson extrapolation;
also retain raw differences and summed absolute per-family changes. Use full
sphere measure $4\pi w_\mu/n_\phi$ when `run.gauss_nodes` weights sum to one;
radial weights include the full panel length and $r^2$. An orthonormal harmonic
projection must include the same full sphere measure.

Proposed diagnostic ceilings for a future preregistration are $10^{-8}$ for
held-out finite-metric/source identities, $10^{-10}$ for scaled angular
projection/support residuals, $10^{-5}$ for independent-stencil PDE and
regular-origin residuals, and $10^{-6}$ for finite-ball flux identities and
lower-order energy controls. Scale a comparison by
$\max(1,|\text{left}|,|\text{right}|)$; retain absolute errors as well. These
distinct diagnostics do not relax the final physical H6 budget below. If the
second-order grids fail a diagnostic, report failure instead of tuning it to
the result. Source, flux and PDE checks also need the expected refinement
behavior; the residual threshold alone is not a convergence claim.

Choose the physical acceptance budget before running:
$B=\mathrm{atol}+\mathrm{rtol}\max(|H_6^{\rm ref}|,|H_6^{\rm prod}|)$, with
no loosening of the declared pilot tolerances. Require separately passed source,
boundary and lower-order gates, then
$$|H_6^{\rm ref}-H_6^{\rm prod}|+E_{\rm ref}+E_{\rm prod}\le B. \tag{13}$$
$E_{\rm ref}$ sums radial/angle/radius indicators, analytic exterior remainder
and source-moment propagation estimates, and a roundoff allowance based on
absolute term sums before cancellation. Report all pieces in physical energy
units. Production local values are already physical, Fourier normalized values
multiply $\pi^{3/2}$, triple normalized values multiply $4\pi$, and canonical
rational factors multiply $\pi^{3/2}/2$. Do not add dimensionless PDE residuals
to an energy error budget. Observed refinements and floating majorants remain
diagnostics, not rigorous interval certificates.

The scalar radial factorizations scale as $O(N_r N_{\rm harm})$ per right-hand
side, while direct angular transforms/source evaluation can cost
$O(N_rN_\Omega N_{\rm harm})$ plus polynomial/jet work. York needs three scalar
component solves and one auxiliary scalar solve; two York responses and three
scalar responses are required, with optional $v$ adding one. Profile source
evaluation, transforms, derivative reconstruction and each solve separately.
Suggested outer ceilings for a subsequently preregistered reference pilot within the applicable standing authorization are
one process, 7,200 seconds wall time, 4 GiB RSS (or a lower declared host budget),
and 25 million point/assembly-stage evaluations over all control and refinement
levels. Count each visit to a spatial point for each source/operator assembly
stage, not only distinct grid points. Preflight harmonic-transform and jet work
separately, and checkpoint cumulative counts. Include a bounded first-grid
profile in the same ceiling; do not silently launch the remaining slices if
its measured work estimate exceeds the declared cap. Store harmonic radial
arrays and stream angular shells rather than retaining every tensor jet on
every grid. These ceilings are design proposals, not resource availability or
runtime predictions, and are not permission to launch. No timing, memory-fit
or twelve-hour completion claim is supported yet. Passing this comparison
establishes one fixed-mode mixed H6 coefficient
under the stated shared-geometry assumptions; it does not establish a complete
H6 tensor, a four-point loop, scattering, or a physical large-rank limit.
