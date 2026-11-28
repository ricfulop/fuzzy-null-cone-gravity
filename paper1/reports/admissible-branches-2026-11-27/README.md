# An exact admissible family with a nonzero angular-shift bracket

27 November 2026. Starting source: `fe44bf54fa99529872b3ab47a5f226573570b861`. The user approved investigation of dynamically admissible branches in `msg_59826d22`. This uses the same [Paper I action](../../sections/nonlinear_action.tex), with no change of regulator. It follows the [higher-preservation obstruction](../higher-order-consistency-2026-11-27/README.md).

**Constructive result.** Static Minkowski spacetime in oblate spatial coordinates lies exactly in the actual coefficient ansatz for every $N\ge3$. For $N=3,4$, a retained polar-quadrupole angular-shift direction has the fixed-momentum Hamiltonian Hessian

$$
\boxed{H''_{WW}[\xi,\xi]
=-\frac{96}{1225}\epsilon^2\int dr\,\frac{w(r)^2}{r^2}
+O(\epsilon^3),\qquad \xi=w(r)D^AP_2\partial_A.}
$$

The parameter $\epsilon$ has dimensions of length squared. On a fixed compact radial shell and for fixed nonzero smooth compact $w$, the leading coefficient is strictly negative. The corresponding primary–secondary bracket is therefore nonzero for all sufficiently small nonzero $\epsilon$ on an **exact stationary retained solution**. Its two-constraint Poisson submatrix has rank two. Higher preservation cannot exclude these data, because the solution exists for every time.

This closes the witness-admissibility gap in the earlier nonzero-bracket argument. It rules out treating every retained temporal constraint as first class on the dynamically admissible set. It does not give the full bracket rank, prove a regular surrounding constraint stratum, or count physical degrees of freedom. The family approaches round **zero-momentum** Minkowski data; it does not approach the previously excluded datum with fixed $\kappa\ne0$.

**Near the excluded datum.** A separate continuity argument shows that a sequence of retained solution slices approaching that datum cannot maintain the specified control of first temporal rates, second temporal rates, and radial derivatives. This is a statement with bounds on derivatives, not an exclusion of an entire neighborhood in coefficient values alone.

The [full release goal](../../GOAL.md), including the complete constraint classification, boundary theory, and simulation matching, remains open.

## 1. Exact solution inside the coefficient band

Let $x=\cos\theta$, $d=1-x^2$, and

$$
\rho=r^2+\epsilon,\qquad\Sigma=r^2+\epsilon x^2.
$$

Work on a shell $0<r_-\le r\le r_+$ with $\rho,\Sigma>0$. The map

$$
X_{\rm Cart}=\sqrt\rho\sqrt d\cos\phi,\quad
Y_{\rm Cart}=\sqrt\rho\sqrt d\sin\phi,\quad
Z_{\rm Cart}=rx
$$

pulls the Euclidean spatial metric back to

$$
\sigma_\epsilon=\frac\Sigma\rho\,dr^2
+\frac\Sigma d\,dx^2+\rho d\,d\phi^2,
\qquad \sqrt{\det\sigma_\epsilon}=\Sigma.
$$

In the actual coefficient reconstruction this is

$$
\lambda=h_{rr}=\frac\Sigma\rho,
\qquad X=r^2+\frac\epsilon2(1+x^2),\qquad V=A_r=0,
$$

$$
C_{xx}=-\frac\epsilon2,\qquad
C_{\phi\phi}=\frac\epsilon2d^2,\qquad C_{x\phi}=0.
$$

The tensor is tracefree with respect to the unit round sphere: $dC_{xx}+C_{\phi\phi}/d=0$. Its polar tensor amplitude is proportional to $P_2''$, hence is a pure $j=2$ mode. The scalar coefficients are explicitly

$$
\lambda=\frac{r^2+\epsilon/3}{\rho}
+\frac{2\epsilon}{3\rho}P_2,
\qquad X=r^2+\frac{2\epsilon}{3}+\frac\epsilon3P_2.
$$

Thus every coefficient has angular degree at most two, even though nonlinear reconstructed operations are unprojected. Set

$$
a=h_{tt}=-1,\qquad b=h_{tr}=0,\qquad W=A_t=0.
$$

The derived lapse and shift are $\alpha=1$, $\beta=0$. The four-metric $-dt^2+\sigma_\epsilon$ is static Minkowski, with $K=0$, retained momentum $p=0$, and $\Lambda=0$. Its unrestricted Einstein equations vanish pointwise. Pulling their variation back to retained coefficients therefore gives an exact solution of all retained Euler–Lagrange equations at every time. Primary momenta are zero and every necessary constraint-preservation equation is satisfied along this trajectory.

The verifier proves the Euclidean pullback using cylindrical-radius squared $\rho d$ and Cartesian axial coordinate $rx$, then independently computes all spatial Ricci components and obtains zero. This avoids depending on a coordinate-name identification alone.

## 2. Fixed-momentum shift Hessian

Use the earlier normalization $\bar I=(16\pi G)/(4\pi)I$ and angular average $d\Omega/(4\pi)$. At fixed spatial coefficients define the spatial DeWitt pairing

$$
G_\sigma(U,V)=\frac1{4\pi}\int dr\,d^2y\,\sqrt\sigma
\left(U_{ij}V^{ij}-\operatorname{tr}_\sigma U\operatorname{tr}_\sigma V\right).
$$

Let $T$ map coefficient velocities to spatial tensor velocities. The restricted spatial Gram matrix $T^*G_\sigma T$ is nondegenerate at the round point and remains so on a sufficiently small neighborhood, by the [canonical Legendre result](../canonical-analysis-2026-11-27/README.md). Its inverse defines the DeWitt-orthogonal tangent projector $\Pi_\sigma$ and normal complement $N_\sigma=1-\Pi_\sigma$. This is an orthogonal decomposition in an indefinite form; it is not a positive-norm minimization.

Vary the actual temporal coefficient $W=s\xi$, holding $q,p=0,a=-1,b=0$ fixed. Then $\alpha=1$, $\beta=s\xi$, and the kinetic Lagrangian is

$$
L_{\rm kin}=\frac14G_\sigma(T\dot q-s\mathcal L_\xi\sigma,
T\dot q-s\mathcal L_\xi\sigma).
$$

The fixed-momentum equation $p=0$ gives $T\dot q=s\Pi_\sigma\mathcal L_\xi\sigma$. The spatial curvature potential is zero and independent of this shift variation. Hence

$$
H(s)=-\frac{s^2}{4}
G_\sigma(N_\sigma\mathcal L_\xi\sigma,N_\sigma\mathcal L_\xi\sigma),
$$

$$
\boxed{H''_{WW}[\xi,\xi]=-\frac{1}{2}
G_\sigma(N_\sigma\mathcal L_\xi\sigma,N_\sigma\mathcal L_\xi\sigma).}
$$

The factor and sign come from the actual finite Legendre transform. No continuum Hamiltonian linearity is assumed after truncation.

## 3. Leading normal component and its sign

Choose the retained polar quadrupole $\xi^x=3w(r)x(1-x^2)$, $\xi^r=\xi^\phi=0$. At every $\epsilon$ its mixed spatial variation is

$$
(\mathcal L_\xi\sigma_\epsilon)_{rx}
=\gamma_{xx}\partial_r\xi^x.
$$

It is exactly represented by the retained tangent coefficient $\delta V^x=\partial_r\xi^x$. Thus all $w'$ terms are included and are tangent. This is not a freeze of radial derivatives.

At $V=0$, the remaining diagonal tangent tensors depend linearly on $(\delta\lambda,\delta X,\delta C)$ and their coefficient-to-tensor map is independent of $\epsilon$. The round shift source is in that tangent space. Subtract it and the exact mixed tangent contribution before expanding the normal projector. At first order,

$$
N_{\sigma_\epsilon}\mathcal L_\xi\sigma_\epsilon
=\frac{\epsilon w(r)}{r^2}Q_0U+O(\epsilon^2),
$$

where $Q_0$ is the discarded round projector for the chosen cutoff. In the earlier component convention

$$
(f,g,\widetilde v,\widetilde c)
=\left(U_{\hat r\hat r},\frac{U_{\hat\theta\hat\theta}+U_{\hat\phi\hat\phi}}2,
\frac{U_{\hat r\hat\theta}}{\sqrt d},
\frac{U_{\hat\theta\hat\theta}-U_{\hat\phi\hat\phi}}{2d}\right),
$$

the dimensionless first-order tensor is

$$
U=\left(6x^2d,\ 3x^2(1-3x^2),\ 0,\ 9x^2\right).
$$

Both $L=2$ and $L=3$ discard the same $j=4$ part:

$$
Q_0U=\left(-\frac{48}{35}P_4,
-\frac{72}{35}P_4,\ 0,\ \frac6{35}P_4''\right).
$$

This tensor is **not tracefree**. Its round trace is $-192P_4/35$, so the negative trace term of the DeWitt form must be retained. Direct integration with measure $dx/2$ gives

$$
\langle Q_0U,Q_0U\rangle_{\rm DW}
=\left\langle-4fg-2g^2+2d\widetilde v^2+2d^2\widetilde c^2\right\rangle
=\frac{192}{1225}>0.
$$

Analyticity of the metric, angular integrals, and inverse finite Gram matrix supplies the expansion on a compact shell in the regular Legendre neighborhood. With fixed $w$, the positive leading coefficient proves the nonzero bracket for sufficiently small nonzero $\epsilon$. Neither exact positivity for arbitrary deformation nor a uniform numerical upper bound on the allowed $\epsilon$ is claimed.

## 4. Canonical consequence on an exact trajectory

Let $\pi_W[\xi]$ be the primary paired with this coefficient direction and $\mathcal C_W[\xi]=\delta H/\delta W[\xi]$ its secondary. On the stationary family they both vanish, but

$$
\{\pi_W[\xi],\mathcal C_W[\xi]\}
=-H''_{WW}[\xi,\xi]
=\frac{96}{1225}\epsilon^2\int\frac{w^2}{r^2}\,dr+O(\epsilon^3)\ne0.
$$

The two-by-two Poisson submatrix is invertible for sufficiently small nonzero $\epsilon$. Thus the constraint Poisson matrix has rank at least two at a dynamically admissible point. The same initial point lies on an exact trajectory for all times, so the earlier failure of second preservation does not recur here.

The result does not identify all second-class pairs or establish a regular full constraint surface through this point. In particular, at $p=0$, zero shift, and zero spatial curvature, the lapse–lapse Hessian vanishes: setting the shift to zero makes the Hamiltonian zero for arbitrary lapse at these fixed spatial data and momentum. Additional degeneracies are therefore unsurprising. Their later preservation and secondary–secondary brackets have not been classified.

The oblate family tends to the round $K=p=0$ solution as $\epsilon\to0$. It provides a controlled admissible reference near that round point, not a family accumulating at the excluded fixed-$\kappa$ rank-one datum. Full diffeomorphism equivalence in continuum Minkowski spacetime does not imply that this particular coefficient restriction retains the corresponding gauge transformations.

## 5. What continuity says near the excluded witness

This conclusion uses the previous obstruction and does not require a constant-rank extension of its null vectors. Fix a nonnegative, nonzero compact bump $w$ and the same retained angular test $\xi$. Write $y=(q,p,z)$, $\mu=\dot z$, and $\nu=\ddot z$. In the regular spatial Legendre domain the twice differentiated, fixed-test secondary can be written

$$
\ddot{\mathcal C}_W[\xi]
=F_\xi(y,\mu)+D_z\mathcal C_W[\xi](y)[\nu].
$$

The functional derivatives include the radial derivatives required by the Hamiltonian. At the excluded witness $y_*$ the entire angular-shift row $D_z\mathcal C_W[\xi](y_*)$ vanishes. For every first-preserving $\mu_*$, the prior calculation gives

$$
F_\xi(y_*,\mu_*)=-c_*,\qquad
c_*=\frac{512}{1575}\kappa^3\int r^2w(r)\,dr\ne0.
$$

Suppose a sequence of actual smooth solution slices has $y_j\to y_*$ strongly in radial norms sufficient for continuity of these jet functionals. Also suppose $\mu_j$ is precompact in the required radial norms; bounds one derivative higher on the compact shell are a sufficient way to obtain that precompactness. Pass to a subsequence $\mu_j\to\mu_*$. First preservation passes to the limit, so $\mu_*$ is in the previously determined affine family. Hence

$$
F_\xi(y_j,\mu_j)\to-c_*,\qquad
D_z\mathcal C_W[\xi](y_j)\to0.
$$

Second preservation therefore requires an unbounded second-rate norm, unless the preceding convergence or precompactness assumptions fail. In norms for which the functional is continuous and the indicated operator norm converges to zero, a necessary estimate for sufficiently large $j$ is

$$
\|\nu_j\|\,\|D_z\mathcal C_W[\xi](y_j)\|_{\rm op}
\ge\frac{|c_*|}{2}.
$$

This excludes an approach with the stated control of temporal rates and radial derivatives. It does **not** prove a solution-free neighborhood in $(q,p,z)$ alone; unbounded rates or short radial scales could evade the assumptions.

There is also a useful partial local reduction. The top two scalar lapse multiplets form a $4L\times4L$ invertible block of $D_a\mathcal C_a$ at $y_*$. The Hamiltonian contains no radial derivatives of $a=h_{tt}$: the derived shift is independent of $a$, and lapse dependence is pointwise. The parameter-dependent implicit function theorem therefore solves those $4L$ secondary equations locally for the top lapse coefficients at each radius, in terms of the remaining fields and their required radial jets. This eliminates 8 coefficients at $N=3$ or 12 at $N=4$ on that partial secondary surface. It does not solve the remaining low-lapse and shift constraints, and the graph itself contains the excluded point.

## 6. Verification and scope of the numerical checks

Run the [local verifier](../../verification/admissible_branches_20261127/oblate_branch.py) from `paper1/`:

```bash
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 \
  python verification/admissible_branches_20261127/oblate_branch.py
```

The [exact and numerical outputs](checks.json), [execution log](verification.log), and [environment/source hashes](reproducibility.json) record one successful run. Symbolic checks cover exact reconstruction, degree-two coefficient membership, the independent Ricci calculation, arbitrary radial profiles, and the rational leading coefficient.

As a separate check, the script computes the deformed tangent Gram projection of the specified shift source, without replacing it by the round projector. The background and source are axisymmetric, so all $m\ne0$ sectors decouple. Off-diagonal tensors are pointwise DeWitt orthogonal to the diagonal source. The remaining block contains every diagonal $m=0$ tangent: both scalar coefficient copies through $L$, and the polar tensor copy for $2\le j\le L$. The mixed source proportional to $w'$ has already been removed as an exact tangent. These facts make the block solve complete for this probe on the regular Legendre domain.

At $r=1$, unit angular profile, and 128-node quadrature, the following block-projection values agree between $N=3$ and $N=4$ to the displayed precision. The 64-node calculations provide an independent integration-resolution check.

| $\epsilon/r^2$ | $H''_{WW}/\epsilon^2$ |
|---:|---:|
| $0.0001$ | $-0.07831961919$ |
| $0.01$ | $-0.07364448578$ |
| $0.1$ | $-0.03527510757$ |
| $0.2$ | $+0.00005262112$ |

The small-deformation limit is $-96/1225\simeq-0.07836734694$. The maximum tangent-orthogonality error is $1.084\times10^{-14}$. The direct normal norm also agrees with the Schur-complement evaluation in the archived checks. Odd diagonal multiplets added at $L=3$ decouple from this even source.

The sign change in this tested diagonal entry at larger deformation is compatible with the indefinite DeWitt form. The moderate-deformation rows are numerical block diagnostics, not a proof that every other spatial Legendre block remains invertible there, a location of an exact zero, or a classification of the full temporal Hessian rank. The analytic nonzero-bracket theorem is confined to sufficiently small nonzero deformation where full spatial Legendre regularity follows by continuity.

Two independent delegated reviews confirmed the exact family, harmonic membership, leading contraction, and canonical interpretation. The nearby-excluded-witness argument is an analytic continuity result; it is not presented as a numerical branch search. No simulation experiment, GPU job, or external compute was launched. Archived calculations and the manuscript are preserved.

## 7. Remaining research decision

This checkpoint provides an exact admissible family and a nonzero canonical bracket on it. The next unresolved task is to analyze the remaining constraints and their preservation near a small nonzero oblate point, and determine whether it belongs to a regular constant-rank branch. Only then can a physical degree-of-freedom count be justified. The boundary reduction and simulation connection remain separate unresolved components of the full objective. Further branch analysis requires the next copilot decision.
