# Canonical analysis of the actual Paper I action: primary–secondary obstruction

27 November 2026. Starting tree: `836b965` (current main incorporated before this work). This note continues the [nonlinear-consistency checkpoint](../nonlinear-consistency-2026-11-27/README.md) and the [full research goal](../../GOAL.md).

**Result.** For the actual metric-coefficient Einstein–Hilbert Galerkin action, the spatial Legendre map is regular near round positive spatial geometry. There are $4N^2-2$ temporal primary constraints in that domain. An explicit vacuum data set satisfies all primaries and secondaries, yet its scalar lapse-coefficient Hessian has rank $4(N-1)$. The ranks are 8 at $N=3$ and 12 at $N=4$. Consequently the Poisson matrix of the primary–secondary constraints has rank at least 16 and 24, respectively, at these data.

**Limit.** This is a theorem on the primary–secondary surface, not on a proved final Dirac surface. Tertiary stabilization could restrict or exclude the witness. This note does not classify every constraint, identify a generic physical degree-of-freedom count, construct the null boundary quantum theory, or connect the simulated action to gravity. The full nonlinear-consistency release gate remains unmet.

The proof below addresses the actual action and its derived lapse. It does not substitute an independently truncated ADM lapse, assume a positive DeWitt form, or rely on an off-constraint spin-two block. The canonical starting formulas are the ordinary ADM decomposition; see Arnowitt, Deser and Misner, Eqs. (3.9)–(3.15). [Original source, pp. 7–8](https://www.alphaxiv.org/pdf/gr-qc/0405109?page=7).

## 1. Exact ADM map and regular Legendre domain

Use the metric in the [current action appendix](../../sections/nonlinear_action.tex), with ordinary unprojected products after angular reconstruction. Put

$$
a=h_{tt},\quad b=h_{tr},\quad \lambda=h_{rr},\quad
V^A=A_r^A,\quad W^A=A_t^A.
$$

Choose a spacelike $t$ foliation on the domain

$$
\lambda>0,\qquad \gamma_{AB}>0,\qquad
\alpha^2=-a+\frac{b^2}{\lambda}>0.
$$

The exact spatial metric and derived lapse and shift are

$$
\sigma_{rr}=\lambda+\gamma_{AB}V^AV^B,\qquad
\sigma_{rA}=\gamma_{AB}V^B,\qquad
\sigma_{AB}=\gamma_{AB},
$$

$$
\boxed{\alpha=\sqrt{-a+b^2/\lambda},\quad
\beta^r=b/\lambda,\quad
\beta^A=W^A-(b/\lambda)V^A},\qquad
\det\sigma=\lambda\det\gamma.
$$

Although $a,b,\lambda,V,W,X,C$ are band limited, these derived ADM fields generally are not. The code checks the complete metric reconstruction algebraically.

Let $q^I=(\lambda,X,V,C)$ denote all spatial coefficient fields and $z^A=(a,b,W)$ all temporal coefficient fields. At each radial position, for $N\ge3$,

$$
d=\dim q=6N^2-10,\qquad m=\dim z=4N^2-2.
$$

The angular truncation leaves a field theory in $(t,r)$; these are coefficient counts per radial point, not a finite total phase-space dimension.

For compact bulk variations, subtract the usual ADM total derivative to use a first-order-in-time representative. This changes the canonical boundary convention, not the bulk equations. It does not supply the still-missing null boundary and joint action. Define

$$
T_I=\frac{\partial\sigma}{\partial q^I},\qquad
D=\mathcal L_\beta\sigma,\qquad E=T_I\dot q^I-D,\qquad K=E/(2\alpha).
$$

At fixed $q,z$, use the symmetric bilinear form

$$
\mathcal B(U,V)=\frac1{32\pi G}\int d^2y\,
\frac{\sqrt{\sigma}}{\alpha}
\left(U_{ij}V^{ij}-\operatorname{tr}_\sigma U\operatorname{tr}_\sigma V\right).
$$

Write $\mathsf W_{IJ}=\mathcal B(T_I,T_J)$, $d_I=\mathcal B(T_I,D)$, and

$$
\mathcal U=\frac1{16\pi G}\int d^2y\,\alpha\sqrt\sigma\,
({}^{(3)}R-2\Lambda).
$$

Then the actual coefficient Lagrangian, momentum map, and Hamiltonian density are

$$
\mathcal L=\tfrac12\dot q^T\mathsf W\dot q-d^T\dot q
+\tfrac12\mathcal B(D,D)+\mathcal U,
\qquad p=\mathsf W\dot q-d,
$$

$$
\boxed{\mathcal H=\tfrac12(p+d)^T\mathsf W^{-1}(p+d)
-\tfrac12\mathcal B(D,D)-\mathcal U.}
$$

This formula includes all spatial velocities. $D$ can depend on radial derivatives of temporal coefficients, so secondary equations and their brackets must be treated as functional differential expressions, not automatically as finite pointwise matrices.

At a round spatial metric with $V=0$ and uniform positive lapse, express a variation in an orthonormal frame as a radial scalar $f$, angular half-trace $g$, mixed vector $v_A$, and angular tracefree tensor $C_{AB}$. Its DeWitt form is

$$
U:U-(\operatorname{tr}U)^2
=-4fg-2g^2+2|v|^2+|C|^2.
$$

For each scalar harmonic the $(f,g)$ block is

$$
\begin{pmatrix}0&-2\\-2&-2\end{pmatrix},
\qquad \det=-4,\quad \text{eigenvalues }-1\pm\sqrt5.
$$

The vector and tracefree tensor blocks are positive. Thus the full spatial Hessian has signature $(5N^2-10,N^2,0)$, up to an overall positive factor. It remains invertible on an open neighborhood of this point. A self-null radial variation is not a kernel vector because it pairs with $g$.

| $N$ | Spatial coefficients / Hessian rank | Positive / negative eigenvalues | Temporal primaries |
|---|---:|---:|---:|
| 3 | 44 | 35 / 9 | 34 |
| 4 | 86 | 70 / 16 | 62 |

These unreduced signs do not establish physical ghosts. Such a statement requires the final constraint reduction.

## 2. First two Dirac stages

The first-order representative has no $\dot z$, giving all temporal primaries

$$
\pi_A(r)=0,\qquad A=1,\ldots,m.
$$

On the regular Legendre domain there are no additional spatial primaries. The total Hamiltonian and secondary equations are

$$
H_T=H[q,p,z]+\int dr\,\mu^A\pi_A,
\qquad \mathcal C_A(r)=\frac{\delta H}{\delta z^A(r)}=0.
$$

The primary–secondary bracket is

$$
\{\pi_A(r),\mathcal C_B(r')\}
=-\frac{\delta^2H}{\delta z^A(r)\delta z^B(r')}.
$$

The secondary functions are weighted variations through the derived $\alpha,\beta$ and the velocity recovered from the finite Legendre map. They are not obtained by declaring an independently truncated continuum Hamiltonian linear in retained lapse and shift coefficients.

The previously established base diffeomorphisms and base-dependent rotations remain exact bulk symmetries. Their five gauge functions do not determine the classification of all the other coefficients. A kernel vector at one background is not, by itself, a first-class generator.

## 3. An admissible vacuum data set

Set $\Lambda=0$ and use Euclidean spatial geometry on a shell away from $r=0$:

$$
\sigma_{ij}dx^idx^j=dr^2+r^2d\Omega^2,\quad
\lambda=1,\quad X=r^2,\quad V=C=0,\quad
a=-1,\quad b=W=0.
$$

In Cartesian coordinates choose constant nonzero $\kappa$ and

$$
K_{ij}=\kappa n_i n_j,\qquad n=(0,0,1).
$$

This satisfies the vacuum Hamiltonian and momentum constraints pointwise:

$$
{}^{(3)}R=0,\quad
(\operatorname{tr}K)^2-K_{ij}K^{ij}=\kappa^2-\kappa^2=0,
\quad D_i(K^i{}_j-\delta^i{}_j\operatorname{tr}K)=0.
$$

The last equality follows from Cartesian constancy and is also checked explicitly in spherical coordinates. With $x=\cos\theta$, the orthonormal spherical components divided by $\kappa$ are

$$
f=x^2,\qquad g=\tfrac12(1-x^2),\qquad
v_\theta=-x\sqrt{1-x^2},\qquad
c=\tfrac12(1-x^2),
$$

where the angular components are $g+c$ and $g-c$. They contain only $j=0,2$ scalar components and $j=2$ vector/tensor components. Consequently $\dot\sigma=2K$ is in the full retained tangent space for every $N\ge3$. Take its momentum using the actual finite Legendre map, and set all $\pi_A=0$.

The envelope identity for the Legendre transform gives $\delta_z H=-\delta_z\mathcal L$ at the recovered velocity. All continuum constraint densities vanish at this data set, so every retained temporal variation vanishes as well. This places the data on the common zero set of **all primaries and secondaries**, not merely a selected constraint. It does not prove that a time evolution within the truncated ansatz preserves the data or that later constraints vanish.

## 4. Exact multiplier Hessian and its rank

Put $L=N-1$. At the round point let $P_L$ be the full spatial tangent harmonic projector and $Q_L=1-P_L$. It retains all scalar, vector and tensor copies at each allowed angular degree. It is self-adjoint for the round DeWitt pairing and commutes with the spatial trace.

Use the normalized action $\bar I=(16\pi G)/(4\pi) I$ and angular average $\langle\cdot\rangle=\int d\Omega\,(\cdot)/(4\pi)$. At $r=1$, the kinetic Lagrangian is

$$
\bar{\mathcal L}_{\rm kin}
=\tfrac14\langle E,\alpha^{-1}E\rangle_{\rm DW}.
$$

Vary the **actual retained coefficient**

$$
a=-(1+\epsilon u),\qquad u\in V_L,
\qquad\alpha=\sqrt{1+\epsilon u},
$$

holding $q,p,b,W$ fixed. Here $V_j$ is the real scalar harmonic space through degree $j$. The potential is zero throughout this variation because the spatial metric remains flat and $\Lambda=0$.

At $\epsilon=0$, $E=2K$. Differentiating the finite inverse Gram operator gives

$$
\bar H''(0)=\langle K,\alpha''(0)K\rangle_{\rm DW}
-2\langle Q_L(\alpha'(0)K),Q_L(\alpha'(0)K)\rangle_{\rm DW}.
$$

The first term vanishes pointwise because $K:K-(\operatorname{tr}K)^2=0$. Using $\alpha'(0)=u/2$ and polarizing yields the exact identity

$$
\boxed{\bar H''[u,v]=-\tfrac12
\langle Q_L(uK),Q_L(vK)\rangle_{\rm DW}.}
$$

The nonlinear relation between $a$ and $\alpha$ has therefore been included, not discarded. The fixed-momentum variation is transverse to the secondary surface; that is the derivative required in the Poisson bracket.

Although the full DeWitt pairing is indefinite, the discarded tensors in this identity are tracefree:

$$
\operatorname{tr}Q_L(uK)=Q_L^{(0)}(\kappa u)=0.
$$

The right-hand pairing is consequently an ordinary positive spatial tensor norm. There is no cancellation from the conformal direction for this witness.

**Kernel theorem.** For every $L\ge2$ and $\kappa\ne0$,

$$
\ker\bar H''\big|_{V_L}=V_{L-2},\qquad
\operatorname{rank}\bar H''\big|_{V_L}=4L.
$$

Proof: $K$ has only $j=0,2$ content, so $u\in V_{L-2}$ produces no discarded harmonics. Conversely, the radial component of $uK$ is $\kappa x^2u$. The coefficient of $Y_{\ell+2,m}$ in $x^2Y_{\ell m}$ is

$$
c_{\ell m}=\sqrt{\frac{((\ell+1)^2-m^2)((\ell+2)^2-m^2)}
{(2\ell+1)(2\ell+3)^2(2\ell+5)}}>0
\quad (|m|\le\ell).
$$

For the highest occupied degree $\ell$ of $u$, lower degrees cannot cancel that contribution. Thus zero leakage forces $u_L=0$ and then $u_{L-1}=0$. The remaining space is exactly $V_{L-2}$. Positivity on discarded tracefree tensors identifies this leakage kernel with the Hessian kernel. Dimension subtraction gives $(L+1)^2-(L-1)^2=4L$.

| $N$ | $L$ | Scalar multiplier dimension | Hessian kernel dimension at witness | Hessian rank | Full primary–secondary Poisson rank, lower bound |
|---|---:|---:|---:|---:|---:|
| 3 | 2 | 9 | 1 | 8 | 16 |
| 4 | 3 | 16 | 4 | 12 | 24 |

The last column follows by selecting lapse primaries and secondaries from the top two scalar multiplets. Their Poisson submatrix has the form $\left(\begin{smallmatrix}0&A\\-A^T&B\end{smallmatrix}\right)$ with $A$ invertible; it is invertible regardless of $B$. Smooth compact radial smearings supported on a shell give the same finite submatrix rank, multiplied by a positive radial integral.

For the unnormalized probe $u=P_2(x)$, both $L=2$ and $L=3$ give

$$
\langle Q_L(uK),Q_L(uK)\rangle_{\rm DW}=\frac{12\kappa^2}{175},
\qquad \bar H''(0)=-\frac{6\kappa^2}{175}.
$$

The corresponding primary–secondary bracket has the opposite sign under the convention $\{\pi_a,\mathcal C_a\}=-\delta_a^2H$. For a radial multiplier bump $w(r)$, restore the physical normalization by multiplying the local unit-$\kappa$ result by

$$
\frac{4\pi}{16\pi G}\,\kappa^2\int dr\,r^2w(r)^2.
$$

The $L=3$ dipoles lie in the lapse-block kernel of this witness because they cannot leak above $L=3$. This is not a proof that they are first class on a general background.

## 5. Verification and independent review

The [verifier](../../verification/canonical_analysis_20261127/actual_adm.py) checks four distinct ingredients:

1. Symbolic ADM reconstruction and determinant identities, plus the full spatial Hessian signature.
2. The vacuum Hamiltonian, momentum, and spatial curvature identities in spherical coordinates.
3. Direct rational polynomial projections for five axisymmetric cases, and separate full-$m$ spin-weighted Gaunt matrices at $N=3,4$. Both spin signs are retained to test off-diagonal cancellations. The reported rank includes every scalar magnetic mode.
4. Independent Gauss–Legendre quadrature and finite differences of the actual inverse Gram Hamiltonian, including its indefinite scalar block and positive vector and tensor blocks. Two quadrature orders and two difference steps give twenty comparisons. Axisymmetry decouples the unused magnetic sectors in this numerical check.

Reproduce from `paper1/`:

```bash
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 \
  python verification/canonical_analysis_20261127/actual_adm.py
```

The exact matrices and numerical diagnostics are written to `checks.json`; environment information and hashes go to `reproducibility.json`. These are local proof-verification calculations, not Monte Carlo measurements or GPU experiments. The logs preserve actual completion status.

The first packaged invocation stopped at a symbolic Ricci-zero assertion. Its residual was

$$
\frac{\sin(2\theta)}{2\tan\theta}-\frac{\cos(2\theta)}2-\frac12=0.
$$

Expanding the trigonometric functions reduces it identically to zero. With the user's approval, the assertion was corrected to normalize this expression before comparison; no equation or tolerance was weakened. The failed invocation is preserved in `verification-attempt1.log` and the corrected invocation in `verification.log`.

The corrected invocation passed every assertion. The maximum absolute difference between the twenty numerical Hamiltonian curvatures and their exact targets was $2.817857613912622\times10^{-10}$. The exact all-$m$ matrices have the ranks stated above. See [checks](checks.json), [environment and source hashes](reproducibility.json), and [verification log](verification.log).

Two independent delegated reviews checked the ADM map, the rank-one secondary-surface witness, the $\alpha''$ cancellation, the axisymmetric coefficients, and the all-$m$ kernel proof. Both confirmed the result and the need to distinguish the primary–secondary surface from the final Dirac surface. Their agreement is a review of the derivation, not an additional physical experiment.

## 6. Next decision

This result prevents treating the entire obtained primary–secondary set as the standard first-class ADM set without further analysis. It neither proves that every retained lapse becomes fixed nor proves that no consistent completion exists.

The next canonical calculation is the stabilization equation

$$
0\approx\dot{\mathcal C}_A
=\{\mathcal C_A,H\}+\int dr'\,
\{\mathcal C_A,\pi_B(r')\}\mu^B(r').
$$

The nonzero Hessian directions determine corresponding primary multipliers when a consistent solution exists. Contracting this equation with its null directions tests for tertiary constraints and for the Noether identities associated with the surviving gauge functions. This must be done with radial derivatives and boundary conditions stated, first at the witness and then on an open nonspherical domain. Only after stabilization can the complete first/second-class classification and physical degree-of-freedom count be claimed.

The proposed next step is to perform that stabilization before changing the regulator or proceeding to boundary quantization. It is a further research decision under the user's copilot setting. Existing manuscript claims, simulation results, and Paper II's distinct regulator have been left intact.
