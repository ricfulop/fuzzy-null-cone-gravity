# First Dirac stabilization and an exact-Einstein evolution obstruction

27 November 2026. Starting source: `2c4abce17ccbb3bee29fb37efed29ae771e98617`. This continues the [actual-action canonical calculation](../canonical-analysis-2026-11-27/README.md), using the same [metric-coefficient action](../../sections/nonlinear_action.tex). The user approved this stabilization step in this chat.

**Result:** all retained secondary constraints can be preserved instantaneously at the rank-one vacuum witness for $N=3,4$. A particular multiplier rate is

$$
\dot h_{tt}=\frac{32\kappa}{9}P_2(\cos\theta),
\qquad \dot\alpha=-\frac{16\kappa}{9}P_2(\cos\theta).
$$

The full temporal Hessian at this witness has ranks 8 and 12. The first preservation step leaves 26 and 50 arbitrary radial functions in the temporal coefficient rates, respectively. These are **not gauge counts**.

**A second, negative result:** an exact continuum vacuum Einstein evolution through this fixed initial datum cannot stay in the current $N=3$ or $N=4$ coefficient ansatz. A nonzero discarded trace in the second time derivative of the spatial metric cannot be removed by any retained lapse or shift rate. This is an exact local obstruction with the hypotheses below; it is not a general impossibility theorem for finite gravity models or for controlled approximation.

**Still open:** existence of a continued solution of the retained equations, the tertiary constraints as functions on an appropriate neighborhood, final first/second-class classification, physical degrees of freedom, boundary quantization, and the simulation connection. The [full release gate](../../GOAL.md) is not met.

## 1. Scope and variables

Retain all coefficient multiplets through $L=N-1$. Use the previous conventions

$$
a=h_{tt},\quad b=h_{tr},\quad \lambda=h_{rr},\quad
V=A_r,\quad W=A_t,
$$

$$
\sigma_{rr}=\lambda+\gamma_{AB}V^AV^B,\quad
\sigma_{rA}=\gamma_{AB}V^B,\quad \sigma_{AB}=\gamma_{AB},
$$

$$
\alpha=\sqrt{-a+b^2/\lambda},\qquad
\beta^r=b/\lambda,\qquad \beta^A=W^A-(b/\lambda)V^A.
$$

All products and inverse metrics are ordinary pointwise operations after reconstruction. The action has not been replaced by a Hamiltonian-first truncation.

On the regular spatial Legendre domain, temporal momenta $\pi_A=0$ are primary and $\mathcal C_A=\delta H/\delta z^A=0$ are secondary, with $z=(a,b,W)$. Their preservation is

$$
\dot{\mathcal C}_A=\{\mathcal C_A,H\}
+\int dr'\,\frac{\delta^2H}{\delta z^A(r)\delta z^B(r')}\mu^B(r')=0,
\qquad\mu^A=\dot z^A.
$$

The calculation uses compact radial smearings on a shell away from $r=0$. Radial derivatives are retained. Normalize the action by $\bar I=(16\pi G)/(4\pi)I$ and use angular average $d\Omega/(4\pi)$ as before. This is a bulk calculation, not a boundary or joint prescription.

## 2. Initial data and the nonlinear embedding acceleration

The witness has $\Lambda=0$,

$$
\sigma=dr^2+r^2d\Omega^2,\quad
a=-1,\quad b=W=V=C=0,\quad \lambda=1,\quad X=r^2,
$$

and constant Cartesian extrinsic curvature $K_{ij}=\kappa n_i n_j$, with $n=(0,0,1)$ and $\kappa\ne0$. Its retained momentum is obtained from $\dot\sigma=2K$ using the actual Legendre map. It satisfies every continuum vacuum constraint and hence all retained primaries and secondaries. The current question is preservation, not a new assumption that the constraints hold.

For an orthonormal spherical tensor write components as

$$
U_{rr}=f,\quad U_{r\theta}=v,\quad
U_{\theta\theta}=g+c,\quad U_{\phi\phi}=g-c.
$$

With $x=\cos\theta$, the entries of $K/\kappa$ are

$$
\left(f,g,v,c\right)=
\left(x^2,\frac{1-x^2}{2},-x\sqrt{1-x^2},\frac{1-x^2}{2}\right).
$$

Let $T_I=\partial\sigma/\partial q^I$ for spatial coefficients $q=(\lambda,X,V,C)$. Even a path with $\ddot q=0$ has spatial metric acceleration

$$
B=T_{IJ}\dot q^I\dot q^J.
$$

Only the products $\gamma V^2$ and $\gamma V$ contribute at this witness. Direct differentiation gives

$$
\frac B{\kappa^2}=
\left(8x^2(1-x^2),\ 0,\ -8x(1-x^2)^{3/2},\ 0\right).
$$

Omitting this term would give incorrect stabilization equations. The verifier checks it by differentiating the full reconstructed spatial metric along a linear coefficient path.

## 3. Normal acceleration and secondary preservation

At the round point let $P_L$ project onto the complete retained spatial tangent space and $Q_L=1-P_L$. The projector is self-adjoint in the round DeWitt pairing. Use $R$ for the difference between the Galerkin spatial acceleration and the continuum Einstein acceleration with the same initial data and temporal rates.

The continuum reference acceleration is

$$
\ddot\sigma_{\rm GR}=2\kappa^2 n\otimes n
+2\dot\alpha K+\mathcal L_{\dot\beta}\sigma.
$$

For example, its zero-rate part follows from the flat continuum metric $-dt^2+dx_1^2+dx_2^2+(1+\kappa t)^2dz^2$. This reference is used to evaluate an initial jet; it is not assumed to remain in the retained coefficient ansatz.

The retained spatial Euler–Lagrange equations set the retained component of the acceleration mismatch to zero. Since $\dot\alpha=-\mu_a/2$ at the witness,

$$
R=Q_L\left(B-2\dot\alpha K-\mathcal L_{\dot\beta}\sigma\right)
=\boxed{Q_LB+Q_L(\mu_aK)}.
$$

Every retained shift rate has a tangent Lie derivative at this round spatial geometry. For a radial shift $\rho(r)Y_{jm}\partial_r$, the spatial components include $2\rho'Y_{jm}$, $\rho D_AY_{jm}$, and $2r\rho Y_{jm}q_{AB}$. For an angular shift $w^A$, they include $r^2q_{AB}\partial_rw^B$ and $r^2(D_Aw_B+D_Bw_A)$. None changes angular degree. This establishes the cancellation while retaining arbitrary radial profiles and their derivatives.

The continuum constraint propagation terms vanish initially because the full continuum constraints are zero. Derivatives of field-dependent smearing weights likewise multiply zero initial constraint densities. The remaining scalar preservation equation is

$$
\boxed{\dot{\mathcal C}_a[f]
=-\frac12\langle Q_L(fK),\,Q_LB+Q_L(\mu_aK)\rangle_{\rm DW}.}
$$

The overall positive radial volume weight is suppressed here. The momentum preservation equations pair $R$ with retained shift Lie derivatives after compact radial integration by parts and therefore vanish. A separate coordinate calculation below verifies those equations without relying solely on this projection argument.

### Full temporal Hessian at the witness

At fixed round spatial geometry, every retained shift is tangent. In the Legendre formula, the shift-quadratic terms cancel and its remaining contribution is linear in the shift. Mixed lapse–shift derivatives vanish. Although the derived lapse depends quadratically on $b$, $\alpha_{,b}=0$ at $b=0$, and its second derivative multiplies the initially zero continuum Hamiltonian constraint density.

Consequently the full temporal Hessian has only the scalar lapse block at this point, with rank exactly $4L$. This is an angular coefficient rank statement, valid for arbitrary compact radial smearings. The full Poisson matrix that also contains secondary–secondary brackets still has only the previously proved lower bound $8L$; its exact rank is not claimed.

## 4. Solution and null compatibility

Define $T(u)=Q_L(uK)$. The previous proof gives $\ker T=V_{L-2}$ and a positive Gram form on its image. Every left-null scalar preservation condition is automatic because it pairs $T(u)=0$ with $Q_LB$. Shift null conditions also vanish by the tangency argument.

For both $L=2,3$, in unit-$\kappa$ angular normalization,

$$
\langle T(P_2),Q_LB\rangle_{\rm DW}=-\frac{128}{525},\qquad
\langle T(P_2),T(P_2)\rangle_{\rm DW}=\frac{12}{175}.
$$

The source has $m=0$ by axial symmetry; parity and the cutoff eliminate the other retained degrees. Thus

$$
\boxed{\mu_a=\frac{32\kappa}{9}P_2+u_{\rm low},
\qquad u_{\rm low}\in V_{L-2}.}
$$

The coefficients of $u_{\rm low}$ are arbitrary functions of $r$ at this step. All $\mu_b$ and $\mu_W$ rates are also free at this step.

| $N$ | Temporal coefficients | Determined rate combinations | Free rate functions at this instant |
|---|---:|---:|---:|
| 3 | 34 | 8 | 26 |
| 4 | 62 | 12 | 50 |

This proves that the first secondary-preservation equations do not immediately exclude the witness. It does not construct a tertiary constraint function on a neighborhood or identify these free rates with gauge freedom.

## 5. Exact-Einstein evolution obstruction

The trace of the normal acceleration mismatch is independent of **every** retained temporal rate, before solving stabilization:

$$
\operatorname{tr}Q_L(\mu_aK)=Q_L^{(0)}(\kappa\mu_a)=0,
\qquad
\boxed{\operatorname{tr}R=-\frac{64\kappa^2}{35}P_4(x)\ne0}
\quad (L=2,3).
$$

Here $\operatorname{tr}B=8\kappa^2(x^2-x^4)$, whose discarded part is the displayed $j=4$ term. Retained shift rates cannot change it because their Lie derivatives are tangent.

**Corollary.** In the prescribed foliation and coefficient parametrization, no sufficiently smooth continuum vacuum Einstein trajectory through the specified $(q,p,z)$ datum can remain in the $N=3$ or $N=4$ ansatz even locally in time. An Einstein trajectory in the ansatz would require $R=0$, which the trace obstruction excludes at the spatial acceleration level.

The hypotheses matter: $\Lambda=0$, nonzero constant $\kappa$, the initial round spatial geometry and temporal coefficients stated above, the actual nonlinear metric reconstruction, the fixed angular cutoff, and regular time derivatives. The initial datum has a continuum vacuum development outside this coefficient ansatz. The result does not exclude approximate Galerkin evolution, a different regulator, different admissible initial data, higher ranks, or a controlled limiting construction.

### Independent constraint-derivative check

For the affine family of rates solving the first preservation step, the entire residual $R$ is unchanged by the remaining free rates. In coordinates $(r,x=\cos\theta,\phi)$ define the undensitized continuum constraints

$$
h=K_{ij}K^{ij}-K^2-{}^{(3)}R,\qquad
m_i=-2\nabla_j(K^j{}_i-\delta^j{}_iK).
$$

Direct differentiation using $\dot\sigma=2K$ and the Galerkin acceleration gives

$$
\dot h=\kappa^3\left(\frac{256}{385}P_4+\frac{128}{231}P_6\right),
$$

$$
\dot m_r=-\frac{64\kappa^2}{21r}P_4,\qquad
\dot m_x=\frac{16\kappa^2}{105}\partial_xP_4,\qquad
\dot m_\phi=0.
$$

Their ADM density versions carry the factor $\sqrt\sigma=r^2$ in these coordinates, apart from the common gravitational normalization. The nonzero $1/r$ term comes from the retained radial derivative calculation. All projections onto the retained scalar and vector harmonics through $L=2$ or 3 vanish. Every nonzero magnetic-number test vanishes by the azimuthal integral, and the axisymmetric axial vector test vanishes because $\dot m_\phi=0$.

Thus the retained secondary constraints are instantaneously preserved while discarded continuum constraint components start to grow. This is consistent with an angular Galerkin approximation; it prevents presenting this check as exact finite-rank Einstein evolution.

## 6. Verification and what remains

The [local verifier](../../verification/dirac_stabilization_20261127/stabilize_witness.py) performs exact symbolic reconstruction, all-$m$ scalar preservation using the previously checked Gaunt matrices, and independent coordinate differentiation of the undensitized continuum constraints. It also checks the source and Gram denominator by numerical quadrature at 24, 48, and 96 nodes for each rank. The largest rate-coefficient discrepancy is $5.329070518200751\times10^{-15}$.

Run from `paper1/`:

```bash
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 \
  python verification/dirac_stabilization_20261127/stabilize_witness.py
```

See [exact outputs and independent checks](checks.json), [verification log](verification.log), and [source hashes and environment](reproducibility.json). No simulation experiment or GPU job was launched. Two independent delegated reviews checked the preservation equation, radial shift argument, multiplier counts, and exact-Einstein corollary. These reviews support the derivation; they are not additional measured experiments.

The full Dirac problem is not resolved by this instantaneous solution. To obtain tertiary constraints as functions that can themselves be preserved, one must determine the rank and kernel of the temporal Hessian along the secondary surface, or analyze different rank strata separately. A pointwise kernel at the round witness is not enough. The next proposed calculation is a neighborhood or time-series consistency analysis with radial dependence, checking whether a continued retained solution exists and which of the presently free rates become fixed. No regulator change is implied or authorized by this note.

The Einstein-lift obstruction is stronger than a numerical residual: it is an exact failure of a particular finite-rank ansatz to contain a specified classical development. It remains narrower than the full nonlinear-consistency objective. The manuscript, archived simulations, and Paper II's distinct regulator are unchanged.
