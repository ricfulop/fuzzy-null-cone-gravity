# Second-order null compatibility near the oblate solution

Preparation checkpoint, 27 November 2026. Scientific source: 157ec0f04866ea36509b325ec65f0a68890f9b17. User approval msg_7c4562eb, relayed by the coordinating chat, authorizes derivation and implementation preparation only. **No scientific checker or experiment has been executed for this checkpoint.**

The calculation uses the [original Paper I coefficient action](../../sections/nonlinear_action.tex), the [certified oblate shift kernel](../oblate-constraint-rank-2026-11-27/README.md), and the [linearized lapse result](../oblate-lapse-consistency-2026-11-27/README.md). No action, regulator, harmonic convention, or canonical normalization changes.

**Derived condition.** For a background null angular descriptor $\xi$ and an admissible first-order retained solution $c_1$, a twice-differentiable family of exact retained solutions must satisfy

$$
\boxed{B_\xi[c_1]
=E_{\rm full}^{(1)}[c_1]\big[Q_c\,DF_\xi(c_0)c_1\big]=0.}
$$

Here $E_{\rm full}^{(1)}$ is the unrestricted continuum coefficient Euler residual, including temporal components; $Q_c$ discards coefficient harmonics above the retained band; and $F_\xi$ is the exact angular diffeomorphism of the coefficient types. This is a necessary second-order preservation condition, not an evaluated obstruction.

For the axisymmetric $N=3$ descriptor $\xi=w(r)x\partial_\phi$, the discarded test has the eight explicit channels in section 4. The condition is independent of the unknown second-order correction. No admissible non-gauge first-order seed has yet been supplied or tested, so no nonzero dynamical source or new final Dirac rank is claimed.

## 1. Which null directions are being examined

Put $L=N-1=2$, $x=\cos\theta$, $d=1-x^2$, $\rho=r^2+\epsilon$, and $\Sigma=r^2+\epsilon x^2$. The background is the exact static oblate solution

$$
g_0=-dt^2+\frac{\Sigma}{\rho}dr^2+\frac{\Sigma}{d}dx^2+\rho d\,d\phi^2,
\qquad p_0=0,\quad a_0=-1,\quad b_0=V_0=W_0=0.
$$

Work on a fixed radial shell away from zero with positive $\rho,\Sigma$. Angular descriptors have smooth compactly supported radial profiles. The normalizations are $\bar I=(16\pi G)/(4\pi)I$ and angular average $d\Omega/(4\pi)$.

At sufficiently small nonzero $\epsilon$, the seven real shift-kernel profiles at $N=3$ split into one radial base direction and six angular directions. With $K_z=\partial_\phi$ and $n_i$ the three Cartesian coordinate functions on the unit sphere, the angular kernel is

$$
\operatorname{span}_{\mathbb R}\{K_x,K_y,K_z\}
\ \oplus\
\operatorname{span}_{\mathbb R}\{n_xK_z,n_yK_z,n_zK_z\}.
$$

All coefficients may depend on radius. The first summand comprises exact rotations. Together with the radial base symmetry it supplies four controls. The three $n_iK_z$ fields are additional null directions at this background; their null status does not make them nonlinear gauge generators.

For example, in the pinned complex harmonic convention,

$$
A_{20}=-3x\partial_\phi,\qquad
P_{11}-\frac{i}{3}A_{21}
=-2i\sqrt d\,e^{i\phi}\partial_\phi.
$$

The minimal prepared test examines the $m=0$ member $\xi=w(r)x\partial_\phi$. The general compatibility identity applies to all three, but the prepared implementation is axisymmetric only.

## 2. Noether derivation and perturbative admissibility

Use unrestricted coefficient types $c=(q,z)$ with $q=(\lambda,X,V,C)$, $z=(a,b,W)$, and $\gamma=Xq_{S^2}+C$. The finite ansatz restricts their harmonic coefficients by a fixed projector $P_c$; set $Q_c=1-P_c$. This is a coefficient projector, not the deformed DeWitt tensor projector used in the Legendre analysis.

For a fixed time-independent angular vector $\xi(r,\Omega)$, the unrestricted transformations are

$$
F_\xi a=\xi(a),\quad F_\xi b=\xi(b),\quad F_\xi\lambda=\xi(\lambda),
$$
$$
F_\xi\gamma=\mathcal L_\xi\gamma,\qquad
F_\xi V=\partial_r\xi+[\xi,V],\qquad
F_\xi W=[\xi,W].
$$

Take the fixed round trace and tracefree parts of $F_\xi\gamma$ to obtain its $X,C$ components. A time-dependent scalar multiplying $\xi$ introduces its time derivative only in the transformation of $W$. With $E_W[\xi]=-C_W[\xi]$, the unrestricted angular Noether identity is

$$
\partial_t C_W[\xi]=-E_{\rm full}(c)[F_\xi(c)].
$$

This is an off-shell identity of the unrestricted action. Along a hypothetical retained solution, the retained Euler equations allow subtraction of the retained, field-dependent test:

$$
\partial_t C_W[\xi]=-E_{\rm full}(c)[Q_cF_\xi(c)].
$$

The subtraction uses the retained equations; the second displayed expression is not asserted off shell.

Let $\eta$ be a perturbation amplitude, distinct from physical time. Suppose the family is twice amplitude differentiable in a topology controlling the time and spatial jets in the Euler equations, so that

$$
c(\eta)=c_0+\eta c_1+\frac{\eta^2}{2}c_2+o(\eta^2).
$$

At this exact continuum background $E_{\rm full}(c_0)=0$. A background null angular descriptor has $Q_cF_\xi(c_0)=0$. Expanding the product consequently gives

$$
\partial_t C_W[\xi](\eta)
=-\eta^2 E_{\rm full}^{(1)}[c_1]
   [Q_cDF_\xi(c_0)c_1]+o(\eta^2).
$$

Thus the coefficient of $\eta^2$ is $-B_\xi$, while the second amplitude derivative is $-2B_\xi$. Neither $c_2$ nor its temporal rates can cancel this necessary condition. The affine term $\partial_r\xi$ is retained in $F_\xi(c_0)$ and has zero amplitude derivative; it has not been dropped by freezing radial profiles.

An admissible input $c_1$ must be a solution germ of all retained linearized equations, including the temporal constraints, the nondegenerate shift solve, the complementary lapse preservation conditions, the associated lapse-determining equations, and Hamilton evolution. A point satisfying only the linearized secondary equations is insufficient. The supplied time and radial derivatives must belong to that same solution.

Second-order data remain free in the other equations. For example the second-order secondary equation has the general form

$$
DC(c_0)c_2+D^2C(c_0)[c_1,c_1]=0
$$

after expressing the necessary jets consistently. It is not legitimate to set $c_2=0$ and interpret the remaining source as an obstruction. The Noether projection above is useful precisely because the full second-order equations eliminate that unknown correction from the selected preservation condition.

If $B_\xi$ is nonzero on a certified first-order solution, it excludes a sufficiently smooth family of retained solutions through $c_0$ with that specified tangent. It leaves the exact stationary solution intact. If it vanishes, the other second-order equations and higher preservation still need to be solved.

## 3. Canonical expression for an independent comparison

Let $y=(q,p,z)$ and $\mu=\dot z$. Hold the background angular shift coefficient test $k$ fixed and define

$$
C_k=D_zH[k],\qquad
\mathcal F_k=\{C_k,H\}_{q,p}+D_zC_k[\mu].
$$

At the stationary solution with $\mu_0=0$, $\mathcal F_k=0$ and $D\mathcal F_k=0$. The latter uses more than the zero temporal-Hessian row: $C_{k,1}=J_k^*p$ and the proved identities $V''J_k=L_fJ_k=0$ kill the remaining linear terms.

Write $H_j,C_{k,j}$ for homogeneous degree-$j$ Taylor terms, including their factorials. The quadratic preservation coefficient is

$$
\mathcal O_k[y_1,\mu_1]
=\{C_{k,1},H_3\}_{q,p}
+\{C_{k,2},H_2\}_{q,p}
+D_zC_{k,2}[\mu_1].
$$

Since $D\mathcal F_k=0$, adjustable second-order phase-space data and second-order temporal rates drop out. On the same admissible first-order solution, the Noether derivation gives $\mathcal O_k=-B_\xi$. This supplies a prospective independent comparison if cubic Hamiltonian derivatives are available. The full quadratic operator export alone does not supply $H_3$ or the discarded continuum residual, so no such comparison has been run or claimed here.

## 4. Explicit axisymmetric coefficient test

Use the unnormalized $m=0$ bases of the existing checker:

$$
Y_\ell=P_\ell(x),\quad
P_\ell^A=(dP_\ell',0),\quad A_\ell^A=(0,-P_\ell'),
$$
$$
(E_\ell)_{xx}=\frac12P_\ell'',\quad
(E_\ell)_{\phi\phi}=-\frac12d^2P_\ell'',\quad
(B_\ell)_{x\phi}=-\frac12dP_\ell''.
$$

Write the first-order angular area coefficient as $X=X_0+X_1P_1+X_2P_2$ and the tensor amplitudes as $e=C_{\mathrm{polar},2}$, $o=C_{\mathrm{axial},2}$. For axisymmetric perturbations,

$$
(DF_\xi c_1)_\gamma=
\begin{pmatrix}
2w\,\gamma_{1,x\phi}&w\,\gamma_{1,\phi\phi}\\
w\,\gamma_{1,\phi\phi}&0
\end{pmatrix},\qquad
(DF_\xi c_1)_V^\phi=-wV_1^x,\qquad
(DF_\xi c_1)_W^\phi=-wW_1^x.
$$

The $a,b,\lambda$ variations vanish in this axisymmetric slice. The $X$ variation does not vanish: it is the round trace of the displayed tensor variation.

The only discarded coefficients at $L=2$ are:

| Discarded coefficient | Amplitude in $Q_cDF_\xi c_1$ |
|---|---|
| $X_4$ | $-12wo/35$ |
| $C_{\mathrm{polar},4}$ | $2wo/35$ |
| $C_{\mathrm{axial},3}$ | $-2wX_1/15$ |
| $C_{\mathrm{axial},4}$ | $-2w(X_2+e)/35$ |
| $V_{\mathrm{axial},3}$ | $-2wV_{\mathrm{polar},1}/15$ |
| $V_{\mathrm{axial},4}$ | $-6wV_{\mathrm{polar},2}/35$ |
| $W_{\mathrm{axial},3}$ | $-2wW_{\mathrm{polar},1}/15$ |
| $W_{\mathrm{axial},4}$ | $-6wW_{\mathrm{polar},2}/35$ |

These coefficients were derived algebraically and independently reviewed by hand. The prepared checker will compare them against exact positive round tensor/vector projections and reconstruct the complete transformed fields. It has not been executed.

Let $\mathscr E_{A,\ell}^{(1)}(t,r)$ denote the full linearized coefficient Euler residual paired angularly with the indicated unit coefficient test, including its radial density. Then $B_\xi$ is the radial integral of the eight table entries times their corresponding $\mathscr E^{(1)}$. In particular the discarded $W$ residuals are required. They are not set to zero merely because the retained temporal constraints vanish.

Every table entry is proportional to $w(r)$, with no derivative of this test profile. Since $w$ is arbitrary compact support, the sum of these eight products divided by $w$ must vanish as a radial distribution on an integrable first-order tangent. Radial derivatives of the perturbation still enter $\mathscr E^{(1)}$.

Reflection $\phi\mapsto-\phi$ makes $\xi$ odd, while the background and its linearized operator are even. The displayed coefficient map flips reflection parity. Consequently a pure-parity seed has zero $B_\xi$ for symmetry reasons. A meaningful search for a nonzero source in this axisymmetric channel must include mixed even/odd admissible seeds. The vanishing of either pure-parity test alone proves no nonlinear symmetry.

## 5. A geometric control that must not be called an obstruction

Consider the exact coordinate change

$$
\phi\mapsto\phi+\eta f(t,r)x
$$

of the static flat four-metric. Its unrestricted coefficient representation is

$$
a=-1,\quad b=0,\quad\lambda=\Sigma/\rho,\quad
V^\phi=\eta x\,\partial_rf,\quad W^\phi=\eta x\,\partial_tf,
$$
$$
\gamma_{xx}=\Sigma/d+\eta^2\rho d f^2,\qquad
\gamma_{x\phi}=\eta\rho d f,\qquad
\gamma_{\phi\phi}=\rho d .
$$

All first-order coefficients are retained at $N=3$. Subtracting the second derivative of the nonlinear reconstruction along the first-order coefficient direction from the coordinate-orbit second derivative leaves only

$$
\Delta\gamma=2\rho d f^2\,dx^2.
$$

The radial and temporal derivative pieces are accounted for by the $V,W$ connection terms and the quadratic metric reconstruction; they have not been frozen. The round trace of the remaining second derivative is

$$
\Delta X=\rho f^2d^2
=\rho f^2\left(\frac8{15}-\frac{16}{21}P_2+\frac8{35}P_4\right).
$$

Thus this finite coordinate orbit leaves the coefficient band at second order. Nevertheless its first-order metric perturbation solves the full linearized Einstein equations, so $E_{\rm full}^{(1)}=0$ and $B_\xi=0$ exactly. It is a necessary zero control for the dynamical test. This example prevents interpreting finite-orbit leakage alone as failed projected evolution.

For the exact rotation $\phi\mapsto\phi+\eta f(t,r)$, the angular metric is unchanged and only retained degree-one $V,W$ are generated. The differentiated coefficient defect vanishes identically. The prepared minimal checker includes this rotation control, not a claimed numerical check of every symmetry generator.

## 6. Full continuum Euler map and the seed contract

The [prepared evaluator](../../verification/oblate_nonlinear_null_20261127/null_compatibility.py) uses coordinates $(t,r,x,\phi)$ and reconstructs the covariant first-order metric by

$$
h_{tt}=a_1,\quad h_{tr}=b_1,\quad h_{tA}=\gamma_{0,AB}W_1^B,
$$
$$
h_{rr}=\lambda_1,\quad h_{rA}=\gamma_{0,AB}V_1^B,\quad
h_{AB}=X_1q_{AB}+C_{1,AB}.
$$

It forms the full coordinate $\delta\Gamma$, $\delta R_{\mu\nu}$, and $\delta G_{\mu\nu}$ without angular projection. For a coefficient test $u$, set $k_{\mu\nu}=Dg(c_0)[u]_{\mu\nu}$. The normalized pairing is

$$
E_{\rm full}^{(1)}[u]
=-\frac1{4\pi}\int dr\,dx\,d\phi\,
\Sigma\,\delta G^{\mu\nu}[h]\,k_{\mu\nu}.
$$

There is one factor $\Sigma=\sqrt{-g_0}$, with no additional $r^2$. Full index contraction counts off-diagonal components twice. Variations of the density and test map multiply $G[g_0]=0$ and do not contribute. The axisymmetric implementation returns the integrand to be integrated with $dr\,dx/2$.

The shared [seed contract](seed-contract.json) preserves named coefficient columns instead of imposing a new canonical order. Required input includes the actual action/source/domain identifiers, $q_1,p_1,z_1$, consistent derivative data, and evidence of the full retained linearized equations and preservation conditions. In particular,

$$
p_1=\frac12T_0^*G_{\sigma_0}
\left(T_0\dot q_1-\mathcal L_{\beta_1}\sigma_0\right),
\qquad \beta_1=(b_1/\lambda_0,W_1)
$$

must hold in the pinned momentum normalization. The independent derivative orders are:

| Block | Needed highest derivatives |
|---|---|
| $q_1$ | two time derivatives; mixed time/spatial derivatives; two spatial derivatives |
| $a_1$ | two spatial derivatives; no time derivative |
| $b_1,W_1$ | mixed time/spatial derivatives and two spatial derivatives; no second time derivative |

Spatial derivatives include radial derivatives. A raw coordinate implementation can carry formal second jets of every metric component before the temporal-acceleration cancellations; it must not introduce those canceled jets as independent physical inputs. Smooth coefficient functions are the initial preferred input representation.

The supplied retained quadratic matrices alone cannot recover discarded Euler residuals. The preprojection metric maps and the unrestricted linear Euler map must be retained. The evaluator's candidate-source function deliberately does not assert that an arbitrary input is on shell. An off-shell nonzero output is not evidence of a compatibility obstruction.

The minimal interface has twelve spatial and ten temporal coefficient profiles in the $m=0$ block. An axisymmetric seed and the axisymmetric background have no $m\ne0$ residual projections by azimuthal symmetry, but every retained $m=0$ equation and preservation condition must still be checked. Restricting this first probe to $m=0$ does not constitute a test of the other two extra angular null directions.

## 7. Prepared minimal run and what it would settle

The [unexecuted check](../../verification/oblate_nonlinear_null_20261127/check_prepared_maps.py) has four groups:

1. Compare the full coordinate Jacobian pullback with the nonlinear coefficient reconstruction, and independently compare first and second Lie derivatives. Keep arbitrary $f(t,r)$.
2. Verify the eight discarded amplitudes by round harmonic integration and complete reconstruction. Check the exact rotation zero and the scalar $P_4$ coefficient with independent 24- and 48-node quadrature.
3. Check the full linear Einstein map on the retained pure-coordinate perturbation, where the answer must be zero including temporal components.
4. Check a nonzero off-shell conformal control against the independent four-dimensional identity $\delta G_{\mu\nu}[2\psi g_0]=-2\nabla_\mu\nabla_\nu\psi+2g_{0,\mu\nu}\Box\psi$. For $\psi(t,r)$, the separate coefficient test $b=1$ must give the Euler pairing density $-4\rho\,\partial_t\partial_r\psi$ before $dr\,dx/2$ integration. This checks the sign, volume factor, index raising, and both off-diagonal entries, as well as catching an implementation that returns zero for every input.

Proposed command, from paper1:

    uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 python verification/oblate_nonlinear_null_20261127/check_prepared_maps.py

This is one local CPU Python process. The proposed initial budget is five minutes and 2 GiB RAM; these are planning estimates, not measured use or a hard limit implemented by the script. No GPU, remote compute, or external API is needed. A future execution should be supervised and stopped if it exceeds that budget. The [run proposal](run-proposal.json) records outputs and stop conditions.

Successful execution would validate the maps and independent controls. It would not evaluate a non-gauge on-shell $B_\xi$, classify the three extra null directions, or establish nonlinear closure. That dynamical evaluation requires the separately prepared operator/domain work to supply an admissible mixed-parity seed. The cubic Hamiltonian comparison is a further independent option, not an output of this minimal run.

Only static source/syntax, file-link, and provenance checks are permitted in this preparation stage. The [preparation record](preparation-provenance.json) identifies the reviewed sources and distinguishes those static checks from scientific execution. No PASS log or measured map results are included. The next decision is whether to authorize the stated local map/control run; the full nonlinear release gate remains open.
