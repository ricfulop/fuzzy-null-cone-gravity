# Contracted gravitational four-point loop: implementation handoff

## Status, provenance, and next deliverable

This is an implementation recommendation, not an evaluated gravitational loop.
The new analytic contractions, adjoint reductions, and heat-kernel formulas
below are **numerically untested**. They were derived from the
[sixth-order prescription](reduced-sixth-order-and-loop-kernel.md), its
[independent audit](sixth-order-independent-audit.md), and source inspected
at `087ab31`; this handoff was written in the primary worktree at `5e04807`.
No experiment, numerical verification, or performance profile was run here.
The [quartic/two-point report](complete-quartic-and-two-point.md) owns the
classical reconstruction, the retained coarse-quadrature failure, and the
subsequent successful two-point calculation checkpointed at `6bd4090`.
That two-point result does not evaluate the higher gravitational vertices here.

The next deliverable is an actual spatial gravitational $H_5,H_6$ oracle,
including all momentum sectors, followed by one contracted four-point number.
Do not replace it with another oscillator-only milestone. Work remains within
the declared reduced canonical candidate and its stated inverse/Gram hypotheses.

## 1. Fixed first target and conventions

Use the ten profiles in `../../verification/quantum_bridge/reduced_gaunt_exact.py`,
with packet width one, $16\pi G=1$, $\hbar=1$, bare midpoint/Weyl ordering,
canonical phase-space measure, and no added order-$\hbar$ counterterm.
Retain interleaved order $(Q_0,P_0,\ldots,Q_9,P_9)$ and all ten internal modes.
The normalization in `../../verification/quantum_bridge/complete_quartic.py` gives
$$
q(x)=\sum_i\sqrt{2/g_i}\,Q_i e_i(x),\qquad
S(x)=\sum_i P_i e_i(x)/\sqrt{2g_i},\qquad
\omega_{0:4}^2=11/2,\quad\omega_{5:9}^2=13/2,
$$
where $g_i$ includes the full spatial integral, including its $\pi^{3/2}$ factor.
Here $q(x)$ is the spatial tensor called $Q$ in the reduction formulas below.

Choose four external axial electric modes: **zero-based normalized mode 1**,
whose raw STF label is $\operatorname{diag}(-1/2,-1/2,1)$, at
$$
(\Omega_1,\Omega_2,\Omega_3,\Omega_4)=(0.2,-0.45,0.8,-0.55),\qquad\sum_r\Omega_r=0.
$$
For the physical source convention $e^{-i\Omega_r\tau}$, the free response is
$$
z_r=G(-\Omega_r)e_{Q_1}
=\frac{e_{Q_1}+\Omega_r e_{P_1}}{\Omega_r^2+11/2}.
$$
Compute the phase-space 1PI coefficient contracted with the amputated directions
$v_r=e_{Q_1}+\Omega_r e_{P_1}$; record the four scalar external denominators
separately. Setting the external momentum components to zero changes the target.
In particular retain $H_5:q^5,q^3p^2,qp^4$ and
$H_6:q^6,q^4p^2,q^2p^4,p^6$.

Axial $SO(2)$ may reduce equivalent internal magnetic-number copies. First derive
its action in the normalized retained basis and verify the transformed Hessians
and contracted traces. Only then use proved equality and multiplicities; do not
drop nonaxial internal modes or assume their contributions vanish. The first
reference implementation should retain all twenty phase-space components.

## 2. Avoid constructing dense fifth- and sixth-order tensors

For $B_S=D_{v_r:r\in S}\nabla^2H_{|S|+2}(0)$, reuse the complete $H_3,H_4$
polynomials and the existing `integrate_vertex` with only $|S|\le2$ supplied.
It gives the quartic bubble, quartic/cubic/cubic triangle, and cubic box.
Add the sixth-order tadpole and fifth/cubic pair below; do not add them twice.

The sixth-order tadpole needs only twenty diagonal contractions:
$$
T_6=\sum_{i=0}^9\left[
\frac{D^6H_6[v_1,v_2,v_3,v_4,e_{Q_i},e_{Q_i}]}{4\omega_i}
+\frac{\omega_iD^6H_6[v_1,v_2,v_3,v_4,e_{P_i},e_{P_i}]}4\right].
$$
Symmetry cancels the odd off-diagonal propagator trace. No off-diagonal $B_{1234}$
entry is needed for this diagram, although momentum-rich diagonal entries are.

All gravitational cubic blocks here are pure configuration. Write
$C^{(r)}_{ij}=D^3H_3[v_r,e_{Q_i},e_{Q_j}]$ and embed $C^{(r)}$ in its QQ block.
For $a=\omega_i,b=\omega_j,\Omega=\Omega_r$, define
$$
I_{ij}(\Omega)=\frac{a+b}{2ab[(a+b)^2+\Omega^2]},\qquad
K_r=\int\frac{d\nu}{2\pi}G(\nu)C^{(r)}G(\nu+\Omega_r).
$$
The required moments of the same two denominators are
$\int\nu=-\Omega a I_{ij}/(a+b)$,
$\int(\nu+\Omega)=\Omega b I_{ij}/(a+b)$, and
$\int\nu(\nu+\Omega)=abI_{ij}$, with measure $d\nu/(2\pi)$ understood.
Consequently the symmetric matrix $K_r$ has blocks
$$
\begin{aligned}
(K_r)_{Q_iQ_j}&=C^{(r)}_{ij}I_{ij},&
(K_r)_{Q_iP_j}&=C^{(r)}_{ij}\frac{\Omega b}{a+b}I_{ij},\\
(K_r)_{P_iQ_j}&=C^{(r)}_{ij}\frac{\Omega a}{a+b}I_{ij},&
(K_r)_{P_iP_j}&=-C^{(r)}_{ij}abI_{ij}.
\end{aligned}
$$
The eight ordered $(1,3)$ partitions give
$$
P_{53}=-\frac12\sum_{r=1}^4 K_r:B_{\{1,2,3,4\}\setminus\{r\}},\qquad
\Gamma^{(1)}[v_1,v_2,v_3,v_4]=\Gamma^{(1)}_{H_3,H_4}+P_{53}+T_6.
$$
For this real target, $K_r=\sum_{\alpha=1}^{20}\lambda_\alpha w_\alpha w_\alpha^T$.
Each pair contraction becomes at most twenty derivatives
$\lambda_\alpha D^5H_5[v_{s\ne r},w_\alpha,w_\alpha]$.
Thus at most 80 fifth-order and 20 sixth-order directional contractions suffice,
before any proved symmetry reduction. This is an operation count, not a timing.

## 3. Spatial oracle: local geometry, York response, and Gram correction

Extend `../../verification/quantum_bridge/reduced_gaunt_vertices.py` to obtain
pointwise $R_2,R_3,R_4$, integrated $R_5,R_6$, and $b_1,b_2,b_3,\mathcal A_1,\mathcal A_2$.
The integrated Gamma--Gamma density can supply the last two curvature integrals;
it cannot replace pointwise curvature in the scalar sources. Christoffels through
order five suffice for the integrated sixth-order Gamma--Gamma coefficient.
Build metric/tangent/divergence series by truncated algebra rather than fitting
finite-amplitude values. Preserve both Fourier parities in Gaussian operations.

Compute $W_1=\mathcal A_0^{-1}b_1$ and
$W_2=\mathcal A_0^{-1}(b_2-\mathcal A_1W_1)$, with the report's $\mathcal N_3,\mathcal N_4$.
The latter still requires $b_3$, despite cancellation of $W_3$.
For $a=\mathsf G_0^{-1}p$ and raw retained basis $e_\alpha$, obtain
$g=\mathsf G_2a$ without forming the whole matrix:
$$
g_\alpha=-\frac1{12}\int\operatorname{tr}([Q,e_\alpha][Q,S])
+2\langle b_1(Q,e_\alpha),W_1(Q,S)\rangle.
$$
This follows by polarizing $\mathcal M_2-\mathcal N_2$ and using inverse symmetry;
it needs ten pairings but no extra York inverse for each basis tensor. Then
$$
I_5=\mathcal N_3[S,S],\qquad
I_6=-\mathcal M_4[S,S]+\mathcal N_4[S,S]+g^T\mathsf G_0^{-1}g.
$$
Use the same $g$ in $Z_2$'s canonical correction $-\sum_\alpha e_\alpha
(\mathsf G_0^{-1}g)^\alpha$. Construct $s_3,s_4$ from the report, including
$\operatorname{tr}(U_1^2)$, not $\operatorname{tr}(U_1^TU_1)$.

## 4. Adjoint reduction removes the production fourth-order scalar solve

Set $\mathcal C=(-\Delta)^{-1}/8$, $B_j=R_j+7s_j$, and
$$
u=\mathcal C(s_2-R_2),\quad v=\mathcal CB_2,\quad
F_3=s_3-R_3+8D_1u,\quad K_3=B_3+8D_1v.
$$
Substituting $u_4$ and moving $\mathcal C,D_1$ between paired factors gives
$$
H_5=I_5-\int R_5-\langle B_3,u\rangle
-\langle v,s_3-R_3\rangle-8\langle v,D_1u\rangle,
$$
$$
\begin{aligned}
H_6={}&I_6-\int R_6-\langle B_4,u\rangle-\langle v,s_4-R_4\rangle
-\langle K_3,\mathcal CF_3\rangle\\
&-8\langle v,D_2u\rangle+\int vB_2u+28\int s_2u^2.
\end{aligned}
$$
Only two quadratic scalar sources and one cubic source need inversion.
Use $\langle v,D_1u\rangle=\int\partial_iv\,Q^{ij}\partial_ju$ and
$\langle v,D_2u\rangle=-\tfrac12\int\partial_iv\,(Q^2)^{ij}\partial_ju$.
The pure-momentum reduction gives $49+28=77$; the configuration reduction
reproduces the report's Schur-complement expressions. These are algebraic
cross-checks; evaluate the original $u_4$ formula independently on a mixed case.

## 5. Exact Gaussian first inverses and numerical nested inverses

For $f_a=e^{-ar^2}$ with positive Gaussian widths, direct heat-kernel integration gives
$$
\mathcal Cf_a=\frac{\sqrt\pi}{32a^{3/2}}
\frac{\operatorname{erf}(\sqrt a\,r)}r
=\frac1{16a}\int_0^1e^{-at^2r^2}\,dt.
$$
Expand polynomial-Gaussian sources in derivatives of $f_a$; derivatives commute
with the inverse and can use the existing Gaussian-derivative routines.
For a constant vector $c$, the analogous York formula is
$$
[\mathcal A_0^{-1}(cf_a)]_i=-\frac1{16a}\int_0^1
[(7+t^2)\delta_{ij}+2at^2(1-t^2)x_ix_j]e^{-at^2r^2}\,dt\,c_j.
$$
Derivation: integrate $e^{-sk^2}[-\delta_{ij}+sk_ik_j/4]$ over $s\ge0$,
apply the Gaussian heat evolution, and set $t=(1+4as)^{-1/2}$.
This represents the complete multiplier $-\delta_{ij}/k^2+k_ik_j/(4k^4)$.
It avoids imposing decay on a separately divergent biharmonic potential.
Two independent Gaussian normalization targets are
$$
\langle f_a,\mathcal Cf_b\rangle=\frac{\pi^{3/2}}{16ab\sqrt{a+b}},\qquad
\int f_a(\mathcal Cf_b)(\mathcal Cf_c)=
\frac{\pi^{3/2}}{256\sqrt a\,(bc)^{3/2}}
\arctan\sqrt{\frac{bc}{a(a+b+c)}}.
$$
Products/derivatives remain Gaussian under these auxiliary integrals. The displayed
sixth-order functionals have at most three inverse factors per expanded connected
spatial term. A wholly auxiliary-integral backend is possible, but joint endpoint
singularities and polynomial growth require care; do not presume it is faster.

Recommended first backend: exact local Gaussian algebra and first inverses,
then all-$m$ radial harmonic coefficients for nested sources. With $A=(-\Delta)^{-1}$,
$$
(Af)_{\ell m}(r)=\frac1{2\ell+1}\left[
r^{-\ell-1}\int_0^r s^{\ell+2}f_{\ell m}(s)ds
+r^\ell\int_r^\infty s^{1-\ell}f_{\ell m}(s)ds\right],\qquad
\mathcal A_0^{-1}b=-Ab-\tfrac14\nabla A^2(\nabla\cdot b).
$$
Generalize the Gaussian radial Green routine to panelwise high-order cumulative
integration with explicit exterior tails. The second scalar potential in the
York formula can contain a bounded nondecaying dipole; only its gradient must
decay. Do not apply a decaying harmonic Robin condition to every intermediate.

The adjoint route needs scalar inverse total $j\le6$ and vector inverse $j\le6$;
Cartesian vector components need orbital $\ell\le7$. Keep larger local support
for $b_3$ and fourth-order sources. Oversample angular products and verify omitted
harmonics before truncation; never project back to the physical free $j=2$ cutoff.
Crucially, $s_4$ contains $(L_0W_1)^2$: a possible $W_1\sim r^{-1}$ tail gives
$s_4\sim r^{-4}$ and an H6 radial remainder of order $R^{-2}$ when paired with
$u\sim r^{-1}$. General mixed tails are not Gaussian. Integrate their multipoles
analytically or use a controlled compactified exterior domain.

## 6. Directional jets, acceptance, and observable limits

For each H5 derivative use $z=\sum_{r=1}^3\eta_rv_r+tw$, with $\eta_r^2=0,t^3=0$.
Twice the coefficient of $\eta_1\eta_2\eta_3t^2$ is the desired derivative:
24 jet slots. H6 uses four external labels and 48 slots. Inverses act coefficientwise.
Do not differentiate a noisy nonlinear constraint solution six times; its error
is amplified by $h^{-6}$. An independent step-free homogeneous polarization is
$$
D^mH_m[w_1,\ldots,w_m]=2^{-m}\sum_{\epsilon_r=\pm1}
\left(\prod_r\epsilon_r\right)H_m\left(\sum_r\epsilon_rw_r\right).
$$
Check cancellation sensitivity with increased arithmetic precision. Use the
existing independent Decimal finite-metric evaluator for local geometry checks.

Required acceptance record before reporting the target:

1. Recover gravitational H3/H4 on held-out mixed electric/magnetic directions.
2. Check inverse residuals, adjointness, $\int s_3=0$, $\int s_4=-a^T\mathsf G_2a$,
   canonical Gram correction, mixed jet/polarization agreement, and momentum reversal.
3. Check the 77 identity, configuration Schur complement, and original-versus-adjoint
   H6 on mixed data. The static square root tests only Q legs at zero momentum/frequency.
4. Independently refine radial panels, angular integration, exterior tails, and loop
   frequency quadrature; use three spatial resolutions and retain failed comparisons.
   Preregister six significant figures as an initial goal, not an achieved accuracy.
5. Export each of the five gravitational diagram families, their sum, leg conventions,
   normalization, and numerical error estimates; check permutations and rotations.

Store inverse fields as radial-harmonic coefficients and stream angular products.
The directional counts above and $O(N_\nu\,75\,20^3)$ dense loop contraction bound
do not predict wall time. Profile one mixed H6 contraction before assigning a budget.

This target is a contracted **phase-space 1PI** contribution. Its physical connected
Q-correlator contribution has the external propagators and the usual Euclidean
minus vertex sign. The full connected one-loop four-point function also needs
corrected three-point exchanges, internal/external self-energy insertions, and a
consistent vacuum/source prescription. Linear free-leg contraction alone is not
the full configuration-only 1PI action after eliminating momenta.
A ten-mode off-shell correlation is not asymptotic four-graviton scattering:
radial-basis completion, physical angular-cutoff limits, asymptotic packets,
Lorentzian continuation, external-state normalization, and a renormalization
prescription remain separate requirements. Auxiliary quadrature refinement is
not any of those physical limits.
