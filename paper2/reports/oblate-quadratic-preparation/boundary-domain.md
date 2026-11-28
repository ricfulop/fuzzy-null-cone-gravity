# Boundary/domain preparation for oblate quadratic dynamics

**Status.** Specification only. No checker, experiment, numerical sweep, remote job, or GPU job is launched here.

This note prepares the radial variational and canonical-domain side of the proposed
oblate quadratic operator export. It is downstream of the Paper I coefficient-Galerkin
Einstein-Hilbert bulk action at commit `157ec0f04866ea36509b325ec65f0a68890f9b17`
and the exact all-mode oblate rank and lapse-preservation archives:

- `/Users/ricfulop/.local/share/openresearch/files/scratch-db17dab4/paper1/oblate-lapse-consistency-2026-11-27/oblate-lapse-consistency.md`
- `/Users/ricfulop/.local/share/openresearch/files/scratch-db17dab4/paper1/oblate-lapse-consistency-2026-11-27/lapse-summary.json`
- `/Users/ricfulop/.local/share/openresearch/files/scratch-db17dab4/paper1/oblate-lapse-consistency-2026-11-27/lapse-leading.json`
- `/Users/ricfulop/.local/share/openresearch/files/scratch-db17dab4/paper1/oblate-lapse-consistency-2026-11-27/sources/sections/nonlinear_action.tex`

The Paper I convention to preserve is

$$
\frac{16\pi G}{4\pi} I,\qquad
q=(\lambda,X,V,C),\qquad
z=(a,b,W),\qquad
\beta^r=b/\lambda,\quad
\beta^A=W-bV/\lambda,\quad
\alpha^2=-a+b^2/\lambda,
$$

with canonical momenta \(p\) conjugate to \(q\). The exact static oblate background is

$$
\rho=r^2+\epsilon,\qquad
\Sigma=r^2+\epsilon x^2,\qquad
\sigma=\frac{\Sigma}{\rho}dr^2+\frac{\Sigma}{1-x^2}dx^2
+\rho(1-x^2)d\phi^2,\qquad
a=-1,\quad b=V=W=0,\quad p=0 .
$$

The immediate target reducer takes a stationary quadratic Hamiltonian

$$
H_2=\frac12 w^T K w+z^T E w+\frac12z^T D z,\qquad w=(q,p),
$$

and stabilizes the finite linear constraints on the canonical phase space. Its inputs
must therefore be the complete \(K,E,D\) blocks or a common radial operator export
from which those blocks, adjoints, and structural null vectors can be assembled; a
rank of \(D\), a reduced bracket, or a frozen-radius angular matrix is not enough.
See `paper2/reports/constraint-stabilization-handoff.md` and
`paper2/verification/constraint_handoff/reduce_quadratic.py`.

The shared canonical assembly target is

$$
S=T^*GT,\qquad
U_s=\mathcal L_{(b/\lambda,W)}\sigma,\qquad
\mathcal J s=S^{-1}T^*GU_s,
$$

and, in canonical variables,

$$
H_2=\frac12\langle q,V''q\rangle
+\langle p,S^{-1}p\rangle
+\langle p,\mathcal J s\rangle
+\langle\delta a,Lq\rangle
+\frac12\langle s,As\rangle,
$$

with

$$
K_{\rm full}=\operatorname{diag}(V'',2S^{-1}),\qquad
E=\begin{pmatrix}L&0\\0&\mathcal J^*\end{pmatrix},\qquad
D=\operatorname{diag}(0,A),\qquad
A=-\frac12(N_\perp U_s)^*G(N_\perp U_s).
$$

Here \(s=(b,W)\). The \(K=2(T^*GT)^{-1}\) used in the Paper I lapse
preservation archive is the kinetic momentum block \(2S^{-1}\), not the full
\(K_{\rm full}\), whose \(qq\) block is the unit-lapse potential Hessian \(V''\).

## 1. Three notions of radial validity

Keep these three layers separate in the export and in every interpretation.

### 1.1 Compact-test bulk identities

The reconstructed Paper I action defines bulk coefficient equations by ordinary
four-dimensional Einstein-Hilbert variation after reconstructing the retained
metric coefficients into \(g_N\). Its current bulk definition is a compact-support
variational statement; it does not choose endpoint data for a compact radial shell
or pull back a null-boundary and joint action. The oblate rank and lapse derivations
initially use smooth compactly supported radial tests, so every radial integration
by parts has zero endpoint pairing.

Consequence: the compact-test identities prove local bulk algebra, not a physical
radial spectrum. They are the right starting point for deriving the differential
expressions and their formal adjoints.

### 1.2 Continuous \(H^1_0\) lapse-bracket extension

The corrected lapse-preservation archive promotes the leading complementary
lapse-tertiary blocks to exact small-deformation weak operators on the fixed
compact shell \(I=[r_-,r_+]\), \(r_->0\). In each complementary channel the leading
term has the form

$$
G(Z_f,Z_g)
=c_{\ell m}\epsilon^2
\int_I r^{-6}(aF+brF')^*(aG+brG')\,dr+O(\epsilon^3),
$$

with nonzero \(c_{\ell m}\) and \((a,b)\in\{(1,-1),(0,1),(1,1)\}\). After signing by
\(\operatorname{sgn}c_{\ell m}\), weighted Poincare control of
\(F'+\eta F/r\) on \(H^1_0(I)\) makes the small-\(\epsilon\) operator
coercive \(H^1_0(I)\to H^{-1}(I)\) on the complementary channels. Its kernel on the
retained lapse space consists of the exact scalar-gradient families identified in
the same archive.

Consequence: the reduced lapse bracket has a legitimate homogeneous weak radial
realization. This is stronger than a compact-support identity and should be reused
as a control in Galerkin checks.

### 1.3 Common canonical generators and physical evolution

The \(H^1_0\) lapse-bracket result does not give a common domain for:

- the original, uneliminated shift block \(A=-\frac12G(n_s,n_t)\),
- the shift tangential map \(\mathcal J\), which contains radial derivatives of the shift
  profiles,
- the full unit-lapse curvature Hessian \(V''\),
- the canonical kinetic map from coefficient velocities to momenta,
- the physical evolution generator after all constraints have been stabilized, or
- boundary traces that are nonzero at \(r_\pm\).

These objects have different differential orders before elimination. A common
domain must make the Green identities used by \(\mathcal J^*\), \(V''\), and
\(L_fK_{\rm eff}L_g^*\) true up to known endpoint bilinears, and then either kill
those bilinears by homogeneous boundary data or cancel them with specified boundary
terms.

Consequence: the homogeneous \(H^1_0\) shell is a good provisional diagnostic, but
nonzero boundary traces and a physical spectrum still require explicit endpoint
data and, if the endpoints are meant to model causal-diamond boundaries, the
appropriate null-boundary/joint variational principle. The shell diagnostic must
not be described as the null-boundary theory.

## 2. Boundary choices that require a decision

### Choice A: homogeneous bulk diagnostic

Use a compact shell with homogeneous radial trial functions:

$$
I=[r_-,r_+],\quad r_->0,\qquad
\mathcal V_h\subset H^1_0(I)
$$

for every retained coefficient profile. Construct the weak forms on
\(\mathcal V_h\) and let radial endpoint traces vanish. This defines a finite
diagnostic using homogeneous tests, not a proved common invariant domain of the
continuum generators. The lapse block is in the archive's form domain; derivative
maps between the different coefficient spaces still require the checks below.

This is **provisional and diagnostic**. It can test algebraic assembly of
\(K,E,D\), Schur elimination, temporal basis invariance, structural nulls whose
radial profiles lie in \(H^1_0(I)\), and convergence under increasing the radial
basis. It deliberately suppresses boundary charges and homogeneous endpoint
solutions.

### Choice B: homogeneous mixed diagnostic

Use homogeneous Dirichlet data for fields whose weak blocks are first order, but
choose Neumann or mixed data for variables whose canonical generator is naturally
second order after the full \(V''\) export is known. This can reduce endpoint
overconstraint if the physical \(K\) block contains Sturm-Liouville suboperators
whose natural self-adjoint realizations are not Dirichlet.

This cannot be chosen abstractly. It needs the full \(V''\) Green identity and a
decision about which boundary traces are fixed, which are free multipliers, and
which are conjugate charges.

### Choice C: physical endpoint data

Allow nonzero endpoint traces and add or pull back the boundary and joint terms
that make the finite-shell action differentiable for those traces. For a
causal-diamond interpretation the data must also specify which radial endpoints
represent which null or joint hypersurfaces and how their null generators are
normalized.

Physical boundary charges and causal-diamond gluing require the corresponding
variational principle. A homogeneous shell can instead describe a chosen cavity
model if its boundary action and domain are justified; it does not thereby describe
a causal diamond. Neither physical interpretation is supplied merely by running a
matrix reducer.

**Recommended immediate choice:** implement Choice A as a homogeneous weak/Galerkin
diagnostic, label its output "bulk compact-shell diagnostic", and reserve a later
decision on endpoint data before assigning physical boundary meaning to the spectrum.
This is a proposal, not a selected boundary model or an authorization to run it.

## 3. Green identities and adjoints to export

Write retained radial profiles as \(u_i(r)\) multiplying real harmonics already
used by `all_modes.py`. For a block differential expression

$$
(\mathcal L u)_\alpha
=A_{\alpha i}(r)u_i''
+B_{\alpha i}(r)u_i'
+C_{\alpha i}(r)u_i,
$$

the weak form must be assembled by integration by parts from the action Hessian:

$$
\mathsf L[v,u]=
\int_I
\left[
v_i'\,P_{ij}(r)\,u_j'
+v_i\,Q_{ij}(r)\,u_j'
+v_i'\,R_{ij}(r)\,u_j
+v_i\,S_{ij}(r)\,u_j
\right]\,dr
+\mathsf B[v,u]\big|_{r_-}^{r_+}.
$$

For the potential weak form, the coefficients correspond to
\(P=A, R=B, Q=B^\dagger, S=C\) in the notation of
[the potential report](potential-hessian.md). These are complete derivative-slot
coefficients; do not add another symmetrization or factor of one half. Export the
strong bulk expression, weak bulk integrand, and Green flux separately. The flux
belongs to the boundary pairing, not to the strong bulk row.

The export should provide:

- the coefficient matrices \(P,Q,R,S\) for each angular \(m\) block;
- the endpoint bilinear \(\mathsf B\) before imposing any boundary condition;
- the formal-adjoint identity \(\mathsf L[v,u]-\mathsf L[u,v]^*\);
- a flag saying whether symmetry requires \(Q=R^*\) after endpoint terms or only
  weak equality on constraints;
- the homogeneous domain on which \(\mathsf B=0\).

The maps that depend on \(r\), especially

$$
T\mathcal J s=\Pi\mathcal L_{v_s}\sigma,\qquad
n_s=N_\perp\mathcal L_{v_s}\sigma,\qquad
\mathcal J_R A_R^{-1}\mathcal J_R^*,
$$

need the same treatment. Even though the proven \(n_s\) and shift Hessian \(A\) do
not contain radial derivatives of shift profiles after exact tangent cancellation,
\(\mathcal J\) still contains those derivatives. Its adjoint cannot be represented
by a transpose until both the coefficient pairing and Green endpoint pairing have
been fixed; neither is the indefinite kinetic Gram by itself.
Explicitly, at \(V=0\),

$$
\begin{aligned}
U_{rr}&=2b'-b\,\partial_r\log\lambda+W^A\partial_A\lambda,\\
U_{rA}&=\partial_A b-b\,\partial_A\log\lambda+\gamma_{AB}(W^B)',\\
U_{AB}&=(b/\lambda)\partial_r\gamma_{AB}+\mathcal L_W\gamma_{AB}.
\end{aligned}
$$

The derivative injections \(b'\mapsto\delta\lambda=2b'\) and
\(W'\mapsto\delta V=W'\) are tangent, so \(n_s\) is algebraic in \(s\); the
tangential map is not. In weak notation

$$
\mathcal J s=J_0s+J_1s',\qquad
(\mathcal J^\dagger p)_{\rm bulk}=J_0^*p-\partial_r(J_1^*p),
$$
$$
\langle p,\mathcal J\delta s\rangle
-\langle(\mathcal J^\dagger p)_{\rm bulk},\delta s\rangle
=[p^*J_1\delta s]_{r_-}^{r_+}.
$$

A formal expression \(A_R^{-1}\mathcal J^*p\) is therefore not automatically an
allowed shift profile. A shift-space inverse \(A_R^{-1}\) may demand exactly the
regularity produced by \(\mathcal J^*p\), or it may not; the export has to prove this on
the chosen radial domain. Eliminating shifts can also introduce \(p'\)-dependence
into the reduced kinetic form \(\mathcal J_R A_R^{-1}\mathcal J_R^*\), so the momentum basis and its
endpoint traces are part of the domain problem, not a post-processing detail.

For the lapse chain, keep the reduced bracket in weak form:

$$
L_fK_{\rm eff}L_g^*
=-\frac12G(Z_f,Z_g).
$$

The leading \(H^1_0\) form controls \(F'\) and \(F/r\) on the compact shell, but the
original \(L_f\), \(K_{\rm eff}\), and \(L_g^*\) must still be assigned domains if a
later export represents them separately.

## 4. Canonical momentum conversion

The reducer assumes

$$
\Omega=\sum_i dp_i\wedge dq_i,\qquad
J=\begin{pmatrix}0&I\\-I&0\end{pmatrix}
$$

in the exported finite vector \(w=(q,p)\). Coefficient velocities do not provide
this automatically because the kinetic term uses the integrated DeWitt Gram and
field-dependent \(r\)-dependent maps.

The projection order must first be named. Canonical radial projection **after**
the angular Legendre transform pulls back \(2S(r)^{-1}\). Restricting the radial
Lagrangian **before** its Legendre transform instead inverts
\(W_M=\frac12\int F^T S F\,dr\). These finite matrices are generally different;
neither order is selected by this preparation.

For canonical projection, write \(u=FQ\), \(p=\Psi c\), and
\(C=\int\Psi^TF\,dr\). The canonical one-form is
\(c^T C\,dQ\), hence the momentum conjugate to \(Q\) is
\(P=C^Tc\). All Hessian blocks must be transformed consistently. In particular,
if both fields use the same real basis \(F\), then
\(M=\int F^TF\,dr\), \(P=Mc\), and \(c=M^{-1}P\).
Equivalently choose the dual momentum basis \(\Psi=FM^{-1}\), so its coefficients
are already \(P\). This pairing mass matrix is distinct from the indefinite
DeWitt kinetic Gram.

Only the Lagrangian-first radial approximation uses
\(P=W_M\dot Q-d_M\), with \(W_M=\frac12\int F^TSF\,dr\) and its derived shift
coupling. The factor one half follows from the normalized kinetic Lagrangian.
The interface must preserve the radial order, pairing matrix, dual basis, and
kinetic Gram as separate objects.

## 5. Provisional radial Galerkin diagnostic

The smallest concrete diagnostic can use exact polynomial "bubble" functions on
the normalized interval \(s=(r-r_-)/(r_+-r_-)\):

$$
\phi_n(r)=s(1-s)\,P_n^{(1,1)}(2s-1),\qquad n=0,\ldots,n_r-1.
$$

Properties:

- \(\phi_n(r_\pm)=0\), so first-order endpoint bilinears vanish;
- \(\operatorname{span}\{\phi_n\}\subset H^1_0(I)\), matching the corrected lapse
  weak-domain theorem;
- derivatives are polynomials in \(s\), so integrals of rational functions of
  \(r,\rho,\Sigma\) against \(\phi_m,\phi_n,\phi'_m,\phi'_n\) can be checked by two
  independent quadratures or by exact symbolic integration at rational
  \(r_\pm,\epsilon\) when feasible;
- the basis does not enforce a false null at a single endpoint by construction
  beyond the chosen homogeneous trace.

Candidate first diagnostic:

$$
r_-=1,\qquad r_+=2,\qquad 0<|\epsilon|\ll1,\qquad N=3,\qquad n_r=3.
$$

Using the same \(\phi_n\) for all coefficient profiles is one inexpensive
diagnostic candidate, not a compatible discretization theorem. A derivative-
compatible alternative should be designed before interpreting any lost kernel as
physical. Neither choice has been selected or launched.

Using one equal radial space can destroy exact tangency and kernel relations: for
example \(s\in\mathcal V_h\) need not imply \(s'\in\mathcal V_h\), and
\(\lambda f\) or \(\partial_r f\) need not lie in the same finite bubble span as
\(f\). The first diagnostic should therefore measure the residual of every known
tangent injection and null template after projection, rather than require those
identities to be exact at finite \(n_r\).

The finite weak matrices should be assembled as bilinear forms, not by collocating
the strong differential expression:

$$
K_{\alpha i,\beta j}
=\int_I\left[
\phi_\alpha' P_{ij}\phi_\beta'
+\phi_\alpha Q_{ij}\phi_\beta'
+\phi_\alpha' R_{ij}\phi_\beta
+\phi_\alpha S_{ij}\phi_\beta
\right]\,dr .
$$

The same rule applies to \(E\) and \(D\). If a temporal variable is a pure
multiplier with no radial derivative in \(D\), its block may be an \(L^2\) mass
matrix. If a shift block couples through \(\mathcal J^*\), any derivative acting on the
shift profile must be moved according to the exported weak form.

## 6. Structural nulls under homogeneous traces

Known structural nulls must be filtered by the radial trace choice.

- The all-mode shift kernel includes \(b=\lambda f\) for \(f\in V_{L-2}\) with
  arbitrary radial profile, and angular \(W\) with \(W(x)\in V_{L-1}\).
- The corrected low lapse families use arbitrary profiles \(F_0(r)\) at \(N=3\)
  and \(F_z,F_\pm\) at \(N=4\), with actual radial shift coefficient
  \(b=\partial_r f\).

On the \(H^1_0\) Galerkin space, a lapse profile \(F\in H^1_0(I)\) has zero trace,
but \(b=\partial_rF\) generally does **not** vanish at the endpoints. If the same
Dirichlet basis is imposed on shift \(b\), these scalar-gradient nulls may be
artificially broken. Three diagnostic options should therefore be tested as
separate boundary choices, not mixed in one matrix:

1. **Strict homogeneous all-fields test.** Enforce \(F,b,W,q,p\in H^1_0\). Expect
   some gradient-generated null identities to be absent; keep only the identities
   whose constructed \(b,W\) lie in the chosen subspace.
2. **Derivative-compatible gradient test.** Use \(F\in H^1_0\) for lapse and a
   shift basis that includes \(\partial_rF\), without requiring \(b(r_\pm)=0\).
   Record the endpoint bilinear that must be cancelled or fixed.
3. **Pure compact-support convergence test.** Use high-order bubbles with profiles
   increasingly localized away from the endpoints to check the compact-test
   identities without asserting endpoint completeness.

The first test exposes a simple finite approximation, with known losses of
structural nulls. The second addresses the endpoint part of the gradient identity;
exact preservation additionally requires closure of every reconstructed product
and derivative in the relevant finite spaces. The third isolates local bulk
algebra from a physical endpoint interpretation. None is automatically a common
domain for the full continuum Hamiltonian.

## 7. Planned verification checklist

Do not execute this checklist until the canonical operator export exists.

1. **Manifest checks.**
   Verify source commit, normalized action, field ordering, real harmonic
   convention, background, \(r_\pm\), \(\epsilon\), \(N\), \(n_r\), and radial basis.

2. **Weak Green identities.**
   For every exported block, assemble the endpoint bilinear \(\mathsf B\). Check
   \(\mathsf B=0\) under the declared homogeneous traces, and check the remaining
   bilinear for any derivative-compatible test.

3. **Tangency and kernel projection defects.**
   Apply the finite basis to \(b'\mapsto\delta\lambda\), \(W'\mapsto\delta V\),
   \(b=\lambda f\), and \(b=\partial_r f\). Record exact residuals caused by
   the radial truncation separately from analytic failures of the continuum
   identities.

4. **Shift-inverse domain check.**
   Given each momentum basis vector \(p_\alpha\), form \(\mathcal J^*p_\alpha\) weakly.
   Test whether it lies in the range of the declared shift block \(A_R\), whether
   \(A_R^{-1}\mathcal J^*p_\alpha\) satisfies the allowed shift traces, and whether the
   reduced kinetic form is symmetric after endpoint terms.

5. **Canonical pairing.**
   Compute the discrete symplectic matrix from
   \(\int p\,\dot q\,dr\). Apply the stated mass/canonical transform and verify
   that the reducer input has \(J=(0,I;-I,0)\).

6. **Temporal block recovery.**
   Recover the archived small-deformation shift-Hessian rank and the equality of
   the normal-map and \(D\)-image ranks in the same angular sectors.

7. **Corrected lapse weak forms.**
   Project the complementary lapse blocks onto the Galerkin basis and verify that
   their leading coefficients match the archived
   \(r^{-6}(aF+brF')^*(aG+brG')\) forms.

8. **Structural nulls.**
   Evaluate the all-mode shift-kernel and lapse-gradient templates that are
   compatible with the chosen radial basis. Report which known nulls are
   intentionally broken by homogeneous endpoint traces.

9. **Direct versus block rank.**
   Build \(D\), \(M=EJE^T\), and the full primary/secondary bracket matrix. Compare
   its rank with
   \[
   2\operatorname{rank}D+\operatorname{rank}(N^TMN)
   \]
   and preserve the null basis used in the comparison.

10. **Schur omission control.**
   Run the reducer once with the shift Schur correction \(E^TQA^{-1}Q^TE\)
   removed. Choose a control known to couple to the eliminated image and predict
   a changed quadratic form, compatibility condition, or spectrum before running
   it. A particular gravitational channel need not change every observable.

11. **Stopped-chain control.**
   Stop the stabilization after \(Bw=0\) and check its residual compatibility and
   Hamiltonian descent. A premature quotient can be undefined, rather than merely
   have a different spectrum. The final pass must include all generated rows.

12. **Boundary sensitivity.**
    Repeat the tiny \(N=3,n_r=3\) assembly under strict homogeneous and
    derivative-compatible gradient traces. Treat differences as boundary-model
    sensitivity, not as a contradiction.

13. **Resolution smoke test.**
    Increase only \(n_r\). Track rank gaps, condition numbers, and low physical
    eigenvalues. No acceptance threshold should be set until the endpoint model is
    chosen and a structural rank gap is visible.

14. **Nonlinear caveat.**
    Record that success of the finite linear reducer does not prove nonlinear
    constant rank, absence of linearization instabilities, a quantum measure, or a
    causal-diamond boundary Hilbert space.

## 8. Current obstruction and smallest missing derivation

A genuinely common domain cannot yet be supplied for the original canonical
generators because the full \(K\) and \(E\) weak forms have not been exported.
The corrected lapse theorem supplies a homogeneous domain for the **reduced**
complementary lapse bracket; it does not specify how endpoint traces in \(\mathcal J\),
\(V''\), and the Legendre transform combine before shift elimination.

The smallest derivation needed from the operator assembly track is:

1. export the weak \(K_{qq}\), \(K_{pp}\), and \(E\) bilinear forms, including their
   endpoint Green bilinears, in the Paper I \(q=(\lambda,X,V,C)\),
   \(z=(a,b,W)\) ordering;
2. prove that \(\mathcal J s=J_0s+J_1s'\), its adjoint, and
   \(A_R^{-1}\mathcal J^*p\) map the proposed shift and momentum spaces into the
   allowed traces, or exhibit the endpoint obstruction;
3. derive \(\int p\,\dot q\,dr\) and the mass matrix converting coefficient
   momenta to reducer-canonical \(p\);
4. evaluate the endpoint bilinears on the \(H^1_0\) bubble domain and on the
   derivative-compatible gradient domain.

If the endpoint bilinears vanish on \(H^1_0\) and the known structural nulls either
survive or are deliberately excluded by traces, the homogeneous realization can be
used as a provisional diagnostic without violating the known constraint analysis.
If the bilinears do not vanish, the domain must be tightened or an explicit
boundary action/trace prescription is needed before any physical spectrum run.
