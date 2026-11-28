# From the oblate constraint calculation to physical linear dynamics

**Status:** downstream derivation and test specification. No new gravitational experiment, remote job, or regulator change is made here. The exact checks below use manufactured finite quadratic systems, not gravitational measurements.

## 1. Reconcile the active work before computing

The Paper I snapshot `71a29e402e0941e2d85f9894c323d5440e0e98cb` establishes an exact static oblate-coordinate Minkowski family inside the actual metric-coefficient ansatz. Unlike the previously excluded nonzero-extrinsic-curvature witness, this family has K=p=0 and is dynamically admissible. A nonzero angular-shift primary–secondary bracket is established on it for sufficiently small nonzero deformation at N=3,4. This gives a rank lower bound on an actual solution, not the full rank or a radiative branch theorem. The source is `paper1/reports/admissible-branches-2026-11-27/README.md`; the saved verification output was read for this handoff.

The Paper I chat is already computing the all-m temporal Hessian, secondary brackets, and remaining preservation conditions. The Paper II snapshot `bab7d132503bfdd05ab0e5fc2ab48575ac0a9aaa` records the distinct Hamiltonian-first mixed-bracket obstruction and develops the uniform positive-area corner quantization. That corner construction has an infinite-dimensional state space; it does not complete finite corner counting.

**Division of work:** reuse the Paper I operator export and Paper II boundary/measure work. This handoff supplies the algebra that converts a stationary quadratic Hamiltonian into stabilized linear constraints and a physical linear spectrum. It does not recreate their tensor or corner calculations.

## 2. What primary–secondary rank can and cannot tell us

Let w=(q,p) be canonical spatial variables and z temporal variables with primary momenta pi_z. Write the Hamiltonian H(w,z), primary constraints pi_z=0, and secondaries C=H_z=0. With J the canonical Poisson matrix on w, define

$$
D=H_{zz},\qquad M=\{C,C\}_w.
$$

On the primary/secondary surface their bracket matrix is

$$
\Omega=\begin{pmatrix}0&-D\\D&M\end{pmatrix}.
$$

For finite real symmetric D of rank r, choose columns Q spanning its image and columns N spanning its kernel. Then A=Q^T D Q is invertible and the temporal change of coordinates P=(Q,N) gives the congruence P^T D P=diag(A,0). Q need not be orthonormal. Eliminating the paired image directions gives the exact rank identity

$$
\boxed{\operatorname{rank}\Omega=2r+\operatorname{rank}(N^TMN).}
$$

To see why off-null secondary brackets do not add another correction, the eliminated pair has matrix [[0,-A],[A,M_xx]] with inverse [[A^{-1}M_xx A^{-1},A^{-1}],[-A^{-1},0]]. The remaining secondary rows couple to this pair only through its secondary columns. Their Schur correction is zero. The remaining primary rows are null, and the remaining secondary block is N^T M N.

If constraint bases depend on phase space, their derivatives contribute terms proportional to constraints; the identity concerns the bracket evaluated weakly on the specified surface. It is not a final nonlinear degree-of-freedom count. For radial differential operators, adjoints, domains, boundary conditions, kernel/range decomposition, and inverses must be justified separately. A numerical matrix rank cannot replace those arguments.

The first preservation condition is

$$
\dot C=F+D\mu=0,\qquad F=\{C,H\}_w,
$$

where mu=dot z. Consequently

$$
N^TF=0
$$

is the remaining compatibility condition. It must be classified as an identity, a combination of existing constraints, or a new constraint and then preserved again. At an exact stationary solution F=0 automatically; its derivatives and nearby functional dependence are essential. Checking its value only at that solution does not close the algorithm.

## 3. A downstream calculation that uses the existing work

At an exact stationary background, with fixed canonical coordinates and specified boundary evolution, write the full quadratic Hamiltonian as

$$
H_2=\tfrac12w^TKw+z^TEw+\tfrac12z^TDz,
\qquad K=K^T,\quad D=D^T.
$$

The matrices are the Hessian of the actual retained Hamiltonian. K includes the spatial-curvature potential and radial derivative contributions. The background is stationary, so no time derivative of the basis is hidden in these formulas. A moving background requires additional terms.

Set z=Qx+Ny. The x equations give x=-A^{-1}Q^TEw. Eliminating x and its primary momentum leaves the canonical w bracket unchanged and yields

$$
K_{\mathrm{eff}}=K-E^TQ A^{-1}Q^TE,\qquad B=N^TE,
$$

$$
H_{2,\mathrm{red}}=\tfrac12w^TK_{\mathrm{eff}}w+y^TBw.
$$

Thus the remaining linear constraints are Bw=0, and

$$
\dot w=J(K_{\mathrm{eff}}w+B^Ty),\qquad
\dot{(Bw)}=BJK_{\mathrm{eff}}w+BJB^Ty=0.
$$

An invertible BJB^T fixes the remaining multipliers; a singular block produces additional compatibility conditions. This is where a Hessian-only analysis can miss constraints or mistake non-gauge directions for gauge freedom.

For an independent current set of constraint rows R, keep the **original** multiplier columns B^T and form

$$
L_R=RJB^T,\qquad
Z_R=\operatorname{basis}\ker L_R^T,\qquad
R_{\mathrm{new}}=\operatorname{rowspan}
\begin{pmatrix}R\\Z_R^TRJK_{\mathrm{eff}}\end{pmatrix}.
$$

Iterate until the row space stops growing. Do not assign new independent multipliers to every generated constraint. At termination check that the compatibility source vanishes on ker R, not merely at w=0.

On this stabilized finite linear system, the allowed perturbations are ker R. Its null symplectic directions are

$$
\operatorname{im}\big[J R^T\ker(RJR^T)\big].
$$

Choose a complement T to these directions inside ker R. The physical linear symplectic form and Hamiltonian are

$$
\omega_{\mathrm{lin}}=T^T(-J)T,\qquad
K_{\mathrm{lin}}=T^TK_{\mathrm{eff}}T,\qquad
\dot u=\omega_{\mathrm{lin}}^{-1}K_{\mathrm{lin}}u.
$$

This supplies a linear spectrum and energy-form test after stabilization, rather than a count inferred from D alone. A regular linear quotient need not extend to a regular nonlinear gravitational branch. In particular, accidental linearized gauge symmetries or linearization instabilities require additional nonlinear checks. The procedure also does not compute the quantum measure or its determinant.

The characteristic polynomial of the physical evolution is invariant under changes of the quotient basis T. Entries of K_lin alone are not canonical oscillator masses or frequencies: the pulled-back symplectic form must also be retained, or explicitly put into Darboux form.

## 4. Executable algebra check

`paper2/verification/constraint_handoff/reduce_quadratic.py` implements the finite exact-rational calculation, including constraint-chain generation and the symplectic quotient. It rejects floating-point/decimal input as exact-rank evidence.

From `paper2/`:

```bash
uv run --isolated --no-project --with sympy==1.14.0 \
  python verification/constraint_handoff/reduce_quadratic.py --self-check
```

Nine manufactured checks cover an ordinary oscillator, a frequency-changing auxiliary Schur term, an indefinite off-diagonal D block, hand elimination with a nonorthogonal image basis, a generated momentum constraint, detection of an artificially stopped constraint chain, a two-step first-class chain, a nonzero secondary bracket when D=0, and invariance under mixing null and non-null temporal coordinates. All passed with SymPy 1.14.0. These are lightweight verification of the adapter, not new results about the gravitational theory.

Each preservation step records the current rows R, multiplier matrix L_R, left-kernel basis Z_R, generated compatibility rows, and the updated independent rows with ranks before and after. This preserves an audit trail when row reduction recombines secondary and higher constraints. An independent read-only review found no blocking algebra issue; the congruence wording, quotient-basis caveat, step diagnostics, and two additional manufactured controls incorporate its recommendations.

The optional JSON input uses K,E,D objects, each with `shape` and row-major `entries`. Canonical order is all q followed by all p. Integers and fraction strings are accepted. The actual deformed gravitational integrals may require symbolic functions or controlled numerical operators; they must not be converted into apparently exact rationals to use this prototype. A numerical version should follow the actual upstream export rather than duplicate it in advance.

## 5. Proposed next decisive test: stabilized oblate quadratic dynamics

**Hypothesis to test:** at a specified small nonzero oblate deformation, the actual finite-rank action admits a stabilized quadratic physical sector with regular evolution and the expected low-energy spin-two behavior. This is a local necessary test, not a proof of the full classical–quantum connection.

**Inputs required from the ongoing Paper I calculation:**

1. The source commit, field ordering, real harmonic conventions, background, lapse/shift, and action normalization.
2. The complete stationary Hessian blocks K,E,D, or radial differential operators from which they are assembled. D or its rank alone is insufficient.
3. The radial inner product, boundary conditions, and canonical normalization. If a radial Galerkin basis is used, momenta must be canonically dual to coordinate coefficients; an unaccounted mass matrix invalidates the use of J above.
4. Exact zero modes or symmetry identities, and independently evaluated matrix elements, to distinguish structural zeros from integration errors.

**Smallest useful launch specification, once inputs exist:** start with N=3 and one small nonzero deformation in the verified spatial-Legendre domain. Include all angular m sectors. Retain radial differential terms through an explicit variational discretization with the same boundary convention. Record the dimension and every approximation. Start with a small radial basis and a second independent integration resolution; extend radial resolution and N only after the initial export and reductions agree. The chosen radial truncation is an additional computational regulator and must be documented as such.

**Required outputs:**

- Spatial-Legendre regularity and temporal block spectra, with precision/integration sensitivity.
- Direct primary–secondary bracket rank versus the block identity above.
- Initial and final constraint row spaces, newly generated rows, and the residual compatibility equations.
- Pulled-back symplectic null directions and physical linear evolution spectrum, with physical normalization.
- Energy-form signature and sensitivity to numerical conditioning, boundary conditions, radial resolution, and small changes of deformation.
- The nonlinear conditions still needed to extend the result beyond the linearized sector.

**Controls:** recover the archived oblate shift-Hessian contraction in its stated scope; independently evaluate the quadratic form from variations of the same action; verify integration by parts with the chosen radial domain; check temporal basis invariance; include a deliberately omitted Schur term and an artificially stopped constraint chain. Use the round background as a separate symmetry-enhanced control, not as evidence of constant rank through the round limit.

**Interpretation rules:**

- A converged physical instability or an incompatible constraint chain is a reason to review that candidate/branch. It is not a theorem against every finite gravitational theory.
- A very large inverse small eigenvalue of D is not by itself a physical instability: gauge restoration or a singular choice of reduced coordinates can produce it. Assess the stabilized, normalized physical operator.
- A successful small-rank linear test justifies broader nonlinear and physical-limit work; it does not close the constraint, loop, or corner gates by itself.
- A local radial principal-symbol calculation, if selected as an earlier inexpensive subtest, answers only the corresponding characteristic question. It cannot be presented as the finite-shell mode spectrum or full radial constraint problem.

No numerical acceptance threshold, runtime, or exact physical mode count is fabricated before the operator export is available. Precision/integration bounds must justify any claimed rank gap. The initial gravitational experiment is therefore **specified but not launch-ready**; the missing input is the full stationary Hessian/operator export, not access to a GPU.

## 6. Efficient resource use and ownership

The user reports Modal B200 access and two DGX Sparks through SSH. The configured environment reports Modal and SSH availability, but this task has not connected to the hosts, verified their identity, changed the compute target, or allocated a device.

- Use local CPU exact arithmetic for the small algebra and formulation checks.
- Use the Sparks for independent numerical evaluation or parallel moderate matrix jobs if profiling demonstrates a benefit and the host/runtime setup is verified.
- Reserve Modal B200s for sufficiently large batched tensor contractions or numerical scans after a representative kernel has a measured performance and precision case. Exact symbolic rank proofs are not automatically accelerated by moving them to a GPU.
- Before launch, specify input commit, experiment hypothesis, configuration, resource estimate, and proposed target. Each experiment and GPU budget still requires the user's approval.

This chat owns the reduction/handoff artifact. The ongoing Paper I chat owns the oblate operator/rank derivation; the ongoing Paper II chat owns its quantum/corner program. No messages assigning new work to those chats and no edits to their branches were made. The next scientific action is to evaluate this test against their exported operators, or report precisely which block is still missing, rather than start a competing implementation.
