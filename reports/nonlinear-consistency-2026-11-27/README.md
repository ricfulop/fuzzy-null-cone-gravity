# Nonlinear consistency: first obstruction and construction checkpoint

**27 November 2026.** Starting point: action/clock repair PR #6, commit `ab303dd3f756c7f3f8d393ef090d4da35a7fdae2`. The governing research objective and completion tests are in [GOAL.md](../../GOAL.md).

This checkpoint narrows the allowed repairs and supplies a small exact boundary construction. It does not establish the full finite-rank gravitational constraint algebra, derive the gravitational boundary Hilbert space, or prove equivalence to the archived simulation. The distinction between a failed candidate and a general impossibility theorem is essential here.

| Question | Established in this checkpoint | Still needed |
|---|---|---|
| Can the angular continuum bracket simply be truncated? | No: an exact L=3 Jacobi counterexample | Canonical constraints of the actual bulk action and a consistent reduction or deformation |
| Can a matrix bracket repair this angular algebra? | Yes for the selected area-preserving subalgebra; it changes the bracket | Full momentum/Hamiltonian constraint algebra and invariance of the same action |
| Can all classical boundary operators act exactly on a finite Hilbert space? | Not if they include an exact canonical pair or nontrivial unitary sl(2,R) representation | Identify the physical reduced boundary algebra before applying the state cutoff |
| Is finite boundary gluing necessarily impossible? | No: an explicit reduced R x SU(2) model has an exact positive gluing projector | Derive that model, its spectrum, and its bulk coupling from gravity |
| Is the archived HMC action directly the real Euclidean Einstein action? | Not under an identification retaining the real conformal direction | A derived reduction, measure, contour, or effective-theory relation |

## 1. Exact failure of the naive angular bracket

Let x, y, z be Cartesian coordinates on the unit sphere, r^2=x^2+y^2+z^2=1, and use the sphere Poisson bracket

$$
\{f,g\}=\boldsymbol{x}\cdot(\nabla f\times\nabla g).
$$

The Hamiltonian vector fields of these functions form the area-preserving subalgebra of angular diffeomorphisms. Let P_L be the orthogonal spherical-harmonic projection onto degrees at most L and propose the finite bracket

$$
[f,g]_L=P_L\{f,g\}.
$$

Take the harmonic polynomials

$$
f=xz,\qquad g=yz,\qquad h=x(5z^2-r^2).
$$

Their degrees are 2, 2, and 3. Exact Fischer decomposition and polynomial differentiation give

$$
[[f,g]_3,h]_3+[[g,h]_3,f]_3+[[h,f]_3,g]_3
=\frac{15}{7}y(5z^2-r^2)\ne0.
$$

At (x,y,z)=(0,1,0), the residual is -15/7. The continuum Jacobiator is exactly zero. Retaining the intermediate degree-four harmonics makes the Jacobiator of this particular triple zero at L=4; that is not a proof of general closure at L=4. Products of other retained modes can generate higher degrees again.

**Consequence:** this naive projected bracket cannot be the claimed Lie algebra of angular gauge generators. This is not yet a calculation of the full ADM constraint bracket of the new Galerkin action. A first-class constraint algebra can also have field-dependent structure functions, and that possibility must be tested from the actual action, on its constraint surface. The established failures of gauge symmetry in other gravitational discretizations are relevant warnings, not proofs about this regulator [1].

**Code and data:** `projected_sphere_algebra.py`, `projected_sphere_algebra.json`. The projection uses exact rational polynomial arithmetic; it is not a numerical quadrature artifact.

## 2. An angular deformation that does close

For Hermitian N x N matrices, the scaled commutator

$$
[F,G]_N=-i\sqrt{j(j+1)}(FG-GF),\qquad j=(N-1)/2,
$$

satisfies Jacobi exactly by associativity. Traceless Hermitian matrices close under this bracket. A symmetrized polynomial quantization based on X_i=J_i/sqrt(j(j+1)) allows a direct comparison with the preceding classical functions.

For the same f,g,h, numerical matrix Jacobi residuals are at most 4.2e-15 over N=4,8,16,32. For f=xz and g=yz, the operator norm of

$$
[Q_N(f),Q_N(g)]_N-Q_N(\{f,g\})
$$

is 0.0172133, 0.00466623, 0.00122789, and 0.000315811 at those ranks. The decrease is consistent with the expected low-polynomial semiclassical behavior; these four values are not a general convergence theorem.

This repairs a specific angular Lie bracket by deforming it. It does not supply the hypersurface-deformation algebra, the Hamiltonian constraint, the mixed brackets, or invariance of the Einstein--Hilbert Galerkin action. One cannot insert this result into the manuscript as full constraint closure.

**Code and data:** `deformed_angular_candidate.py`, `deformed_angular_candidate.json`.

## 3. Boundary operators: what finite dimension excludes

If a chosen unreduced corner polarization contains an area/boost-angle canonical pair, write its nonzero canonical constant as kappa. An exact representation

$$
[\widehat A,\widehat\eta]=i\kappa\mathbf1
$$

on a d-dimensional space is impossible: the trace of the left side is zero while the right side has trace i kappa d. This elementary trace obstruction is independent of numerical accuracy. It applies only if both members of this canonical pair are required to be operators on that same finite space.

There is also no nontrivial finite-dimensional unitary representation of sl(2,R). For Hermitian generators with

$$
[J,K_1]=iK_2,\quad[J,K_2]=-iK_1,\quad[K_1,K_2]=-iJ,
$$

trace cyclicity gives

$$
i\operatorname{Tr}K_2^2
=\operatorname{Tr}([J,K_1]K_2)
=\operatorname{Tr}(J[K_1,K_2])
=-i\operatorname{Tr}J^2.
$$

Positivity forces J=K_2=0, and then K_1=0. This does not exclude finite representations of the Abelian boost group R alone, compact rotations, gauge-invariant reduced observables, or a finite physical sector obtained from a larger kinematical space. Gravity's edge-mode and corner reductions must determine which operators actually survive [2,3].

### Quantified compression instead of an exact canonical claim

For a projection P and Q=1-P,

$$
[PAP,PBP]=P[A,B]P-PAQBP+PBQAP.
$$

For a normalized state psi in ran(P), the magnitude of the expectation of the defect is at most

$$
2\|QA\psi\|\,\|QB\psi\|.
$$

Thus a compression can be controlled on specified states when leakage is small. It is not an exact representation merely because the retained dimension is large. For the d-level oscillator truncation,

$$
[a_d,a_d^\dagger]=\mathbf1-d|d-1\rangle\langle d-1|.
$$

The canonical defect vanishes on the vacuum but has operator norm d. Exact tests at d=2,4,8,16 and general matrix checks of the compression identity are supplied. Statewise accuracy and uniform operator-norm accuracy are different claims.

**Code and data:** `boundary_compression.py`, `boundary_compression.json`.

## 4. A constructive reduced boundary gluing model

A consistent finite boundary model need not contain a boost-angle operator. Choose a reduced boundary algebra with area labels and compact spins. Treat the boost angle eta as an external gauge parameter acting by area characters. For two regions to be glued, use dual orientations:

$$
U_L(\eta)\otimes U_R(\eta)^*,\qquad
C_A=A_L\otimes\mathbf1-\mathbf1\otimes A_R^T.
$$

Area matching follows from the invariant subspace of this action. For a finite set of distinct area characters, the normalized long-time group average projects onto equal-area blocks. This uses a representation of the Abelian group R, which is compatible with finite dimension, and does not impose a canonical angle operator.

For a matching spin-j block, compact rotation invariance is implemented in V_j tensor V_j* by

$$
|\Omega_j\rangle=\frac1{\sqrt{2j+1}}
\sum_{m=-j}^j|jm\rangle\otimes\overline{|jm\rangle}.
$$

The supplied example chooses two distinct area labels and one spin-1/2 and one spin-1 block. Each boundary has dimension five. The tensor product has dimension 25, and

$$
P_{\rm glue}=|\Omega_{1/2}\rangle\langle\Omega_{1/2}|
+|\Omega_1\rangle\langle\Omega_1|
$$

is an exact Hermitian rank-two projector. It is the complete common kernel of C_A and the three dual-orientation rotation constraints. The same-sign area sum has no invariant state for these positive area labels, demonstrating why orientation matters. An explicit Hermitian Hamiltonian on the two-dimensional reduced space has unitary rational Cayley evolution, and embedded kernels sew associatively. All these checks use exact symbolic arithmetic.

This is a usable finite gluing construction for the selected R x SU(2) reduced algebra. Its area labels and spin content are test inputs; they are not a derivation of the manuscript's corner spectrum. The full gravitational corner algebra includes additional structure [2,3]. Deriving this reduction, proving that omitted generators are gauge or controlled approximations, and coupling it to the bulk remain required. A complete null-boundary variational principle also needs the null-generator and joint prescriptions [4].

**Code and data:** `reduced_boundary_gluing.py`, `reduced_boundary_gluing.json`.

## 5. Direct real Euclidean action equivalence fails a sign test

Use the standard convention

$$
S_E[g]=-\frac1{16\pi G}\int\sqrt g\,R[g]
$$

at zero cosmological constant. For a conformal rescaling of a Ricci-flat four-metric, g=A gbar with A>0, the conformal derivative contribution, after a compactly supported integration by parts, is

$$
S_E[A\bar g]-S_E[\bar g]
=-\frac3{32\pi G}\int\sqrt{\bar g}\,
\frac{|\bar\nabla A|^2}{A}.
$$

It is negative for a nonconstant compactly supported conformal deformation. The general curvature routine independently checks the flat example g=A(u) delta, obtaining

$$
\sqrt g R=-3A''+\frac32\frac{(A')^2}{A}.
$$

For A=1+epsilon cos(u) on a periodic interval of length 2pi, per unit transverse volume and G=1, the quadratic Euclidean coefficient is -3/32. The compactly supported version establishes the local sign issue without changing boundary data; the periodic example is only a convenient independent calculation.

In contrast, set Z=0 and choose a spatially constant matrix quadrupole Phi=epsilon T_20 in the archived HMC action. With Tr(T_20^2)=N, mu^2=0.25, and lambda_4=0.1,

$$
S_{\rm HMC}=\frac12(6+0.25)N\epsilon^2
+\frac{0.1}{4N}\operatorname{Tr}(T_{20}^4)\epsilon^4.
$$

The scalar quadratic coefficients at N=3,4,6,8 are 9.375, 12.5, 18.75, and 25. The DeWitt term is linear and vanishes on this traceless probe; the shear contribution is zero. The check extracts the original scalar action statements and harmonic-basis builder from `modal_fsg_hmc.py` and executes them with a NumPy tensor adapter. It does not import Modal or launch a sampler. Agreement with the analytic scalar expression is better than 4e-16 in the tested evaluations.

These two probes are not asserted to be an already established field map between the theories. They identify the necessary sign problem for any proposed direct real map retaining the Einstein conformal direction. More generally, a smooth real invertible change of variables at corresponding stationary points transforms the Hessian by congruence and cannot remove negative directions. The HMC model's positive real scalar stiffness therefore cannot simply be called the unreduced real Euclidean Einstein conformal action.

A reduced effective theory, a complex integration contour, or a universality relation is not excluded by this test; each needs its own derivation and observable map. In particular, imposing gravitational constraints is not a general proof that a real Euclidean measure becomes bounded: explicit nonlinear counterexamples are known [5]. Lorentzian contour approaches provide another route to examine [6]. Neither reference establishes an equivalence to this HMC model.

**Code and data:** `euclidean_signature.py`, `euclidean_signature.json`. The archived simulation driver is unchanged; its hash is recorded.

## Reproduction and next decision

```bash
python verification/closure_audit_20261127/run_all.py
```

Dependencies: Python, NumPy, SciPy, SymPy. Five checks run locally on CPU. `reproducibility.json` records versions, run status, timing, and source hashes.

The working direction is to derive and reduce the actual gravitational constraints before applying the finite state restriction, derive the boundary algebra on that reduced space, and then derive a computational measure. The reduced gluing example shows one way a finite physical boundary sector can avoid an impossible canonical-operator requirement. The angular deformation shows a possible finite algebra but does not determine whether it is compatible with the chosen action. The existing simulations remain results for their stated matrix model until a matching calculation succeeds.

## Primary sources

1. B. Bahr and B. Dittrich, *(Broken) Gauge Symmetries and Constraints in Regge Calculus*, Class. Quantum Grav. 26, 225011 (2009), [arXiv:0905.1670](https://arxiv.org/abs/0905.1670). Relevant discretization/constraint distinction; not a theorem about the present angular regulator.
2. W. Donnelly and L. Freidel, *Local subsystems in gauge theory and gravity*, JHEP 09 (2016) 102, [arXiv:1601.04744](https://arxiv.org/abs/1601.04744). Boundary edge modes, surface symmetries, and gluing framework.
3. L. Freidel, M. Geiller, and D. Pranzetti, *Edge modes of gravity. Part III. Corner simplicity constraints*, [arXiv:2007.12635](https://arxiv.org/abs/2007.12635). Corner reduction and its noncompact algebra; the present finite example is not a quantization of their complete construction.
4. L. Lehner, R. C. Myers, E. Poisson, and R. D. Sorkin, *Gravitational action with null boundaries*, Phys. Rev. D 94, 084046 (2016), [arXiv:1609.00207](https://arxiv.org/abs/1609.00207).
5. G. T. Horowitz, D. Marolf, and J. E. Santos, *Constraints are not enough*, [arXiv:2505.13600](https://arxiv.org/abs/2505.13600). Nonlinear counterexamples to a general positivity claim after imposing constraints.
6. D. Marolf, *Gravitational thermodynamics without the conformal factor problem: Partition functions and Euclidean saddles from Lorentzian Path Integrals*, [arXiv:2203.07421](https://arxiv.org/abs/2203.07421).
