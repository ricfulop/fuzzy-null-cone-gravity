# A common classical and quantum gravitational construction

## User's objective

Calculate finite-N one-loop four-graviton amplitudes with Gaunt interactions
and their controlled physical large-N limit; derive the gravitational corner
states and their dynamics; construct the Euclidean--Lorentzian relation from
the same dynamics. These are three parts of one theory, not independent toy
examples. **The objective is active and has not been achieved.**

The [common-regulator proposal](reports/common-regulator-bridge.md) removes an
operation-order mismatch by evaluating the complete angular contractions of
one metric-coefficient Galerkin action. This is a proposed revision of Paper
II's repeated binary-product prescription, not an equivalence already proved
between the historical models. It does not, by itself, define the quantum
measure, constraints, or gravitational boundary Hilbert space.

## What completion requires

1. **A defined common theory.** Fix the coefficient-to-metric map, complete
   Gaunt tensors, canonical constraints and their classification, boundary
   action, physical state space, quantum measure, and the meaning of any
   occupation restriction. If a regulator changes, describe the change and
   redo the affected calculations. Gauge fixing cannot silently substitute a
   different theory.
2. **The requested amplitudes.** Define the finite-volume/rank observable and
   its relation to asymptotic scattering; derive all required vertices,
   determinants and intermediate states from item 1; evaluate the one-loop
   four-graviton observable, including rational contributions where relevant;
   check unitarity cuts, factorization, the applicable symmetry identities and
   the matching Einstein benchmark. Separate infrared subtraction,
   ultraviolet regularization and renormalization. Demonstrate convergence at
   fixed physical kinematics, not merely fixed harmonic number on a sphere
   whose radius is changing. A bare trace-log formula or a scalar loop is not
   this result.
3. **Dynamical corner states.** Pull back the chosen gravitational boundary
   and joint action, derive its symplectic structure and charges, impose the
   actual reduction/gluing constraints, quantize with a specified positive
   inner product, and derive the spectrum, multiplicities and Hamiltonian.
   Explain the origin of any finite dimension and relation to area. Chosen
   spin labels, entropy-matching degeneracies or a two-state gluing example
   do not complete this task.
4. **Euclidean--Lorentzian reconstruction.** Derive the Euclidean measure or
   transfer operator from that physical Hamiltonian and prove the needed
   positivity, semigroup/sewing, analytic continuation and unitary real-time
   evolution. Establish the observable map to any claimed simulation.
   Hermiticity of a separately chosen finite matrix is not enough.
5. **Controlled limits.** State which cutoffs are fundamental and which are
   regulators; control angular, radial/time, occupation, volume and infrared
   limits separately. Recover the classical Einstein regime and the
   appropriate quantum effective-theory predictions without assuming the
   desired microscopic completion.

## Current evidence and next decision

The archived v0.2 results establish continuum weak-field/tree benchmarks and
several regulator diagnostics. They do not establish items 1--5. Paper I's
[nonlinear-consistency report](../paper1/reports/nonlinear-consistency-2026-11-27/README.md)
also records unresolved canonical and Euclidean issues.

The first tracked calculation has established a nonzero lapse primary--
secondary bracket on continuum vacuum initial data in the common action.
Exact harmonic arithmetic and an independent finite-amplitude Legendre
calculation agree. The [result and analytic cutoff-shell limit](reports/quantum-bridge/constraint-result-and-next-calculation.md)
show why a full Dirac completion or a declared reduction-before-cutoff
revision is required before importing continuum lapse/ghost assumptions.
The [kinetic-completion result](reports/quantum-bridge/kinetic-completion-result.md)
verifies the leading Schur correction and its associated angular determinant.
The second lapse derivative cancels in the five tested cases; higher-order
Hamiltonian terms still survive at finite auxiliary cutoff. This is a
component of a repair, not a complete repair. The
[mixed spatial/lapse constraint bracket](reports/quantum-bridge/canonical-mixed-bracket.md)
of a separately declared Hamiltonian-first candidate is also nonzero on
the proposed constraint surface: 22 exact cases and ten independent full
Hamiltonian-flow checks confirm the prediction. Thus lapse linearity alone
does not restore the gravitational algebra. Equivalence to the original
Lagrangian regulator is not assumed. A valid loop calculation still needs
either its actual completed constraints and measure or a fully specified
revision that derives the corresponding physical observables.

The [corner and reconstruction derivations](reports/quantum-bridge/corners-and-reconstruction.md)
identify the continuum area--angle charge, the missing microscopic state
count, and the exact transfer/reflection-positivity construction conditional
on a physical Hamiltonian. They do not fill the missing finite-rank
gravitational Hilbert space or establish equivalence to the simulations.
The [explicit uniform-corner quantization](reports/quantum-bridge/corner-quantization-and-transfer.md)
constructs positive-area states, solves uniform boost matching, and gives
Euclidean and Lorentzian transfer kernels for the same boundary charge.
Its Hilbert space remains infinite inside an area interval and has no
finite Gibbs trace. The undeformed normal-frame algebra also cannot have
a nontrivial finite-dimensional unitary representation. These facts
identify requirements for a finite physical reduction; they are not a
derivation of the postulated microscopic multiplicities.

All three principal deliverables remain pending. No claimed percentage of
completion or guarantee of success is assigned to an unresolved theory.

The [continuum all-plus one-loop reference](reports/quantum-bridge/continuum-one-loop-reference.md)
provides an explicit quantum comparison target, including the rational term
lost by a strictly four-dimensional cut calculation. Its scalar tree and
state-count inputs belong to continuum Einstein gravity. The reference has
no finite angular rank and is not evidence that those inputs survive the
Gaunt constraint obstruction.

The [reduction-first canonical candidate](reports/quantum-bridge/reduced-canonical-gaunt.md)
now specifies a distinct repair route: restrict free TT data, reconstruct
the dependent geometry through the continuum constraints, and use a dual
Gram matrix to make the induced symplectic form canonical. Its local
hypotheses, change of regulator and unresolved interacting limit are
explicit. The same quadratic Hamiltonian induces a leading uniform
far-sphere area spectrum and positive Euclidean/Lorentzian evolution.
The first tracked nonlinear calculation passes its specified checks:
ten exact TT modes, independently verified momentum-direction interaction
coefficients, and separate convergence tests for dependent angular and
radial geometry. The free rank stays $N=3$; the full sampled constraint
residual is about $2.03\times10^{-6}$ even when the separate angular defect
is about $2.34\times10^{-14}$. Neither that limited boundary sector nor
this numerical result completes the microscopic corner or interacting
reconstruction goals. The next calculation derives the full cubic and
quartic Hamiltonian, including the constraint-generated nonlocal terms.
The [interaction derivation](reports/quantum-bridge/reduced-vertices.md)
now gives those full classical functionals for the declared reduced chart.
Its seven exact coefficient directions, held-out rotational invariants and
32 independent finite-metric evaluations passed their specified checks.
The largest scaled geometric-coefficient error is about
$3.04\times10^{-14}$. These are verified classical ingredients. A
finite-basis four-point loop can also require fifth/sixth-order vertices and
quantization terms, so this step does not complete the requested amplitude.

The [sixth-order reduction and loop prescription](reports/quantum-bridge/reduced-sixth-order-and-loop-kernel.md)
now derives the full higher-order bookkeeping. Scalar responses through
fourth order and two York vector inverses suffice for the required
constraint contributions; pointwise sources and the finite-basis inverse
Gram corrections remain essential. An ordered-partition formula retains
all one-loop four-point vertex families. These are derived evaluation
identities, not evaluated higher gravitational couplings or an amplitude.
The next implementation must supply those couplings, evaluate the loop,
assemble the connected observable and check its physical limits.
