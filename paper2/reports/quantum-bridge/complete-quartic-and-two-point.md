# Complete quartic tensors and a finite-mode gravitational loop

## Status and scope

This calculation, preregistered before its first run, uses the same
[reduced canonical Gaunt candidate](reduced-canonical-gaunt.md) and its
[verified cubic/quartic functionals](reduced-vertices.md).
It reconstructs the full quartic polynomial in the twenty phase-space
variables and evaluates the complete **one-loop two-point function** for
the ten retained physical modes. The four-point loop integrator is tested
separately on controls. Its complete gravitational fifth/sixth-order
inputs remain unevaluated, so this is not the requested four-graviton
amplitude.

The quantum prescription is explicitly bare midpoint/Weyl quantization
with canonical phase-space measure and no order-$\hbar$ counterterm.
Results are the coefficient of $\hbar$, not a claim that a loop expansion
is numerically small when all scales are set to one. A constant vacuum
energy subtraction does not remove these field-dependent contributions.
Neither a global interacting Hamiltonian nor its reflection positivity
is established by a formal loop coefficient.

The physical free rank is $N=3$: five electric and five magnetic total-$j=2$
profiles, with one selected radial profile per parity. The packet width
is one and $16\pi G=1$. Neither the angular rank nor the radial basis is
extrapolated in this run.

The first run completed the exact classical reconstruction but failed its
frequency-quadrature convergence control before reaching the gravitational
two-point stage. Section 5 records this partial result and the resolution
change for a separate child; no gravitational loop result is available yet.

## 1. Complete finite-mode interaction polynomial

Let $E_q,B_q$ be the STF configuration labels and $E_s,B_s$ the
STF labels of the dual-Gram momentum field $S$.
The five-invariant basis for either pure configuration or pure momentum
quartics is

$$
\mathcal I(A,B)=
\left((\operatorname{tr}A^2)^2,\,
(\operatorname{tr}B^2)^2,\,
\operatorname{tr}A^2\operatorname{tr}B^2,\,
(\operatorname{tr}AB)^2,\,
\operatorname{tr}(A^2B^2)\right).
$$

Each pure species pair $A_q^2C_s^2$ contributes three invariants,

$$
\mathcal J(A,C)=
\left(\operatorname{tr}A^2\operatorname{tr}C^2,\,
(\operatorname{tr}AC)^2,\,
\operatorname{tr}(A^2C^2)\right).
$$

There are four such electric/magnetic species choices.
The remaining even-parity sector is linear in each of
$A=E_q,B=B_q,C=E_s,D=B_s$, with five invariants

$$
\mathcal X(A,B,C,D)=
\left(
\operatorname{tr}AB\operatorname{tr}CD,\,
\operatorname{tr}AC\operatorname{tr}BD,\,
\operatorname{tr}AD\operatorname{tr}BC,\,
\operatorname{tr}ABCD,\,
\operatorname{tr}ACBD
\right).
$$

The third cyclic trace is redundant for symmetric trace-free
$3\times3$ matrices:

$$
\operatorname{tr}ABCD+\operatorname{tr}ABDC+\operatorname{tr}ACBD
=\frac12\left(
\operatorname{tr}AB\operatorname{tr}CD+
\operatorname{tr}AC\operatorname{tr}BD+
\operatorname{tr}AD\operatorname{tr}BC\right).
$$

These give $5+5+4\cdot3+5=27$ quartic coefficients. The representation
count provides an independent completeness check: the symmetric square
of spin two has channels $0,2,4$, whereas two distinct spin-two fields
have channels $0,1,2,3,4$. Spatial parity and momentum reversal exclude
the other species sectors.

The pure-configuration coefficients and the two cubic coefficients are
read from the frozen, previously verified source result. The other
coefficients are computed from the spatial Hamiltonian functionals:
twelve pure-species mixed directions, five cross-species directions
selected by a fixed seed and exact rank test, and five pure-momentum
directions. The resulting polynomial must reproduce every saved
direction from the preceding run, including its held-out mixed and
rotated cases.

The seven historical cases have cross-invariant row rank three. Two further
held-out directions, selected with a separate fixed seed and excluded from
the fit, raise the combined validation rank to five. Their direct spatial
mixed functional must equal the reconstructed polynomial exactly. This
adds two contractions, for 24 new directions in total. These checks are
essential: the entire cross-species sector drops out of the one-loop
two-point tadpole, so the loop comparison cannot validate its five
coefficients.

All quartic coefficients multiply $\pi^{3/2}/\sqrt2$.
The cubic coefficients multiply $(2\pi/3)^{3/2}$.
No loop output is used in fitting the classical coefficients.

## 2. Canonical normalization and derivatives

For each orthogonal profile with Gram norm $g_\alpha$, the normalized
canonical variables obey

$$
q^\alpha=\sqrt{\frac2{g_\alpha}}\,Q_\alpha,\qquad
S^\alpha=\frac{P_\alpha}{\sqrt{2g_\alpha}}.
$$

Consequently $H_2=\frac12\sum_\alpha(P_\alpha^2+\omega_\alpha^2Q_\alpha^2)$,
with $\omega_E=\sqrt{11/2}$ and $\omega_B=\sqrt{13/2}$.
These are basis frequencies, not graviton masses.
The complete polynomial is stored sparsely in interleaved raw variables
$(q^0,S^0,\ldots,q^9,S^9)$, preserving rational coefficients.
The declared factors map its derivatives to normalized
$(Q_0,P_0,\ldots,Q_9,P_9)$ variables before loop contraction.
Raw rational coefficients and normalization factors are exported separately.

Hessian derivatives supply all loop blocks:
$C_{Aij}=\partial_A\partial_{Q_i}\partial_{Q_j}H_3$ and
$D_{ABij}=\partial_A\partial_B\partial_i\partial_jH_4$.
Repeated indices carry their ordinary derivative factorials.
The absence of cubic momenta is verified in the reconstructed polynomial,
not assumed for a different regulator.

## 3. Complete one-loop two-point coefficient

The normalized phase-space free inverse Hessian is

$$
\mathcal G_i(\nu)=
\frac1{\nu^2+\omega_i^2}
\begin{pmatrix}1&\nu\\-\nu&\omega_i^2\end{pmatrix}.
$$

For external phase-space indices $A,B$ and external frequencies
$\Omega,-\Omega$, the one-loop two-point tensor is

$$
\Pi_{AB}(\Omega)=
\frac12\int\frac{d\nu}{2\pi}
\operatorname{tr}\big(\mathcal G(\nu)D_{AB}\big)
-\frac12\int\frac{d\nu}{2\pi}
\operatorname{tr}\big(\mathcal G(\nu)C_A
\mathcal G(\nu+\Omega)C_B\big).
$$

This includes the momentum-dependent quartic tadpoles.
Since all cubic blocks are purely configurational, an independent
evaluation reduces it to

$$
\begin{aligned}
\Pi_{AB}(\Omega)
={}&\sum_i\left[
\frac{D_{AB,Q_iQ_i}}{4\omega_i}
+\frac{\omega_iD_{AB,P_iP_i}}4\right]\\
&-\frac12\sum_{ij}C_{Aij}C_{Bji}
\frac{\omega_i+\omega_j}
{2\omega_i\omega_j[(\omega_i+\omega_j)^2+\Omega^2]}.
\end{aligned}
$$

The bubble integral follows independently by Fourier transforming the
product
$e^{-\omega_i|\tau|}/(2\omega_i)$ times
$e^{-\omega_j|\tau|}/(2\omega_j)$.
This analytic expression is compared with direct real-frequency
quadrature of the full phase-space blocks.

The mixed $\Pi_{QP}$ component vanishes at this order: the cubic
momentum block is zero, and the symmetric quartic mixed block contracts
with an antisymmetric off-diagonal free propagator. The run tests this
cancellation directly. Only after that check is the configuration-space
inverse propagator correction formed as

$$
\delta\mathcal D_Q(\Omega)=
\Pi_{QQ}(\Omega)-\Omega^2\Pi_{PP}(\Omega).
$$

The coefficient $\Pi_{PP}$ changes the kinetic term. Keeping only
$\Pi_{QQ}$ would not describe the complete two-point correction of this
canonical theory.

All ten internal modes participate. Three Euclidean frequencies are
used for each representative parity, with both $Q$ and $P$ external
legs; additional magnetic components test rotation-equivariance.
The total is twenty diagonal cases and two mixed checks.
The scalar frequency quadrature uses two resolutions, independently
of the exact spatial Gaussian/Gaunt contractions.

## 4. Four-point kernel controls and limits of the result

The shared loop integrator also evaluates the ordered-partition formula
in the [sixth-order derivation](reduced-sixth-order-and-loop-kernel.md).
Controls must include fifth- and sixth-order interactions, nonzero
external frequencies, permutations and canonical coordinate changes.
Such controls test combinatorics and frequency integration; their
vertices are not gravitational evidence.

For the reduced gravitational theory, the two-point loop requires only
$H_3,H_4$ and is therefore computable from the present complete tensor.
The four-point loop additionally needs $H_5,H_6$. A static trace-square-root
check tests only four configuration legs at zero frequency and $p=0$;
it cannot validate the missing general phase-space components.

The retained oscillator spectrum has discrete two-particle poles.
It is not a continuum scattering spectrum. A four-graviton result still
requires the full connected observable, external normalization,
analytic continuation, renormalization and controlled physical
angular/radial/asymptotic limits.
No cutoff-dependent frequency shift is identified here with a
graviton mass or a Lorentz-invariant prediction.

The calculation uses a declared perturbative quantization choice.
It does not exponentiate the unbounded quartic Taylor polynomial as a
global interacting Hamiltonian. Complete corner quantization,
microscopic state counting and interacting Euclidean--Lorentzian
reconstruction remain open.

## 5. First run: classical tensor verified, loop convergence failed

Experiment 0f758de2-3e35-4b9d-9e66-7916c13fc06e used source commit
087ab31c06d22046ef2dbf7d6da43e61dd098f47.
Run 196cf70e-cad9-4208-9819-2aaaf78359de terminated with an explicit
control failure. Its effective quadrature orders were 32 and 64,
with unchanged relative tolerance $10^{-9}$ and control absolute
tolerance $10^{-12}$. The complete raw log was read before this review.

The earlier stages completed:

- All 22 new fit contractions and both additional held-out contractions
  were evaluated with rational arithmetic. Combined held-out
  cross-invariant coverage has rank five.
- The assembled polynomial contains 24 cubic and 515 quartic monomials.
  It reproduces all seven historical full-polynomial values exactly.
- The normalized one-point tadpole contractions vanish in all twenty
  phase-space directions.

The two new held-out mixed coefficients were, in the declared common
quartic unit,

$$
-\frac{22606581439103}{220200960},\qquad
-\frac{600659397549413}{18790481920}.
$$

Their direct spatial and reconstructed invariant values coincide exactly.
The complete rational invariant coefficients and the effective execution
configuration are preserved in
[the partial-result record](data/complete-quartic-partial-results.json).
This record is extracted from the failed run's log; it is not a fabricated
successful final export.

At order 32, the unequal-frequency quartic-bubble controls have absolute
errors up to $1.0372\times10^{-9}$, approximately
$7.90\times10^{-9}$ relative. Mixed four-point controls also fail the
32-versus-64 convergence criterion. At order 64, every specified
reference/control comparison passes; the maximum absolute comparison
error is $5.55\times10^{-17}$. Nevertheless, the required agreement
between resolutions fails, so the control suite correctly returns false
and the runner stops.

The gravitational two-point stage and final export were not reached.
There is therefore no gravitational loop number from this run. The next
child uses 64 and 128 quadrature points while retaining the tolerance,
Hamiltonian, cases, and fixed run recipe. The failed source node and log
remain immutable evidence of the resolution requirement.
