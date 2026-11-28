# What can be derived about corners and Euclidean evolution

**Status:** continuum corner dynamics under stated boundary conventions and
a conditional reconstruction theorem. Neither derives Paper I's finite
corner multiplicities, nor identifies the archived Monte Carlo measure
with the common gravitational action. The requested complete construction
remains open.

## 1. The gravitational joint action fixes an area--angle charge

For a spacelike joint C of two null hypersurfaces, choose the null-normal
normalizations and orientations. In the additive convention the joint term is

$$
S_C=\frac{\varepsilon_C}{8\pi G}\int_C d^2y\sqrt\sigma\,\eta,
\qquad \eta=\log\left| -\frac12 k\cdot\bar k\right|,
$$

with the orientation sign epsilon_C fixed by the adjoining segments. The
factor 1/(8 pi G) follows from converting Lehner et al.'s S=16 pi G I and
their joint contribution 2 integral a dS. Their null-normal prescription and
its remaining normalization dependence are on
[pp. 24--25](https://www.alphaxiv.org/pdf/1609.00207?page=25).
The classical additivity result is not a quantum gluing theorem; their
discussion explicitly identifies obstacles beyond the on-shell action on
[p. 5](https://www.alphaxiv.org/pdf/1609.00207?page=5).

The action term alone should not be mistaken for an already complete corner
symplectic form: which data vary, boundary polarization, exact variations,
and the bulk constraints must also be specified. In the extended
area--angle polarization, where the relative boundary angle is allowed to
vary, the established conjugate pair is the area density and that angle.
An explicit derivation for black-hole topology appears in Carlip and
Teitelboim's variation of the action and change of polarization,
[pp. 3--5](https://www.alphaxiv.org/pdf/gr-qc/9312002v3?page=4).
It concerns that boundary setup; importing it into every finite causal
diamond requires deriving the corresponding symplectic reduction.

For the spatially constant corner mode, with an orientation convention
chosen so that the symplectic form is positive in this ordering, write

$$
a=\frac{A}{8\pi G},\qquad
\Theta_C=a\,\delta\eta,\qquad
\Omega_C=\delta a\wedge\delta\eta,\qquad
\{\eta,a\}=1.
$$

Changing polarization to -eta delta a changes Theta by an exact variation
and leaves Omega unchanged. The charge generating eta -> eta+lambda is a.
For a prescribed boost rate nu(t), a choice of oriented corner evolution is

$$
H_C=\nu(t)a,\qquad \dot a=0,\qquad\dot\eta=\nu(t).
$$

This gives a genuine boundary charge/evolution relation within the stated
extended boundary phase space. It is not yet the physical Hamiltonian of a
whole diamond: normal-frame transformations may be gauge after gluing, and
bulk constraints and the boundary observer determine what survives.

## 2. This dynamics does not derive a finite microscopic count

The Lorentzian boost angle is noncompact. For any area interval of positive
width, the unreduced phase-space volume

$$
\int_{a_1}^{a_2}da\int_{-\infty}^{\infty}d\eta
$$

is infinite. Finite angular rank changes the number of coefficient pairs,
not their continuous ranges. A finite occupation/area sector is an
additional reduction or microscopic law until it is derived from constraints.

If one tries to retain an exact canonical commutator on a nonzero finite
Hilbert space, taking a trace of [eta,a]=i hbar I yields 0=i hbar dim(H).
This rules out that particular finite representation, not a reduced or
deformed corner algebra. On a positive-area half-line, an unrestricted
self-adjoint canonical angle translation also cannot shift the area
spectrum through all real values; the representation and domains must be
chosen with care. An affine representation or a constrained boundary
quantization is possible, but neither automatically supplies the postulated
finite degeneracy.

Wick-rotating a smooth horizon sets the *geometric* Euclidean opening angle
to 2 pi. This is not, by itself, a compact identification of the Lorentzian
boost coordinate and does not prove a discrete spectrum. Carlip and
Teitelboim recover the semiclassical area entropy but state: “We still lack
a ‘microscopic’ explanation for the exponential weight in the integration
measure for the surface degrees of freedom”.
[p. 5](https://www.alphaxiv.org/pdf/gr-qc/9312002v3?page=5)

Thus A/(4G hbar) in a saddle-point entropy and dim H_corner(A) are different
claims. The current bulk-plus-joint data do not derive Paper I's 2^M or 3^M
corner multiplicities, its area spacing, or its rank/area law. A successful
completion must derive the reduced corner algebra, its positive inner
product and multiplicities, its coupling to bulk dynamics, and its gluing
measure from one specified construction.

## 3. Conditional Euclidean--Lorentzian reconstruction is exact

Suppose that the missing constraint reduction has supplied a positive
physical Hilbert space and a self-adjoint physical H bounded below. Subtract
its lower bound E_0. The spectral theorem then defines

$$
T(\tau)=e^{-\tau(H-E_0)/\hbar},\qquad \tau\ge0,
\qquad U(t)=e^{-it(H-E_0)/\hbar}.
$$

T is a positive operator and a strongly continuous contraction semigroup;
U is unitary. They obey T(tau+sigma)=T(tau)T(sigma), and the spectral
matrix elements admit continuation from positive Euclidean time to real
time with the usual domain qualifications. In a finite-dimensional physical
sector all these statements are elementary finite spectral sums.

Reflection positivity follows from a norm identity. For a finite-dimensional
sector, let beta=2T and use any sequence of bounded observables at
0<t_1<...<t_n<T. Define

$$
B_F=e^{-(T-t_n)H/\hbar}O_n e^{-(t_n-t_{n-1})H/\hbar}
\cdots O_1 e^{-t_1H/\hbar}.
$$

Extend F -> B_F linearly to linear combinations of such histories. Reflection
reverses ordering and takes adjoints, so the reflected quadratic form is

$$
\langle\theta F\,F\rangle_\beta
=\frac{\operatorname{Tr}(B_F^\dagger B_F)}
{\operatorname{Tr}(e^{-\beta H/\hbar})}\ge0.
$$

For infinite dimension the stated traces must exist; alternatively a
ground-state version uses a Hilbert-space norm instead of a thermal trace.
This is a proof given a physical H, not a derivation of H from an arbitrary
Euclidean sampling action. Positive operator T does not imply that every
matrix element in a chosen basis is nonnegative, so a positive local
Boltzmann measure for Monte Carlo is an additional property.

## 4. Projection and reconstruction do not commute automatically

For a projector P, the two candidates

$$
T_P(\tau)=e^{-\tau PHP/\hbar},\qquad
\widetilde T_P(\tau)=P e^{-\tau H/\hbar}P
$$

generally differ. The first is a semigroup for a changed projected
Hamiltonian. The second retains excursions through Q=1-P and obeys

$$
\widetilde T_P(\tau)\widetilde T_P(\sigma)
-\widetilde T_P(\tau+\sigma)
=-P e^{-\tau H/\hbar}Qe^{-\sigma H/\hbar}P.
$$

Its leading difference is -tau sigma PHQHP/hbar^2, when the expansion is
defined. Consequently an exact integration of discarded dynamical states
usually yields memory or energy-dependent effective interactions, rather
than simply the original budget-compressed H. The matching resolvent is

$$
P(E-H)^{-1}P
=\left[E-PHP-PHQ(E-QHQ)^{-1}QHP\right]^{-1}
$$

on the stated resolvent domains. This is precisely the distinction that
also matters for virtual intermediate states in a loop amplitude.

## Current conclusion

We have a continuum area--angle boundary relation and an exact conditional
transfer construction. The finite-rank gravitational physical Hilbert space,
its corner multiplicities and its Hamiltonian are still needed to join them.
Completing a formal one-loop effective action or a classical boundary charge
does not by itself provide these missing nonperturbative ingredients.
