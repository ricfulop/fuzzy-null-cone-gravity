# A continuum four-graviton quantum reference

## Scope and status

This report derives an evaluation recipe for the one-loop **all-plus**
four-graviton amplitude of continuum pure Einstein gravity. It supplies a
quantum reference for the ongoing finite-rank investigation. **It is not a
finite-$N$ Gaunt amplitude, a large-$N$ extrapolation, or a constraint repair.**

The numerical verification described below is preregistered here. Its
tracked run has not yet been reviewed in this revision. An analytic formula
below is distinguished from a measured numerical residual.

The inputs are published scalar Compton tree amplitudes, their KLT relation,
and a continuum state-count identity. We evaluate their dimensionally
regulated loop integrals independently of the published final answer.
This is a reproduction of known Einstein effective-theory physics, not a
claim of a new amplitude or a derivation from Paper I's action.

The earlier [mixed-bracket result](canonical-mixed-bracket.md) explains why
the continuum gauge structure and physical intermediate-state sum cannot
simply be assigned to the finite-rank candidate.

## Conventions and admitted inputs

Use $\hbar=1$, $\kappa^2=32\pi G$, metric signature $(+---)$ and four
outgoing null momenta satisfying $\sum_i k_i=0$. Define

$$
s=(k_1+k_2)^2,\qquad t=(k_2+k_3)^2,\qquad
u=(k_1+k_3)^2,\qquad s+t+u=0.
$$

Our spinor convention is $\langle ij\rangle[ji]=2k_i\cdot k_j$.
The loop dimension is $D=4-2\epsilon$; the external states and the count of
physical internal polarizations retain their four-dimensional values.
Split a loop momentum into a four-dimensional part and an extra-dimensional
part, $L=\ell+\mu$, with the notation $\mu^2\geq0$ for the squared effective
four-dimensional mass on a cut.

In this scheme the continuum supersymmetry identities express the all-plus
loop contribution as the signed number of physical states times the
single-real-scalar contribution. Pure Einstein gravity has two such
bosonic states. This identity is an **input**, not a result established
for an angular truncation. See
[Bern, Dixon, Perelstein and Rozowsky, PDF p. 17, Eq. (4.2)](https://www.alphaxiv.org/pdf/hep-th/9811140?page=17).

The published, coupling-stripped scalar/gluon trees needed on a cut are

$$
A_1=i\frac{\mu^2[12]}{\langle12\rangle D_1},\qquad
A_2=-i\frac{\mu^2[12]}{\langle12\rangle}
\left(\frac1{D_1}+\frac1{D_2}\right),
$$

where $D_i=(\ell-k_i)^2-\mu^2$. Their KLT product gives

$$
-iD_1 A_1 A_2
=-i\left(\frac{\mu^2[12]}{\langle12\rangle}\right)^2
\left(\frac1{D_1}+\frac1{D_2}\right).
$$

The gravitational coupling multiplying this four-point tree is
$(\kappa/2)^2$.
These trees and their relative signs are taken from
[the same paper, PDF p. 20, Eqs. (4.9)–(4.10)](https://www.alphaxiv.org/pdf/hep-th/9811140?page=20).
The code checks their algebraic KLT relation. It does not describe that
check as a fresh derivation of the input trees.

## Reconstructing the cuts without double counting

On the $12|34$ cut, the two denominator choices on either side produce four
terms. Identifying box orderings related by a cyclic permutation or reversal
leaves two box classes, each with coefficient two. The other two channel
cuts provide the third class and repeat already identified ones.

The implementation enumerates these choices. A box seen by two different
cuts must have the same coefficient in both; the coefficient is recovered
once, not summed over the cuts. The resulting three box classes give

$$
\mathcal M_4^{(1)}(++++)
=2\left(\frac{\kappa}{2}\right)^4
\left(\frac{[12][34]}{\langle12\rangle\langle34\rangle}\right)^2
\left[I_4(s,t;\mu^8)+I_4(s,u;\mu^8)+I_4(t,u;\mu^8)\right].
$$

In each two-particle cut the two physical states cancel the
identical-particle factor of one half. The coefficient two of each box
then comes from equivalent denominator pairings. No additional state-count
factor should be applied to this expression.
The assumed dimensional cut construction and the above box representation
are the continuum Einstein prescription of
[Bern et al., PDF pp. 19–20, Eqs. (4.8)–(4.13)](https://www.alphaxiv.org/pdf/hep-th/9811140?page=19).
In this prescription scaleless integrals vanish; no extra higher-curvature
interaction has been introduced. A finite-volume or finite-rank regulator
can change those premises.

## Evaluating the rational quantum term

With the loop normalization $d^D L/(2\pi)^D$ and Feynman boundary values,
dimension shifting gives

$$
I_4(s,t;\mu^8)
=-\epsilon(1-\epsilon)(2-\epsilon)(3-\epsilon)(4\pi)^4
I_4^{D=12-2\epsilon}(s,t).
$$

After the loop momentum integral,

$$
I_4(s,t;\mu^8)
=-\frac{i}{(4\pi)^{2-\epsilon}}\,
\epsilon(1-\epsilon)(2-\epsilon)(3-\epsilon)\Gamma(\epsilon-2)
\int_{\Delta_3}[-s x_1x_3-t x_2x_4-i0]^{2-\epsilon}.
$$

Here $\Delta_3$ is $x_i\geq0$, $\sum_i x_i=1$, with the delta-function
simplex measure. The integral and dimension-shift conventions are recorded
in [Bern et al., PDF pp. 47 and 49, Eqs. (D.1)–(D.3), (D.9)](https://www.alphaxiv.org/pdf/hep-th/9811140?page=49).

The apparent zero from the prefactor cancels the ultraviolet pole of the
dimension-shifted integral. The recurrence for the Gamma function yields
the exact simplification

$$
-\epsilon(1-\epsilon)(2-\epsilon)(3-\epsilon)\Gamma(\epsilon-2)
=-(3-\epsilon)\Gamma(1+\epsilon)\longrightarrow-3.
$$

For nonnegative integer exponents the simplex moments are

$$
\int_{\Delta_3}\prod_{i=1}^4 x_i^{a_i}
=\frac{\prod_i a_i!}{(\sum_i a_i+3)!}.
$$

Consequently,

$$
\int_{\Delta_3}x_1^2x_3^2=\frac1{1260},\qquad
\int_{\Delta_3}x_1x_2x_3x_4=\frac1{5040}.
$$

Expanding the numerator and integrating each monomial therefore gives

$$
\boxed{I_4(s,t;\mu^8)
=-\frac{i}{(4\pi)^2}\frac{2s^2+st+2t^2}{840}
+O(\epsilon).}
$$

The program obtains these coefficients from factorial moments and the
Gamma recurrence. The number $840$ is not an input to its integration
routine. It then sums the box classes reconstructed from the cuts and
uses $s+t+u=0$.

The derived amplitude is

$$
\boxed{
\mathcal M_4^{(1)}(++++)
=-\frac{i}{(4\pi)^2}
\left(\frac{\kappa}{2}\right)^4
\left(\frac{s t}
{\langle12\rangle\langle23\rangle\langle34\rangle\langle41\rangle}\right)^2
\frac{s^2+t^2+u^2}{120}.
}
$$

It agrees analytically with the known reference,
[Bern et al., PDF p. 20, Eq. (4.14)](https://www.alphaxiv.org/pdf/hep-th/9811140?page=20).
The amplitude has mass dimension zero: the fourth power of the
dimensionful coupling cancels the fourth power of momentum.

Taking $\mu=0$ in the trees before integrating makes both cut trees vanish
and loses this nonzero answer. Thus a test using strictly four-dimensional
cuts can miss the entire rational contribution. This is one reason that
matching the earlier tree amplitude does not determine the quantum result.

## Independent numerical checks to run

The implementation is
[one_loop_all_plus.py](../../verification/quantum_bridge/one_loop_all_plus.py).
The fixed experiment entry point dispatches to it through the committed
configuration. No third-party Python packages, GPU or remote compute are
required.

The planned checks are:

1. Exact rational simplex moments, Gamma-polynomial cancellation and the
   full on-shell polynomial identity, plus box multiplicities reconstructed
   independently in all three channels.
2. Six rational $(s,t)$ pairs, evaluated with 60-digit arithmetic by
   independent Gauss cubature. A cube-to-simplex transformation carries
   Jacobian $(1-u)^2(1-v)$. Four and six nodes per axis both integrate the
   degree-six polynomial exactly apart from arithmetic rounding. They must
   agree with the moment calculation within $10^{-45}$.
3. Twelve nonexceptional real scattering configurations, spanning three
   energies, four scattering angles, changing azimuths and three Lorentz
   frames. The loop result is compared with the published spinor expression
   in both magnitude and phase, with relative tolerance $10^{-10}$.
4. All 24 external permutations in every configuration, spin-two
   little-group weights, and the fourth-power energy scaling at fixed
   coupling.
5. Controlled failures: deleting the extra-dimensional numerator loses the
   entire amplitude; retaining only one real scalar state gives half the
   pure-gravity result.

These are finite, deterministic coverage checks. They do not constitute
an off-shell Ward-identity proof or evaluations of the other helicity
channels. The numerical result will be reported from the tracked output,
including its immutable source commit.

## What a finite-$N$ calculation still must supply

The reference has no angular rank parameter. It does not tell us that an
angular cutoff commutes with the continuum state-count identity or the
dimensional limit. A valid finite-$N$ comparison must instead:

- derive the physical propagating and constrained states from its own
  completed constraints and measure;
- derive the Gaunt vertices, including induced interactions or determinants
  required by any declared completion;
- retain the full internal angular sums and the remaining radial/time loop
  integrations, with an explicit ultraviolet and infrared prescription;
- compute any rational terms in that same prescription;
- relate its finite-region observable to scattering at fixed physical
  kinematics, then control the limit and its symmetry identities.

Projecting only the external states of the boxed continuum answer would
not perform those steps. Likewise, the two-state scalar identity cannot
silently replace a finite-$N$ graviton/constraint calculation when that
identity has not been established for the regulator.

The corner construction and its positive transfer kernels remain the
separate, explicitly bounded-scope results in
[corner-quantization-and-transfer.md](corner-quantization-and-transfer.md).
The three principal quantum deliverables in [GOAL.md](../../GOAL.md)
remain pending.
