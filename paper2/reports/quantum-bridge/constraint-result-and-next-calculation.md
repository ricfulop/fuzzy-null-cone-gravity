# Constraint result and the first constructive repair

## Verified result

The tracked baseline ran source commit
`15493f83615d64308b9e8e81ad7ea32f8fc658a6` as run
`82f9a0a7-3fd7-432b-9913-27738446ad38`. Its full log identifies Python
3.13.7, the effective configuration, exact results for 11 cutoffs, independent
numerical results for nine cases, and a final success summary. The measured
calculation time was 1.186 seconds, excluding process setup. The preserved
[result data](data/lapse-constraint-results.json) give the exact arithmetic
and both quadrature sequences.

The [analytical prediction](lapse-constraint-on-vacuum-data.md) was confirmed:
on flat spatial data with constant rank-one Cartesian extrinsic curvature,
both continuum vacuum constraints vanish, whereas the finite-rank lapse
Hamiltonian Hessian has a nonzero axisymmetric block. At N=3 its two nonzero
entries for u=P_1 and u=P_2 are respectively -4/15 and -24/175, in the stated
normalization. The largest absolute discrepancy of the independent
finite-amplitude calculation from its exact target was
1.1403568037217492e-14. Changing quadrature from 96 to 128 nodes changed the
extrapolated curvatures by less than 4.7e-52 at 60-digit precision.

**Interpretation:** the primary--secondary constraint structure of the
finite-rank pullback is not the one obtained by simply restricting continuum
Einstein constraints. A complete Dirac analysis or an explicit change of
regulated construction is necessary. This result does not count final
physical degrees of freedom and is not a no-go theorem for all finite-state
gravity proposals.

## An analytic cutoff-shell limit

The axisymmetric results admit a closed form, derived from harmonic
coefficients rather than fitted to the measurements. For u=P_ell and the
rank-one K used in the calculation, define

$$
c_\ell=\frac{(\ell+1)(\ell+2)}{(2\ell+1)(2\ell+3)}.
$$

The degree ell+2 part of uK has components, in the scalar/vector/tensor
polynomial representation of the derivation,

$$
\left(c_\ell P_{\ell+2},\,-c_\ell P_{\ell+2},\,
-\frac{c_\ell}{\ell+2}P'_{\ell+2},\,
\frac{c_\ell}{2(\ell+1)(\ell+2)}P''_{\ell+2}\right).
$$

The usual factors sin(theta) and sin(theta)^2 multiply the last two
components. These coefficients follow by comparing the highest powers of x
in P_ell and P_(ell+2). Lower degrees are retained for ell=L or L-1.
Using the associated-Legendre norms and the full DeWitt contraction gives

$$
\|Q(P_\ell K)\|^2
=\frac{(\ell+1)(\ell+2)}{(2\ell+1)^2(2\ell+3)},
\qquad\ell=L-1,L.
$$

After unit-L2 normalization of the lapse perturbation, its curvature is

$$
\widehat{\mathcal H}''_{\ell}
=-\frac{2(\ell+1)(\ell+2)}{(2\ell+1)(2\ell+3)}
\longrightarrow -\frac12.
$$

For ell<=L-2 it is zero at this order. The two omitted-shell contributions
are mutually orthogonal. This proves the rank-two diagonal axisymmetric
block at this particular phase-space point for every L>=2. The numerical
table through L=12 is a check, not the proof of the limit.

This is a nonuniform cutoff-shell statement. Fixed low lapse harmonics
eventually have no second-order defect, whereas a unit-normalized harmonic
whose degree grows with L retains one. This is not yet a statement about
fixed physical external momentum, a scattering amplitude, or the effect
of renormalization. The defect is proportional to k^2 and hence vanishes at
K=0; this makes importing the gauge structure from the flat-background
quadratic action particularly unsafe.

## The missing term can be identified

Let A be the full velocity Hessian, with retained and discarded blocks P,Q.
The Galerkin pullback sets discarded velocities to zero and uses A_PP.
Eliminating those variables in a Gaussian calculation instead gives

$$
A_{\rm eff}=A_{PP}-A_{PQ}A_{QQ}^{-1}A_{QP},\qquad
P A^{-1}P=A_{\rm eff}^{-1}.
$$

This identity assumes the specified blocks are invertible. The null/gauge
directions of the full gravitational problem require their own reduction;
they cannot be removed by this identity alone.

At the spherical spatial point used here, the kinetic A=P alpha^-1 D P
extends to all angular channels with the same DeWitt operator. Its first
off-diagonal term is -epsilon P u D Q. The missing second-order Schur term
therefore supplies, in the Hamiltonian, the correction

$$
\delta\mathcal H_N^{(4)}
=+\epsilon^2\langle Q(uK),D Q(uK)\rangle.
$$

The superscript counts two lapse perturbations and two momentum fields.
It cancels the verified -epsilon^2 leakage term in this kinetic restriction.
Its angular coefficient is fixed by the complete Gaunt contractions; no
fitting is involved. Restoring channels through L+2 is enough for this
second lapse derivative on these data. It is not enough for arbitrary
finite-amplitude nonlinearities or all orders of gravity.

This is a change of operation order and hence of the regulated construction,
not a claim that the original pullback already contains the correction.
General discarded metric modes have dynamics and spatial derivative terms.
Integrating them out also produces nonlocal interactions; introducing
instantaneous auxiliary velocities alone is not the full gravity repair.

## Why the determinant matters for the requested loop

The same block elimination obeys

$$
\det A=\det A_{QQ}\det A_{\rm eff}.
$$

Keeping a corrected classical interaction but omitting the discarded-block
determinant changes the quantum theory. Ghosts and the gravitational
constraint measure need analogous treatment once they are actually derived.
For a correctly defined invertible bosonic Hessian

$$
K[\varphi]=K_0+K_1[\varphi]+K_2[\varphi^2]
+K_3[\varphi^3]+K_4[\varphi^4]+\cdots,
\qquad X_n=K_0^{-1}K_n,
$$

the homogeneous fourth-order part of the one-loop logarithm is

$$
\frac{\hbar}{2}\operatorname{Tr}
\left[X_4-X_1X_3-\frac12X_2^2+X_1^2X_2-\frac14X_1^4\right].
$$

The Lorentzian convention has the corresponding factor of i, and complex
Grassmann ghost determinants have the opposite sign and their own
multiplicity. The formula is an algebraic expansion, not an evaluated
four-graviton amplitude. K_4 contains sixth-order action vertices; the
full scattering amplitude also needs one-particle-reducible contributions,
external-state normalization and the chosen regulator/counterterms.

The next tracked calculation will check the Schur correction and the
fourth-order determinant identity on the actual angular kinetic matrices.
This is a necessary component of a repair, not completion of the loop.

## A possible route beyond the kinetic example

Integrating out a gauge-complete field sector must preserve the gauge and
boundary identities. BV pushforward provides a formal framework for that
operation when its hypotheses hold. Cattaneo, Mnev and Reshetikhin prove the
finite-dimensional pushforward statement and the corresponding quantum
master-equation property on [pp. 12--13](https://www.alphaxiv.org/pdf/1507.01221v2?page=12).
They expressly limit the field-theory extension: “theorems that are proven
for measure-theoretic integrals have to be checked, model by model, on the
level of Feynman diagrams.”
[p. 11](https://www.alphaxiv.org/pdf/1507.01221v2?page=11)

Their explicit gluing theorem concerns the listed BF-like theories, not a
finished four-dimensional Einstein quantization or a positive finite
gravitational Hilbert space. Applying this construction here is a research
direction with remaining derivations. Its advantage is that the necessary
bulk interactions, determinants and boundary operator would be derived
together, rather than selected separately.
