# Kinetic completion: a verified cancellation with remaining constraints

## Result and provenance

The Schur correction cancels the previously measured second lapse derivative
in all five stated angular kinetic cases. The determinant factorization
also holds exactly through fourth order in the lapse perturbation. These
results validate this component of a proposed repair; they do not establish
a complete gravitational constraint algebra or a one-loop amplitude.

- Experiment: 62f737a6-3a58-458e-9320-5763aa298690.
- Run: d29c46d3-f648-426b-ab3f-b83a014c45c1, local CPU, exit code zero.
- Immutable source: fac157a2ec8f8ea29cbcca67cbbe9034635700a1.
- Source: [kinetic_completion.py](../../verification/quantum_bridge/kinetic_completion.py).
- Machine-readable result: [kinetic-completion-results.json](data/kinetic-completion-results.json).
- Raw log: lines 4 (configuration), 11--35 (exact blocks), 39--99
  (independent finite-amplitude checks), and 102 (final summary).
- Calculation runtime reported by the script: 1.1221 seconds.

The starting data and normalization are those of the
[lapse-constraint calculation](lapse-constraint-on-vacuum-data.md):
flat spatial geometry, constant Cartesian \(K_{ij}=k\delta_{iz}\delta_{jz}\),
and \(h_{tt}=-1-2\epsilon P_\ell(\cos\theta)\).
Both continuum vacuum constraints vanish pointwise. The tables set \(k=1\)
and suppress the common positive radial and gravitational normalization.

## What was changed

Write the angular kinetic matrix at auxiliary cutoff \(J\) in retained
\(P=P_L\) and auxiliary \(Q=P_J-P_L\) blocks. Instead of using only
\(A_{PP}\), form

\[
 A_{\mathrm{eff}}=A_{PP}-A_{PQ}A_{QQ}^{-1}A_{QP}.
\]

The exact block inverse identity gives
\((A^{-1})_{PP}=A_{\mathrm{eff}}^{-1}\). Thus the corrected retained-momentum
Hamiltonian includes the virtual kinetic response of the auxiliary block.
This modifies the original cutoff Lagrangian. It is not an identity
equating the original \(P_L\) Legendre transform with the corrected one.

For these data the missing second-order Hamiltonian contribution is
\(\epsilon^2\langle Q(uK),D Q(uK)\rangle\), where
\(u=P_\ell\) and \(D\) is the DeWitt bilinear form. Keeping \(J\geq\ell+2\)
captures that particular leakage. It does not capture every order in
the nonpolynomial lapse.

| \(L=N-1\) | \(\ell\) | \(J\) | Original \(H''(0)\) | Corrected \(H''(0)\) | Corrected coefficient of \(\epsilon^4\) |
|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 4 | \(-4/15\) | 0 | 0 |
| 2 | 2 | 4 | \(-24/175\) | 0 | \(-6/2695\) |
| 3 | 3 | 5 | \(-40/441\) | 0 | \(890/3006003\) |
| 4 | 4 | 6 | \(-20/297\) | 0 | \(116200/478690641\) |
| 6 | 6 | 8 | \(-112/2535\) | 0 | \(5155290/1128581148409\) |

These are exact rational coefficients, not fitted floating-point zeros.
The nonzero fourth-order entries are a material limit on the cancellation.
Since the Hamiltonian is quadratic in the momentum, a term
\(p^2\epsilon^4\) is sixth order in weak fields. Such orders cannot be
discarded without justification in a background-Hessian calculation of a
one-loop four-point quantity.

At finite amplitudes \(\epsilon=0.1,0.2,0.3\), the script separately
evaluated the full enlarged kinetic matrix at \(J=2,4,\ldots,14\), using
60-digit arithmetic and two Gauss--Legendre quadratures (96 and 128 nodes).
For \(\epsilon=0.3\), the remaining \(|H|\) decreased from
\(5.6813176568\times10^{-3}\) at \(J=2\) to
\(2.7367517922\times10^{-14}\) at \(J=14\); the two quadratures agreed to
within \(1.5\times10^{-59}\) along that sequence. The continuum
phase-space restriction has \(H=0\) on the specified data.

This is convergence of one angular kinetic observable as an auxiliary
space is enlarged. It is not the requested physical large-\(N\) scattering
limit and does not establish a uniform bound over gravitational fields.

## Determinants and the missing amplitude

The calculation also verified, for each of the five cases and three
axisymmetric spin blocks,

\[
 \log\det A=\log\det A_{QQ}+\log\det A_{\mathrm{eff}}
\]

as an exact series through \(\epsilon^4\). Dropping either the auxiliary
determinant or the mixed Schur terms changes the fourth-order coefficient.
The result therefore identifies terms that a matched quantum calculation
must keep; there is no license to replace the full determinant by the
corrected low block alone.

An independently organized series calculation checked the cyclic identity

\[
 [\epsilon^4]\operatorname{Tr}\log
 \bigl(1+\epsilon X_1+\epsilon^2X_2+
       \epsilon^3X_3+\epsilon^4X_4+\cdots\bigr)
 =
 \operatorname{Tr}
 \left(X_4-X_1X_3-\tfrac12X_2^2+
 X_1^2X_2-\tfrac14X_1^4\right).
\]

Only the \(m=0\) angular kinetic blocks were evaluated. A background with
axial symmetry still has fluctuating modes with other \(m\); they cannot
be omitted from a gravitational loop. The calculation includes neither
frequency/radial propagators nor the potential, ghosts, the physical
constraint measure, boundary contributions, renormalization, external
state reduction, or rational amplitude terms. A kinetic log determinant
is therefore not a one-loop four-graviton amplitude.

## Consequence for the next calculation

The cancellation motivates a separately declared Hamiltonian-first
candidate: restrict canonical metric and momentum variables, and restrict
the continuum Hamiltonian and spatial-diffeomorphism functionals using
complete Gaunt contractions. Lapse linearity is built into this candidate.
It changes the regulator and is not assumed equivalent to Paper I's
metric-coefficient Lagrangian restriction.

The decisive next local question is whether the candidate's constraint
brackets vanish on its constraint surface. In particular, evaluate the
mixed bracket of a radial spatial deformation and a lapse deformation
on the same vacuum data. A vanishing lapse Hessian alone does not answer
that question. The test should compare an exact projection calculation
with a direct Hamiltonian-flow evaluation and keep all potential terms.

For the full project, the [corner and reconstruction analysis](corners-and-reconstruction.md)
still supplies conditional continuum relations, not a finite-rank physical
Hilbert space. Finite-\(N\) one-loop scattering, a dynamical microscopic
corner-state count, and reconstruction from one physical measure all
remain pending.
