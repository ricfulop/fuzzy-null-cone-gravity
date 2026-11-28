# Higher-order geometry and first Green responses: verified ingredient checks

Run `6ffbc441-ac18-4c69-90e3-e3d3827feb85` completed all the specified checks
at source `2a59f1b60f630bd2b4e4d922dc693a7810ab2608` in 39.56 seconds.
The complete [unaltered run JSON](data/higher-geometry-results.json) records
the configuration, comparisons and tolerances. This supplies required
ingredients for the [contracted four-point calculation](contracted-four-point-implementation.md)
after the complete quartic and finite-mode two-point checkpoint `6bd4090`.
The nonlinear regulator remains the declared reduction-first canonical candidate.

## Observed result

- Both exact geometry directions matched their archived low-order formulas.
  Independent finite-metric extraction used 17 distinct amplitudes at each of
  four points, with maximum scaled coefficient error $1.33144\times10^{-15}$
  and maximum extrapolation change $8.38866\times10^{-14}$.
- All 27 Green checks passed, including four actual scalar sources and two
  vector sources. Their exact PDE residuals vanished; the largest separately
  evaluated numerical PDE error was $7.94068\times10^{-16}$.
- Ordinary heat moments agreed with independent Decimal quadrature to relative
  error $7.73258\times10^{-16}$. The deliberately large-degree logarithmic
  stress cases had a larger error, $1.16609\times10^{-11}$, within their
  separate preregistered $2\times10^{-10}$ tolerance. They are not claimed
  to have ordinary-degree machine precision.
- The pure-momentum sixth-order control was $32618.557080168644$ in the
  stated raw axial-profile units. Its 96/128-point integral refinement changed
  by scaled error $2.67740\times10^{-16}$; the full pointwise radial-inverse
  comparison reached $3.08826\times10^{-14}$. The conservative omitted
  sixth-order tail bound was $2.40415\times10^{-47}$.

These results validate the selected ingredients and the previously derived
pure-momentum control. They do not evaluate the full mixed $H_5,H_6$ or the
gravitational four-point loop. The protocol and pre-launch review follow.

## Independent local geometry check

The exact sparse-polynomial module constructs the ordinary metric exponential,
Christoffels through fifth order, pointwise curvature through fourth order,
integrated Gamma--Gamma densities through sixth order, and the covariant York
source through cubic order in the configuration field. Every generated angular
component is retained. The input directions are the axial electric case and a
noncommuting mixed-parity case in the fixed ten-mode basis.

The archived second-to-fourth-order formulas provide exact rational comparisons.
A separate finite matrix-exponential evaluator then constructs the full metric,
inverse, Ricci tensor and tangent divergence at two Cartesian points per
direction using 90-digit Decimal arithmetic. Nine centered amplitude values
extract coefficients through degree eight; three step sizes and Richardson
extrapolation test the new curvature, York and raw-Gram fourth-order terms.
The extraction never subtracts the expected Taylor coefficients. Both the
scaled error and the change between extrapolations must be below $10^{-12}$.

The raw-Gram check includes the pointwise coefficient
$\operatorname{tr}([Q,[Q,S]]^2)/360$, which enters the sixth-order Hamiltonian.
Pointwise curvature is kept distinct from its integrated Gamma--Gamma density.

## First Green responses

The scalar and York inverses act on complete polynomial-Gaussian sources on
$\mathbb R^3$. Exact Hermite reconstruction and analytic derivatives retain
their algebraic tails. Geometry metadata uses $e^{-nr^2/2}$, whereas the Green
API uses $e^{-ar^2}$; the protocol explicitly sets $a=n/2$.

Actual second-order curvature, momentum-square and first York sources must
pass exact recurrence-reduced PDE identities. Evaluated residuals, independent
moment quadrature, regular values near the origin and Gaussian monopole tails
provide numerical checks of the evaluator. Numerical tolerances are recorded
with each check; exact residuals must vanish exactly.

The axial momentum source is also evaluated through the independent archived
radial Legendre inverse. Its Coulomb integral is compared with the exact Fourier
normalization, and $77\int s_2(\mathcal C s_2)^2$ is evaluated at two radial
quadratures. This rechecks the already derived pure-momentum sixth-order sector.
The final integral is cut off at radius twelve while each potential retains
its full-space tail; positivity bounds the omitted sixth-order contribution
by $77 B^2\int_{r>12}s_2$, using the established uniform potential bound $B$.
This control supplies no mixed or configuration sixth-order coefficient.

## Execution and scope

The fixed tracked command remains `bash paper2/verification/quantum_bridge/run.sh`.
The calculation selector is `higher_geometry_and_green`; the committed config
declares a 1,800-second local wall-clock limit and no GPU or remote work.
Completed stage records and the final JSON preserve the source commit and run ID.

The pass validates these selected higher-order ingredients. Nested nonlinear
Green sources, the full $H_5,H_6$ functional, the contracted gravitational
four-point integral and connected scattering conversion would remain pending.
No physical large-rank limit or microscopic corner multiplicity follows from
these checks.

## Pre-launch review

A separate read-only model audit checked the local geometry index conventions,
Gaussian rates, interpolation factorials, Richardson powers, Green constants,
quadrature normalization and omitted-tail inequalities. It identified two useful
acceptance gates that were added before launch: enforce the finite metric's
determinant/connection/jet-tail diagnostics, and explicitly compare the axial
Cartesian source with the independent positive-frame expression at each test
point. This is an internal mathematical/code review, not external peer review
or a substitute for the tracked calculation recorded above.
