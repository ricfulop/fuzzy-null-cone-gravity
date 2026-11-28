# Full fifth-order Hamiltonian: evaluated directions and validation

This is a calculation protocol for the full classical $H_5$ functional
of the declared reduction-first candidate. It builds on the verified higher
geometry and first Green responses, source `2a59f1b`, reviewed at `ef48f21`.
After a preserved spatial-refinement failure, the finer child passed all
specified checks for two directions of the full functional. These are
directional coefficients, not a reconstructed vertex tensor or loop amplitude.

## Reviewed result

Run `116088d3-0ed7-4fb7-95c0-fdeafcdc9072`, source
`502f75d96fe8c838469570eb7ae7f41eebf3e742`, completed in 400.89 seconds.
Its [complete output](data/fifth-order-results.json) was copied byte for byte
from the isolated run. The effective configuration and final result are
also present in the retained raw log.

| Raw field direction | Full $H_5$ value | Maximum scaled spatial/refinement diagnostic |
|---|---:|---:|
| `electric_axis` | $-1641.1412606219328$ | $1.82585\times10^{-10}$ |
| `noncommuting_labels` | $-76826.87759611206$ | $7.57660\times10^{-10}$ |

These use width one, $16\pi G=1$, and the declared raw TT directions;
$S=\mathsf G_0^{-1}p$ is not the normalized canonical momentum coordinate.
The last column is a checked diagnostic, not a rigorous numerical error bar.
Both directions satisfy the unchanged $10^{-6}$ spatial acceptance gate,
including every family refinement, independent spatial comparison, York
and scalar integration-by-parts comparison, reference-radius tail bound,
and cancellation margin. All 42 recorded `passed`/`converged` flags are true.

The shared 20-point electric-axis family integrals reproduce the failed
parent exactly, as do both collapsed total values. The improvement is in
the required 20-to-24-point refinement. No tolerance was relaxed. The two
integration methods still share the first HEAT fields and moment evaluator;
their independent operations and limitations are specified below.

## Functional and integration method

The [independent pairing audit](h5-pairing-audit.md) derives
$$
\begin{aligned}
H_5={}&-\int\mathcal L_5-4\int b_2\cdot W
-2\int T_0:T_1+2\int\operatorname{tr}(QT_0^2)\\
&+\int R_3(v-u)+2\int(S:T_0)(7u+v)
-8\int(\nabla v)^TQ\nabla u.
\end{aligned}
$$
Here $W=\mathcal A_0^{-1}b_1$, $T_n=L_nW$,
$u=\mathcal C(s_2-R_2)$, $v=\mathcal C(R_2+7s_2)$ and $S=\mathsf G_0^{-1}p$.
All sectors $q^5,q^3p^2,qp^4$ are retained. Each nonlocal term has a
polynomial-Gaussian multiplier and at most two first-inverse HEAT factors.
This property is specific to $H_5$ and does not eliminate the nested inverses
or possible algebraic integrand tails at sixth order.

The primary integration exactly contracts Cartesian angular monomials using
$$
\langle n_x^{a}n_y^{b}n_z^{c}\rangle_{S^2}
=\frac{(a-1)!!(b-1)!!(c-1)!!}{(a+b+c+1)!!}
$$
for even powers and zero otherwise. It merges signed rational coefficients
before radial integration. The Gaussian multiplier exponent remains separate
from the rates inside each HEAT factor. Three radial Gauss orders, termwise
refinement, absolute-bin sums and analytic Gaussian tail inequalities are
recorded. Those diagnostics are not rigorous quadrature remainder certificates.

## Independent acceptance checks

1. Small rational dictionaries test anisotropic HEAT contractions, rescuing
   parity in triple products, both convolution branches, equal-rate index
   exchange, exact negative cancellation and split-versus-direct merges.
   Four numerical controls then compare analytic Gaussian integrals with
   zero, one and two Green factors, including an exactly odd integral.
2. Separate $H_3,H_4$ families must reproduce the previous exact Fourier
   calculation in physical spatial units. The complete quartic result is
   not obtained by adding the cubic family.
3. The York operator coefficients $\mathcal A_1,\mathcal A_2$ are compared
   with independent 90-digit finite metric evaluations. The contravariant
   test-vector jet is held fixed while the metric varies. Flat recovery,
   conformal tracelessness, determinant, connection trace and exponential
   jet convergence also gate acceptance.
4. The axial electric and a noncommuting mixed-parity direction each supply
   six $H_5$ families and their full sum. A separate full-sphere spatial
   quadrature integrates the strong York form $2W\cdot\mathcal A_1W$.
   It also compares $-8(\nabla v)^TQ\nabla u$ with $8vQ:\nabla^2u$;
   flat transversality makes the latter the strong $-8vD_1u$ form.
5. Every family and the full reference integral must refine. Bounds are
   evaluated at the reference's actual radius, separately from the primary
   integrator's automatically chosen radius. A factor of 100 multiplies
   the machine-epsilon times absolute-bin-sum indicator; its scaled value
   must remain below one tenth of the $10^{-6}$ spatial comparison tolerance.
   A large cancellation indicator requires additional precision work rather
   than an unqualified accuracy claim.

Every retained free input is total angular momentum $j=2$. Rotation-equivariant
differentiation, inverse operators and contractions therefore make each
homogeneous degree-five scalar family have $j\le10$. For this specific basis,
$n_\phi>10$ and $2n_\mu-1\ge10$ suffice for exact angular integration before
floating evaluation errors. All three reference grids satisfy this bound;
the larger grids additionally check those numerical errors. This is not a
generic angular bound for arbitrary fields accepted by the assembly API.

The two integrations share the exact source and HEAT-field construction and
the elementary heat-moment evaluator. Their angular contractions, organization
of cancellation and strong/weak operator pairings differ. Prior first-inverse
PDE and independent moment tests support those shared inputs, but the spatial
comparison is not a second independent solution of all constraint equations.
The factor-100 roundoff indicator is a diagnostic margin, not a certified
error bound. A separate read-only model audit checked the signs, rates,
normalizations, API alignment and asserted acceptance gates before launch;
it did not run the calculation or supply an external peer review.

## Execution and limitations

The fixed tracked command is unchanged; the committed selector is
`full_fifth_order`. The run is one local process with a 1,800-second wall cap.
All tolerances, work limits and quadrature grids are recorded before launch.
Results preserve the exact source commit and failed comparisons if any.

### First attempt: failed spatial refinement

Run `b7057ab6-0465-4dd9-8e1d-3e01db6ca97a`, source `bf88c08`, used
spatial radial orders 12/16/20 and angular grids 8x16, 10x24, 12x32.
It completed the exact/numerical backend controls, both directions' H3/H4
calibrations and York coefficient checks, and both collapsed H5 integrals.
It stopped after the first direction's spatial reference, before testing
the second reference or writing an overall passing result.

The failed quantity was the electric-axis scalar-gradient strong-form
16-to-20-point refinement: its scaled change was
$1.4682497276062995\times10^{-6}$, exceeding the declared $10^{-6}$ gate.
The weak scalar-gradient refinement also exceeded that gate at
$1.1111768757696966\times10^{-6}$. At the final 20-point level, the total
agreed with the collapsed calculation to $8.93\times10^{-13}$ scaled,
and the individual families agreed within $1.31\times10^{-10}$ scaled.
Those agreements do not waive the failed refinement requirement.

The [partial JSON records](data/full-fifth-order-first-attempt.json) are
copied from the immutable raw log and explicitly marked unaccepted. The
original source and failed node remain unchanged. The integrator's
`accuracy_certified=false` qualification means that its convergence
diagnostics are not rigorous error certificates; it is distinct from
the actual failed assertion. All H3/H4 calibration acceptance checks passed.

The successful child retained the same fields, Gaussian pairing options,
radius 24, panel width 3, wall cap and acceptance tolerances. It changed only the
spatial grids to radial 16/20/24 and angular 10x24, 12x32, 14x40, and
reports the failure metric before asserting acceptance. Both directions
were rerun. The shared 20-point result supplied a reproducibility
cross-check, and the new 20-to-24-point refinement passed the original
threshold. No revised acceptance threshold was selected from the failed data.

Two directional evaluations do not reconstruct the rank-five vertex tensor
or complete its contracted loop contribution. Complete mixed $H_6$, the
four-point diagram sum, the connected observable, physical angular/radial
limits, microscopic corner dynamics and the relation to the historical
simulation remain open. No equivalence of the two archived regulators is
assumed.
