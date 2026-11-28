# Full fifth-order Hamiltonian: assembly and validation protocol

This is a prospective calculation of the full classical $H_5$ functional
of the declared reduction-first candidate. It builds on the verified higher
geometry and first Green responses, source `2a59f1b`, reviewed at `ef48f21`.
No fifth-order numerical result is claimed here before the tracked run.

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

Two directional evaluations do not reconstruct the rank-five vertex tensor
or complete its contracted loop contribution. Complete mixed $H_6$, the
four-point diagram sum, the connected observable, physical angular/radial
limits, microscopic corner dynamics and the relation to the historical
simulation remain open. No equivalence of the two archived regulators is
assumed.
