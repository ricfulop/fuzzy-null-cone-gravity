# Stationary germ to nonlinear obstruction linkage

**Preparation only.** This note records the conditional proof chain that would turn a future positive direct-state stationary germ certificate into a local Paper I nonlinear obstruction. No physical-radius minor, regular interval, discarded response, or solution germ is certified here; the live direct-state checker was still in its interval search when this note was written.

## Conditional manuscript claim

The current `canonical_oblate.tex` statement `prop:oblate_rank` is only a primary/secondary rank theorem on the exact static oblate Minkowski family, and the adjacent higher-preservation paragraph concerns the older round witness with $K_{ij}=\kappa n_i n_j$. A responding stationary germ should therefore be a separate $N=3$ local Lagrangian obstruction proposition, not an extension of the rank-count proposition and not a final Dirac or degree-of-freedom claim.

A safe conditional claim is:

> At $N=3$ and $\epsilon=1/100$, suppose the stationary direct-state certificate exhibits an interval-local reflection-even polar retained linearized germ $e$ with
> $$
> D[e](r_0)\ne0,\qquad
> D[e]=6\mathscr E_{X,4}^{(1)}[e]-\mathscr E_{C_{\rm p},4}^{(1)}[e],
> $$
> at a certified regular radius. Then there are compactly supported profiles $f,w$ such that
> $$
> c_1=e+g_f,\qquad
> g_f:\ C_{{\rm a},2}=-\frac{2\rho f}{3},\quad
> V_{{\rm a},2}=-\frac{f'}3,\quad W=0,
> $$
> solves the retained linearized equations on the shell but cannot be the spacetime tangent of any retained bulk solution family on a common spacetime neighborhood that is $C^2$ in amplitude in a topology controlling the Euler jets.

This is an interval-local Lagrangian nonintegrability statement for a specified mixed tangent. It would not certify a global radial boundary domain, identify the accidental null shift as a nonlinear gauge generator, cover all representatives after adding linearized gauge directions, prove an $N=4$ analogue, or produce a final Dirac classification.

## Already established inputs

The nonlinear null preparation gives the required Noether reduction for the descriptor $\xi=w(r)x\partial_\phi$, which is retained at the oblate background. For
$$
c(\eta)=c_0+\eta c_1+\frac{\eta^2}{2}c_2+o(\eta^2)
$$
near the exact continuum oblate background,
$$
\partial_t C_W[\xi](\eta)
=-\eta^2E_{\rm full}^{(1)}[c_1]
  \left[Q_cDF_\xi(c_0)c_1\right]+o(\eta^2).
$$
Every occurrence of the second amplitude correction $c_2$, including time-dependent second-order temporal data, is multiplied away by $E_{\rm full}(c_0)=0$ or $Q_cF_\xi(c_0)=0$. Thus a nonzero
$$
B_\xi=E_{\rm full}^{(1)}[c_1]\left[Q_cDF_\xi(c_0)c_1\right]
$$
cannot be cancelled by choosing a time-dependent second-order correction. The preparation also records the necessary scope: the first-order input must be a retained linearized solution germ, not just a point on the linearized secondary surface. See [the Noether derivation](../oblate-nonlinear-null-preparation-2026-11-27/README.md#2-noether-derivation-and-perturbative-admissibility).

The passed map/control check certified the descriptor and Euler machinery used in $B_\xi$. Its `expected_discarded_generator()` has the eight discarded channels, in particular
$$
X_4=-\frac{12}{35}w C_{{\rm a},2},\qquad
C_{{\rm p},4}=\frac{2}{35}w C_{{\rm a},2},
$$
with the additional axial tensor, radial-vector, and temporal-vector channels. Its `euler_pair_density()` uses exactly one factor
$$
\Sigma=r^2+\epsilon x^2
$$
and the full symmetric index contraction, with no extra $r^2$ radial measure. The archived `map-checks.json` explicitly leaves `admissible_non_gauge_seed_supplied`, `quadratic_source_evaluated_on_admissible_non_gauge_seed`, `nonlinear_obstruction_established`, and `full_nonlinear_Dirac_classification_established` false. See the [map-check report](../oblate-null-map-check-2026-11-27/README.md) and its [result JSON](../oblate-null-map-check-2026-11-27/map-checks.json).

For the odd coordinate leg
$$
g_f:\quad C_{{\rm a},2}=-\frac{2\rho f}{3},\qquad
V_{{\rm a},2}=-\frac{f'}3,\qquad W=0,
$$
the checked map gives
$$
Q_cDF_\xi g_f
=\frac{8\rho wf}{35}e_{X,4}
-\frac{4\rho wf}{105}e_{C_{\rm p},4}.
$$
For a stationary polar germ $e$, reflection kills pure-parity pairings and the odd coordinate leg has zero unrestricted linearized Einstein residual. The compact-test pairing is therefore
$$
B_\xi[e+g_f]
=\frac4{105}\int dr\,\rho\,wf\,D[e],
\qquad
D[e]=6\mathscr E_{X,4}^{(1)}[e]
-\mathscr E_{C_{\rm p},4}^{(1)}[e],
$$
with no additional radial volume factor beyond the one already built into the coefficient Euler densities. See the [stationary mixed-parity seed preparation](../oblate-stationary-seed-preparation-2026-11-27/README.md#2-a-known-odd-coordinate-mode-supplies-the-other-parity).

## What the future certificate must supply

The direct-state certificate must supply a real stationary polar germ in the exact original $N=3$, $m=0$, reflection-even, $X_0=0$, $\epsilon=1/100$, zero-remaining-control problem. It must bind the pinned original `stationary-operators.json` to the terminal state snapshot with
$$
C\in\mathbb Q(r,T)^{6\times14},\quad
K\in\mathbb Q(r,T)^{14\times14},\quad
B\in\mathbb Q(r,T)^{14\times0},\quad
M\in\mathbb Q(r,T)^{14\times8}.
$$
Here $T(r)=\arctan(1/(10r))$ and $\mathcal D=\partial_r-10(100r^2+1)^{-1}\partial_T$ is the total radial derivative on coefficient matrices. The core exact identities are
$$
CM=0,\qquad JM=I_8,\qquad (\mathcal DC+CK)M=0.
$$
At the same positive physical radius as the response witness, the certificate must prove all proof denominators nonzero and exhibit a nonzero rank-six minor of the original $6\times14$ constraint matrix $C$.

For the physical lift back to the stationary coefficient equations, it must also check
$$
(P_vK-P_w)M=0
$$
and all 22 restored retained Euler rows, including the old $X_0$ equation:
$$
\left(E_0P_v+E_1P_w+E_2P_wK\right)M=0.
$$
Finally, it must lift the discarded row independently to a fourteen-state row $d$ and certify
$$
(dM)_j(r_0)\ne0
$$
at that same regular radius for some unit initial vector $e_j$. The [direct-state preparation](../oblate-stationary-state-germ-preparation-2026-11-28/README.md) explains why these obligations produce an analytic constrained ODE germ: $JM=I_8$ and rank six give $\operatorname{im}M=\ker C$, while $S=(\mathcal DC+CK)$ annihilates $\ker C$ and so satisfies $S=LC$ locally for an analytic multiplier $L$.

## How a positive response is used

Once a certificate gives $D[e](r_0)\ne0$, continuity lets us shrink to a shell where $D[e]$ has fixed sign. Choosing real smooth compactly supported radial profiles $f,w$ inside that shell with $fw\ge0$ and positive overlap gives
$$
\frac4{105}\int dr\,\rho\,wf\,D[e]\ne0,
$$
because $\rho=r^2+\epsilon>0$. The mixed tangent $e+g_f$ then violates the Noether condition above, excluding sufficiently smooth exact retained spacetime families with that specified tangent.

If the certified row $dM$ is instead identically zero, that proves only that this forward ODE subfamily is silent. A universal hidden-identity claim would still need a converse showing that every stationary germ in the stated sector and gauge enters the certified parametrization.
