# Local oblate obstruction integrated into Paper I

The new Proposition 4 gives a local Lagrangian nonintegrability result for the original coefficient-Galerkin bulk action at $N=3$, $\epsilon=1/100$. A certified stationary polar germ has nonzero discarded Euler response at $r_0=1$. Adding a retained odd continuum coordinate leg produces a retained linearized solution that cannot be the specified spacetime tangent of any exact retained solution family that is twice amplitude differentiable in a topology controlling the Euler jets. Arbitrary time-dependent second-order corrections cannot cancel the selected Noether compatibility functional.

[Reading PDF](paper_reading.pdf) · [Technical PDF](paper.pdf) · [Source-to-claim audit](source-to-claim-audit.json) · [Build validation](validation.json)

## Evidence and proof

The [accepted run archive](../oblate-stationary-state-germ-pass-2026-11-28/README.md) contains run `c925d4a2-1531-4fb7-96e0-a8a0f4baa5e1`, executed source `0034162b4410ecfc229c48daaa47fcf06da4eafe`, and certificate SHA-256 `f07a00573347ed4dd037c4395124675b0dca02b5ab355909fe09e8d787030ca8`. It supplies all 22 retained rows, the exact kernel/preservation/profile identities, the physical original rank-six minor, all 57 proof denominators and the response numerator at the same radius and arctangent bracket. The initial state is $M(1,\arctan(1/10))(0,1,0,0,0,0,0,0)^t$. Twenty manufactured controls passed. Independent certificate review `msg_aceca875` found no material proof gap; the [input-domain audit](../oblate-stationary-input-domain-audit-2026-11-28/README.md) supplies the analytic interpretation of the original coefficients and the forward candidate matrices.

The manuscript combines that new certificate with the previously reviewed [Noether derivation](../oblate-nonlinear-null-preparation-2026-11-27/README.md), [passed coefficient map](../oblate-null-map-check-2026-11-27/README.md), and [mixed-parity normalization](../oblate-stationary-seed-preparation-2026-11-27/README.md). With $D[e]=6\mathscr E_{X,4}^{(1)}[e]-\mathscr E_{C_{\mathrm p},4}^{(1)}[e]$, it obtains

$$B_\xi[e+g_f]=\frac4{105}\int_I dr\,\rho w_*fD[e]\ne0.$$

Continuity supplies a fixed-sign subinterval and smooth compact profiles with positive overlap. The nonzero functional contradicts the second-order Noether condition for an exact retained family. This is an analytic inference from the accepted germ and established maps, not a new nonlinear numerical invocation. The run certificate itself retains `nonlinear_consistency_established: false` without modification. Independent manuscript review `msg_b9f3143f` found no material mathematical or scope gap. Its notation and compact-support suggestions were applied.

## Manuscript and builds

The new subsection follows the existing canonical oblate analysis. The primary/secondary rank proposition is unchanged. The abstract, action overview and discussion state the local result with its fixed-sector scope. Both layouts share the same source. The technical PDF has 27 pages; the reading PDF has 45. Both were built twice with `pdflatex`, all four processes exited zero, and the final logs have no undefined references, multiply-defined labels, fatal errors or overfull boxes. Nonfatal package/PDF-string warnings and underfull lines remain recorded. The new proposition pages were visually checked in both layouts.

The result concerns the specified mixed tangent. It does not establish a physical global boundary domain, an obstruction for the polar leg alone or for all gauge-related representatives, an $N=4$ result, a complete stationary parametrization, final Dirac classification, a physical mode count, or a completed continuum quantum theory. The eight free initial coordinates are algebraic coordinates of the certified forward family. The broader objective in `GOAL.md` remains active. Historical failed runs and preparation reports are preserved unchanged.
