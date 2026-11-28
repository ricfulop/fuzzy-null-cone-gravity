# Independent inverse primitives: refined result

All declared required controls passed in run
`b054b86a-be75-4406-9d77-43d01dccc3e6`, source
`e803109b79af5e7fdd3853953661260b51574680`. The result contains 26 angular
control groups, nine scalar controls and twelve York case-levels across
independent refinement slices. The final callback ledger has 903 required
comparisons, zero required failures and nine coarse diagnostic failures.
Those coarse failures are retained, as are the nested coarse scalar failures
which the callback diagnostic counter does not count.

The [complete result](data/independent-inverse-results.json) is byte-identical
to the actual run output. The [provenance record](data/independent-inverse-provenance.json)
records SHA-256 hashes, raw-log locations, source identity, and the complete
shared-parent comparison. All 1,255 complete checkpoint objects agree exactly
with their raw-log events. Configuration is at raw-log line 4 and the terminal
summary at line 4260. The run-card collection initially lagged the terminal
output; its final inspection shows the correct 903 comparisons and passed flag.
An independent read-only audit reproduced the hashes, event reconciliation,
acceptance-count replay, shared-parent checks and numerical pass flags. It
found no remaining concrete acceptance, unit or provenance defect.

## What changed from the failed parent

The [initial run](independent-inverse-initial-result.md) passed its angular
controls but failed two high-angular-order scalar regular third derivatives.
The [refinement](independent-inverse-refinement.md) changed only scalar radial
levels and their driver pin. All accuracy gates, the fixed physical origin
patch, solver, fixtures, angular/York schedules and resource caps remained
unchanged. Newly incorporated main-branch oblate archives do not alter this
verification directory; they describe a distinct regulator track.

At R=4, the scaled regular-jet errors were:

| Scalar mode | 1024 intervals | 2048 intervals | 4096 intervals | Gate |
|---|---:|---:|---:|---:|
| ell=8, m=-3 | 0.00861877925 | 0.00210472752 | 0.00051995154 | 0.002 |
| ell=10, m=7 | 0.02257143544 | 0.00542949694 | 0.00133096473 | 0.002 |

Thus the measured finest levels meet the unchanged gate. This agrees with
the prior second-order prediction but is not an interval error certificate.
All 26 angular control records and all seven shared R=4,N=1024 scalar level
records reproduce the frozen parent exactly. The parent failure remains
archived and its experiment branch is unchanged.

## York evidence and limits

Both compact fixtures were evaluated by an actual two-stage scalar solve:
V=P f, h=div V, chi=P h, W=-V-grad(chi)/4. The second source uses computed
first-response origin and exterior coefficients. Analytic targets do not set
those computed coefficients. The localized auxiliary monopole and bounded
dipole are retained, along with the traced exterior r^-3 contribution.

Required structural gates passed at all radial levels. Physical field,
derivative, pairing, longitudinal norm and exterior flux targets passed at
the finest radial level and at each single-level angular/radius slice. Nine
coarse physical rows missed final tolerances and remain diagnostic under the
preregistered acceptance roles.

The largest cross-slice scaled changes were 2.9643e-14 for the angular slice,
2.2457e-5 for R=2 versus R=4, and 2.2501e-6 for R=6 versus R=4. These compare
common physical points and global pairings/norms. Exterior norm/flux values
depend on R and were instead tested against their own exact targets.
Exterior radial quadrature order remained 12 while sphere quadrature changed.
The largest required York error-to-its-own-tolerance ratio was 0.324565,
from the final ddW refinement (0.0009736934 against 0.003). This ratio
describes the declared diagnostic; it is not a universal error bound.

The reported wall time was 294.7094285 seconds, with sampled peak RSS
422,969,344 bytes (403.375 MiB). The cumulative vector count was
1,606,392 point evaluations under the five-million cap. Angular counted work
was 31,146,657 under its fifty-million cap. These counters have their declared
operation meanings; they are not hardware FLOP counts or memory reservations.

This validates the specified synthetic angular, scalar and York primitives
at the declared resolutions and tolerances. It does not yet validate the
noncompact gravitational source adapters, propagate their Gaussian exterior
remainders through nested solves, or evaluate the original mixed H6 functional.
The [unreduced reference](unreduced-sixth-order-reference.md) still requires
actual gravitational W2, u3, u4 responses, source anchors, lower-order
normalizations, and finite-ball flux and energy comparisons. It must retain
the exterior-source correction inside the u4 solve.

No four-point loop, physical large-N limit, microscopic corner-state
multiplicity or reconstruction of the original microscopic theory follows
from this run. Numerical accuracy certification remains false.
