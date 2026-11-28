# Complete sixth-order energy: pilot protocol and first partial result

The protocol below was committed before execution. The first run then
completed the pure-momentum case and seventeen mixed contributions but
stopped at its declared parameter-evaluation work cap. **No full mixed H6
value or four-point quantum observable was obtained.**

## Reviewed first attempt

- Experiment: `84e19d09-b756-4cfd-b07c-2c22d83172e8`.
- Run: `2624f0fc-47a2-499d-af2f-ec8b923998b3`.
- Immutable source: `afac83580b7eafd5b329aff57d66bfd93f72211b`.
- Recorded calculation time: 312.326 seconds; runner duration 312.706 seconds.
- [Complete byte-preserved partial output](data/sixth-energy-first-attempt.json).

All nine triple-integral controls and the indicator-conversion controls
passed. The driver preserved thirteen completed control/radius records.
Pointwise source and T1/T2 checks passed for both directions, together with
87 pure-momentum and 76 mixed exact canonical checks.

Pure momentum completed all twenty entries, with only the two required
families nonzero. Their exact bin ratio passed. The total was
32618.557080168644, versus the independent archived radial value
32618.55708016862, a scaled difference of 7.81e-16. The summed physical
numerical indicator was 1.74e-8. The separate 49I and 28I comparisons also
passed; no cancellation between them was needed to pass the reference.

The mixed direction completed fifteen local and two canonical entries,
each satisfying its own numerical gate. Both nested kernels compiled:

| Kernel | Compiled terms | Two-parameter terms | Compilation seconds | Planned term evaluations |
|---|---:|---:|---:|---:|
| Scalar | 18,065 | 14,753 | 16.092 | 76,805,028 |
| York | 42,435 | 36,391 | 147.856 | 189,253,918 |

The planned counts include all 16/24/32 grids, independent final-axis
refinements and the 96-digit repeat. The scalar preflight rejected work
above the declared 50-million limit before integrating either nested term.
The triple term and total energy gate were consequently not reached.
This is a computational work-limit failure, not a measured failure of
the gravitational formula, and it leaves the energy question unanswered.

The partial JSON, complete raw log and exact source archive are retained
as durable project artifacts. The proposed implementation repair is exact
algebraic consolidation of the redundant Laurent kernel representation;
the original source/run remain preserved. No tolerance or work-limit
relaxation has been selected, and the repaired evaluator has not yet run.

## Exact representation repair, pending execution

The same provisional experiment now adds
[exact kernel consolidation](gaussian-kernel-consolidation.md) before the
existing numerical integrator. Each heat weight is replaced by the rational
identity t^(2m)=(A-beta)^m/a^m in the original variable t. Signed Laurent
coefficients merge exactly, and a sparse Horner construction forms a common
denominator. The compiler, fields, physical functional, all quadrature grids,
64/96-digit precision comparison, numerical tolerances and 50-million
parameter-evaluation cap remain unchanged. This is an implementation repair
of the unanswered pilot, not a new physical regulator.

Ten tracked controls precede the original suite: seven exact representation
fixtures with input immutability and output idempotence, a deterministic
rewrite-work cap boundary, and scalar/York two-parameter derivative identities
against the spatial two-HEAT backend. Completed successes and failures are
checkpointed. Both complete raw and compact production coefficient maps are
saved to run output files with SHA-256 hashes before numerical integration.
Five exact rational point checks accompany each production kernel. They
check symmetrized stripped integrands for equal endpoints because the existing
canonical orientation can swap individual terms after heat-weight elimination.
These are sampled implementation checks supporting the algebraic derivation,
not a proof of functional equality from a finite sample.

The rewrite has its own preregistered 50-million-operation cap, a 200,000
live-map support cap, numerator degree 128 and coefficient bit-length 4096.
Exact representation checks have a one-million term-visit cap. No nonzero
coefficient is pruned. The original raw compilation statistics are retained.
Parameter convergence and cancellation diagnostics now refer to the compact
Laurent basis; they need not equal or dominate the raw-basis diagnostics.
Common-denominator expansion can worsen conditioning. A new work or
convergence failure remains possible and must be preserved.

The [independent unreduced reference design](unreduced-sixth-order-reference.md)
specifies an additional spatial W2/u3/u4 calculation, with exterior-source
data and finite-ball flux checks. It has received algebraic review only;
its implementation and numerical evaluation remain separate pending work.

## Preregistered question and parent evidence

The parent experiment `2a5669d5-55a0-445f-8538-0941b3615cc9` established the
source, Fourier-kernel and exterior-tail ingredients, as documented in
[its reviewed report](sixth-order-backend-validation.md). The new child
asks whether the complete integrated H6 expression can be assembled with
all twenty signed terms and fixed numerical acceptance gates. It preserves
the parent's source and result snapshots.

The two raw directions, in order, are:

1. Pure momentum: Q=0, S=E(A), with A=diag(-1/2,-1/2,1).
2. Mixed monopole: Q=E(A), S=B(XZ), with XZ_xz=XZ_zx=1.

The second has a nonzero York monopole and therefore tests the genuine
algebraic exterior contribution. No canonical-coordinate amplitude factors
are applied to these raw profiles. The retained modes remain the five
electric plus five magnetic TT profiles, including all ten components in
the canonical Gram calculation.

## Assembly and normalization

The [audited reduction](sixth-order-bilinear-and-tail-reduction.md) is
implemented by `reduced_sixth_energy.py` and `sixth_energy_check.py`:

- Fifteen signed local Gaussian/first-response families, each integrated
  separately with the existing angular-collapse backend.
- The exact canonical Gram square and canonical cross term. Their rational
  factors multiply pi^(3/2)/2.
- The nested scalar term -<K,C F> and York term -2<JL,Y JR>. The scalar
  compiler already includes C=(-Delta)^(-1)/8. Normalized Decimal values
  multiply pi^(3/2); the outer -1/-2 are applied exactly once.
- The algebraic-tail term -integral w T:T, including both off-diagonal
  tensor multiplicities. Its new Decimal integrator returns physical/(4*pi).

Every value and error diagnostic is converted to physical spatial-integral
units before entering the twenty-row ledger. Positive indicators are summed
and scaled in Decimal before upward binary64 conversion; an otherwise
underflowing positive indicator is retained as the least positive subnormal.
This conversion safeguard is tested in the tracked driver.

## Acceptance gates

The new full triple-integral backend first checks analytic negative and
zero integrals at both radii 10 and 12, without changing tolerances. Independent
Decimal AGM pi targets are distinct from the production Machin calculation.
Positive Gaussian quadrature checks the production Kummer moments at the
origin and unequal moment/argument pairs, including m=90 and z=800.
Series-cap rejection is mandatory. Every completed control and radius row
is checkpointed before proceeding, including failed rows.

For each gravitational direction, the source checks are retained and the
stripped T1/T2 fields are compared with an independent pointwise York-series
assembly using the original W jets and metric coefficients. T2 symmetry
and its trace identity are gated. The strong and weak A2 densities are
recorded separately: they need not be pointwise equal, and this pilot does
not claim an integrated boundary-flux check.

All ten exact g/c components are checked against the archived quartic
invariant polynomial and old scalar/York Fourier formulas. Full vectors,
Gram matrix, normalized exact factors, compiler statistics and control
rows are retained. Canonical controls run separately from assembly.

For pure momentum, only scalar_B2 and scalar_s2_u2 may contribute. Their
exact collapsed bins must obey 4*scalar_B2=7*scalar_s2_u2. Their separate
physical values are checked against 49 and 28 times the archived independent
radial integral, as is their total 77 times that integral. The archived
HEAT result, radial refinement and tail bound are also recorded. The
reference is pinned to higher-geometry source
`2a59f1b60f630bd2b4e4d922dc693a7810ab2608`; the quartic data are pinned to
`a078aca003f0bb3878c0edf992489424866b0532`.

Each family must satisfy its own refinement/precision/tail gates. In
addition, the **sum of all physical absolute error diagnostics**, including
an arithmetic allowance, must meet 1e-6 + 1e-7*abs(total H6). Thus separately
small relative errors cannot conceal an unresolved cancellation in the
total. No estimated error here is an interval certificate.

## Fixed numerical settings and work budget

`verification/quantum_bridge/config.json` is the authoritative prospective
configuration. Local radial, Fourier-parameter and triple radial orders are
16/24/32. Fourier and triple arithmetic use 64 digits with a 96-digit check.
Fourier absolute tolerance 1e-10 is in units of physical/pi^(3/2); triple
absolute tolerance 1e-10 is in units of physical/(4*pi). Both relative
tolerances are 1e-9. Local absolute/relative tolerances are 1e-8/1e-9 in
physical units. The pure-momentum comparison tolerance is 1e-8 scaled.

Production nested compilation and parameter evaluation each have a
50-million work cap; Fourier and final kernel terms are capped at 200,000.
The canonical/local pair cap remains 20 million and the exterior expansion
cap remains one million. The triple radius is fixed at 12 with panel width
2. Local Gaussian radii use the existing deterministic analytic-tail search
(target absolute tail at most one quarter of the local absolute tolerance);
the chosen radii and search records are retained. These settings are declared
before any production result. They are not
automatic retry parameters.

The calculation uses one local standard-library Python process, at most
1800 seconds and never past the standing authorization deadline. This is
a wall-clock cap, not an empirically calibrated runtime estimate. Production
nested compilation was measured in the partial attempt above; the repaired
rewrite and complete integration costs are unmeasured. A work-limit failure will be retained as
such; it cannot be reported as a physical result. No GPU or remote budget
is used. The fixed entry point remains `bash verification/quantum_bridge/run.sh`
from the Paper II folder through the tracked runner.

## What success would and would not establish

Success would supply two numerical evaluations of the complete reduced H6
formula, with an independent pure-momentum comparison and detailed checks
of its mixed ingredients. The mixed result would still require an independent
unreduced W2/u4 calculation or a separately derived full spatial reference
before claiming independent full-functional validation.

The next layers are directional jets, higher-vertex loop contractions,
connected observables and physical cutoff/asymptotic limits. This pilot
does not supply those layers, microscopic corner multiplicities, equivalence
of the original nonlinear regulators, or a uniform-rank quantum completion.
