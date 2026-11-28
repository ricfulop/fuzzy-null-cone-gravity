# Complete sixth-order energy: prospective bounded pilot

This protocol is written **before execution**. It specifies a numerical
assembly of the reduction-first, fixed ten-mode Hamiltonian. No mixed H6
value or full four-point quantum observable has yet been obtained by this
pilot. The code has independent read-only algebra and interface review;
AST/import checks do not constitute numerical validation.

## Question and parent evidence

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
nested-kernel cost is unmeasured. A work-limit failure will be retained as
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
