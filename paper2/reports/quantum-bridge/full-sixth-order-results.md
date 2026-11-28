# Complete reduced sixth-order energy: first numerical result

The repaired pilot completed all twenty contributions for the mixed raw
configuration/momentum direction in the ten-mode reduction-first candidate.
Every declared numerical gate passed. The recorded mixed coefficient is

    H6 = 14796721.40459518,
    summed physical absolute error diagnostic = 5.785261212270265e-5,
    preregistered total acceptance budget = 1.4796731404595178.

This is a homogeneous sixth-order Taylor coefficient under common scaling of
the specified configuration and canonical momentum direction. It is not the
full Hamiltonian evaluated at unit perturbation amplitude. The raw profiles
have the fixed Gaussian width and normalization specified in the
[pilot protocol](full-sixth-order-pilot.md). The number supplies a classical
input to the pending higher-vertex loop contractions. Independent full mixed
functional validation, the four-point observable and physical cutoff limits
remain incomplete.

## Provenance and reconciliation

- Experiment: `84e19d09-b756-4cfd-b07c-2c22d83172e8`.
- Successful run: `dbf69d65-7e3c-440c-a6c7-9c31c5ff8550`.
- Immutable source: `f5a64e3b1a198407507fa7209a3e26da1224ec4b`.
- Recorded calculation time: 337.620 seconds; runner duration 337.989 seconds.
- [Complete byte-preserved result](data/sixth-energy-results.json).
- [Result/log/kernel hash manifest](data/sixth-energy-provenance.json).
- [Prior incomplete attempt](data/sixth-energy-first-attempt.json), run
  `2624f0fc-47a2-499d-af2f-ec8b923998b3`, source `afac83580b7eafd5b329aff57d66bfd93f72211b`.

The successful result matches the two complete directional records emitted
in the raw log. All eight exact raw/compact kernel files match their recorded
SHA-256 hashes and are retained under `data/sixth-energy-kernels/`. The raw
log and complete source archive are durable project artifacts. All 29 original
configuration keys are unchanged. All 37 contribution values shared with the
incomplete attempt are binary64-identical, and the exact canonical factors
and vectors are unchanged.

A separate read-only code/data review independently reconciled these hashes,
shared values, unit conversions and acceptance gates. It found no numerical
or normalization discrepancy. This review did not rerun the computation or
provide the still-required independent gravitational functional evaluation.

The final card now includes the twentieth contribution and final mixed value;
the early completion notification had only collected nineteen entries.
`nested_consolidation[*].status` in the immutable JSON describes the stage
before integration and still says that its integral was pending. The later
complete energy ledgers, final status and raw final record establish that
integration completed; the archived stage message has not been edited.

This experiment has answered its pilot question and its branch is frozen.
Further physical or numerical hypotheses require a child experiment.

## Direction, units and complete signed ledger

The mixed input is Q=E(A), S=B(XZ), with A=diag(-1/2,-1/2,1) and
XZ_xz=XZ_zx=1. Both raw fields include exp(-r^2/2). All ten retained raw
TT basis modes enter the canonical correction. This is a selected radial
and angular basis, not a complete finite-angular-rank scattering problem.

All ledger entries include their full spatial integral factors. Local
outputs already have these units. Fourier normalized values multiply
pi^(3/2), triple normalized values multiply 4*pi, and exact canonical factors
multiply pi^(3/2)/2. The nested scalar outer factor -1, York outer factor -2,
and triple minus sign/tensor multiplicities were applied once in assembly.

| Contribution | Recorded physical integral |
|---|---:|
| curvature | 92.965689850761 |
| raw_gram_fourth | -496.42966290724 |
| york_source_b3 | 20766.875458607 |
| york_A2_T1_square | -358777.91593834 |
| york_A2_T_T2 | -21096.745134086 |
| york_A2_Q_T_T1 | 41187.466863446 |
| york_A2_Q2_T_square | -22870.247804856 |
| york_A2_QTQT | 503.01170236459 |
| scalar_curvature | -19545.217016247 |
| weighted_SQT | 101535.86211014 |
| weighted_S_T1 | 7079582.1329153 |
| weighted_commutator | -114035.16797369 |
| scalar_gradient | 40407.894596142 |
| scalar_B2 | 9787489.0175751 |
| scalar_s2_u2 | 6907363.2761458 |
| gram_square | 281229.56141844 |
| canonical_cross | -6590184.6063258 |
| scalar_nested | -22269.374401325 |
| york_nested | -123925.11971849 |
| triple_tail | -2190235.8359043 |
| **Sum** | **14796721.40459518** |

The sum of absolute family values is 33723594.72435526. Each backend's
convergence result and its refinement, tail, precision and arithmetic
indicators are retained in the full result. The final absolute diagnostic
sums all families in common units, plus a separate arithmetic allowance.
These diagnostics are not rigorous quadrature or interval error certificates.

## Exact consolidation removed the evaluation-cap obstruction

The raw kernels were retained. Exact rational heat-weight elimination and
common-denominator construction reduced the evaluator workload as follows:

| Nested kernel | Raw terms | Compact terms | Raw planned term evaluations | Compact planned term evaluations | Rewrite seconds |
|---|---:|---:|---:|---:|---:|
| Scalar | 18,065 | 207 | 76,805,028 | 971,386 | 4.263 |
| York | 42,435 | 265 | 189,253,918 | 1,256,818 | 13.043 |

Both fit the unchanged 50-million evaluation cap. Rewrite work was 2,108,985
and 6,430,804 declared operations, respectively, under its separate
50-million cap. No coefficient was pruned and the gravitational fields and
functional were unchanged.

Each production kernel passed fifteen exact sampled stripped-integrand
comparisons across its three endpoint blocks. Equal-endpoint comparisons
were symmetrized to respect termwise canonical swaps. These finite samples
support the independently audited algebraic rewrite; they are not a proof
of a rational identity solely by sampling.

The compact scalar/York parameter termwise changes were 3.24e-14 and
6.38e-13 in units of physical/pi^(3/2). Their 64-to-96-digit termwise
changes were 1.48e-54 and 1.54e-53. Absolute compact-term integrals remain
large (5.69e8 and 4.39e9 in those normalized units), so the high-precision
and separate-axis checks remain material. All diagnostics now refer to
the compact Laurent basis; they are not claimed to preserve or dominate
the expensive raw-basis termwise diagnostics.

## Independent controls and their scope

Ten new consolidation cases passed: seven exact fixtures, including raw-map
immutability and compact idempotence; exact work-cap success/rejection; and
two unequal-rate, two-parameter scalar/York derivative identities against
the spatial two-HEAT backend. The independent spatial scaled discrepancies
were 3.47e-18 and 2.78e-17. Nine Decimal triple-integral/moment controls also
passed, with full completed-control and radius records retained.

The two gravitational directions passed 340 source-component comparisons,
350 local T1/T2 comparisons and 163 exact canonical checks. The maximum
source scaled discrepancy was 3.24e-14; the maximum local discrepancy was
5.96e-15. The first attempt's successful shared measurements were reproduced.

Pure momentum returned H6=32618.557080168644 against the independent archived
radial reference 32618.55708016862 (scaled discrepancy 7.81e-16). Only the
49I and 28I families contributed, with both separate numerical targets and
the exact 4*scalar_B2=7*scalar_s2_u2 bin identity passing. Its physical summed
diagnostic was 1.74e-8.

The mixed triple contribution was -2190235.8359043417. Its production
radius remained 12; the recorded normalized omitted-Gaussian remainder
bound was 1.69e-52. The algebraic exterior was retained explicitly.
A small omitted-Gaussian bound does not certify the inner quadrature or the
full gravitational functional.

The next independent check is the [unreduced spatial reference](unreduced-sixth-order-reference.md),
with explicit W2/u3/u4 responses, exterior-source data and finite-ball fluxes.
Its [analytic inverse fixtures](spherical-green-control-design.md) include
forced powers, resonance, bounded dipoles and localized auxiliary monopoles.
That implementation and its gravitational comparison are pending. Subsequent
work still needs directional higher-vertex derivatives, all loop families,
connected external observables, radial/asymptotic and physical large-rank
limits. Microscopic corner dynamics and a quantum reconstruction of the
same physical theory remain separate unfinished parts of the full goal.
