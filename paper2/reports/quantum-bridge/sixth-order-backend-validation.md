# Sixth-order sources, inverse bilinears and algebraic exterior: validation protocol

This protocol is recorded before the first numerical evaluation of the new
backend. It builds on the validated full fifth-order functional at source
`502f75d`, reviewed in `75f3200`. It tests ingredients for the integrated
sixth-order Hamiltonian of the declared reduction-first candidate. It does
not evaluate that full Hamiltonian, its directional jet contractions, a
four-point loop, or a physical large-rank limit.

## Reviewed completed suite

Run `0fe39553-1539-459e-8fa4-411496c79808`, source
`7ad6a54c6e58f97c49cfdb7ea26fb089be4455e8`, completed in 94.96 seconds.
The [complete result](data/sixth-backend-results.json) is copied byte for
byte from the isolated run. Its effective configuration, intermediate
measurements and final summary were checked against the retained raw log.

All nine Fourier cases, 72 synthetic tail checks, 680 source-component
checks, exact momentum parity and both gravitational-tail directions passed.
The maximum scaled source discrepancy was $8.08\times10^{-14}$; the maximum
sphere/collapsed tail discrepancy was $1.48\times10^{-12}$. These are
numerical validation diagnostics, not certified error bars.

The nonzero-monopole direction reproduced the entire preserved raw exterior
coefficient dictionary exactly, including its leading coefficient. It used
235,291 counted expansion/scatter operations. The second direction completed
with 348,341 operations, below the unchanged one-million cap. No original
numerical tolerance, field direction or quadrature was changed.

This validates the specified source and integration ingredients. Production
nested gravitational pairings, the integral of the complete mixed H6
Hamiltonian, higher-vertex contractions and the four-point loop remain
unevaluated at this checkpoint.

## Calculation and fixed conventions

The [bilinear reduction](sixth-order-bilinear-and-tail-reduction.md) leaves
two nested inverse pairings, with localized endpoints,
$$
-\langle K,\mathcal CF\rangle-2\langle r+2t,\mathcal Yr\rangle,
\qquad \mathcal C=(-\Delta)^{-1}/8,\quad \mathcal Y=\Delta_L^{-1}.
$$
The new [Fourier backend](gaussian-fourier-bilinear-backend.md) integrates
spatial/Fourier variables analytically before zero, one or two finite heat
parameter integrals. Its absolute tolerances and Decimal values omit the
common factor $\pi^{3/2}$; displayed binary64 values include that factor.

The separate local density $-wT:T$, with $w=7u+v$ and $T=L_0W$, has no
Gaussian multiplier. Exact angular collapse and exterior recurrence isolate
its algebraic integral and a bound on the Gaussian remainder. A positive
Gaussian tail bound is not substituted for this algebraic term.

## Acceptance specified before execution

1. Nine Fourier controls cover scalar/York constant Gaussians, nonradial
   pure-Gaussian derivative identities, one-heat inverse adjointness, two
   genuinely variable heat parameters with unequal widths/moments, and the
   negative York spectral bounds on a real vector self-pairing. Reciprocity
   swaps both endpoints including their heat weights. Signed merging must
   cancel exactly. Pure Gaussian references use 96-digit direct moments;
   heat references use the independently organized existing spatial pairing.
2. Parameter orders 8/16/24, independent final-axis refinements and 64-to-96
   digit comparisons must meet the integrator's termwise acceptance gates.
   The normalized tolerances are absolute $10^{-12}$ and relative $10^{-10}$.
   Analytic zero-parameter controls require scaled agreement $10^{-40}$;
   independent spatial controls require $10^{-9}$ and must themselves
   converge on 16/24/32-point radial quadrature.
3. Synthetic triple-HEAT controls check the exact recurrence, anisotropic
   parity contractions, the identical-pair off-diagonal factor, permutations,
   merging, work-cap failures, rejected nonintegrable tails, and range guards.
   Closed primitive controls test both a canceled algebraic tail and the
   nonzero algebraic-term Gaussian remainder. A high-moment, large-radius
   control uses an independent factorial expression and Decimal Machin pi
   to detect false zero from intermediate underflow.
4. Four required raw TT directions are reconstructed: the historical
   `noncommuting_labels`, a nonzero-monopole mixed case, pure configuration,
   and pure momentum. Source values at five declared Cartesian points are
   compared with separately assembled physical metric/York jets, retaining
   both connection terms, complete Gaussian envelopes, and all nine tensor
   entries. Every component, symmetry and trace comparison must satisfy
   $|a-b|/\max(1,|b|)\le10^{-8}$.
5. Pure configuration must give $F=-K$ and zero vector sources. Pure momentum
   must give $F=K=r=J_R=0$ and $J_L=2t$; the test explicitly permits and
   verifies nonzero $t$. Recomputing from $(Q,-S)$ must leave $F,K$ even and
   the vector/tensor sources odd, with exact rational source comparison.
6. The two mixed directions undergo full-sphere checks of $-wT:T$ at radii
   $1/2,1,2,4,8,16$ on $8\times16$ and $12\times32$ angular grids. Their
   $j\le12$ premise is inherited from six retained $j=2$ inputs and
   rotation-equivariant operations; both grids satisfy the angular floors.
   All radial/collapsed and grid comparisons require $10^{-8}$ scaled
   agreement. Exterior summaries at radii 8/12/16 remain distinct from
   pointwise comparisons.

The new mixed direction is $Q=E(A)$ and $S=B(XZ)$, with
$A=\operatorname{diag}(-1/2,-1/2,1)$ and $XZ_{xz}=XZ_{zx}=1$.
The [independent analytic audit](h6-bilinear-independent-audit.md) predicts
$$
\frac{M}{\pi^{3/2}}=(0,2079/32,0),\qquad
\frac{M_h}{\pi^{3/2}}=251559/64,
$$
where $M=\int b_1$ and $M_h=\int(14s_2-6R_2)$.
Its exact exterior coefficient must be
$$
A_{3,-3}=-\frac{15}{4096}\frac{M_h}{\pi^{3/2}}
\sum_i\left(\frac{M_i}{\pi^{3/2}}\right)^2.
$$
This is a required nonzero test. The historical mixed direction has zero
monopole and cannot replace it. Terms with slower algebraic falloff must
cancel exactly; a residual slow coefficient fails rather than being deleted.

## Execution, provenance and limits

The fixed command remains `bash paper2/verification/quantum_bridge/run.sh`;
the selector is `sixth_backend_validation`. One local standard-library
process has a 1,800-second cap, further limited by the recorded overnight
authorization deadline. No GPU or remote service is used. Partial results
are checkpointed and failed gates retained. The complete effective config,
source commit and run id are emitted with the output.

### First attempt and bounded implementation repair

Run `c9bfec9d-1fb5-476f-a96a-54e2dccdddb4`, source
`4dd4c03a2b498833581bdb9d23032fc56f450a46`, stopped after 95.91 seconds
at the one-million-operation exterior-expansion cap in the second mixed
direction. This was a computation-budget exception, not a failed physical
identity. The overall suite remains incomplete. The [partial output](data/sixth-backend-first-attempt.json)
was preserved byte for byte, along with its source and raw log.

Before that exception, all nine Fourier controls, 67 synthetic tail checks,
680 source-component comparisons across four directions, and exact momentum
parity passed. The largest source discrepancy was $8.08\times10^{-14}$
scaled. The nonzero-monopole direction passed all 29 tail checks, including
the exact leading coefficient $-16309479355785/268435456$; its maximum
sphere/collapsed discrepancy was $9.79\times10^{-13}$ scaled. These are
ingredient results, not a complete mixed sixth-order energy.

The repair expands each distinct sorted moment triple once per call and
then shifts its powers for each radial degree. It preserves every exact
coefficient, all cancellations, the million-operation cap, and every
numerical tolerance. Work counts include convolution and final scatter
entries. The prior uncached implementation remains as an exact reference
path. New controls compare signed repeated triples, exact work boundaries
and deterministic per-call counts. The first gravitational direction must
also reproduce the entire preserved raw exterior coefficient dictionary.
The same provisional experiment is retried because its requested full
suite was interrupted by an implementation work limit. No original run
snapshot or partial measurement is altered.

Independent code review has inspected normalization, signs, phases,
Gaussian rates, source contractions and the asserted checks. It supplies
no measured result. Shared first-HEAT constructors and evaluators limit the
independence of the source and sphere checks. Fourier parameter evaluation
uses separate rational/Decimal arithmetic, while its spatial references
share the established Green inputs. Quadrature refinement and binary64
tail-bound evaluation are diagnostics, not interval certificates.

The local sixth-order assembly module is included for review but its full
15-family energy is not computed here. Canonical Gram terms, production
nested bilinears, full mixed energy, jet contractions, connected observables,
physical cutoff limits and microscopic corner dynamics remain pending.
