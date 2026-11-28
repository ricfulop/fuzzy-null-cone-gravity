# Exact interval-search cache and instrumentation preparation

**Reviewed preparation only.** The original direct-state run passed the exact algebraic identities and reported a formally nonzero response, then reached its 900-second wall limit during physical-radius certification. It exported no germ certificate. This patch removes repeated polynomial preparation and exposes the work inside that search. Its runtime and new manufactured controls have not been measured. [Original failed run](../oblate-stationary-state-germ-interval-timeout-2026-11-28/README.md)

The checker keeps canonical exact polynomials, extracts the minor's numerator/denominator pairs and response numerators once, and caches one coefficient tuple per unique polynomial at each candidate radius. The cache is recreated at each radius. Alternating arctangent sums advance incrementally through the same rational brackets. Rational Horner evaluation, interval elimination, radius/refinement/response-column order, and the requirement that every denominator, the original rank-six minor and the chosen response numerator be certified at the same radius are preserved.

Only `EnclosureUndecided` triggers interval refinement. Parsing, structural or serialization errors propagate as failures. Decimal evidence formatting occurs after selecting the result. Nested `PROGRESS` events identify polynomial preparation, radius specialization, refinements, denominator batches, minor pivots, response columns and selected-result serialization. Events report endpoint bit lengths without printing the endpoints. These events remain inside the existing combined interval phase; they do not represent separate completed proof stages.

## Validation and unchanged obligations

Independent source review `msg_686ca093` found no blocker. Root's static comparison against `0ae1c15e11b2b8c9a2f82ca86b6774b68f7e0e79` found 23 critical input, field, proof-product and interval-arithmetic functions unchanged. Python AST, JSON and view XML parsing passed. The stage plan retains its 17 phases. No project module was imported and no mathematical checker was executed during preparation. [Static validation](static-validation.json)

The prepared manufactured path now has 19 cases. The new case compares cached and original bracket/refinement verdicts, exact regularity evidence, and response enclosures at two radii; the row-swap control now requires determinant interval `[-1,-1]` as well as one swap. Existing positive, silent, unresolved, singular-radius, physical-lift and all-22-row defect controls remain required. These are prepared checks, not new passing results.

The mathematical input pins, 900-second wall limit, sampled 2 GiB process-group RSS threshold, 256-term bound and approval cutoff are unchanged. The configuration changes only explanatory purpose text. The view now displays the already reported formal algebraic response before radius certification and keeps it separate from a final certified witness. [Pinned configuration and source hashes](preparation-provenance.json)

## Before another launch

The source copies in this preparation include the original supervisor implementation and its historical 20-test validation. A separate infrastructure review identified future-run supervisor corrections; that reviewed supervisor and its test evidence must be integrated before the next invocation. This report does not certify those pending changes or authorize a larger budget.

A successful future run must still export a valid, input-matched certificate after all original proof obligations. A formal nonzero response alone is insufficient. Exact silence is scoped to the forward ODE family; an unresolved physical response remains `null`. No physical interval, responding germ, nonlinear obstruction, global radial boundary realization or Dirac classification is established by this preparation.
