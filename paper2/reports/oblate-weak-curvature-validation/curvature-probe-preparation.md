# Targeted exact-curvature diagnostic

Preparation only. Two incomplete weak-validation attempts are preserved. This diagnostic tests the evaluation strategy before spending another full validation budget. A diagnostic pass is not a pass of the complete weak-form validator.

## Evidence motivating the change

The first source completed the three representative sampled sectors in 19.924 seconds, then reached the monitored memory budget during the first symbolic gauge control. The second source introduced rational cancellation throughout the common curvature path. It needed 35.885 seconds for the first sampled sector and timed out in the second at 900 seconds, without reaching any gauge control. Its peak sampled group RSS was 221,315,072 bytes and cleanup reached an empty group. These observations motivate restoring the previously fast sampled path and examining a selective symbolic path. They do not locate the exact expensive operation within the second sector.

## Exact representation and contraction

The symbolic gauge path continues to form the metric perturbation by the coordinate Lie derivative, compute its covariant jet, and vary the Levi-Civita connection. It does not replace the calculation with the known pure-gauge curvature identity.

For the connection variation $C^a{}_{ij}$, the direct contraction is

$$
\delta R_{ij}=\nabla_a C^a{}_{ij}-\nabla_j\tau_i,
\qquad \tau_i=\sum_a C^a{}_{ia}.
$$

It retains all nine Ricci entries, without mirroring a triangular block. Contracting the trace before differentiating avoids constructing unused entries of the full three-index covariant derivative. Scalar curvature and Einstein tensor are formed from the same Ricci matrix.

The chart has $d=1-x^2>0$ for $-1<x<1$. Each relevant expression is represented exactly in the one-root extension with $s=\sqrt d$ and $s^2=d$. Only existing half-powers with polynomially identical base $d$ may be rewritten into this field. Numerator and denominator polynomial remainders retain every coefficient. Unsupported radicals or an unhandled denominator structure cause a diagnostic failure rather than an assumed identity. Rational coefficients are canceled before restoring the original chart root. This neither selects a physical radial boundary condition nor asserts a formula at the chart endpoints.

The sampled weak/Green sector path is restored to its earlier raw expression construction, followed by exact rational-point substitution. Its comparison with the prepared physics modules is unchanged. The final zero predicate keeps its existing branch discipline.

## Five diagnostic batches

1. Exact one-root arithmetic controls, including nonzero and unsupported-expression controls.
2. Agreement of both raw and normalized contracted Ricci with the original dense raw varied-Ricci calculation on an off-shell complex-harmonic perturbation at a recorded rational chart point. At least one curvature entry must be nonzero, with that witness recorded.
3. The original global symbolic pure-gauge control at $m=1$.
4. The original global symbolic pure-gauge control at $m=0$.
5. The original global symbolic pure-gauge control at $m=2$.

Fine start/done events identify expensive operations without pretending to complete a batch. The diagnostic uses its own output and success metric. The full validator retains its eight batches and cannot acquire a success metric from this diagnostic.

## Proposed execution

One immutable local CPU invocation, maximum 300 seconds and 2 GiB sampled process-group RSS, with the existing fixed runner and pinned dependencies. The supervisor explicitly records `curvature-probe` mode. No GPU or remote compute is involved. The user's standing overnight approval ends at 2026-11-28 16:29 UTC.

This run cannot establish finite radial `K,E,D`, a common physical operator domain, an admissible non-gauge seed, a spectrum, or nonlinear consistency. Its purpose is to decide whether the exact evaluator is suitable for completing the previously specified validation.

## Source review before execution

Review caught and corrected an intermediate draft in which the new chart algebra was attached to the final zero predicate while the costly curvature path still used generic cancellation. The final design keeps the original zero predicate, uses the strict chart normalizer only in the selected curvature path, and does not apply generic cancellation after recombining the already reduced rational and root coefficients. Additional guards reject nonfinite expressions and nonrational coefficient functions. These were preparation fixes; no diagnostic was executed on the rejected draft.
