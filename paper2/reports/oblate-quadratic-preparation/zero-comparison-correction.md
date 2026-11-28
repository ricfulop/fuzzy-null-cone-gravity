# Exact-zero comparison correction and approved single rerun

User approval `msg_da46f19d` covers a narrow equality-predicate correction, exact zero/nonzero regression controls, and one local rerun with the original kinematic coverage. The [first failed invocation](first-kinematic-check/README.md) and all its raw output remain unchanged.

## Correction

The archived residual contains `sqrt(-(x-1)*(x+1))` and `sqrt(1-x**2)`. Their radicands are the same polynomial, but the old `exact_zero` path did not recognize the resulting cancellation.

The corrected predicate expands the polynomial base of half-integer powers while retaining the exponent. This changes the expression representation, not its mathematical value. It does not split roots of products or assume `sqrt(x**2) == x`. The subsequent `powsimp` call uses `force=False`. Multiplication by the original chart weights remains valid on `-1 < x < 1`, where `sqrt(1-x**2)` is nonzero. The zero criterion remains exact; no floating-point tolerance or random numerical probe is used.

The same invocation first tests three zero controls: the archived residual, unevaluated `0*I`, and reciprocals of the equivalent chart roots. It must also reject five nonzero expressions: offsets of `1`, the exact rational `10^-30`, and `1-x**2`; the sum instead of the difference of the roots; and `sqrt(x**2)-x`, which is not identically zero on the chart.

Only `verify_kinematics.py` changes. The action, harmonic coefficients, kinematic tensors, density assembly, canonical block formulas, realification, and pinned Paper I reference remain untouched. The original scientific check coverage follows the new predicate controls without modification.

## Approved command and evidence

Run once from `paper2/`, after recording the correction commit:

```bash
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 \
  python verification/oblate_quadratic_preparation/verify_kinematics.py \
  --output-dir reports/oblate-quadratic-preparation/kinematic-recheck
```

Use a fresh output directory, preserve the entire stdout/stderr and exit status, and record the source commit and hashes. This is one local CPU invocation under the original resource estimate and three-minute inspection threshold. No remote or GPU compute is involved. Stop after a failure; further correction or execution needs the user's decision.

This document records the correction and authorization, not a successful result. A pass would validate only the prescribed kinematic/reference, selected density/integration, normalization, and realification checks. It would not validate the potential/lapse weak forms or Green identities, select a radial domain/projection order, generate a physical spectrum or an admissible seed, or establish nonlinear consistency.
