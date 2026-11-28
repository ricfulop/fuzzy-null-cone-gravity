# Single Ricci-entry diagnostic preparation

This is preparation, not a mathematical result. The previous five-batch diagnostic stopped at its 300-second cap while constructing the global m=1 gauge Ricci entry (1,2). Its two earlier batches passed: strict chart-arithmetic controls and all nine contracted Ricci entries compared with the dense raw oracle at a recorded nonzero offshell sample. The complete failure is preserved in [curvature-probe-attempt/README.md](curvature-probe-attempt/README.md), result commit `e280107`, source `0031aca219129d89026422278b521f1e87a3160e`.

The source normalizes the divergence A and covector trace derivative B separately, then normalizes A-B again. This redundant work is a plausible explanation for the remaining symbolic cost; the old entry timer cannot isolate the actual operation. The proposed repair retains chart-normalized varied Christoffel components and their trace, forms A and B raw, and normalizes only the combined Ricci entry. The varied-metric formula and all nine independently computed entries remain unchanged. The strict chart normalizer, original exact-zero predicate, raw dense oracle, background metric, weak densities and harmonic conventions are unchanged.

The new fixed `ricci-entry-probe` mode calls the same single-entry helper as the complete Ricci builder. It reports start/done elapsed timings for divergence construction, trace derivative construction, and normalization of their difference. There is no curvature-Lie identity shortcut, numerical tolerance, symmetry mirroring, or assignment of a known zero.

## Bounded invocation

One tracked local CPU invocation uses the unchanged fixed shell entry point and isolated Python 3.12 / SymPy 1.14.0 / NumPy 2.3.5 environment. Its cap remains 300 seconds and sampled process-group RSS remains 2 GiB. Cleanup and source/environment provenance remain active. Authorization is the user's overnight approval `msg_ad636446` and subsequent all-other-approvals message, ending conservatively at 2026-11-28 16:29 UTC. The new scope's runtime is uncalibrated; 300 seconds is a stop limit, not a predicted completion time.

The run owns a three-batch diagnostic plan:

1. The existing strict half-root arithmetic and rejection controls.
2. The existing nine-entry offshell m=1 raw/normalized contracted-versus-dense comparison, including an exact nonzero Ricci witness at the recorded rational point.
3. Only the same global symbolic gauge m=1 Ricci entry (1,2), constructed through the coordinate Lie metric perturbation, tensor jet, independently varied connection, and shared single-entry contraction, then tested with the original exact-zero predicate.

The result filename is `ricci-entry-probe-checks.json`. `diagnostic_passed` may be reported only after successful export. No full-validation success metric is emitted. The five-batch curvature mode and eight-batch weak-form mode remain available with their original scope, but this invocation does not execute them. A successful single-entry check only supports a subsequent decision on the complete gauge diagnostic; it does not complete the three symbolic gauge sectors, weak potential/lapse Green tests, endpoint tests, finite radial assembly, physical boundary choice, spectrum, or nonlinear consistency theorem.

On another time stop, use the new internal timers to choose the next representation change. Do not raise the limit or repeat the full validator without evidence that its bottleneck has been resolved. All prior attempts remain intact.
