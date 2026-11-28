# Single Ricci-entry diagnostic passed

Run `9b6ab336-ed3c-460b-9308-b2dc21c470e4` from source `3445f140ac8352e7e83a07edefb9ecdcaad83e73` passed its three diagnostic batches. The supervisor recorded exit 0 in 36.74804387497716 seconds, with 163,774,464 bytes peak sampled process-group RSS across 125 samples and an empty final group. The checker reported 35.677620332920924 seconds. Source copies match the execution hashes and pinned Git commit; the result's four scientific source hashes agree with the supervisor manifest.

The passed scope is deliberately narrow:

1. Four positive, two exact nonzero and six unsupported-expression chart controls.
2. All nine raw and chart-normalized contracted Ricci entries agree exactly with the dense raw oracle at the same nonzero offshell m=1 sample. The recorded nonzero (0,0) witness remains `104182865513061*sqrt(6)/6237887788090` at r=23/14, x=1/5, epsilon=1/4.
3. The single previously stalled global m=1 pure-gauge Ricci entry (1,2) passes the original symbolic exact-zero predicate, through the same varied-metric connection and shared Ricci-entry helper used by the complete evaluator.

For the gauge component, normalized Christoffel/trace construction took 1.062515 seconds, raw divergence 0.025211 seconds, raw trace derivative 0.003403 seconds, final chart normalization of their difference 2.029145 seconds, and the exact-zero test 0.000056 seconds. Timings are nested operation durations, not independent whole-run sums. The prior attempt stopped during entry (1,2) without timing its internal operations; no isolated speedup factor can be inferred. This pass supports the combined-normalization strategy for the checked case.

The offshell comparison took 32.393543 seconds cumulatively and was the dominant workload. The third checkpoint arrived at 35.667204 seconds. Successful JSON export preceded the `diagnostic_passed` event. The runtime plan and progress both declared three batches.

The watcher flagged a missing complete gauge result. That is an expected scope exclusion: the configured invocation intentionally checks only R[1,2], and its JSON explicitly sets `all_gauge_entries_performed`, `curvature_probe_performed`, and `full_validation_performed` to false. No complete three-sector gauge result or full weak-form pass is claimed.

Next decision: run the existing five-batch curvature diagnostic with the identical checker and physics, restoring global gauge m=1,0,2 and all nine Einstein components plus scalar in each sector. Its 300-second and sampled 2-GiB limits remain unchanged. Any later eight-batch weak/Green validation is separate. This result does not assemble finite radial K/E/D, choose a physical boundary domain, construct an admissible seed, compute a spectrum, or establish nonlinear consistency.
