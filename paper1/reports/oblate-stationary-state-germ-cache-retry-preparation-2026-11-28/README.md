# Integrated direct-state cache retry preparation

This is a reviewed source integration and mocked supervisor validation, with no scientific execution of the revised checker. It combines the exact interval-cache preparation from `da24c65` with the two-file supervisor repair from `698d7fcd65c3e7eb9f25f6d311d96af27a73e069`. The latter arrived by incorporating main at `2b0c34947c5d8f8b5f5ab98e416aae0b8f3db1d7`; the only merge conflict was documentation, resolved to retain both the cached-method description and the distinction between historical and current tests.

## Reproducible source and validation

[validation.json](validation.json) captures the exact commands, combined outputs, exit codes, and eight before/after source hashes. The 27 mocked infrastructure tests passed in 0.395 seconds; shell syntax, stdlib AST/JSON/XML parsing, and both independently pinned scientific input hashes also passed. No checker module or scientific dependency was imported, and no mathematical control was executed. The supervisor and its test file are byte-identical to the reviewed two-file commit. The checker, run configuration, shell entry point, and view configuration are byte-identical to the reviewed cache preparation.

The updated supervisor preserves the raw child log independently of best-effort console forwarding, bounds process-table sampling by the remaining wall budget, records absent RSS measurements as unknown, and freezes cancellation before terminal evidence publication. A reader must verify the hashes of `execution.json` and `metrics.json` against `terminal-evidence.json` before accepting a run pass. The old 20-test preparation archive and all prior run snapshots remain unchanged.

## Planned invocation

The same direct-state experiment and fixed shell entry point will run one local CPU child, under the same 900-second wall limit, sampled 2 GiB process-group RSS threshold, maximum 256 arctangent terms, and 2026-11-28 16:29 UTC cutoff. The original operator and terminal candidate input pins are unchanged. All 19 prepared mathematical controls run inside that tracked child before proof calculations. The new cache retains exact rational interval arithmetic and the radius, refinement, and response-column order. Detailed nested timings report actual preparation, denominator batches, minor pivots, response candidates and selected serialization.

The prior run `09c32343-1dee-4d53-a1da-c6d1bb2e1fe3` is preserved in the [interval-timeout archive](../oblate-stationary-state-germ-interval-timeout-2026-11-28/README.md). That run reported exact identities and formal nonzero response but produced no physical-radius certificate. Neither the present preparation nor a runtime estimate supplies a physical rank minor, regular interval, response witness, solution germ, global boundary realization or nonlinear conclusion. A bounded timeout remains failed evidence; a formally nonzero response is not a physical witness.

The separate rational-specialization/transcendence method is under preparation in another directory and is excluded from this launch source set.
