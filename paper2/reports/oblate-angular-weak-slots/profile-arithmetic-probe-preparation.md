# Actual-profile arithmetic comparison preparation

Status: independently source-reviewed and statically checked; no scientific execution of this diagnostic yet.

The previous [dense attempt](first-attempt/README.md) and [pairwise attempt](profile-pair-attempt/README.md) both stopped at the unchanged 300-second wall limit before completing an angular snapshot. Both frozen source sets, logs and runtime records remain preserved. The pairwise attempt completed three controls and at least 120 of the 144 terms in its first potential-density family block; no direct reconstruction comparison completed. Its component profile jets were built in 0.522 seconds, while the first density/flux slot exports took 121.717/13.770 seconds. Row 10 began but did not finish. The row timers do not identify which inner symbolic operation dominated.

This follow-on tests a smaller arithmetic hypothesis before another full-case retry: exact component compaction after differentiation and fixed-radius substitution may reduce repeated expansion during the weak-density contraction. It does not add the compaction method to either production snapshot mode.

## Fixed scope

- Original metric-coefficient Galerkin weak form; L=2, m=0, r=3/2, epsilon=1/4.
- Left and right profile families both 0, left spatial index 10, all 12 right spatial indices. The result must record the actual angular labels.
- Existing three control batches, augmented by explicit compaction support/rejection controls within the first batch.
- Build each actual coefficient-times-profile-times-basis tensor and differentiate it through `tensor_jet` before substituting exact r and epsilon.
- Compare raw jets to exact rational component compaction over QQ(I)(x). Individual component coordinate poles are permitted; assembled scalar densities still pass the existing strict rationality and endpoint/interior pole guards.
- Require exact component equality and exact equality of each assembled raw/compact density. Integrate both versions at 24, 48 and 96 Gauss-Legendre nodes with the existing dx/2 measure and refinement tolerance 2e-10. Require a nonzero density witness.
- Record each completed pair separately and a distinct final arithmetic result. No slot matrices, lapse/flux comparison, realification, Hermiticity or angular snapshot are claimed by this row diagnostic.

The direct values continue to use the actual weak-density formula. Exported derivative-slot matrices do not enter this diagnostic. The earlier four-family reconstruction obligation and all six angular snapshots remain pending.

## Timing and result contract

Optional instrumentation records component compaction, density construction, substitution, scalar expansion, the combined cancellation/QQ-coefficient/denominator-root guard step, lambdification and quadrature. These operation groups emit start/done records; cancellation and its guards are timed together. Raw and compact evaluator caches are isolated. Timings share one Python process and may share SymPy internals; they are diagnostic measurements, not a controlled speedup benchmark or a performance acceptance gate.

Plan: setup 1 environment; controls 3 batches; arithmetic 12 pairs; export 1 result. `profile_arithmetic_probe_passed` is emitted only after the atomic write of `profile-arithmetic-probe-checks.json`. Completed pair files are partial diagnostic evidence until the entire configured probe passes. Neither `angular_export_passed` nor `completed_snapshots` is emitted by this mode.

## Authorization and limits

The user renewed approval in `msg_c7ac8e6c`, chat `chat_5351fc4a-9bc3-4992-9e39-0a9505c63137`: “ok I am going back to sleep, you have my approval for any of the upcoming experiments that come up on the chat or other requests to keep going till I wake up”. This directly read authorization covers reviewed bounded local follow-ons; it supersedes routine repeated approval prompts. The owning coordinator has been informed of the two failures and this smaller diagnostic proposal.

The fixed launch command, isolated Python 3.12 / NumPy 2.3.5 / SymPy 1.14.0 environment, one local scientific process, 300-second wall cap, 2-GiB sampled process-group RSS threshold and 2026-11-28 16:29 UTC cutoff remain unchanged. No GPU or remote compute is authorized here. Runtime is unmeasured; the cap is not a duration estimate. Source review, static checks and source authentication are required before the single launch.

A pass establishes only the specified row arithmetic comparison. It does not complete any angular snapshot, radial projection/operator, endpoint-domain choice, finite K/E/D assembly, kernel-preservation check, spectrum, seed or nonlinear theorem. The reduction-first Gaunt bridge is a separate regulator and research track.

## Final source review

Independent design review `msg_7544de5e`, initial source review `msg_33b5ea9b`, hardening review `msg_3dd741a7` and final review `msg_53ecdd46` found no remaining blocker for this bounded row diagnostic. The final checks cover operation-level start/done records, fixed labels and dimensions, exact component equality, independently constructed raw/compact densities, nonzero witnesses, isolated probe-only controls, and self-describing pair partials.

AST comparison against executed source `8eab4afd1a153331c2252669e3cfde2b6017b13f` confirms the pre-existing formula and evaluator bodies are unchanged. The only changed pre-existing definitions are the optional control entry point and top-level argument parsing/dispatch. The physics modules, launch shim and earlier weak-form verifier match that source byte-for-byte. All seven selected Python/shell sources are hashed in the [preflight record](profile-arithmetic-probe-preflight.json); AST/in-memory compilation, JSON/XML checks and nine stdlib-only mocked supervisor tests passed. No scientific module was imported or executed by this preparation.
