# First stationary DAE closure attempt — incomplete

One local invocation from source `868dfb9f98875a39de26210fef3dbbc5919fd5b4`
ended with framework exit code 1. It produced no terminal `dae-closure.json`,
solution germ, or nonlinear compatibility result. Run:
`ef17ca1d-2663-4c0d-acc8-2f8607ae7881`; experiment:
`9032485c-47b6-4b0b-8c59-b6aebe68afd5`.

The [complete framework log](engine.log) records both manufactured controls,
the fixed-field conversion, the 23-row radial gauge check, and the principal
rank-four probe before entering closure. Step 0 completed in 68.611125 seconds;
its [partial checkpoint](outputs/dae-step-00.json) records four state constraints
and three remaining acceleration controls. Step 1 began, but no completion was
logged. These intermediate data do not certify the stationary DAE.

The wrapper then raised `PermissionError` from `os.killpg(pgid, 0)` while
checking termination. Its `finally` block repeated that probe and raised before
saving `execution.json` or `resource-samples.json`. Both files are absent. The
configured bound was 900 seconds and the framework ended at approximately that
duration, consistent with a wall-budget stop. **The primary stop reason, child
exit code, and maximum sampled RSS were not persisted and are not reconstructed
as measurements.** No scientific assertion failure appears in the log.

The [startup manifest](launch.json) was written before child launch, so its
`scientific_invocations: 0` is a startup value. The actual child invocation is
evidenced by its log and [invocation identity](outputs/invocation.json). The
[framework exit](framework-exit.json) describes the outer wrapper's exit;
its `verified` outcome is process bookkeeping, not scientific validation.

Post-run inspection found no member of the original process group 52076 and no
surviving known child 52078. The timestamped inspection and the independent
comparison of all six frozen source files with the recorded Git commit and
startup hashes are in [archive-status.json](archive-status.json). This later
source check is distinguished from the missing supervisor completion check.

The archive preserves the original logs, frozen sources, startup identity and
partial checkpoints. A repair will use a fresh invocation, retain all original
22 Euler rows plus the discarded response, add detailed elimination timings,
and correct cleanup so failure evidence is saved even if process-group probing
raises. The separately prepared fraction-free exact-arithmetic option changes
row-operation implementation; its equivalence controls and existing terminal
identities must pass before any result is accepted. Performance is unmeasured.

The completed principal probe and stationary operator export remain separate
earlier results. This incomplete closure attempt establishes neither a
responding seed nor a no-go theorem.
