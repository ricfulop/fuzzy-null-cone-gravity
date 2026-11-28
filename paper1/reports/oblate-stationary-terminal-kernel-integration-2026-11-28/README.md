# Full stationary reduction with the measured terminal kernel

Preparation for one fresh local CPU invocation of the original stationary DAE reducer, with a 900-second wall cap and a 2-GiB sampled process-group RSS threshold. The existing approval window ends at 2026-11-28 16:29 UTC; the supervisor checks that the full budget and cleanup allowance fit before launching. The fixed command remains `bash paper1/verification/oblate_stationary_reduction_20261128/run.sh`.

The earlier full invocation reached terminal embedding after 360.670 seconds and stopped at its 900-second budget without closure. The separate [kernel diagnostic](../oblate-terminal-embedding-probe-2026-11-28/README.md) then passed in 117.565 seconds. It verified a kernel of the saved constraint matrix only. Neither result supplies a radial solution germ.

## Narrow implementation change

`terminal_kernel.py` uses fraction-free `rref_den(method="CD", keep_domain=False)` on the final accumulated constraint matrix C alone. It normalizes the returned scalar denominator, builds the free-coordinate embedding M, verifies the free rows are the identity, and checks the original C times M is exactly zero over `QQ(r,T)`. It returns the same tuple consumed by the full reducer.

The historical control-only discovery pivot loop, original residual basis selection, differential derivation, and original operators are unchanged. New timings separate the six kernel substages and the subsequent constraint-preservation, original-22-row, reduced-dynamics, reconstruction, and discarded-response checks. None of these downstream obligations is replaced by the kernel diagnostic. The original X0 equation remains among the retained rows after gauge removal.

The full run starts from the completed, independently pinned stationary operator export. It does not load any partial elimination checkpoint. Five new terminal-kernel control groups exercise dependent rational rows with a swap, empty dimensions, a zero matrix, a full-column-rank empty kernel, and Float rejection. They precede the original calculation together with six exact-row controls and two manufactured DAE controls. They have been prepared but not executed from this integrated source.

## Validation and limits

All Python sources passed stdlib AST parsing; both JSON configs and the XML view parsed; the restored shell passed `sh -n`. The discovery helper, its tests, DAE control source, supervisor, supervisor tests, progress helper and shell dispatcher are byte-identical to the previous control-only source `e928d6c89324465018a6b5b209152dffde1a454b`. The existing source-hash globs include both new Python files. No scientific module was imported and no checker was executed for this preparation.

The implementation author recorded its static comparison in `msg_00b9b416`. Independent static review `msg_683c54e2` found no blocker, confirmed MemoryRecorder compatibility and unchanged downstream obligations, and inspected the restored full-run dispatch. [Validation](preparation-validation.json) and [source/input hashes](preparation-provenance.json) record the concrete source.

The diagnostic's elapsed time is not a promise for full closure: downstream verification costs remain unmeasured. A resource stop will preserve incomplete evidence. A passing full reduction would establish formal differential identities at epsilon=1/100, not interval regularity, a responding non-gauge germ, nonlinear integrability, a canonical Dirac classification, or a physical mode count. The independent terminal germ checker remains a separate stage requiring an actual completed closure.
