# Terminal preservation over the exact coefficient field

The full reduction `157032da-1468-473e-8cd7-37763ff86822` completed its kernel, then reached the 900-second wall budget while constructing DC+CK. The [failure archive](../oblate-stationary-preservation-timeout-2026-11-28/README.md) preserves the unchanged source and partial evidence. It establishes neither preservation failure nor closure. The timer did not distinguish differentiation, multiplication and expression simplification.

## Prepared change

Only the terminal preservation calculation moves to sparse matrices over QQ(r,T). Each coefficient uses the derivation D=partial_r-10/(100r^2+1) partial_T, implemented through exact fraction-field differentiation. Nonzero derivative results alone enter the sparse matrix, so differentiating constants cannot leave explicit zero entries that confuse a sparse zero predicate.

The verifier separately times conversion, DC, CK, the sum S=DC+CK, CB, SM, and the two zero checks. It requires **CB=0 and SM=0**. S itself need not be zero. The original discovery preservation function, pivot/residual basis, 22-row terminal equations, reduced-dynamics reconstruction and response evaluation are unchanged.

Eight manufactured control groups are prepared: an independently specified total derivative, nonzero S=C with vanishing projection, independent unit SM and CB defects that must reject, empty dimensions, Float rejection in every input, exact snapshot shape/string encoding, and nonfinite snapshot rejection. These run alongside the six row, five kernel and two existing DAE control groups in the next tracked invocation. No scientific test from this source has been executed during preparation.

After kernel verification, the real recorder writes `terminal-state-before-preservation.json` atomically with no overwrite. It contains exact C,K,B,M, explicit state-coordinate labels for the seven profiles and first derivatives, shapes, free rows, remaining-control labels and full source/input/attempt identity. Floats and nonfinite entries are rejected before serialization. The file explicitly leaves preservation, original retained rows, reconstructed dynamics, discarded response, closure and germ uncertified. There is no reader or resume path.

## Envelope and review

The fixed command and supervisor remain unchanged: one fresh local CPU child, SymPy 1.14.0, 900 seconds, and a 2-GiB sampled process-group RSS threshold. The full original completed operator export is independently pinned. The user-approved window ends at 2026-11-28 16:29 UTC; startup and immediate pre-spawn guards enforce the complete budget plus cleanup allowance.

AST/JSON/XML parsing, shell syntax, and unchanged-source comparisons passed. All new Python files are automatically covered by the existing child and supervisor hash globs. The new source has no measured runtime yet. The unchanged supervisor tests were not repeated. [Provenance](preparation-provenance.json) and [validation](preparation-validation.json) record the source and review status.

A passing full reduction would still establish only formal fixed-deformation differential identities. An independent terminal certificate must then verify the actual local analytic solution germ and any nonzero discarded response before a Lagrangian nonintegrability conclusion. No physical boundary domain, canonical Dirac rank or degree-of-freedom count follows from this repair.

Independent review `msg_a8352e9f` found no launch blocker after committing the two new source files. Follow-up `msg_40ee982c` confirmed the explicit coordinate-label addition without changing that verdict. Both reviews were static only.
