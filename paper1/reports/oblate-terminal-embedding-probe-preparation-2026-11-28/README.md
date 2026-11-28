# Exact terminal-kernel diagnostic preparation

The latest full closure attempt completed three elimination stages, then
reached its wall budget during terminal embedding. Its unchanged evidence is
archived in [the terminal-timeout report](../oblate-stationary-closure-terminal-timeout-2026-11-28/README.md).

This separate diagnostic consumes only the independently pinned six-by-fourteen
state-constraint matrix in that run's step-02 checkpoint. It checks the hash
before parsing, performs fraction-free row reduction over `QQ(r,T)`, normalizes
the returned denominator, constructs a kernel with identity free-coordinate
rows, and verifies the exact product against the original matrix. Each phase
is timed. The calculation does not repeat the earlier elimination or change
the historical constraint rows differentiated there.

Six manufactured algebraic controls exercise a known rational kernel, a
deliberately wrong kernel, empty dimensions, zero matrices, full column rank,
and inexact-literal rejection. They remain unexecuted in this preparation and
will run inside the diagnostic. Static source/API review found no blocking
issue. Final view files are present and pass static parsing.

The adapter reuses the previously reviewed and tested reduction supervisor's
resource and cleanup machinery. Four additional mocked dispatch tests passed
in 0.076 seconds with no scientific process or real signal. They test pinned
copies, changed source and input rejection, and the approval cutoff. The
[validation record](runner-validation.json) preserves their output and hashes.

The bounded invocation is one local CPU child, 300 seconds, and a 2-GiB sampled
process-group RSS threshold. The same provisional experiment and fixed shell
command are retained, with the entry point routed to this diagnostic for the
next immutable snapshot. No remote or GPU resources are used. The user's
standing approval window ends at 2026-11-28 16:29 UTC.

A pass establishes only the formal right kernel of this pinned intermediate
matrix. The export explicitly leaves constraint preservation, radial dynamics,
all 22 retained Euler equations, the discarded response, closure, and a smooth
solution germ unchecked. This is diagnostic use of partial evidence, not a
resume or certified completion of the original reduction. Performance of the
alternative remains unmeasured until the bounded invocation.

Source hashes, review references, and the exact input pin are in
[preparation-provenance.json](preparation-provenance.json). The source-level
specification is [here](../../verification/oblate_terminal_embedding_probe_20261128/README.md).
