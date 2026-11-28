# Bounded independent terminal-germ verification: preparation

This package prepares the independent check of a completed stationary DAE
closure against the original Paper I operator export. It contains no executed
germ check, no terminal closure input, and no claimed solution or responding
seed. The active closure calculation remains separate.

The scientific checker retains the original 22 Euler rows, including the
equation associated with the fixed `X0` profile. Its independent profile and
derivative identities, augmented constant-control check, rational interval
regularity certificate, and same-radius response test are described in the
[checker specification](../../verification/oblate_stationary_germ_20261128/README.md).

The new runner snapshots independently pinned original-operator and closure
JSON inputs, verifies materialized source against the immutable run commit,
and invokes the checker on those frozen copies. The supervisor preserves the
raw log, exit status, resource samples, process-group cleanup evidence, and
original/snapshot hashes. A clean zero exit is accepted only with an
input-matched terminal certificate. The checker records its actual runtime
environment as well as the pinned input identities.

The prepared limits are one local CPU child, 900 seconds, and 2 GiB sampled
process-group RSS, fitting the user's standing approval window ending at
2026-11-28 16:29 UTC. The memory threshold is sampled, not a hard OS cap.
The committed input configuration will be supplied only after a completed
closure has been inspected and pinned. The existing tracked command will be
preserved on the eventual child experiment by routing its shell entry point
to this runner. No live reduction source is changed by this preparation.

## Verification and interpretation

The runner's 19 mocked stdlib tests passed with exit code zero in 0.256 seconds.
They cover input and source mutation, immutable dispatch, cutoff, interruption,
cleanup, raw-log relay, and silent/unresolved/responding certificate states.
No scientific child or actual process-group signal is used by those tests.
The captured command, output, and source hashes are in
[runner-validation.json](runner-validation.json).

The instrumentation adds stage timing and cheap numeric progress values;
it does not repeat symbolic calculations or change the terminal identities.
Python syntax, JSON, view XML, and whitespace checks are static preparation.
The mathematical manufactured controls have not been executed in this stage.

A formal nonzero response without interval separation stays unresolved.
An exact zero establishes silence only on the verified forward subfamily.
A responding witness requires every relevant denominator and a response
numerator to be separated from zero at the same radius. Neither the ODE's
coordinate count nor a successful terminal check is a canonical Dirac rank or
a physical degree-of-freedom count.
