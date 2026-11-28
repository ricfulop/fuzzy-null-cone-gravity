# Control-field closure reached terminal embedding, then timed out

Run `8fa984ea-3c9b-4e85-b984-c23956fc95d8` used source
`e928d6c89324465018a6b5b209152dffde1a454b`. The supervisor stopped its single
scientific invocation at the 900-second wall budget. It produced three
elimination checkpoints but no `dae-closure.json`, no final retained-row
certificate, and no solution germ.

The [full log](engine.log) records the six exact-row controls, two manufactured
DAE controls, field conversion, radial gauge identity, and principal rank-four
probe. The historical control-only field arithmetic completed step zero,
including differentiation, at 35.177 seconds of checker time. Step one
finished at 198.593 seconds; step two finished at 360.670 seconds.
The [last checkpoint](outputs/dae-step-02.json) records six accumulated state
constraints, no remaining acceleration controls, and no new constraint at
that step. These are intermediate reduction data, not a terminal proof.

The `terminal_embedding` stage began at 360.670 seconds and did not emit an
end event. That stage includes state-constraint row reduction, construction
of a kernel embedding, checking the original constraint product, and final
normalization. The log cannot identify which of those operations exhausted
the remaining time. No assertion failure or traceback is present.

The [supervisor record](execution.json) reports 900.348 seconds elapsed,
launcher exit code 143, `wall budget reached`, and an empty final process
group. Cleanup and log relay both completed without errors. The framework
separately reports exit code 1. The largest of 3032 process-group RSS samples
was 472268800 bytes (450.4 MiB); the instantaneous peak was not measured and
the 2-GiB sampled threshold is not a hard operating-system limit.

All nine frozen source files match their before/after digests and the bytes
at the recorded commit. An independent post-run process-table check found the
original group empty. The first checkpoint's constraints and both control
substitution matrices match the original implementation's serialized
equations exactly. Source checks, checkpoint hashes, and process evidence
are recorded in [archive-status.json](archive-status.json). Stage timings are
nested and listed in [timing-summary.json](timing-summary.json).

## Next diagnostic and scope

The next bounded diagnostic may consume the independently pinned final
checkpoint's constraint matrix to measure and verify an exact algebraic
kernel calculation. Such use does not resume the scientific closure or give
a partial checkpoint the status of a certified export. In particular, it
does not supply the final dynamics, the 22 retained Euler identities, an
interval regularity proof, or a responding seed.

A fraction-free row reduction of the final state-constraint matrix and a
direct rational-field check of its kernel are a narrow implementation
candidate. Their performance remains unmeasured. The historical discovery
constraint basis, radial derivation, original operator export, and all final
solution obligations remain unchanged.
