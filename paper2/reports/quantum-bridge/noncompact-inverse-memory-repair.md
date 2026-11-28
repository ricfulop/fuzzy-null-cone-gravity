# Noncompact inverse validation: preserved memory interruption

Run `99b9bc29-8627-4993-aa6c-c7cea2a0d341`, source
`162383cf313d60bf50b4c92e15834834d189c6ff`, stopped at its declared sampled
memory guard after 339.055 seconds. The peak recorded resident memory was
2,176,417,792 bytes (2.026947 GiB), exceeding the 2,147,483,648-byte ceiling.
This was an explicit guard interruption, not an operating-system OOM or a
failed numerical tolerance.

The [unaltered result](data/noncompact-inverse-memory-limit-results.json) and
[provenance](data/noncompact-inverse-memory-limit-provenance.json) preserve
the configuration, complete checkpoints, source hashes and failure. The
result and raw log have SHA-256 hashes
`e801207d7f6909d5657191c237aabc09b1b31fa40812c1259872ab247f1da2bc`
and `ea883eacae6a08bb6ab0b4ed1d8a4d2a59eea349afcc2a0e70da53af11b2d27c`.
Durable artifacts also retain the raw log and exact source archive.

The completed suites report passing 18 original-recursion algebra groups,
10 Gaussian exterior groups, and nine corrected synthetic inverse levels.
The actual mixed-source radial levels at N=256, 512 and 1024 all completed
with their declared gates passing. The final N=1024 angular slice was
unfinished. The only saved failed comparison row is the `angular_1024`
MemoryError record; no earlier comparison row failed. The driver's 2,157
required-row count includes this exception row and is not a count of
independent scalar identities.

The exception handler attempted to checkpoint a named failed row while the
sampled high-water memory was still over the ceiling, so the raw traceback
also records a second MemoryError during exception handling. The outer driver
nevertheless saved the final failed status and all preceding checkpoints.
The final suite conjunction was not reached. This attempt therefore does
not establish complete first-response acceptance.

## Memory-only repeat

The user requested a larger memory cap; the local machine reports 128 GiB of
physical memory. The same provisional experiment was repeated with an 8-GiB
(8,589,934,592-byte) sampled ceiling. Scientific source, fields, grids,
tolerances, all work counters and the 1800-second wall limit remain unchanged.
No lower-resolution shortcut or omitted angular test is introduced.

The repeated calculation remains the [declared first-response protocol](noncompact-first-response-protocol.md).
The [completed repeat](noncompact-first-response-results.md) passed the final
angular slice and every required comparison, with 2,178,908,160 bytes peak RSS.
Its scientific source and all shared completed measurements reproduce the
initial attempt. Provider moment enclosures do not certify the combined finite-grid
field error. Neither this attempt nor its memory repair evaluates W2, u3, u4,
the independent mixed H6 functional, or the four-point loop.
