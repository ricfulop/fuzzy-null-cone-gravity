# Closure retry stopped after formula growth

Run `7fbed460-e7b2-4f37-b328-1b4932d91867` from source
`8335010cf170654c7e53b2371b94761d997c7efa` was deliberately cancelled while
differentiating the first set of state constraints. It did not produce a
terminal closure or solution germ. The same closure question remains open.

The [full log](engine.log) records successful manufactured controls, exact
field conversion, radial gauge verification and the principal rank-four probe.
The new fraction-free reduction of the full augmented matrix completed the
first control solve in 21.043 seconds and saved its checkpoint at 35.934
seconds of checker time. Constraint differentiation then ran for more than
five minutes without a completion event. The first historical implementation
had completed that entire first step, including differentiation, in 68.611
seconds. No scientific assertion failed in this attempt.

The [saved checkpoint comparison](archive-status.json) shows the source of the
performance concern: the four constraint rows grew from 748, 971, 670 and 8451
expression characters to 48891, 46031, 42022 and 42574. The substitution state
grew from 8621 to 172200 characters; the checkpoint file grew from 21882 to
354139 bytes. These are serialized formula-size measurements, not ranks or
physical results. Full augmented row reduction is algebraically equivalent,
but it chooses a different state-constraint basis that is expensive to
differentiate here. Restoring the historical control-only pivot order is a
justified implementation change; its future speed remains unmeasured.

The [supervisor record](execution.json) preserves one invocation, 363.362
seconds elapsed, a signal-triggered stop, launcher exit code 143, and a final
empty process group. The largest of 1220 RSS samples was 802406400 bytes
(765.2 MiB); this is a sampled group maximum, not an instantaneous peak or hard
memory cap. Cleanup encountered `EPERM`, classified the group using the
process table, and completed with no cleanup or relay error. Both final
evidence files were saved. The [framework exit](framework-exit.json) separately
records cancellation by signal 15.

All nine before/after source hashes match the frozen copies and recorded Git
commit. An independent post-run process-table check found the original group
empty. Partial checkpoints remain diagnostic evidence only and will not be
loaded as a certified closure. The next prepared variant preserves the
original control-only elimination basis while performing its arithmetic
directly over `QQ(r,T)`; all original 22-row terminal checks remain required.
