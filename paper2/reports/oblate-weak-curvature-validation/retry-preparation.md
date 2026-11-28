# Bounded weak-validation retry preparation

This is preparation for one local retry, not a successful scientific result. The original invocation `85aa6197-97a5-4a4a-9239-4119a350d7ef` and source `42200f8c135e9ae5837a25515ba33b3c4055e1fe` remain unchanged. Its failure archive was committed as `b189395`.

## Reason for the repair

The original log records the three representative weak-density/Green sector passes, then entry into symbolic pure-gauge `m=0`. The supervisor observed process-group RSS above its sampled 2-GiB budget and began termination. A subsequent process-group permission error interrupted cleanup accounting. The child's return code was not recovered. Later process-table inspections found neither the original group nor its recorded children; those observations are separately recorded in the failure archive.

Expression growth in the symbolic curvature calculation is the source-level explanation being tested. The first run did not instrument the exact operation at the memory stop, so this is an inference, not a measured attribution to one SymPy call. The failed attempt supplies no final validation result or physics contradiction.

## Algebra-preserving checker changes

The independent varied-connection path now cancels rational expressions component by component after constructing the connection variation, its covariant derivatives, and the Ricci contraction. This reduces expression size before later differentiation. It retains the same tensor contractions and exact arithmetic. Existing radicals are kept in their original branch; the repair does not split radicals, force power identities, replace symbolic zeros by floating tolerances, or discard any tensor component. Identities are interpreted on the original nonsingular oblate chart.

Each pure-gauge sector computes its varied Ricci tensor once and derives both the scalar curvature and Einstein tensor from it. The previous implementation built that tensor twice. All three global symbolic gauge controls remain enabled, as do the representative sampled Green checks, Cartesian TT normalization, and nonzero ablations. Flushed substage timings separate connection, coordinate Lie derivative, tensor jet, varied Ricci, scalar/Einstein contraction, and final exact-zero work.

The action, harmonic conventions, `weak_forms.py`, `potential_hessian.py`, and `kinematic_blocks.py` are unchanged. The repaired evaluation order is not a change to the gravitational regulator.

## Execution and evidence

The fixed command, isolated dependency pins, one local CPU process, 900-second wall budget, and 2-GiB sampled process-group RSS budget are retained. The user's overnight authorization ends at 2026-11-28 16:29 UTC. The 12-hour authorization is not a requested run duration.

The lifecycle repair records signal and process-table observations separately, preserves the primary stop reason, and records unknown cleanup or child status without substituting a successful exit. Infrastructure validation uses mocked process and signal calls, not real scientific invocations. Exact source hashes and the final static/mock validation results accompany the retry preparation record.

Setup and export now each report their actual single-unit total. The scientific total remains eight completed batches; substage timing messages do not advance that count. The success metric is emitted only after the result JSON is written.

Even a complete pass would establish only the stated bounded validation: representative `L=2`, `m=0,1,2` density/flux checks at recorded rational chart points, three continuum pure-gauge controls, and the stated normalization/negative controls. It would not assemble finite radial `K,E,D`, select a physical boundary domain, certify a retained nonlinear solution, compute a spectrum, or prove full nonlinear consistency.
