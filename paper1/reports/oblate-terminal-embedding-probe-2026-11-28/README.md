# Terminal constraint-kernel diagnostic

One local CPU invocation passed at source `bdb5d4cb4a5d719b36bb0e0c9a154bedcf48b071`, run `6128c4ef-4b35-4ddf-9763-1b40369bd6f6`. This is an algebraic diagnostic on an authenticated intermediate checkpoint from the earlier timed-out closure attempt. It does not certify DAE closure or a radial solution germ.

The saved 6-by-14 constraint matrix has formal rank six over `QQ(r,T)`. The checker constructed a 14-by-8 embedding, verified that its free-coordinate rows form the identity, and checked the product of the **original** constraint matrix with that embedding is exactly zero. Pivot columns are `[0,1,2,3,4,6]`; free columns are `[5,7,8,9,10,11,12,13]`. These are algebraic coordinates, not physical degrees of freedom.

All six manufactured controls passed: a known rational kernel, rejection of a perturbed non-kernel, empty dimensions, a zero matrix, full column rank, and rejection of inexact literals. The diagnostic used SymPy 1.14.0 and Python 3.12.11. The supervisor separately records Python 3.13.7.

The supervised elapsed time was **117.565 seconds**, exit 0. The exact constraint RREF took 23.310 seconds; its denominator normalization took 1.170 seconds, kernel construction 0.003 seconds, and the original-matrix product check 90.064 seconds. Serialization took 2.075 seconds. Timings come from the flushed raw log; the result JSON ends its timing table before serialization.

The largest of 397 process-group RSS samples was **132,972,544 bytes (126.8 MiB)**. The 2-GiB setting was a sampled stop threshold, not an OS hard cap; an inter-sample peak was unmeasured. The 300-second wall limit was not reached. The final process group was empty and the supervisor recorded no cleanup or log-relay error.

## Evidence

- [Complete engine log](engine.log) and [child log](verification.log).
- [Execution/environment record](execution.json), [launch identity](launch.json), [RSS samples](resource-samples.json).
- [Kernel result and exact embedding](outputs/embedding-probe.json).
- [Authenticated partial input](inputs/dae-step-02.json), SHA-256 `ca630b313848bc4230f72f05642d0821297e50322ed2f0dc5ccbaa43eb104ccc`.
- [Archive checks](archive-status.json), frozen sources under `sources/`.
- [Prepared diagnostic and scope](../oblate-terminal-embedding-probe-preparation-2026-11-28/README.md).

All nine committed source snapshots match their recorded digests and the bytes at the source commit. Original and copied inputs match the independent pin and before/after records. Archival used file copying, JSON inspection, and hashing only; it did not rerun scientific code.

## Interpretation and next step

The result supports using this C-only fraction-free kernel calculation in the full reducer's terminal embedding stage. The historical control-only discovery pivot basis can remain unchanged. It does not identify which suboperation consumed the earlier combined `terminal_embedding` timer, so no speedup ratio for that old stage is claimed.

Constraint preservation, the reduced radial dynamics, all 22 original retained Euler rows, the discarded response, and interval regularity remain unchecked by this diagnostic. All six corresponding output flags are false. The partial input has not become a certified closure or a resumable final solution. A fresh full reduction must reconstruct and verify the dynamics before an independent germ certifier can be used.

Independent read-only result/provenance review (`msg_b3a096e1`) found no material issue. The final archive status records that review; no scientific calculation was repeated.
