# Direct-state stationary germ check reached the wall limit during interval search

Run `09c32343-1dee-4d53-a1da-c6d1bb2e1fe3` completed the exact field identities for the candidate at $\epsilon=1/100$, including all 22 retained equations. The discarded response row had a nonzero coefficient over $\mathbb Q(r,T)$, where $T=\arctan(1/(10r))$. The checker then collected 57 denominator polynomials and began the combined physical regularity, rank-minor and response search. Its final child event was the start of that search at **264.654623 seconds**. The wall limit expired with the physical rank-six minor, regular interval, solution germ and same-radius response witness all uncertified. [Child log](verification.log) · [Supervisor execution record](execution.json)

## Scientific evidence

The child authenticated the original stationary operators and the candidate fourteen-state snapshot from run `c2a1482f`. It verified their linkage and passed 18 manufactured cases. For the six constraints, eight free state coordinates and zero remaining controls, the following stages ended successfully. Durations are elapsed child seconds. [All 12 completed stages, their raw metrics, and the unfinished stage](timing-summary.json)

| Completed stage | Recorded result | Duration (s) |
|---|---|---:|
| Kernel identities | $CM=0$ and $JM=I_8$ | 88.358759 |
| Preservation | $(DC+CK)M=0$ for all six rows | 70.809932 |
| Physical profile lift | $(P_vK-P_w)M=0$ for seven profiles | 0.000161 |
| Retained Euler equations | All 22 rows vanish on the candidate embedding | 43.248327 |
| Discarded response | Formal nonzero among eight columns | 8.572675 |
| Denominator collection | 57 distinct polynomials | 8.742826 |

No radius, refinement count, denominator sign table, determinant interval or responding initial datum was emitted by the physical search. The run directory contained no final certificate. The supervisor's failure fallback therefore records certificate acceptance as zero and response values as `null`. The earlier child log preserves the positive *formal* response metric. The physical response remains unknown. [Final fallback metrics](metrics.json) · [Framework status and kept flag](framework-status.json)

## Resource and sampling evidence

The supervisor measured 900.2542698751204 seconds for one invocation and recorded a wall-budget stop, child exit 143, and framework exit 1. Across 3,080 process-group RSS samples, the largest observed value was 245,121,024 bytes (233.77 MiB), below the sampled 2 GiB threshold. SIGTERM cleared process group 63736. The final probe returned ESRCH, and a later process-table check found zero group members. The cleanup and log-relay error lists are empty. [Resource samples](resource-samples.json) · [Framework exit](framework-exit.json) · [Archive status](archive-status.json)

At 08:09:45.836 UTC, one read-only macOS sampling observation targeted Python PID 63738 for one second at 1 ms intervals. It recorded 790 main-thread samples. The call graph includes a 245-sample `builtin_eval` branch. The separate top-of-stack counts include 19 for `k_mul` and eight for `_PyLong_GCD`. These counts describe a single observation. It supplies neither a Python source line nor a timing breakdown for the entire interval search. [Original macOS sample](macos-stack-sample.txt)

## Provenance

All eight frozen source files match Git commit `0ae1c15e11b2b8c9a2f82ca86b6774b68f7e0e79`, their materialized originals and the supervisor's preflight and postflight hashes. The frozen checker SHA-256 is `94af8bfe83084b691d55920c001cb64d0d07bebb553a75d1a3384a09a70b5d50`. Both input copies match the independently pinned source files, the child-reported digests and the postflight digests. The recorded child used Python 3.12.11 and SymPy 1.14.0. The 38 child progress events agree between the verification log and the engine log. [Pinned configuration](sources/verification/oblate_stationary_state_germ_20261128/run-config.json) · [Child launch record](launch.json) · [Engine log](engine.log)

The archive preserves the executed implementation. The separately prepared cache and instrumentation patch contributed no code or measurements to this run. Archival verification used file hashes, Git objects, JSON, logs, and process metadata. [Archive manifest](archive-manifest.json)
