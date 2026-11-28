# Stationary reconstruction reached the 1800-second wall limit

Run `874b4e9c-d962-4595-8e76-ece483c13e07` completed construction of the reduced dynamics at 1363.891531 seconds of child time. Its last event opened the full state reconstruction check at 1363.891817 seconds. The supervisor stopped the child at the 1800-second wall budget before it returned a verdict. The discarded response remains unknown. [Child log](verification.log) · [Nested timings](timing-summary.json) · [Execution record](execution.json)

## Completed checks and stopping point

At fixed $\epsilon=1/100$, the three discovery steps added four, two and zero constraints, leaving three, two and zero controls respectively. The principal rank-four probe and radial gauge check completed. The saved control artifacts record two manufactured DAE cases, six row-operation tests, five kernel cases, eight preservation cases and eight terminal-dynamics cases. [Principal probe](outputs/principal-block.json) · [Checks and checkpoint digests](archive-status.json)

| Stage | Duration (s) | End on child clock (s) |
|---|---:|---:|
| Terminal kernel embedding | 114.582787 | 484.486543 |
| Snapshot publication | 2.276816 | 486.763429 |
| Constraint preservation | 75.123611 | 561.887093 |
| All 22 retained Euler rows | 50.028805 | 611.915939 |
| $KM$ | 234.080508 | 845.996489 |
| Total radial derivative $DM$ | 517.413152 | 1363.409768 |
| $KM-DM$ | 0.481470 | 1363.891347 |
| Free-row extraction | 0.000056 | 1363.891503 |
| Reduced dynamics, including the four preceding operations | 751.975571 | 1363.891531 |

The checker formed $W=KM-DM$ and extracted $K_{\mathrm{red}}$ from its free rows. The next check asks whether $M K_{\mathrm{red}}=W$ on **all** rows. The journal has a begin event for this product and no matching end or error event. The control reconstruction and discarded-response stages have no begin events. The child published neither a final closure export nor a summary or failure JSON. This leaves the reconstruction identity undecided and supplies no stationary germ or nonlinear obstruction. [Frozen reconstruction routine](sources/verification/oblate_stationary_reduction_20261128/terminal_dynamics.py) · [Progress journal](outputs/progress.jsonl)

## Snapshot and provenance

All three DAE checkpoint files are byte-identical to those of run `c2a1482f`. A JSON comparison of the two terminal snapshots found identical matrix entries ($C$ of shape $6\times14$, $K$ of shape $14\times14$, $B$ of shape $14\times0$, $M$ of shape $14\times8$), coefficient and derivative labels, free rows, and publication flags. The only differing fields are the run ID, attempt ID, identity object and identity digest. The current snapshot is 1,004,745 bytes with SHA-256 `7036a5b6bb2f4c3680a5b5877512f0ea845d19775d95ede06bcbf6f000411b5e`. Its flags describe publication before preservation. The subsequent checks are recorded in the logs. [Original snapshot](outputs/terminal-state-before-preservation.json) · [Comparison and review](archive-status.json)

The run executed source `c763d3a4624526b5d9b3ef995993cc1fe000554a` once. All 15 frozen source files match the executed Git revision, the materialized files, and their launch and exit digests. All Python files match the previous attempt. The completion, stationary-operator and reproducibility input hashes match the run's input directory and the earlier operator archive. The child reported Python 3.12.11 and SymPy 1.14.0. [Launch](launch.json) · [Child invocation](outputs/invocation.json)

The supervisor recorded 1800.2224154998548 seconds elapsed, child exit 143, and an empty child process group after SIGTERM. The framework recorded supervisor exit 1. A later process-table check also found group 26630 empty. The largest of 6,140 sampled group RSS values was 470,630,400 bytes (448.83 MiB), below the 2 GiB sampled threshold. The execution record lists zero cleanup and relay errors. The 100 progress events match across the child log, engine log and journal. [Resource samples](resource-samples.json) · [Engine log](engine.log) · [Framework exit](framework-exit.json) · [Framework status](framework-status.json)

Independent review `msg_fefbe5ce` corroborated the failure boundary, completed checks, input and source provenance, and resource totals. Archival verification used byte comparisons, Git objects, JSON, logs and process metadata. The separately prepared direct-state checker uses the earlier authenticated candidate and must establish its own identities and physical-radius certificate. [Archive manifest](archive-manifest.json)
