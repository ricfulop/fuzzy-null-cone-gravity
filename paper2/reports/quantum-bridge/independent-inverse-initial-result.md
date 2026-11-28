# Independent inverse controls: initial scalar resolution failure

Run `226eb930-9667-4bdf-91b6-ae8a4e0552d4`, source
`c771712679e729c9e60c7beb26a17cf71418d138`, completed all 26 angular
controls and all nine scalar controls, then stopped at the scalar acceptance
gate. Two scalar cases failed. The York stage was not reached. The script
recorded 3.3839 seconds and 45,776,896 bytes sampled peak RSS.

The complete immutable output is
[independent-inverse-initial-results.json](data/independent-inverse-initial-results.json);
[provenance](data/independent-inverse-initial-provenance.json) records the
source, result/log SHA-256 hashes, and failed cases. The original experiment
branch remains frozen at this source. The
[initial protocol](independent-inverse-primitive-protocol.md) is preserved.

All declared angular, exact polynomial, Cartesian-jet, Gaussian-profile and
Kelvin-exterior controls passed, including full magnetic-mode checks through
l=10. The suite counted 31,146,657 operations under its 50-million cap.
These are implementation checks on their stated finite fixtures, not an
independent gravitational functional or physical cutoff limit.

The scalar failures were both third derivatives of the regular radial
amplitude G(s)=u_l(r)/r^l, s=r^2, evaluated immediately to the right of the
fixed r=1/2 origin patch. The exact interior reference has G_sss=0.

| Radial intervals at R=4 | l=8 maximum regular-jet error | l=10 maximum regular-jet error |
|---:|---:|---:|
| 256 | 0.1572746427 | 0.4467857992 |
| 512 | 0.03608143949 | 0.09724646751 |
| 1024 | 0.008618779249 | 0.02257143544 |

The unchanged final derivative gate is 0.002. Both cases failed this gate
also in the fixed-h outer-radius slice. Every other finest-level gate in
these cases passed: raw radial jets, origin coefficient, independent
fourth-order raw-grid PDE diagnostic, core moments, source-bound checks and
outer moments. The other five scalar fixtures, exterior-source omission
ablation and domain/join guards passed. In particular, omitting the bounded
dipole exterior source produced a visible error, rather than accidentally
passing the physical reference.

The errors decrease by approximately four per radial doubling, consistent
with the declared second-order solver and derivative reconstruction. This
supports a resolution refinement; it does not establish the next result in
advance. The planned child uses 1024/2048/4096 intervals at R=4 and retains
R=2/4/6 at fixed h=1/1024. It preserves the physical origin patch, every
tolerance, all fixtures, the York schedule, and existing work/memory/wall
caps. The planned scalar work is 140,331 unknowns, below the unchanged
200,000 limit. The shared 1024-point case supplies a reproducibility check.

This is a negative result for the initial requested accuracy, so it is
preserved as an answered experiment rather than overwritten by a retry.
No full mixed-H6 reference, four-point loop, microscopic corner-state
derivation, or controlled physical large-N limit follows from this run.
