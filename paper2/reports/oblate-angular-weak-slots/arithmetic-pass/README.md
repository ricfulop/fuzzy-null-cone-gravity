# Passed actual-profile row arithmetic diagnostic

Run `b173211f-b492-4d14-9866-84342c279986`, frozen source `aeafd47aae703908dc72ba3dde534cc5afce86b9`, completed its configured `profile-arithmetic-probe` with exit zero. All three control batches and all twelve arithmetic pairs passed. This is a row diagnostic: no angular snapshot or complete profile block was assembled.

## Verified checks

The fixed target was L=2, m=0, r=3/2, epsilon=1/4, left/right profile families0/0. Left row10 is C_polar(2,0) with coefficient11/13 and profile r^2+11r/7+3/5; every one of the twelve right columns was included with the original coefficients and profiles.

Actual coefficient-times-profile-times-basis fields were differentiated before fixed-radius substitution. The run checked559 exact raw/compacted component equalities across13 jets, then twelve exact equalities of the independently constructed weak-potential densities. Nine right columns supplied exact nonzero witnesses: [0, 1, 2, 3, 4, 6, 7, 8, 10]. Both raw and compact scalars retained the strict rationality and endpoint/interior-pole guards and dx/2 quadrature at24/48/96 nodes. Maximum scaled48/96 refinement indicator was 5.4550544749133603e-15, below the unchanged2e-10 gate. Every order also passed the raw-versus-compact comparison. Refinement is numerical evidence, not a certified error bound.

Component compaction allows coordinate poles in individual tensor components; final assembled scalar pole exclusion is unchanged. Existing Legendre/quadrature, exact round weak/flux, rationality rejection and manufactured realification controls also ran, with the additional component-compaction support controls.

## Timing and resources

- Checker: 153.358560000s; supervisor: 154.875650709s, under300s.
- Peak sampled group RSS: 245514240 bytes = 234.140625MiB, across524samples, under2GiB.
- Component compaction plus559 equality checks: 1.686086542s.
- Sum of recorded normalization/integration operation times over all pairs/orders: raw99.509191326s; compact3.137111375s.
- Separate exact raw-minus-compact density checks: 48.273113081s.
- Final process group empty; a later PID/PGID/RSS-only process check also found no members.

For first-order24 evaluation across this row, raw expansion took13.155s and combined cancellation/QQ/Sturm checking80.501s; compact counterparts took0.423s and2.474s. Density construction itself took about0.055s raw and0.048s compact. These measurements localize the work to expansion and the grouped symbolic normalization/guard path, without distinguishing cancellation from its guards internally. They support testing component compaction in the full first-case profile oracle. They do not establish a general speedup: arms shared one Python process and possible SymPy caches, despite separate evaluator caches and denominator-cache resets.

## Provenance and interpretation

All twelve pair JSONs match the aggregate result exactly. The raw log, five runtime records, twelve pair files and six executed sources are preserved without modification. `timing-events.json` is a derived parse; `archive-status.json` records root reconciliation. `manifest.json` hashes every payload except itself. No environment-launch records are included.

Explicit per-run approval was user message `msg_37a2dc0f`. The subsequent standing instruction was “you dont need to ask me for approval till 10am”, recorded literally with America/New_York and15:00UTC. The message-ref presentation may display the assistant reply to that interrupted turn; this does not replace the captured user instruction or the independent explicit launch approval.

The original metric-coefficient Galerkin regulator remains distinct from the reduction-first Gaunt bridge. This pass does not complete the first angular snapshot or six-snapshot export, radial projection/operator, physical endpoint domain, finite K/E/D, finite-basis kernel preservation, a spectrum, an admissible seed, or nonlinear consistency. Production compaction remains unintegrated in this frozen source.
