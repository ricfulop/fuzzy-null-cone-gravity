# Component-compacted first angular snapshot

Status: independently source-reviewed and statically checked; unrun preparation.

This child builds on the [accepted row arithmetic diagnostic](arithmetic-pass/README.md), run b173211f-b492-4d14-9866-84342c279986 from frozen source aeafd47aae703908dc72ba3dde534cc5afce86b9. Independent result audit msg_09e12055 and archive audit msg_f102f3e5 verified all three control batches, 559 exact component equalities, twelve exact density equalities, nine nonzero witnesses and the quadrature gates. The row diagnostic produced no angular snapshot. Its recorded raw/compact normalization and integration times were 99.509/3.137 seconds, with shared-process cache caveats; component compaction plus equality checks took 1.686 seconds. This motivates integration, not a guaranteed full-case runtime.

The parent experiment remains frozen. This child tests the complete first angular snapshot only: L=2, m=0, r=3/2, epsilon=1/4. It retains all twelve spatial and three lapse labels, all four radial derivative slots and separately exported flux slots, all four profile-family combinations, exact profile-span checks, potential bulk Hermiticity and existing quadrature/refinement/reconstruction gates.

## Narrow change

The raw density/flux slot export remains unchanged. The candidate compacts only the independently built actual-profile objects, after symbolic differentiation and exact radius/epsilon substitution. Spatial objects begin as coefficient times actual radial profile times the original tensor basis, then pass through tensor_jet. All h, nabla, trace, trace-gradient and divergence components are compacted over QQ(I)(x) and checked exactly against their raw counterparts. The lapse scalar and its already computed radial derivative receive the same exact transformation and equality check.

All four family combinations then use the compacted actual-profile objects in the unchanged potential/lapse density and flux formulas. Their pairwise integrated sums must match the independent raw slot contractions. No exported-slot data produces the direct answer. Component coordinate poles remain permissible; final scalar rationality and endpoint/interior-root exclusion remain strict. Per-term raw-density repetition from the arithmetic comparison is not repeated: exact component equality and the original end-to-end slot reconstruction are the integration gates.

The new compact-profile-probe mode is default-off and rejects nonzero m. Existing full, direct-profile-probe and profile-arithmetic-probe behavior is retained. No half-root field generalization or change to the physics formulas is included.

Independent source review msg_3ee20a8f found no blocker. Static AST compilation (without scientific imports), JSON/XML parsing and unchanged-file byte comparisons passed; nine mocked supervisor tests passed during preparation. The preflight records the seven selected source hashes. Expected equality accounting is 2,064 spatial components plus 12 lapse components.

## Reporting and execution

The plan is setup1, controls3, angular_export1 snapshot and export1. The first control batch opts into the already tested component-compaction controls. Timings identify compaction by side/family, equality counts and the existing direct pair blocks/rows. The result is compact-profile-probe-checks.json in its own output directory; compact_profile_probe_passed is emitted only after the final result write. A completed first-case snapshot is persisted separately, with exact labels, slot/measure conventions, numerical matrices, controls and source hashes.

User msg_908bb84e directly states “you dont need to ask me for approval till 10am”. The local timezone is America/New_York; the authorization cutoff is November 28, 2026 at 10:00-05:00, equivalent to 15:00 UTC. Limits remain one local scientific process, 300 seconds wall and 2 GiB sampled process-group RSS. The latter is not an OS hard limit. The fixed launch command and isolated Python 3.12 / NumPy 2.3.5 / SymPy 1.14.0 environment remain unchanged. Runtime is unknown; 300 seconds is the cap, not a completion estimate. No GPU/remote budget is used.

This could establish only the first fixed-radius angular snapshot. The remaining five snapshots, any radial operator/projection, physical endpoint domain, finite K/E/D, kernel preservation, spectrum, seed or nonlinear consistency claim remain outside scope. The original metric-coefficient regulator remains distinct from the reduction-first Gaunt bridge.
