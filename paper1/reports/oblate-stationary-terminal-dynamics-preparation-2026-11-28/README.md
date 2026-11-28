# Exact-field terminal Euler and dynamics checks

The preceding [attempt](../oblate-stationary-dynamics-timeout-2026-11-28/README.md) completed preservation and all 22 original-row checks, but reached its 900-second limit during reduced-dynamics construction. Its partial checks are retained as evidence, with no closure or germ claim.

## Narrow implementation change

`terminal_dynamics.py` evaluates the two original Euler residuals over sparse QQ(r,T):

- ([E0|E1]+E2 Ka)M=0 for the embedded state;
- E2 Ba=0 for every remaining acceleration control.

E0,E1,E2 include all 22 original equations after removing only the X0 profile column. The separate discarded row remains a response and need not vanish. The production path converts K,B,M once for this terminal stage, and retains the original operator labels and derivative convention.

For reconstruction, compute W=KM-DM once with D=partial_r-10/(100r^2+1) partial_T. The ordered free-row selector J gives Kred=JW and Bred=JB. Both full identities MKred-W=0 and MBred-B=0 must hold; the check is not confined to free rows. Conversion, products, differentiation, extraction, full reconstruction and conversion back to SymPy expressions have separate timers.

Discovery, the kernel, preservation, the write-only intermediate snapshot, and the discarded-response calculation remain unchanged. No checkpoint loader or restart path is added.

## Controls and evidence

Eight new manufactured groups cover a valid nonzero discarded control response, a state defect in the X0 equation, a control-only defect outside the eight discovery rows, a variable embedding with valid nonzero control, independent non-free-row K and B defects, empty widths, and invalid coefficient rejection. They will run with all existing row, kernel, preservation and DAE controls before the actual reduction. No scientific code was imported or executed during this preparation; these new controls remain unexecuted.

Static AST/source comparisons, JSON/XML parsing and shell syntax checks are recorded in preparation-validation.json; source hashes are in preparation-provenance.json. The unchanged supervisor suite is not rerun.

## Next bounded invocation

One fresh local CPU process uses the fixed original completed operator export, the same 900-second wall limit and 2-GiB sampled process-group RSS threshold, and the approved cutoff 2026-11-28 16:29 UTC. Source, dependency and input identities remain recorded. The remaining runtime is unknown; a timeout remains possible.

Only completion of every terminal check can yield a formal closure export. The independent terminal certifier must subsequently verify all original equations and actual regularity on an open radial interval before an admissible responding germ or a Lagrangian nonintegrability statement. This preparation establishes no physical boundary domain, canonical Dirac rank or physical degree-of-freedom count.

Independent source review `msg_40088be0` found no blocker for one fresh bounded invocation after freezing both new files. It confirmed separate retained state/control checks, the full-row reconstruction identities, and unchanged discovery/preservation/response. Runtime and the new controls remain unverified until execution.
