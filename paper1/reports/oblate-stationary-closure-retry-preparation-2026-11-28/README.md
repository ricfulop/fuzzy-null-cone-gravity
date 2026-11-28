# Exact stationary closure retry preparation

The first closure attempt is preserved at result commit
`e9b18f89da73047d2fb79705cbfe95da65c2d777`. Its log ends during the second
elimination step, followed by a process-group cleanup permission error. No
closure export was produced. A later process-table inspection confirmed that
the original child group was empty. The missing final resource and stop records
are disclosed in that attempt's archive, not reconstructed.

This preparation retains the same scientific question and upstream input:
close the stationary radial DAE at `epsilon=1/100` over `QQ(r,T)`, preserve all
22 original retained equations and the separate discarded response, and leave
every unsolved radial control explicit. The original stationary export and
completion-marker pins are unchanged.

The arithmetic repair uses SymPy 1.14's fraction-free row reduction after
clearing denominators. In the augmented system `[A | C]`, only pivots in `A`
solve controls; every state-only row is retained. Changes to the upper
substitution by residual state constraints preserve the same solution set on
the constraint surface. Selection of independent constraints uses one
transpose-pivot computation and returns original candidate rows in order.
The radial derivation and all original-row terminal identities are unchanged.
Five new manufactured algebra tests precede the two existing differential
constraint controls. No performance improvement has yet been measured.

The supervisor repair distinguishes a permission error on an empty process
group from a live or unknown group. Bounded cleanup, final group disposition,
and evidence writing survive cleanup exceptions. Nine mocked stdlib lifecycle
tests cover these cases, including the timeout regression; their validation
record is stored alongside this note. Those tests use no real child launches,
process signals, or scientific modules.

The fixed entry command remains:

```
bash paper1/verification/oblate_stationary_reduction_20261128/run.sh
```

The next invocation uses one local scientific CPU child, a 900-second wall
budget, and a 2-GiB sampled process-group RSS stop threshold. It uses the
standing overnight approval through 2026-11-28 16:29 UTC and must finish within
that window. No GPU or remote job is planned. Fresh output directories preserve
the failed invocation. The run display has an interim principal-rank headline;
fine stage timings distinguish control solving, substitution, independent-row
selection, checkpointing, and constraint differentiation.

All symbolic equivalence controls and the actual DAE retry run within the
tracked invocation, after freezing this source. At preparation time only static
source/AST review and mocked stdlib supervisor tests have executed. The
independent terminal germ certifier is prepared separately and is not invoked
by this run. A formal closure would still need its regularity and original-row
certificate before a solution germ or mixed nonlinear obstruction is claimed.
