# Control-only field arithmetic for the stationary closure

This repair keeps the historical control-only Gauss–Jordan pivot sequence and
performs its arithmetic on exact `QQ(r,T)` field elements. The raw operator
export, radial derivation, selected equations, all 22 terminal retained rows,
discarded response, and explicit remaining controls are unchanged.

The preceding full-augmented variant retained equivalent constraints, but its
new state basis enlarged the expressions differentiated next. Its controlled
stop and complete evidence are preserved in
`../oblate-stationary-closure-basis-growth-2026-11-28/`. That implementation
performance issue does not establish failure or closure of the mathematical
DAE.

The repair converts each input entry once, rejects inexact floats, scans only
control columns, and applies paired swaps/scaling/elimination to the complete
control-plus-state row. The same input order and first-nonzero pivot rule give
the historical rational functions by induction over the operations. It does
not normalize residual state rows by pivoting them. Batched independent-row
selection remains and returns original candidate rows in order.

Six manufactured row-operation tests now include entrywise historical residual
rows, paired row swaps, unchanged state rows when controls are absent, empty
dimensions, independent-row selection, and float rejection in either block.
The two existing differentiated-constraint controls remain. They run inside
the next tracked invocation before actual reduction. Preparation itself used
only static AST/source inspection. The unchanged supervisor retains the nine
previously passed mocked lifecycle tests and preserved final evidence during
the controlled stop of the full-augmented attempt.

Independent source review found no material issue in the two-file arithmetic
repair. Its design and review references, final source/configuration hashes,
and the unchanged input pin are in [preparation-provenance.json](preparation-provenance.json).
Performance is not assumed and will be measured.

The same provisional closure experiment and fixed command are used:

```
bash paper1/verification/oblate_stationary_reduction_20261128/run.sh
```

The next invocation has one local scientific CPU child, a 900-second wall
budget, a 2-GiB sampled process-group RSS threshold, fresh evidence output, and
the standing approval cutoff of 2026-11-28 16:29 UTC. No GPU or remote job is
planned. A formal closure still does not certify an actual interval germ; the
separate prepared terminal checker must verify that further claim.
