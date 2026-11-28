# Maintained supervisor and provenance review fixes

These are infrastructure changes prepared for PR #10. Historical reports,
frozen source snapshots, exported operators, failed-run evidence, and scientific
experiment refs are preserved. No scientific checker or manufactured symbolic
control was executed as part of this preparation.

The reducer authenticates each complete input byte buffer before parsing it and
records that buffer's digest. Its supervisor executes the authenticated source
copy, verifies originals and copied sources after cleanup, checks the wall limit
after observing child exit, and limits process-table sampling by the remaining
budget. The shared diagnostic wrapper's explicitly labelled input files remain
separate from copied source files in that validation.

The maintained reduction and germ supervisors archive child output before
nonblocking console forwarding. The reducer's stderr progress mirror can fail
without hiding the original failure record; failure of its authoritative journal
still propagates. Missing RSS samples are unknown, not a measured zero.

Both germ supervisors define a terminal cancellation boundary after certificate
validation. Signals received before that boundary determine the same terminal
status, metrics and exit decision. A new `terminal-evidence.json` marker hashes
both `execution.json` and `metrics.json` after both are written; a standalone
metrics file is provisional. Passing telemetry is emitted only afterward.
The views distinguish "not passed" from a pending state and distinguish a valid
germ from the separate question of a responding coefficient. The germ README
also documents the checker's mandatory denominator-regularity test at `r=1`.

The older stationary-extraction supervisor now catches SIGTERM/SIGINT during
dispatch and rejects source changes in its success decision. Source deletion,
unreadable supervisor bytes and cleanup exceptions produce failed evidence and
restore signal handlers. Its tests replace
all subprocess and signal operations; they launch no extractor. The standalone
reducer control entry point now refuses optimized Python so its scientific
assertions cannot be silently disabled. This guard was inspected statically;
the symbolic controls were not run.

The following commands were run from `paper1/` with Python's stdlib only. Each
returned exit code zero. Counts record the maintained suites at this review;
the earlier preparation archives keep their original test counts and hashes.

```sh
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_reduction_20261128 -p test_run_bounded.py -q
# 17 tests passed
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_reduction_20261128 -p test_review_infrastructure.py -q
# 6 tests passed
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_reduction_20261128 -p test_supervisor_lifecycle.py -q
# 2 tests passed, each exercising all three maintained reduction/germ supervisors
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_germ_20261128 -p test_run_bounded.py -q
# 26 tests passed
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_state_germ_20261128 -p test_run_bounded.py -q
# 27 tests passed
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_seed_20261127 -p test_progress_io.py -q
# 6 tests passed
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_stationary_seed_20261127 -p test_supervise_stationary_check.py -q
# 6 tests passed
python3 -B -W error::ResourceWarning -m unittest discover -s verification/oblate_terminal_embedding_probe_20261128 -p test_run_bounded.py -q
# 4 tests passed
```

The lifecycle tests start only synthetic Python sleep/pipe processes. They check
that a full console pipe cannot block forwarding and independently observe that
a timed-out descendant has exited. All other listed tests use mocked dispatch
or temporary-file fixtures. No test supplies a physical closure, radial germ,
discarded response, spectrum, or degree-of-freedom result.

Syntax validation uses AST parsing without project imports, JSON/XML parsing,
and `git diff --check`. Published archived review comments remain attached to
their original snapshots; their maintained counterparts carry these repairs.
