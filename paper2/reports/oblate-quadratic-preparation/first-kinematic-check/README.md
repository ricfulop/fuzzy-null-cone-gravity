# First oblate kinematic check: failed reference comparison

The user authorized one local invocation in `msg_78c9db37`, together with an independent Paper I map/control check. This report covers only the Paper II invocation. Input source was the clean preparation commit `ee7d0a359a7b67a6b4444653b04fd68248e6b58f`; the original Paper I physics reference remains `157ec0f04866ea36509b325ec65f0a68890f9b17`.

## Observed outcome

The checker exited **1** after **2.02 seconds** wall time, as recorded by the system timing wrapper. The reported maximum resident set size was **79,265,792 bytes**. The full [verification log](verification.log) is preserved, including the traceback and timing output.

The log records a pass for geometry and pinned-reference comparison in `L=2,m=0`. Execution then stopped in the `r=1,m=1` reference comparison at tensor entry `(2,2)`. It did not complete the all-sector check or reach the later DeWitt-density, integrated Schur, manufactured Legendre, and realification controls. No success JSON was produced.

The failed residual was

```text
3*I*sqrt(-(x - 1)*(x + 1))*(x - 1)*(x + 1)
- 3*I*sqrt(1 - x**2)*(x**2 - 1)
```

## Interpretation from the recorded expression

The polynomial identities `-(x-1)*(x+1) = 1-x**2` and `(x-1)*(x+1) = x**2-1` make the two terms equal. The displayed residual is therefore algebraically zero. This points to a symbolic representation/cancellation problem in the verifier's `exact_zero` predicate, rather than evidence of a different tensor coefficient. This interpretation uses the printed expression only: no additional symbolic diagnostic, code correction, or rerun was performed.

The observed outcome remains a failed checker invocation. It does not validate the unvisited checks or establish a physical spectrum, an admissible linearized seed, a common evolution domain, or nonlinear consistency.

## Preserved execution contract

The command was run from `paper2/`, with complete stdout/stderr redirected to `verification.log` and `/usr/bin/time -l` recording resource use:

```bash
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 \
  python verification/oblate_quadratic_preparation/verify_kinematics.py \
  --output-dir reports/oblate-quadratic-preparation/first-kinematic-check
```

[Provenance](failure-provenance.json) records the exact input commit, source hashes, command, approval, and execution boundaries. The dependency versions above are the requested pins; the checker's success-only environment manifest was not reached. The original preparation archive remains unchanged. There was one scientific invocation, no retry, and no remote or GPU compute.

## Proposed next action, awaiting approval

Inspect and narrowly correct exact-zero normalization so equivalent polynomial radicands are compared canonically on the stated chart. Include an exact-zero regression and a deliberately nonzero control, preserve this attempt, and make one separately authorized rerun of the same kinematic coverage into a new output directory. No action or harmonic coefficient change is proposed. A successful rerun would still leave curvature/Green-identity validation and radial-domain choices as subsequent gates.
