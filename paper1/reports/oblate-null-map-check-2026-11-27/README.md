# Oblate nonlinear-null map and control check

**Result: passed in one authorized local invocation.** The prepared axisymmetric $N=3$ geometry, eight-channel coefficient map, and full linearized Einstein evaluator passed their independent controls. This is an implementation checkpoint for the second-order compatibility test. No admissible non-gauge first-order solution, evaluated nonlinear compatibility source, final Dirac classification, or physical degree-of-freedom count is supplied.

The invocation used the unchanged scientific sources from commit `53b06316547d0a5a47754f0de71f8c91c7a42a40`. Execution was separately approved by `msg_aa8cfa1c` and by the coordinated approval `msg_78c9db37`. The original [preparation checkpoint](../oblate-nonlinear-null-preparation-2026-11-27/README.md), [proposal](approved-run-proposal.json), and [preparation provenance](preparation-provenance.json) describe the earlier unexecuted stage; their unrun flags are historical.

## Recorded execution

The [full merged stdout/stderr](verification.log), [supervisor record](execution.json), [resource samples](resource-samples.json), [check results](map-checks.json), and [checker provenance](map-reproducibility.json) preserve the evidence.

| Item | Recorded value |
|---|---|
| Local date | 27 November 2026, America/New_York |
| UTC interval | 28 November 2026, 03:29:37.873785–03:29:40.798896 |
| Exit status | 0 |
| Supervised elapsed time | 2.924 seconds |
| Internal checker time | 2.108 seconds |
| Sampled peak process-group RSS | 116,703,232 bytes, approximately 111.3 MiB |
| Approved envelope | 300 seconds and 2 GiB; one invocation |
| Scientific runtime | Python 3.12.11, NumPy 2.3.5, SymPy 1.14.0 |
| Remote/GPU jobs; tracked experiment launches | 0; 0 |

The [stdlib supervisor](supervise.py) sampled process-group RSS approximately every 0.25 seconds and would terminate on an observed resource-budget violation. This is sampled memory monitoring, not a hard OS memory cap or a measurement of the true instantaneous peak. No budget stop occurred. The supervisor used Python 3.13.7, separately from the checker's isolated Python 3.12.11 environment. Numerical-library thread limits were set to one; their exact environment overrides are recorded.

All four pinned source/action/basis hashes matched before and after execution. The corresponding files are frozen under [sources](sources/). Neither scientific implementation was changed for this invocation, and no automatic retry was made.

## What passed

1. **Nonlinear reconstruction and changing tangent.** The direct coordinate pullback for $\phi\mapsto\phi+\eta f(t,r)x$ agreed with the coefficient reconstruction and with independent first and second Lie derivatives. Arbitrary time and radial profile derivatives were retained. Subtracting the embedding's second derivative left the expected angular defect, with scalar degree-four coefficient $8(r^2+\epsilon)f^2/35$. The base-dependent rotation control stayed exactly in the ansatz.
2. **Eight-channel coefficient map.** Exact round harmonic projections and complete reconstruction agreed with every predicted discarded scalar, tensor, radial-vector, and temporal-vector coefficient for $\xi=w(r)x\partial_\phi$. The rotation-generator defect was zero. Separate 24- and 48-node scalar quadratures reproduced $8/35$ with absolute errors $2.50\times10^{-16}$ and $1.19\times10^{-15}$, respectively.
3. **Full linearized Einstein and Euler-density controls.** The full retained pure-coordinate perturbation gave $\delta G_{\mu\nu}=0$. The off-shell conformal perturbation obeyed the independent four-dimensional identity
   $$
   \delta G_{\mu\nu}[2\psi g_0]
   =-2\nabla_\mu\nabla_\nu\psi+2g_{0,\mu\nu}\Box\psi.
   $$
   The nonzero $b=1$ test gave Euler-pairing density $-4(r^2+\epsilon)\partial_t\partial_r\psi$ before $dr\,dx/2$ integration. This checks the sign, volume density, index raising, and both symmetric off-diagonal entries. Discarded temporal residuals were retained.

The three PASS lines in the log cover these groups; the last log line explicitly states that no admissible non-gauge seed or nonlinear obstruction was evaluated. The conformal perturbation is an off-shell control, not an admissible physical seed. The coordinate-orbit harmonic leak alone is not a dynamical obstruction.

## Remaining nonlinear question

The prepared necessary condition is
$$
B_\xi[c_1]=E_{\rm full}^{(1)}[c_1]\!\left[Q_cDF_\xi(c_0)c_1\right]=0
$$
for a first-order retained solution tangent to a sufficiently smooth exact solution family. This invocation validates maps used in that condition on the stated tests. It does not evaluate $B_\xi$ on an admissible non-gauge solution.

The next dynamical calculation still requires a compatible radial domain or local solution germ, the complete retained linearized equations and preservation conditions, and a mixed-parity non-gauge seed with consistent time/radial jets. Pure-parity seeds have a symmetry-enforced zero for this probe. The present implementation examines only the axisymmetric member of the three extra background angular nulls. Background null directions have not been established as nonlinear gauge generators.

The separate quadratic-operator/domain work is a dependency, not part of this invocation. Its independently reported check failure does not invalidate these results or supply the missing admissible seed. Further scientific execution requires a new scoped decision.

## Reproduction and archive

The approved command, run from `paper1/`, was:

```sh
uv run --isolated --no-project --with numpy==2.3.5 --with sympy==1.14.0 python verification/oblate_nonlinear_null_20261127/check_prepared_maps.py
```

The frozen checker writes its two JSON outputs to the preparation directory. Their byte-identical copies are archived here and their original output locations are preserved. The historical preparation artifact is unchanged. The execution archive includes its own source snapshots, proposal, selected environment information, hashes, complete log, exit status, and sampled resource record. No prior scientific calculation was rerun.
