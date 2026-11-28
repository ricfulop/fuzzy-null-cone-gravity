# Selected oblate lapse endpoint identities: passed

Run `95e69683-94dd-433d-b081-aba5b9e32fec`, source `cae4d6104c1d7e3b65693331787a2c3d75eae674`, passed all 11 named cases and eight control records. The single local invocation took 7.896592 seconds of supervisor time and peaked at 106,921,984 bytes (101.96875 MiB) of sampled process-group RSS across 27 samples. The child exited 0, the final process group was empty, and no signal, relay or cleanup error was recorded. The limits were 300 seconds and 2 GiB sampled RSS. [Execution](execution.json), [certificate](outputs/lapse-endpoint-checks.json), [raw log](verification.log).

The calculation fixes the original Paper I coefficient ansatz at $N=3$, $\epsilon=1/100$, $m=0$, $\ell=2$, on the spatial shell $[1,2]$. It uses $Y=P_2(x)$, the retained variation $U_{AB}=Yq_{AB}$, and independent radial profiles $f=F(r)Y$, $g=G(r)Y$. Write $\rho=r^2+\epsilon$ and $\Sigma=r^2+\epsilon x^2$. The checker verifies $\det\sigma=\Sigma^2$ and zero background Ricci curvature, differentiates the metric inverse and connections to construct the scalar-curvature variation, and compares the two separately constructed bulk densities with their coordinate flux divergences.

After imposing simultaneous Dirichlet endpoint jets, the radial coordinate fluxes are exactly

$$
\rho F'Y^2\left(\frac1\Sigma+\frac1\rho\right),
\qquad
-\frac{\rho^2 r}{\Sigma}\left(\frac1\rho+\frac1\Sigma\right)F'G'Y^2.
$$

These are the scalar-curvature Green flux and two-Hessian flux respectively, before the common $1/(8\pi)$ and $1/(4\pi)$ prefactors in the reviewed derivation. The inner endpoint has the opposite outward orientation. Both angular pole fluxes vanish, both clamped endpoint controls vanish, and the wrong-sign scalar and wrong-factor Hessian controls remain nonzero. The [frozen checker](sources/verification/oblate_lapse_endpoint_20261128/check_lapse_endpoints.py) shows the independent construction; the [published preparation](sources/reports/oblate-common-domain-preparation-2026-11-28/README.md) states the analytical Green identities and their domain qualifications.

The exact sign controls use $x=0$ and unit nonzero slopes. They are point controls, not a separate global angular sign certificate:

| Flux before its prefactor | Outer endpoint, $r=2$ | Inner endpoint with outward orientation, $r=1$ |
|---|---:|---:|
| Scalar curvature, $F'=1$ | $801/1600$ | $-201/400$ |
| Two Hessians, $F'=G'=1$ | $-321201/320000$ | $20301/40000$ |

All nine checker stages completed; the [timing summary](timing-summary.json) preserves their raw durations. The [engine log](engine.log) and child log agree on all 28 raw events and 27 child progress reports. The engine adds source setup and the final acceptance report. The initial notification omitted the final counts while the log was being collected; the archived [final view](framework-view.txt) and [status](framework-status.json) show pass, 11 cases and eight controls, with no open flag.

The certificate binds checker SHA-256 `24f46a839f99a6c511fd539999a8a5ff618c49c6da2dcb1de7f9d1e1e44ad168` and configuration SHA-256 `41c2f8e898685cdfd34d3968d4b053874b5a4f605174292d23f4f8533056f58e`. Its actual runtime is Python 3.12.11 and SymPy 1.14.0. All eight checker-directory files and the published report match their Git blobs and their pre/post-run snapshots. The report is separately pinned to `ae32314d27e20ac682c5ed1b71b4d7f11bee7133`; its derivation sources were reviewed at `773ccc4e3e8a1a89e06086163042fbf45e171c76`. The [terminal marker](terminal-evidence.json) authenticates execution and metrics, and the execution record authenticates certificate SHA-256 `5f8201fd2362b6fd4e360599f57bda4e1017b99efb3bacab3cc61fa530fdd269`.

[Independent review](review-notes.json) and the [archive status](archive-status.json) distinguish the accepted mathematical scope from the remaining work. This result validates the selected original lapse/adjoint endpoint identities and trace controls. Dirichlet lapse values alone leave a normal-derivative/area-trace term; the two-Hessian pairing has a separate mean-curvature boundary term. Specifying a common function-space domain for all original generators, proving its preservation under evolution and stabilization, and relating it to physical null-boundary data remain necessary. No final Dirac classification or physical degree-of-freedom count follows. The previously accepted stationary-germ and specified-tangent nonintegrability result is unchanged.

This archive uses an explicit whitelist of scientific output, the supervisor's own execution records, source snapshots and safe framework summaries. It contains no framework environment-launch record. Archival work only read, copied, parsed and hashed evidence; it ran no scientific code. The [manifest](archive-manifest.json) covers every payload other than itself, and the durable artifact copy is byte-identical.
