# A Finite-State Model of Quantum Gravity in a Causal Diamond

**Author:** Ric Fulop, Center for Bits and Atoms, Massachusetts Institute of Technology

This is a research manuscript and reproducibility archive for a finite-resolution model motivated by Venkatesa Chandrasekaran's null-cone construction. It contains analytic examples, symbolic and numerical checks, selected Lean lemmas, and Euclidean matrix simulations. It does **not** establish a continuum quantum theory of general relativity.

This is **Paper I**, contained in `paper1/`. Run all commands below from this folder; paths in the manuscript and archived reports are relative to this folder unless a historical commit is specified. The companion manuscript has a separate [Paper II folder](../paper2/README.md).

## Read the current draft

- [Single-column reading copy](paper_reading.pdf), built from [paper_reading.tex](paper_reading.tex).
- [Two-column technical copy](paper.pdf), built from [paper.tex](paper.tex).
- Both layouts use the same scientific source, including the files under `sections/`. The reading copy uses REVTeX in one-column mode; it is not an official Science submission template.
- [Targeted B200 follow-up report](reports/modal-b200-followup-report.md) and [archived follow-up analysis](verification/modal_hmc_followup_analysis.json).

The review incorporates the Open Research follow-up at `d288b95712c735c2ef1c3df293d1c87abef76849` and preserves its simulation driver, analysis, and raw result files. It retains the corrections to pair-state counting, benchmark precision, and numerical scope, and replaces the figure presentation with a self-contained narrative and regenerated vector artwork.

The [canonical manuscript integration](reports/canonical-manuscript-integration-2026-11-28/README.md) brings the certified oblate primary/secondary ranks and linearized lapse results into both current layouts. The [direct stationary-germ certificate](reports/oblate-stationary-state-germ-pass-2026-11-28/README.md) now supplies a responding local linearized solution at $N=3$, $\epsilon=1/100$. Its separate Noether linkage gives a local nonlinear obstruction for a specified mixed tangent, included in both layouts. A global radial boundary domain and final Dirac classification remain open.

## Implemented action and clock repairs

The current draft includes citation PR #5 and replaces the failed Appendix A action with an explicit Einstein--Hilbert angular Galerkin **bulk** construction. All metric components are varied before null gauge fixing; exact spherical and quadratic checks and a finite-amplitude nonspherical check are supplied. This reconstruction uses pointwise products of harmonic symbols and is distinct from the archived fuzzy-product simulation model. Finite-rank constraint closure, null-boundary quantization, and equivalence of the models remain open.

The clock now has a proved uniform oscillator lower bound and uses the actual total-occupation compression before exponentiation. With the stated area, rank, and classical-clock assumptions, a new proposition bounds the number of area sectors. Old unrestricted-occupation exponential estimates and the unsupported scalar-QEI attribution are retired. The recomputed vacuum second-moment bound at s=5000 is at least 0.9999989079 for a clock threshold twice its classical value.

- [Repair calculations and provenance](reports/repair-audit-2026-11-27/README.md)
- [Calculation supplement PDF](reports/repair-audit-2026-11-27/repair-calculations.pdf)
- Reproduce all seven checks: `python verification/repair_audit_20261127/run_all.py`

## What the current evidence supports

The next research objective is tracked in [GOAL.md](GOAL.md). The [first nonlinear-consistency checkpoint](reports/nonlinear-consistency-2026-11-27/README.md) supplies exact obstruction tests, an angular deformation candidate, and a reduced boundary-gluing example. The full completion tests remain open.

The [actual-action canonical analysis](reports/canonical-analysis-2026-11-27/README.md) now gives the derived ADM lapse/shift, a regular spatial Legendre map, and an exact nonzero primary–secondary bracket on explicit vacuum data. The scalar lapse Hessian has ranks 8 and 12 at $N=3,4$, with a general-rank proof. This is a result on the primary–secondary surface; tertiary stabilization and a physical degree-of-freedom count remain open. Reproduce its local checks with `python verification/canonical_analysis_20261127/actual_adm.py` using NumPy and SymPy.

The [first Dirac preservation test](reports/dirac-stabilization-2026-11-27/README.md) finds an instantaneous solution for all retained secondary constraints, including radial terms. It also proves that this fixed vacuum datum has no exact continuum Einstein development staying in the $N=3,4$ coefficient ansatz: a discarded trace appears in the spatial acceleration for every retained temporal rate. That calculation alone does not settle continued projected evolution or final constraint classification. Reproduce with `python verification/dirac_stabilization_20261127/stabilize_witness.py` (NumPy and SymPy).

The [second preservation test](reports/higher-order-consistency-2026-11-27/README.md) resolves the continuation question for that particular datum: no sufficiently smooth projected solution passes through it at $N=3,4$. A retained angular-shift constraint acquires a nonzero second time derivative independent of the remaining temporal rates. The earlier bracket ranks therefore remain primary–secondary-surface results and cannot be promoted using this witness to ranks on a final dynamically admissible surface. Other initial data and the final Dirac classification remain open. Reproduce with `python verification/higher_order_consistency_20261127/second_preservation.py` (NumPy and SymPy).

The [admissible-family calculation](reports/admissible-branches-2026-11-27/README.md) supplies an exact stationary oblate-Minkowski solution in the $N=3,4$ ansatz. Its angular-shift primary–secondary bracket is nonzero for sufficiently small nonzero oblateness, giving a Poisson rank lower bound of two on an actual retained trajectory. The full constraint rank and physical degree-of-freedom count remain open. This family approaches round zero-momentum data; it does not continue the excluded rank-one witness. Reproduce with `python verification/admissible_branches_20261127/oblate_branch.py` (NumPy and SymPy).

The [all-mode oblate calculation](reports/oblate-constraint-rank-2026-11-27/README.md) sharpens the result to exact temporal-Hessian ranks 18/26 and primary-plus-secondary Poisson ranks 36/52 for sufficiently small nonzero deformation. The rank proof includes all magnetic modes and the mixed secondary brackets. Its original failed lapse follow-up remains archived.

The [corrected lapse-preservation check](reports/oblate-lapse-consistency-2026-11-27/README.md) now passes. The defect was an unevaluated symbolic zero in a comparison; no harmonic coefficient or descriptor changed. Exact retained-gradient families and complementary radial forms give one/four redundant profile families and eight/twelve active linearized lapse channels at $N=3/4$, with a small-deformation proof on the stated continuous radial form extension. The common canonical boundary domain, surrounding nonlinear constraint stratum, final Dirac classification, and degree-of-freedom count remain open. Reproduce with `python verification/oblate_constraint_rank_20261127/lapse_prolongation.py` using NumPy and SymPy; the report gives the pinned-environment command.

The [second-order null preparation](reports/oblate-nonlinear-null-preparation-2026-11-27/README.md) derives a necessary nonlinear preservation test near the oblate solution, starting with one extra $N=3$ angular null direction. Its separately approved [single map/control check](reports/oblate-null-map-check-2026-11-27/README.md) now passes: nonlinear reconstruction, all eight discarded channels, and the full linearized Einstein evaluator agree with independent zero and nonzero controls. The test retains discarded temporal residuals. A finite coordinate-orbit leak alone is not a dynamical obstruction. The later responding-germ certificate supplies the admissible polar leg for a local mixed-parity test; it does not select a global canonical boundary domain.

The [stationary seed preparation](reports/oblate-stationary-seed-preparation-2026-11-27/README.md) reduces one possible mixed-parity test to eight stationary radial profiles and the response $6\mathscr E_{X,4}-\mathscr E_{C_{\mathrm p},4}$. A known odd coordinate mode supplies the second parity; the even profile must still solve every retained linearized equation throughout an interval. The single approved [extraction attempt](reports/oblate-stationary-system-check-2026-11-27/README.md) stopped at its five-minute wall limit after five column-completion messages. It produced no operator export and did not reach the final controls; no assertion failure was reported. The separately approved [instrumented rerun](reports/oblate-stationary-system-recheck-2026-11-27/README.md) subsequently passed all eight columns and final controls in 463.842 seconds, exporting the exact three 23-by-8 radial coefficient matrices. All 22 retained rows and the separate discarded response are preserved; 772 independent quadrature comparisons pass. That operator export alone supplies no responding seed. The later independent direct-state certificate constructs a forward solution germ without completing the reduced-coordinate reconstruction.

The total occupation budget makes the specified radiative sector finite; the corner spectrum and budget remain model assumptions. Symbolic identities and known black-hole benchmarks test particular parts of the construction. They do not prove that the finite-state Hamiltonian, full constrained gravitational action, and sampled Euclidean model are equivalent.

The [fixed-deformation principal probe](reports/oblate-stationary-principal-2026-11-28/README.md) certifies rank four for the stationary second-radial-derivative block after the verified radial gauge removal. It does not solve the full differential-algebraic system. The [incomplete closure attempt](reports/oblate-stationary-closure-attempt-2026-11-28/README.md) and [reviewed retry preparation](reports/oblate-stationary-closure-retry-preparation-2026-11-28/README.md) preserve that distinction. A separately prepared terminal checker must verify all original equations and actual local regularity before any radial solution germ is claimed.

The [terminal-embedding timeout](reports/oblate-stationary-closure-terminal-timeout-2026-11-28/README.md) preserved three elimination checkpoints and clean shutdown at 900 seconds, without a closure export. A separate [kernel diagnostic](reports/oblate-terminal-embedding-probe-2026-11-28/README.md) passed in 117.565 seconds: the saved 6-by-14 constraint matrix has rank six and a verified eight-column algebraic kernel. This supports a terminal-algebra implementation repair; it does not establish constraint preservation, all original Euler equations, a solution germ, or a physical mode count.

The subsequent [full retry](reports/oblate-stationary-preservation-timeout-2026-11-28/README.md) completed the repaired kernel in 114.775 seconds but reached the 900-second limit during terminal constraint preservation. Its shutdown was clean and all partial evidence is preserved. The next repair changes only terminal preservation arithmetic; all original equations and the independent germ certificate remain required.


The next [terminal verification attempt](reports/oblate-stationary-dynamics-timeout-2026-11-28/README.md) passed exact constraint preservation and all 22 original retained-row checks, then stopped at 900 seconds while constructing reduced dynamics. Its authenticated intermediate state is preserved; full dynamics reconstruction, the discarded response, and interval/germ certification remain unfinished.


The [exact-field dynamics attempt](reports/oblate-stationary-field-dynamics-timeout-2026-11-28/README.md) reduced the original-row check to 50.614 seconds and completed KM in 240.467 seconds, but stopped at the 900-second limit during differentiation of the embedding. The next invocation extends the wall budget with scientific code and checks unchanged. Reduced dynamics, response and interval/germ certification remain open.

The [shared local weak-form validation](reports/shared-oblate-weak-validation-2026-11-28/README.md) now passes all eight batches for the original Paper I action: representative exact sampled potential/lapse/Green checks, three fixed symbolic gauge vectors, and nonzero sign/ablation controls. This is separate supporting evidence; it does not certify the stationary DAE, a physical radial domain, an admissible seed, or nonlinear consistency.

The [independent direct-state check](reports/oblate-stationary-state-germ-interval-timeout-2026-11-28/README.md) subsequently passed exact kernel, preservation, physical-lift and all-22 retained-row identities, and found a formally nonzero discarded response. It reached its 900-second limit during the physical-radius interval search, so that invocation certified no physical minor or responding germ. The [conditional linkage note](reports/oblate-stationary-nonlinear-linkage-preparation-2026-11-28/README.md) records the local nonintegrability implication subsequently applied to the completed germ certificate.

With matched measured zero-coupling controls at fixed scalar self-interaction, the strongest tested coupling reduces the matrix model's traceless-field variance Q2 by **47.68–48.02%** across four ranks. A larger lattice at rank 24 gives **47.621%**, a change of **−0.060 ± 0.047 percentage points**. This supports persistence over the tested volume change, not an exact universal percentage or a continuum extrapolation.

Extended local collective moves improve local-area diagnostics at ranks 32 and 48: rank-normalized split-R-hat is **1.0098** and **1.0021**, respectively. The larger-lattice rank-24 run has well-sampled roughness but poorly mixed local area (**R-hat 1.3268; bulk ESS 20**). The main four-rank comparison still uses archived summaries for the coupled ensembles. The matrix exponential is positive by construction; positivity alone does not demonstrate the absence of spacetime caustics.

The six Lean files certify their encoded algebraic identities and conditional implications. Their physical hypotheses and the identification of abstract variables with gravitational operators require separate arguments. In particular, the conditional Appendix E.2 lemma assumes the relevant implication rather than formalizing the Hilbert-algebra reconstruction.

## Build the manuscript

Use a TeX distribution with `revtex4-2`, AMS packages, `graphicx`, `bm`, `mathrsfs`, `microtype`, and `hyperref`:

```bash
pdflatex paper.tex
pdflatex paper.tex
pdflatex paper_reading.tex
pdflatex paper_reading.tex
```

In Overleaf, select `paper_reading.tex` or `paper.tex` as the main document. A separate concise Science submission can be prepared once its scientific claims and main/supplement split are settled.

## Reproduce the manuscript figures

After reproducing the repair calculations above, the canonical builder produces all 11 figures as vector PDF, editable SVG, and PNG previews. It reads archived results and evaluates the displayed deterministic formulas in memory; it does not launch GPU jobs or overwrite numerical result files.

```bash
python verification/plot_manuscript_figures.py
```

Install `numpy scipy sympy matplotlib`. The pair-amplitude panel explicitly enforces the two-particle cost `2*omega <= B`; the old archived single-particle-cutoff JSON is preserved for provenance. The Schwarzschild spectrum is recomputed on the documented grid. Legacy benchmark plotting files supply some panel data to the canonical builder; run the canonical builder to reproduce the manuscript layout.

- [Figure captions and source mapping](figs/README.md).
- [Captioned figure review](figures_review.pdf), built from `figures_review.tex`.
- The main text now explains the original proposal, the continuum gaps, and the added model assumptions without requiring the source paper beside it. Detailed source-section and equation locations are in an appendix.
- Six main figures emphasize the construction and evidence; five additional benchmark figures appear in the appendix. Sampling diagnostics distinguish roughness from local area, and captions separate model tests from continuum claims.

## Simulation and formal sources

- `verification/modal_fsg_hmc.py`: Modal simulation driver, including the follow-up local collective moves.
- `verification/analyze_followup_hmc.py`: follow-up statistics and diagnostic calculation.
- `verification/modal_hmc_followup_batch1_results.json` and `verification/modal_hmc_followup_b2_both_results.json`: raw follow-up results.
- `lean4/`: six Lean modules, with 50 selected statements cataloged in the manuscript.
- `figs/`: vector PDFs and editable SVGs.
- `reports/`: historical reviews and the targeted follow-up report. The manuscript and this README reflect the newer numerical results.

To check the Lean files using the archived toolchain version (4.36.0):

```bash
for f in lean4/*.lean; do
  lean "$f"
done
```

GPU runs incur compute charges and are separate from figure generation. Consult the driver and follow-up report for the specific run configurations.

## License

MIT License — see [LICENSE](../LICENSE).

The [cached direct-state attempt](reports/oblate-stationary-state-germ-serialization-failure-2026-11-28/README.md) passed 19 controls and all exact identities, then selected a responding interval candidate at radius 1 with 25 arctangent terms. It failed while converting a determinant interval endpoint beyond Python's 4,300-digit limit and exported no physical certificate. The narrow serialization repair preserves the exact equations, interval arithmetic, search order and resource limits.

The [completed direct-state retry](reports/oblate-stationary-state-germ-pass-2026-11-28/README.md) exports an authenticated physical certificate at $r_0=1$. Exact kernel, preservation, profile-lift and all-22-row identities pass; rational interval arithmetic separates all 57 proof denominators, an original rank-six constraint minor and the selected response numerator using the same 25-term arctangent enclosure. The second free initial coordinate supplies the responding germ. The run passed 20 mathematical controls in 262.778 seconds of supervisor time, with 229.875 MiB peak sampled process-group RSS and clean shutdown. The [input-domain audit](reports/oblate-stationary-input-domain-audit-2026-11-28/README.md) supplies the analytic interpretation of the pinned inputs. The new manuscript proposition combines this certificate with the existing odd coordinate leg and Noether identity to exclude exact retained families with that specified mixed tangent, including arbitrary time-dependent second-order corrections. It makes no global-boundary, $N=4$, final Dirac or physical-mode-count claim.
