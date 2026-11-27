# A Finite-State Model of Quantum Gravity in a Causal Diamond

**Author:** Ric Fulop, Center for Bits and Atoms, Massachusetts Institute of Technology

This is a research manuscript and reproducibility archive for a finite-resolution model motivated by Venkatesa Chandrasekaran's null-cone construction. It contains analytic examples, symbolic and numerical checks, selected Lean lemmas, and Euclidean matrix simulations. It does **not** establish a continuum quantum theory of general relativity.

## Read the current draft

- [Single-column reading copy](paper_reading.pdf), built from [paper_reading.tex](paper_reading.tex).
- [Two-column technical copy](paper.pdf), built from [paper.tex](paper.tex).
- Both layouts use the same scientific source. The reading copy uses REVTeX in one-column mode; it is not an official Science submission template.
- [Targeted B200 follow-up report](reports/modal-b200-followup-report.md) and [archived follow-up analysis](verification/modal_hmc_followup_analysis.json).

The review incorporates the Open Research follow-up at `d288b95712c735c2ef1c3df293d1c87abef76849` and preserves its simulation driver, analysis, and raw result files. It retains the corrections to pair-state counting, benchmark precision, and numerical scope, and replaces the figure presentation with a self-contained narrative and regenerated vector artwork.

## What the current evidence supports

The total occupation budget makes the specified radiative sector finite; the corner spectrum and budget remain model assumptions. Symbolic identities and known black-hole benchmarks test particular parts of the construction. They do not prove that the finite-state Hamiltonian, full constrained gravitational action, and sampled Euclidean model are equivalent.

With matched measured zero-coupling controls at fixed scalar self-interaction, the strongest tested coupling reduces the matrix model's roughness observable by **47.68–48.02%** across four ranks. A larger lattice at rank 24 gives **47.621%**, a change of **−0.060 ± 0.047 percentage points**. This supports persistence over the tested volume change, not an exact universal percentage or a continuum extrapolation.

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

The canonical builder produces all 11 figures as vector PDF, editable SVG, and PNG previews. It reads archived results and evaluates the displayed deterministic formulas in memory; it does not launch GPU jobs or overwrite numerical result files.

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

MIT License — see [LICENSE](LICENSE).
