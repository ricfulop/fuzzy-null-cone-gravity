# Reviewed manuscript revision

Based on main at 589a99b9a5fa5695ebaf20cd3082b76f94a02346. The title and plain-language abstract had already been adopted upstream. This branch extends that revision using the completed production data. Additional Modal calculations were running during editing.

## Numerical interpretation

- The abstract reports the directly measured 40–44% roughness decrease when coupling rises from 0.5 to 2 with the quartic scalar term held fixed. The body also retains the 47.69–48.10% comparison with the analytic Gaussian baseline. The new plot propagates both independent ensemble standard errors. These intervals describe Monte Carlo sampling uncertainty and exclude finite-volume systematics.
- At N = 48, there are 331,392 coupled variables per independent chain. Four parallel chains have 1,325,568 variables in total.
- Ordinary archived R-hat is observable-specific. Roughness is close to one at strong coupling, while local mean area reaches 1.181 at N = 48 on the scaled branch and 1.416 at N = 32 at fixed coupling. The present plots reproduce archived diagnostics. A rerun with saved chain histories will support rank-normalized split R-hat and bulk/tail ESS.
- Positive area eigenvalues follow from exponentiating a Hermitian matrix. The reported extrema are means of configuration extrema sampled at one lattice site. The manuscript describes their narrowing with the convergence qualification attached.
- Mass-fit intercepts vary in sign with rank. The revised manuscript reports them as estimates, with uncertainties and continuum extrapolation pending.
- The invariant roughness Binder statistic alone cannot determine spontaneous breaking of SU(2). The revised text states this limitation.
- Early measure comparisons have matched lattices at N = 8, 16 and 24. The N = 32 pair has different temporal lengths. The manuscript records these conditions and avoids an unconditional phase-transition conclusion.

## Pair-state budget and checks

The same-mode two-particle intermediate state costs 2 omega. Equation (12), its verification script, regenerated JSON and corresponding figure now use 2 omega ≤ B. At B = 4 the pair amplitude vanishes. The 31-state radiative example has 30 one-particle states and one vacuum. The O(s^-2) tail argument uses L = B/2 and continues to preserve the leading asymptotic terms. The numerical scan gives s² times the reference difference about 0.987 to 0.939, rather than the earlier single-mode-cutoff coefficient about 0.33.

The original four two-loop numbers are retained as single-mode-regulated diagnostic values. Their reproduction has not inserted a total-occupation projector into every many-particle intermediate state. The manuscript now makes that scope explicit.

The Page calculation gives alpha = 3.8366580527e-6 when reproduced on the archived grid. Page’s 1976 Table I reports the two-helicity total as 0.038 × 10^-4. The percentage 0.017% used the more precise reference stored in the archived code. The manuscript now compares with the printed precision. [Page, Physical Review D 13, 198, Table I](https://doi.org/10.1103/PhysRevD.13.198).

The one-loop comparison plot displays both the five-radius fit and the final doubling. The manuscript table identifies its final-doubling estimator. Earlier numerical outputs remain available in the Git history.

## Artwork and manuscript format

The artwork uses 18.4 cm full-width layouts, sans-serif lettering, embedded fonts, uppercase 10 pt bold panel labels, at least 7 pt ordinary text, a colorblind-oriented palette with distinct symbols, and vector PDF plus editable SVG exports. Titles are placed in captions. Grid lines and minor ticks are removed. Scientific curves are calculated or taken from archived data. The introductory light paths are explicitly schematic.

Science’s live instructions returned HTTP 403 during this revision. The figure settings follow the accessible transcription of its archived initial/revised-manuscript guidance. The manuscript remains in the repository’s REVTeX format. Submission category, manuscript length, display-item limits and other full-journal requirements need a separate Science submission conversion. [Science instructions](https://www.science.org/content/page/instructions-preparing-initial-manuscript) and [accessible archived-guidance transcription](https://github.com/surPoudel/make-my-figure/blob/main/journal_preset_research/official_guidelines_notes.md).

## Pending simulation updates

The matched zero-coupling controls, larger-lattice pair and longer area chains are pending. Their results should be incorporated after checking the exact configurations and chain diagnostics. The archived scalar j = 2 correlator uses a variable overwritten with a matrix norm in modal_fsg_hmc.py. The manuscript omits the blanket propagator-validation claim affected by that measurement issue. The GPU driver and raw Monte Carlo production datasets are preserved on this branch so the ongoing run can be reconciled independently.

## Revision validation

All ten figure generators completed with clean layout audits. The exported PDFs were checked for 18.4 cm width and embedded fonts, and the manuscript figure captions remain below 200 words each. The revised pair-budget calculation and the CPU graviton-spectrum calculation were rerun. The 19-page manuscript compiles without undefined references or overfull boxes, and the rendered figures and manuscript were visually inspected. These checks validate the revision artifacts and selected calculations; they do not establish a continuum theory of quantum gravity.

The generated spectrum and two-loop CSVs, HMC summary JSON, and regenerated pair-budget JSON remain local pending approval to publish numerical datasets. The repository retains its archived pair-budget JSON. Run the revised pair-budget verification before rebuilding the figures, as listed in the README. The other outputs are regenerated by the included figure scripts.
