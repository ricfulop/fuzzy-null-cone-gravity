# Latest draft reconciliation and template assessment

Review date: 27 November 2026.

Inputs: Chandrasekaran, *An algebraic quantization of gravity from the null cone*, arXiv:2611.34924v1 (83 pages, including appendices); the available conversation history and retrieved earlier decisions; Open Research commit `d288b95712c735c2ef1c3df293d1c87abef76849`; and review branch commit `9d6bbbf0b00d378a73a8ad2ceb25da854478db69`.

The repository's newer main branch and the earlier review branch diverged from `589a99b9a5fa5695ebaf20cd3082b76f94a02346`. The new work did not delete the review branch. This revision brings their relevant changes together, preserving the new simulation inputs and data.

This is a manuscript and evidence review, not a complete independent derivation of every equation, a rerun of all GPU simulations, or a formal verification of quantum gravity. Earlier conversation retrieval supplied summaries rather than a guaranteed verbatim export of every prior turn.

## 1. Template recommendation

Keep REVTeX 4.2 for the complete technical manuscript. Add a single-column 11-point reading copy built from the same source, with no duplicate scientific text. This is the most useful immediate change for a long derivation-heavy paper.

A Science submission should be a separate editorial condensation: a short central argument, a small set of main figures, and detailed derivations and diagnostics in supplementary materials. The Science-family template available through Overleaf uses the standard 12-point `article` class, Times-like fonts, letter paper, one-inch margins, and expanded line spacing. It does not reproduce the journal's published multi-column pages. Changing the document class does not establish scientific suitability or compliance with an article type's current limits.

Primary template sources checked:
- https://journals.aps.org/revtex
- https://www.overleaf.com/latex/templates/science-and-science-family-journals-template/ptkzkvxkznbh (Siegert, Smith, and Arseneau; identifies Science's preparation page as its source)

Recommended main narrative: the physical question and causal diamond; explicit assumptions and finite-state construction; what the latest simulations actually show; tests needed to connect the model to classical and quantum gravity. A potential condensed figure sequence is geometry, finite-state budget, matched-control roughness, and sampling/benchmark validation. This is an editorial recommendation, not a journal-imposed figure limit.

## 2. What the new computation resolves

| Previous concern | New evidence | Current assessment |
|---|---|---|
| Gaussian reference omitted the scalar quartic interaction | Measured g=0 controls keep lambda4=0.1 and mu2=0.25 | The roughly 48% effect survives matched controls; quartic-only shifts are at most 0.17% |
| Larger-lattice robustness untested | Rank 24: 8x4 to 12x8, both directions increased | Suppression changes by -0.060 ± 0.047 percentage points; one useful volume check, not continuum convergence |
| Poor local-area mixing at ranks 32 and 48 | New local Fourier-wave collective moves; 16 chains, 2000 retained draws each | Rank-normalized split-R-hat 1.0098 and 1.0021; bulk ESS 2576 and 4020; these two follow-ups meet the stated targets |
| No saved histories for stronger diagnostics | Histories in the follow-up; rank-normalized split-R-hat and bulk/tail ESS | Applies to the new runs. Archived coupled summaries do not acquire new diagnostics retroactively |

The primary matched-control suppression estimates are 48.016%, 47.758%, 47.681%, and 47.772% for ranks 8, 16, 24, and 32. Rank 32 uses a 12x4 lattice; the other primary ranks use 8x4. The larger 12x8 rank-24 pair gives 47.621%. The extended rank-32 coupled run gives a 47.789% cross-check using the new control.

The **rank-24 larger-lattice local-area observable still fails the targets**: rank-normalized split-R-hat 1.3268, bulk ESS 20, tail ESS 42. Roughness on that same ensemble has R-hat 1.0041 and bulk ESS 896. The revised figure shows both. An all-observables convergence claim would be incorrect.

The archived ordinary unsplit R-hat values 1.416 and 1.181 are not the same estimator as the follow-up rank-normalized split values. Compare Batch 1 and Batch 2 using the same estimator (1.1082 to 1.0098, and 1.0325 to 1.0021).

Intervals in the roughness figure come from the archived propagated-error analysis and describe Monte Carlo uncertainty conditional on the estimated means and errors. They do not include regulator, lattice, or modeling uncertainty. At g=0 the plotted ratio of a mean to itself is exactly one; assigning independent numerator and denominator errors there would be incorrect.

## 3. What reading the original changes

The original contains a substantial classical construction and an explicitly conditional axiomatic correspondence. The strongest trace concern is specific: the Cauchy–Schwarz estimate for open amplitudes in Eqs. (4.29)–(4.30) does not supply convergence of the closed traces invoked in Sec. 5.2. With a uniform positive lower bound on the trace operator, a finite ordinary trace already bounds the Hilbert-space dimension. An independent justification of the closed contraction is therefore essential.

The original clock discussion needs more care than the blanket assertion that its Hilbert space is empty. Appendix C.2 allows a limiting form domain whose closure is M (Eq. C.13). Its fixed-resolution estimates and the explicitly heuristic noncommuting estimate in Eq. (5.16), footnote 18, leave a genuine continuum-domain question. Our quadratic oscillator witnesses do not, without an operator-identification argument, settle the full nonlinear constrained clock. Appendix D.2 explicitly subtracts the local vacuum singularity in Eq. (D.21), so saying the original never discusses subtraction is also inaccurate.

Likewise, JT plus a scalar is not automatically a counterexample to every hypothesis of the four-dimensional proof: its clock/matter relation differs. The qubit construction demonstrates a missing mathematical implication about traces, not an entire gravitational countermodel.

The revised manuscript replaces the blanket 'Continuum Obstructions' theorem with a scoped assessment. It retains the constructive motivation and credit to Chandrasekaran while separating justified gaps from model-dependent illustrations.

## 4. Scientific limits that the new data do not remove

1. **Equivalence of constructions.** The sampled positive Euclidean Phi/Z action, the full constrained null-cone gravity action, and the small projected Hamiltonian examples have not been shown to define one equivalent theory. Matching a cubic vertex is a useful check, not such a proof.
2. **Classical and continuum limits.** Known Schwarzschild, Regge–Wheeler, heat-kernel, or entropy benchmarks do not establish recovery of full Einstein dynamics, constraints and symmetries for general states, or a predictive continuum quantum limit.
3. **Model assumptions.** The finite occupation budget and corner spectrum are imposed ingredients. Finite-state counting under those assumptions is meaningful but does not derive those assumptions from gravity.
4. **Masslessness and caustics.** The finite-rank mass-fit intercepts are not a controlled massless extrapolation. exp(g Phi) is positive for Hermitian Phi by construction. Positive sampled eigenvalues do not establish the absence of physical spacetime caustics. The reported lower/upper area values are means of extreme eigenvalues at one site, not bounds on the full spectrum.
5. **Formal verification scope.** Lean certifies the encoded statements. For example, `appendix_E2_valid_conditionally` applies an implication assumed inside `SoundRules`, and the clock incompatibility lemma assumes the growth bound it uses. Exact polynomial identities and `(1-P)PHP=0` do not certify an infinite-dimensional gravitational construction. The paper now says this explicitly.
6. **Further derivations.** The gravitational application of the scalar QEI, all-mode nonlinear clock bounds, detailed black-hole action/entropy claims, and full gluing consistency require deeper mathematical review. The retained calculations should be assessed as scoped examples and benchmarks until those derivations are complete.

## 5. Corrections retained from the earlier review

- A pair intermediate state must satisfy **2 omega <= B**. The B=4 sector has the vacuum plus 30 one-particle states and no allowed pair; the projected pair-creation amplitude there is zero.
- The corrected pair-tail calculation uses the effective single-particle ceiling B/2 and retains an O(s^-2) bound. The earlier 0.33 coefficient is not carried into the corrected plot.
- The archived two-loop numbers are a mode-cutoff diagnostic, not a calculation applying the total many-body projector to every internal state of a full gravitational diagram.
- Page's printed graviton luminosity table gives 3.8e-6; the claimed 0.017% comparison used extra stored reference digits not printed in that table. The numerical coefficient 3.83666e-6 remains a benchmark, with appropriate precision language.
- The largest original chain has 331,392 coupled variables. 1,325,568 is the total across four independent replicas.
- The geometry, 31-state example, two-loop budget, and graviton-spectrum figures are restored. The classical caustic illustration ends at zero rather than continuing through a negative radius.
- The one-loop coefficient figure distinguishes signed residuals, a multi-radius fit, and a final-doubling estimate. These are numerical resolution effects, distinct from sampling drift.

The historical review report is retained with a notice that its pending-controls and pending-area-convergence statements are superseded where documented above.

## 6. Figure and layout changes

The two simulation figures now have separate jobs. The roughness figure shows the matched-control coupling response and the strongest-coupling size/volume checks. The diagnostic figure compares roughness and area R-hat, displays area bulk/tail ESS, labels the one-site eigenvalue extrema correctly, and identifies the archived propagator data. The poor larger-lattice area result remains visible.

The manuscript defines sampling drift, rank-normalized split-R-hat, ESS, and numerical resolution dependence. Bars and intervals are identified in captions; finite sampling and model uncertainty are kept distinct.

Both PDFs are built from `paper.tex`; `paper_reading.tex` selects the reading layout. The preserved GPU driver and raw simulation results are not changed by this review. No new GPU jobs are launched. Earlier unpublished generated numerical CSV/JSON artifacts are not included; their regeneration code remains available.
