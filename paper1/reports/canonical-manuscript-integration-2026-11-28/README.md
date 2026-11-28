# Canonical results in the two Paper I manuscript layouts

Both layouts compiled twice from the Paper I directory with pdfLaTeX. The technical layout has 26 pages and the reading layout has 43. Each layout used a separate fresh temporary output directory. The durable project artifact contains the resulting paper.pdf and paper_reading.pdf alongside a copy of this report.

## Scientific scope and review

The appendix now derives the ADM variables of the actual angular Galerkin action and gives local angular primary-plus-secondary Poisson ranks 36 and 52 at cutoffs 3 and 4 on the exact static oblate family. The result assumes zero cosmological constant, compact radial tests on a shell away from the origin, and sufficiently small nonzero deformation. After eliminating the nondegenerate shift pairs, the linearized calculation gives one and four redundant lapse tertiary profile families and eight and twelve complementary angular channels, respectively. The reduced form extends to radial profiles in H¹₀. A common physical boundary domain and the final nonlinear Dirac classification remain open.

The independent source review (msg_2397292f) checked the rank and kernel argument, radial derivatives in the brackets, the indefinite DeWitt pairing, the reduced lapse domain, and the second-preservation exclusion of the separate round rank-one witness. It reported no blocking proof or scope issue and requested that the lapse paragraph explicitly name the eliminated shift pairs. The final source includes that clarification. The reviewer performed a prose and proof-scope review. This archive records the subsequent builds. No scientific checker was run for this manuscript integration.

## Build evidence

[Build provenance](build-provenance.json) records the base commit, working-tree status of the manuscript, engine, commands, exit codes, and SHA-256 hashes of the project inputs identified by TeX’s recorder, the outputs, and the build logs. Each layout has a log and captured stdout/stderr for each pass, named by layout and pass number.

Neither final pass reported unresolved references or overfull boxes. The technical layout reported a stuck float at manuscript input line 301 and a deferred float during page clearing. Both final logs contain PDF-string and label-package warnings. The builder initially left the checked-in PDFs untouched. Root subsequently verified every recorded source, log, and PDF hash and copied these exact outputs to the current checked-in PDF paths.
