# Quantum Gravity repository

- Paper I lives in `paper1/`; Paper II lives in `paper2/`. Each folder owns its manuscript, figures, calculation reports, and verification scripts.
- Read the relevant paper's README before changing it. Run its build and verification commands from that paper's folder.
- Fetch and incorporate `origin/main` before continuing an older branch. Paper I files previously at the repository root have moved to `paper1/`.
- Keep a paper-specific change within that paper's folder. Cross-paper edits must be intentional and described in the pull request.
- Preserve archived calculations and provenance. A historical reference to a pinned commit describes that version, even if the current paths differ.
- Do not equate the papers' distinct nonlinear regulators or promote a limited check into a proof of full gravitational consistency.
- Use a paper-specific branch and PR title, targeting `main`. Root files are for repository-wide navigation and configuration.
- GPU runs require a separately authorized compute budget; local build and verification commands should not launch them.
