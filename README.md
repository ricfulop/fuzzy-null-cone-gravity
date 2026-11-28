# Quantum Gravity

**Author:** Ric Fulop, Center for Bits and Atoms, Massachusetts Institute of Technology under Prof. Neil Gershenfeld supervision. 

This repository contains two companion research manuscripts. Each paper has its own source, PDFs, figures, reports, and verification scripts.

| Paper | Working folder | Read the draft | Build and verification instructions |
|---|---|---|---|
| I — A Finite-State Model of Quantum Gravity in a Causal Diamond | [`paper1/`](paper1/) | [Reading PDF](paper1/paper_reading.pdf) · [Technical PDF](paper1/paper.pdf) | [Paper I README](paper1/README.md) |
| II — Toward Interacting Einstein Gravity from Finite Causal-Diamond States | [`paper2/`](paper2/) | [Working draft v0.2 PDF](paper2/paper2_einstein_gravity.pdf) | [Paper II README](paper2/README.md) |

These are research proposals with stated limits. Consult each manuscript and its verification reports for the claims established by its calculations.

## Working on either paper

Run each paper's build and verification commands from its own folder. For example:

```bash
cd paper1
pdflatex paper.tex
pdflatex paper.tex
pdflatex paper_reading.tex
pdflatex paper_reading.tex
```

Use `cd paper2` for the companion manuscript's commands. Open the appropriate folder in your IDE or select that folder's main TeX file in Overleaf.

Before continuing an existing IDE branch, fetch and merge `origin/main` so it includes the folder organization. Paper I paths formerly at the repository root now start with `paper1/`; Paper II remains in `paper2/`. Keep changes and generated artifacts within the paper being edited. Use `paper1/…` and `paper2/…` branch names and identify the paper in pull-request titles; both target `main`.

## License

MIT License — see [LICENSE](LICENSE).
