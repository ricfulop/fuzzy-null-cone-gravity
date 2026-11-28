#!/usr/bin/env python3
"""
Generate publication-grade 2-panel diagnostic figure from real Test 1 and Test 2 outputs:
  (a) Test 1A: Divergence of Z_P[X] in opposite directions under product vs cycle exhaustions
      on H = C^2 (x) l^2(N) (refuting Section 5.2 eq. 5.38 without closed trace-class convergence).
  (b) Test 2B: Exact Williamson zero-point focusing lower bound E_{0,min}(N, J) vs mode cutoff
      N = J for a macroscopic diamond (b = 10^5 l_P, l(D) = 10 u_cl), crossing the collapse
      threshold ln(10) above which P_D = 0 identically.
"""

import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from orx_figstyle import BASELINE, MUTED, PALETTE, WIDE, figure_grid, panel_labels, save, use_style

use_style()

with open(os.path.join(HERE, "test1_results.json")) as f:
    t1 = json.load(f)
with open(os.path.join(HERE, "test2_results.json")) as f:
    t2 = json.load(f)

fig, (ax1, ax2) = figure_grid(1, 2, width=WIDE, ratio=0.42)

# --- Panel (a): Test 1A counterexample ---
prod_M = np.array([r["M"] for r in t1["test_1a_qubit_clock"]["product_rows"]])
prod_Z = np.array([r["Z_P_X"] for r in t1["test_1a_qubit_clock"]["product_rows"]])
cyc_K = np.array([r["K"] for r in t1["test_1a_qubit_clock"]["cycle_rows"]])
cyc_Z = np.array([r["Z_P_X_max"] for r in t1["test_1a_qubit_clock"]["cycle_rows"]])

ax1.axhline(0.0, color=BASELINE, linestyle="--", linewidth=0.8, zorder=1)
ax1.plot(prod_M, prod_Z, marker="o", color=PALETTE["blue"], label="Product exhaustion $P_M$")
ax1.plot(cyc_K, cyc_Z, marker="s", color=PALETTE["red"], label="Cycle exhaustion $P_{\\mathrm{cyc}}(K, r=200)$")
ax1.set_yscale("symlog", linthresh=1.0, linscale=0.6)
ax1.set_xscale("log")
ax1.set_ylim(-5e4, 6e2)
ax1.set_yticks([-1e4, -1e2, -1, 0, 1, 1e2])
ax1.set_yticklabels(["$-10^4$", "$-10^2$", "$-1$", "$0$", "$1$", "$10^2$"])
ax1.minorticks_off()
ax1.set_xlabel("Included blocks ($M$) or cycles ($K$)")
ax1.set_ylabel("Partial amplitude $Z_P[X]$ (dimensionless)")
ax1.annotate(
    "$Z_{P_M}[X] \\to +\\infty$",
    xy=(prod_M[-1], prod_Z[-1]),
    xytext=(55, 2.2),
    textcoords="data",
    fontsize=7,
    color=PALETTE["blue"],
    backgroundcolor="white",
)
ax1.annotate(
    "$Z_P[X] \\to -\\infty$",
    xy=(cyc_K[-1], cyc_Z[-1]),
    xytext=(3.0, -8000),
    textcoords="data",
    fontsize=7,
    color=PALETTE["red"],
    backgroundcolor="white",
)
ax1.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0))

# --- Panel (b): Test 2B continuum clock zero-point collapse ---
sweep = t2["continuum_sweep"]
Ns = np.array([r["N"] for r in sweep], dtype=float)
E0s = np.array([r["E0_total"] for r in sweep], dtype=float)
log_cap = float(np.log(t2["l_over_ucl"]))

ax2.axhspan(log_cap, 1e6, color=MUTED, alpha=0.35, zorder=0)
ax2.axhline(log_cap, color=BASELINE, linestyle="--", linewidth=0.8, zorder=1)
ax2.loglog(Ns, E0s, marker="o", color=PALETTE["blue"], label="Williamson floor $E_{0,\\min}(N,J)$")

# Mark Planck crossover where E0_min = log(10): c_rad * N^4 / s = log(10)
c_rad = t2["asymptotic_c_rad"]
s_macro = t2["s_macro"]
N_cross = (log_cap * s_macro / c_rad) ** 0.25
ax2.axvline(N_cross, color=PALETTE["orange"], linestyle=":", linewidth=1.0, zorder=2)

ax2.set_ylim(1e-6, 1e6)
ax2.set_xlim(8, 7000)
ax2.set_yticks([1e-6, 1e-3, 1e0, 1e3, 1e6])
ax2.set_xticks([1e1, 1e2, 1e3])
ax2.minorticks_off()
ax2.set_xlabel("Radial & angular cutoff $N = J$ (modes)")
ax2.set_ylabel("Focusing floor $E_{0,\\min}$ (nats)")
ax2.text(
    11,
    log_cap * 12.0,
    "Collapse region:\n$P_{\\mathscr{D}} = 0$ identically",
    fontsize=6.8,
    color="#333333",
)
ax2.text(
    N_cross * 1.2,
    8e-5,
    f"Planck cutoff\n$N_\\star \\approx {int(round(N_cross))}$",
    fontsize=6.8,
    color=PALETTE["orange"],
    backgroundcolor="white",
)
ax2.legend(loc="lower right")

panel_labels([ax1, ax2])
out_stem = os.path.join(HERE, "verification_summary")
save(fig, out_stem)
print("Saved verification_summary.svg and .pdf")
