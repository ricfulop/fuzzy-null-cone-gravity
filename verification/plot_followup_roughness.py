#!/usr/bin/env python3
"""
Generate publication-grade 2-panel follow-up roughness figure (fig_followup_roughness.pdf):
  (a) Horizon roughness ratio Q_2(g_N) / Q_{2,free} across g_N in {0, 0.5, 1.0, 1.5, 2.0}
      for N in {8, 16, 24, 32} (mu^2 = 0.25, lam4 = 0.1), plotting MEASURED zero-coupling
      controls (with error bars) rather than a pinned analytic point, plus the N=24 12x8 pair.
  (b) Primary spin-2 smoothing percentage S(%) = 100 * [1 - Q_2(g_N=2) / Q_2(g_N=0)]
      with explicit measured g_N=0 denominator and 95% MC intervals propagating uncertainty
      from both independent ensembles, compared against the analytic Q_{2,free} denominator,
      the N=24 12x8 larger-lattice matched pair, and the extended 16-chain N=32 cross-check.
"""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from orx_figstyle import BASELINE, PALETTE, WIDE, figure_grid, panel_labels, save, use_style

use_style()

ana = json.load(open(os.path.join(HERE, "modal_hmc_followup_analysis.json")))
grid = ana["crosscheck_gN_grid"]
rows = ana["primary_4rank_table"]
lat = ana["lattice_12x8_test"]
ext32_sc = ana["extended_n32_smoothing_crosscheck"]

fig, axes = figure_grid(1, 2, width=WIDE, ratio=0.44)
ax1, ax2 = axes[0], axes[1]

# =========================================================================
# Panel (a): Measured Q_2(g_N) / Q_{2,free} including measured g_N=0 controls
# =========================================================================
specs = [
    (8, PALETTE["red"], "o", "--", r"$N=8$ ($8\times 4$)"),
    (16, PALETTE["orange"], "s", "-.", r"$N=16$ ($8\times 4$)"),
    (24, PALETTE["green"], "^", ":", r"$N=24$ ($8\times 4$)"),
    (32, PALETTE["blue"], "D", "-", r"$N=32$ ($12\times 4$)"),
]
for N_v, col, mk, ls, lab in specs:
    sub = sorted([r for r in grid if r["N"] == N_v], key=lambda r: r["g_N"])
    xs = np.array([r["g_N"] for r in sub], dtype=float)
    ys = np.array([r["ratio_vs_analytic_free"] for r in sub], dtype=float)
    # 95% MC error bars (+/- 1.96 SE)
    yerr = 1.96 * np.array([r["ratio_vs_analytic_free_se"] for r in sub], dtype=float)
    ax1.errorbar(
        xs, ys, yerr=yerr, fmt=f"{mk}{ls}", color=col,
        markersize=3.6, linewidth=1.2, elinewidth=0.9, capsize=2.0, label=lab,
    )

# Overlay N=24 12x8 matched pair at g_N = 0 and g_N = 2
q2f_12x8 = lat["lat_12x8_g0"]["params"]["Q2_free"]
y_lat = np.array([
    lat["lat_12x8_g0"]["diagnostics"]["Q2"]["mean"] / q2f_12x8,
    lat["lat_12x8_g2"]["diagnostics"]["Q2"]["mean"] / q2f_12x8,
])
yerr_lat = 1.96 * np.array([
    lat["lat_12x8_g0"]["diagnostics"]["Q2"]["se_chains"] / q2f_12x8,
    lat["lat_12x8_g2"]["diagnostics"]["Q2"]["se_chains"] / q2f_12x8,
])
ax1.errorbar(
    [0.04, 2.04], y_lat, yerr=yerr_lat, fmt="*", color=PALETTE["purple"],
    markersize=5.8, markerfacecolor="white", markeredgewidth=1.0, elinewidth=0.9, capsize=2.0,
    label=r"$N=24$ ($12\times 8$, measured)",
)

ax1.axhline(1.0, color=BASELINE, linestyle=":", linewidth=0.85, label=r"Analytic free ($\lambda_4=0$)")
ax1.set_xlim(-0.12, 2.18)
ax1.set_ylim(0.46, 1.055)
ax1.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax1.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax1.set_xlabel(r"Fixed microscopic coupling $g_N$ ($\mu^2 = 0.25,\ \lambda_4 = 0.1$)")
ax1.set_ylabel(r"Measured roughness ratio $Q_2(g_N) / Q_{2,\mathrm{free}}$")
ax1.legend(loc="lower left", fontsize=5.1, labelspacing=0.22, frameon=True, facecolor="white", framealpha=0.92, edgecolor="none")

# =========================================================================
# Panel (b): Primary Smoothing % = 100 * [1 - Q_2(g_N=2) / Q_2(g_N=0)]
#            with explicit denominator comparison and 95% MC intervals
# =========================================================================
x_pos = np.arange(4, dtype=float)
Ns = [r["N"] for r in rows]

s_meas = np.array([r["smoothing_vs_measured_g0"]["smoothing_pct"] for r in rows])
lo_meas = s_meas - np.array([r["smoothing_vs_measured_g0"]["ci95_mc_lo"] for r in rows])
hi_meas = np.array([r["smoothing_vs_measured_g0"]["ci95_mc_hi"] for r in rows]) - s_meas

s_free = np.array([r["smoothing_vs_analytic_free"]["smoothing_pct"] for r in rows])
lo_free = s_free - np.array([r["smoothing_vs_analytic_free"]["ci95_mc_lo"] for r in rows])
hi_free = np.array([r["smoothing_vs_analytic_free"]["ci95_mc_hi"] for r in rows]) - s_free

ax2.errorbar(
    x_pos - 0.10, s_meas, yerr=[lo_meas, hi_meas], fmt="o", color=PALETTE["blue"],
    markersize=4.5, elinewidth=1.1, capsize=2.6,
    label=r"Measured control denom. $\overline{Q}_2(g_N=0)$ (95% MC)",
)
ax2.errorbar(
    x_pos + 0.10, s_free, yerr=[lo_free, hi_free], fmt="s", color=PALETTE["orange"],
    markersize=4.0, markerfacecolor="white", markeredgewidth=1.0, elinewidth=0.95, capsize=2.2,
    label=r"Analytic free denom. $Q_{2,\mathrm{free}}$ ($\lambda_4=0$)",
)

# N=24 12x8 larger-lattice matched pair
s_12x8 = lat["smoothing_12x8_vs_measured_g0"]["smoothing_pct"]
lo_12x8 = s_12x8 - lat["smoothing_12x8_vs_measured_g0"]["ci95_mc_lo"]
hi_12x8 = lat["smoothing_12x8_vs_measured_g0"]["ci95_mc_hi"] - s_12x8
ax2.errorbar(
    [2.0 - 0.28], [s_12x8], yerr=[[lo_12x8], [hi_12x8]], fmt="*", color=PALETTE["purple"],
    markersize=6.2, elinewidth=1.1, capsize=2.6,
    label=r"$N=24,\ 12\times 8$ pair ($\overline{Q}_2(g_N=0)$ denom.)",
)

# Extended 16-chain N=32 cross-check
s_ext32 = ext32_sc["smoothing_pct"]
lo_ext32 = s_ext32 - ext32_sc["ci95_mc_lo"]
hi_ext32 = ext32_sc["ci95_mc_hi"] - s_ext32
ax2.errorbar(
    [3.0 - 0.28], [s_ext32], yerr=[[lo_ext32], [hi_ext32]], fmt="D", color=PALETTE["green"],
    markersize=4.0, elinewidth=1.1, capsize=2.4,
    label=r"$N=32$ ext. 16-chain ($\overline{Q}_2(g_N=0)$ denom.)",
)

ax2.set_xlim(-0.52, 3.45)
ax2.set_ylim(47.42, 48.24)
ax2.set_xticks(x_pos)
ax2.set_xticklabels([r"$N=8$" + "\n" + r"($8\times 4$)", r"$N=16$" + "\n" + r"($8\times 4$)", r"$N=24$" + "\n" + r"($8\times 4$)", r"$N=32$" + "\n" + r"($12\times 4$)"])
ax2.set_yticks([47.5, 47.7, 47.9, 48.1])
ax2.set_xlabel(r"Fuzzy-sphere matrix rank $N$ and lattice $N_\tau \times N_r$")
ax2.set_ylabel(r"Roughness smoothing $100[1 - \overline{Q}_2(2)/\mathrm{denom.}]$ (%)")
ax2.legend(loc="upper right", fontsize=5.1, labelspacing=0.20, frameon=True, facecolor="white", framealpha=0.94, edgecolor="none")

panel_labels([ax1, ax2])
out_stem = os.path.join(HERE, "fig_followup_roughness")
save(fig, out_stem)

# Also copy PDF directly into worktree figs/
wt_figs = "/Users/ricfulop/.local/share/openresearch/worktrees/6439af1c-583a-47ac-bbe7-dba0c8996921/chat_c6a471e6-1d45-4fa3-83ca-b1f73d9c3af3/figs"
if os.path.isdir(wt_figs):
    import shutil
    shutil.copy2(out_stem + ".pdf", os.path.join(wt_figs, "fig_followup_roughness.pdf"))
print("Saved fig_followup_roughness.pdf and .svg")
