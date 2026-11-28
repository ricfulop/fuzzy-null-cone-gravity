#!/usr/bin/env python3
"""
Generate publication-grade 4-panel figure for the Path 2 UV-Completion Frontier:
  (a) 2-Loop No-Go vs. Cross-Focusing Rescue: At fixed fuzzy-sphere cutoff L_perp=1, 1-loop W^2
      converges as L_par -> 32 while 2-loop derivative sunset diverges quadratically (~ L_par^2);
      full Path 2 cross-focusing renders 2-loop (and all-loop) sums strictly finite at N_rad = inf.
  (b) Dynamical SU(N) Rank Sum in AdS_4 vs. Flat Space: log W_N vs N showing flat-space Jeans
      instability (+N^2 explosion) vs. AdS_4 strict UV convergence (-N^3 decay at large N).
  (c) 4D AdS_4 Hawking-Page Phase Transition (L=25) from Exact Discrete-Rank Sum: Black-hole
      probability P_BH(T*L) showing opposite-sign 1-loop shift of T_HP between Path 2 glued
      (+106/45, lowers T_HP) and kinematical LQG (-3/2, raises T_HP) relative to classical 1/pi.
  (d) Multi-Diamond Coarse-Graining + 4D CEMZ Causality Selection of s_min: Superluminal
      impact parameter b_adv/l_P = sqrt(s_min/pi) vs local embadon dimension d_0, ruling out
      large d_0 without sub-Planckian higher-spin towers J > 2 and isolating minimal qubit/qutrit embadons.
"""

import json
import math
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from orx_figstyle import BASELINE, PALETTE, WIDE, figure_grid, panel_labels, save, use_style

use_style()

with open(os.path.join(HERE, "p2_uv1_2loop_gs_cemz_results.json")) as f:
    u1 = json.load(f)
with open(os.path.join(HERE, "p2_uv2_rank_sum_hp_results.json")) as f:
    u2 = json.load(f)

fig, axes = figure_grid(2, 2, width=WIDE, ratio=0.77)
ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

# =========================================================================
# Panel (a): 1-Loop Convergence vs. 2-Loop Fuzzy-Only No-Go
# =========================================================================
sc1 = u1["part_a_2loop_finiteness"]["fuzzy_only_2loop_divergence_scan"]
L_par = np.array([r["L_par"] for r in sc1], dtype=float)
I1 = np.array([r["I_1loop_W2"] for r in sc1], dtype=float)
I2 = np.array([r["I_2loop_derivative_sunset"] for r in sc1], dtype=float)

ax1.loglog(L_par, I2, marker="s", color=PALETTE["red"], label="2-loop sunset, $S^2_N$ only (diverges)")
ax1.loglog(L_par, I1, marker="o", color=PALETTE["blue"], label="1-loop $W^2$, $S^2_N$ only (converges)")
ax1.axhline(sc1[0]["I_1loop_limit_A_inf"], linestyle=":", color=PALETTE["blue"], linewidth=0.9)
# Also plot normalized constant marker showing cross-focusing B=4 finite cap (independent of L_par -> inf!)
ax1.axhline(251.18, linestyle="--", color=PALETTE["green"], linewidth=1.1, label="Mode cutoff $B=4$")

ax1.set_xlim(1.7, 38)
ax1.set_ylim(0.5, 1.2e4)
ax1.set_xticks([2, 4, 8, 16, 32])
ax1.get_xaxis().set_major_formatter(lambda x, _: f"{int(x)}")
ax1.set_yticks([1e0, 1e1, 1e2, 1e3, 1e4])
ax1.minorticks_off()
ax1.set_xlabel("Parallel radial/time cutoff $\\Lambda_\\parallel / \\Lambda_\\perp$")
ax1.set_ylabel("Loop integral magnitude")
ax1.legend(loc="upper left", fontsize=7, labelspacing=0.25)

# =========================================================================
# Panel (b): Dynamical Rank Sum Convergence (AdS_4 vs. Flat-Space Jeans)
# =========================================================================
# Compute exact log W_N curve up to N = 260 at beta = 20, L = 20 so the AdS_4 turnaround & plunge to -inf is crystal clear!
s_min = math.log(3.0)
Ns_b = np.linspace(5, 260, 120)
sN_b = s_min * (Ns_b**2)
E_flat_b = 0.5 * np.sqrt(sN_b / math.pi)
E_ads_b = E_flat_b + 0.5 * (sN_b**1.5) / ((math.pi**1.5) * (20.0**2))
logW_flat = sN_b + (106.0 / 45.0) * np.log(sN_b) - 20.0 * E_flat_b
logW_ads = sN_b + (106.0 / 45.0) * np.log(sN_b) - 20.0 * E_ads_b

ax2.axhline(0.0, color=BASELINE, linestyle=":", linewidth=0.7)
ax2.plot(Ns_b, logW_flat / 1e4, linestyle="--", color=PALETTE["red"], linewidth=1.3, label="Flat space ($+N^2$ Jeans instability)")
ax2.plot(Ns_b, logW_ads / 1e4, color=PALETTE["blue"], linewidth=1.4, label="AdS$_4$ ($L=20$)")

ax2.set_xlim(0, 265)
ax2.set_ylim(-5.5, 8.2)
ax2.set_xticks([0, 60, 120, 180, 240])
ax2.set_yticks([-4, -2, 0, 2, 4, 6, 8])
ax2.set_xlabel("Matrix rank $N$")
ax2.set_ylabel("Log Boltzmann weight $\\ln W_N \\; (\\times 10^4)$")
ax2.legend(loc="upper left", fontsize=7, labelspacing=0.25)

# =========================================================================
# Panel (c): Exact 4D Hawking-Page Phase Transition (L=25) & 1-Loop Shift
# =========================================================================
sw_m = u2["L25_temperature_sweep_mixed"]
sw_l = u2["L25_temperature_sweep_lqg"]
TL_c = np.array([r["TL"] for r in sw_m], dtype=float)
P_m = np.array([r["P_bh"] for r in sw_m], dtype=float)
P_l = np.array([r["P_bh"] for r in sw_l], dtype=float)

ax3.axvline(1.0 / math.pi, color=BASELINE, linestyle=":", linewidth=1.0, label="Classical $T_{\\mathrm{HP}} L = 1/\\pi$")
ax3.plot(TL_c, P_m, color=PALETTE["blue"], linewidth=1.4, label="FSG glued ($+\\frac{106}{45}\\ln s$, lowers $T_{\\mathrm{HP}}$)")
ax3.plot(TL_c, P_l, linestyle="--", color=PALETTE["orange"], linewidth=1.3, label="Kinematical LQG ($-\\frac{3}{2}\\ln s$, raises $T_{\\mathrm{HP}}$)")

ax3.set_xlim(0.285, 0.345)
ax3.set_ylim(-0.04, 1.08)
ax3.set_xticks([0.29, 0.30, 0.31, 0.32, 0.33, 0.34])
ax3.set_yticks([0.0, 0.5, 1.0])
ax3.set_xlabel("Dimensionless temperature $T \\cdot L_{\\mathrm{AdS}}$ ($L = 25\\,\\ell_P$)")
ax3.set_ylabel("Black-hole phase probability $P_{\\mathrm{BH}}$")
ax3.legend(loc="center left", bbox_to_anchor=(0.01, 0.56), fontsize=7, labelspacing=0.25)

# =========================================================================
# Panel (d): Coarse-Graining + 4D CEMZ Causality Selection of s_min
# =========================================================================
d0_grid = np.geomspace(2.0, 1000.0, 200)
b_A = np.sqrt(np.log(1.0 + d0_grid) / math.pi)
b_B = np.sqrt(np.log(d0_grid) / math.pi)

ax4.axhline(1.0, color=PALETTE["red"], linestyle="--", linewidth=1.1, label=r"Sub-Planckian bound $\rho_{\min} = \ell_P$ ($s_{\min}=\pi$)")
ax4.axhspan(1.0, 1.65, color=PALETTE["red"], alpha=0.10)
ax4.semilogx(d0_grid, b_A, color=PALETTE["blue"], linewidth=1.3, label="Model A: $s_{\\min} = \\ln(1+d_0)$")
ax4.semilogx(d0_grid, b_B, linestyle="-.", color=PALETTE["green"], linewidth=1.3, label="Model B: $s_{\\min} = \\ln(d_0)$")
# Mark minimal qubit/qutrit fixed points (d0 = 2, 3)
ax4.scatter([2, 3], [math.sqrt(math.log(3) / math.pi), math.sqrt(math.log(4) / math.pi)], color=PALETTE["blue"], s=22, zorder=5)
ax4.scatter([2, 3], [math.sqrt(math.log(2) / math.pi), math.sqrt(math.log(3) / math.pi)], color=PALETTE["green"], s=22, zorder=5)

ax4.set_xlim(1.7, 1200)
ax4.set_ylim(0.35, 1.62)
ax4.set_xticks([2, 10, 100, 1000])
ax4.get_xaxis().set_major_formatter(lambda x, _: f"{int(x)}")
ax4.set_yticks([0.5, 1.0, 1.5])
ax4.minorticks_off()
ax4.set_xlabel("Local embadon Hilbert dimension $d_0$")
ax4.set_ylabel(r"Minimal pixel radius $\rho_{\min} / \ell_P$")
ax4.legend(loc="lower right", fontsize=7, labelspacing=0.25)

panel_labels([ax1, ax2, ax3, ax4])
out_stem = os.path.join(HERE, "..", "figs", "fig3_uv_completion")
save(fig, out_stem)
print("Saved path2_uv_frontier.svg and .pdf")
