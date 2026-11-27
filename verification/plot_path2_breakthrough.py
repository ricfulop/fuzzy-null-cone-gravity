#!/usr/bin/env python3
"""
Generate publication-grade 4-panel summary figure for the Path 2 Breakthrough:
  (a) Transverse-coupled QEI clock survival on fuzzy sphere S^2_J (C^infty_c bump, N_rad -> inf):
      Ultralocal vs 4D transverse-coupled variance Var_0(Xi_avg^{NO}) & QEI floor |Xi_QEI| vs area s.
  (b) Flat-space spin-2 graviton entanglement entropy on fuzzy sphere S^2_J:
      Extracted c_log vs fuzzy cutoff alpha = J_max / R_eff for Sharp Step vs Smooth Casimir
      exp(-j(j+1)/(alpha R_eff)^2) hitting Benedetti-Casini's -61/45 to 0.02% at alpha = 2.0.
  (c) Curved 4D Schwarzschild discrete-harmonic gluing (N_rad = inf, j <= J_max(r)):
      Subleading entropy S(s) - s showing how Lichnerowicz Weyl tidal coupling (+3 ln s) flips
      kinematical embadon corners (c_log = 0, -3/2) straight onto Sen's +106/45 and +77/90!
  (d) Unified 4D logarithmic-correction scorecard across Flat Sphere, Schwarzschild, and dS S^4.
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

with open(os.path.join(HERE, "p2_test1_transverse_qei_results.json")) as f:
    t1 = json.load(f)
with open(os.path.join(HERE, "p2_test2_fuzzy_graviton_ee_results.json")) as f:
    t2 = json.load(f)
with open(os.path.join(HERE, "p2_test5_schwarzschild_s4_results.json")) as f:
    t5 = json.load(f)

fig, axes = figure_grid(2, 2, width=WIDE, ratio=0.76)
ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

# =========================================================================
# Panel (a): Transverse-Coupled QEI Clock Floor & Variance vs Area s
# =========================================================================
rows1 = t1["fuzzy_sphere_coupled_survival"]
s1 = np.array([r["s"] for r in rows1], dtype=float)
fl_u = np.array([abs(r["Xi_floor_ultralocal"]) for r in rows1], dtype=float)
fl_c = np.array([abs(r["Xi_floor_coupled"]) for r in rows1], dtype=float)
vr_u = np.array([r["Var_Xi_avg_ultralocal"] for r in rows1], dtype=float)
vr_c = np.array([r["Var_Xi_avg_coupled"] for r in rows1], dtype=float)

ax1.loglog(s1, fl_u, marker="o", linestyle="--", color=PALETTE["orange"], label="Ultralocal QEI floor $|\\Xi_{\\mathrm{QEI}}|$")
ax1.loglog(s1, fl_c, marker="s", color=PALETTE["red"], label="4D coupled QEI floor $|\\Xi_{\\mathrm{QEI}}|$")
ax1.loglog(s1, vr_u, marker="^", linestyle="--", color=PALETTE["green"], label="Ultralocal $\\mathrm{Var}_0(\\Xi_{\\mathrm{avg}}^{\\mathrm{NO}})$")
ax1.loglog(s1, vr_c, marker="D", color=PALETTE["blue"], label="4D coupled $\\mathrm{Var}_0(\\Xi_{\\mathrm{avg}}^{\\mathrm{NO}})$")

ax1.set_xlim(35, 7000)
ax1.set_ylim(2e-7, 1.2e-1)
ax1.set_yticks([1e-6, 1e-4, 1e-2])
ax1.set_xticks([50, 200, 1000, 5000])
ax1.get_xaxis().set_major_formatter(lambda x, _: f"{int(x)}")
ax1.minorticks_off()
ax1.set_xlabel("Dimensionless horizon area $s$")
ax1.set_ylabel("Clock operator scale (dimensionless)")
ax1.legend(loc="center right", bbox_to_anchor=(0.99, 0.48), fontsize=5.3, labelspacing=0.25)

# =========================================================================
# Panel (b): Flat-Space Spin-2 Graviton EE on Fuzzy Sphere S^2_J
# =========================================================================
sharp = t2["path2_fuzzy_sharp_interp"]
smooth = [r for r in t2["path2_fuzzy_smooth_casimir"] if r["alpha"] <= 4.0]

a_sh = np.array([r["alpha"] for r in sharp], dtype=float)
c_sh = np.array([r["c_log"] for r in sharp], dtype=float)
a_sm = np.array([r["alpha"] for r in smooth], dtype=float)
c_sm_gr = np.array([r["c_log_graviton"] for r in smooth], dtype=float)
c_sm_mx = np.array([r["c_log_maxwell"] for r in smooth], dtype=float)

ax2.axhline(-61.0 / 45.0, color=PALETTE["blue"], linestyle=":", linewidth=1.0, label="Exact graviton $-\\frac{61}{45}$")
ax2.axhline(-16.0 / 45.0, color=PALETTE["green"], linestyle=":", linewidth=1.0, label="Exact Maxwell $-\\frac{16}{45}$")
ax2.plot(a_sh, c_sh, marker="o", linestyle="--", color=PALETTE["orange"], label="Sharp step $j \\leq \\alpha R_{\\mathrm{eff}}$")
ax2.plot(a_sm, c_sm_gr, marker="s", color=PALETTE["blue"], label="Smooth Casimir graviton")
ax2.plot(a_sm, c_sm_mx, marker="^", color=PALETTE["green"], label="Smooth Casimir Maxwell")

ax2.set_xlim(1.5, 12.6)
ax2.set_ylim(-1.75, 1.15)
ax2.set_xticks([2, 4, 6, 8, 10, 12])
ax2.set_yticks([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0])
ax2.set_xlabel("Fuzzy-sphere cutoff ratio $\\alpha = J_{\\max}/R_{\\mathrm{eff}}$")
ax2.set_ylabel("Extracted log coeff. $c_{\\log}$ (in $\\ln R$)")
ax2.legend(loc="upper right", fontsize=5.8)

# =========================================================================
# Panel (c): Curved Schwarzschild Discrete-Harmonic Gluing vs Sen
# =========================================================================
dyn = t5["part_a_schwarzschild"]["dynamical_Jmax_r_rows"]
s5 = np.array([r["s"] for r in dyn], dtype=float)
dS_bulk = np.array([r["Delta_S_1loop"] for r in dyn], dtype=float)
# Shift by common O(1) constant at s_0 = s5[0] to compare pure log slopes cleanly!
dS_mixed_rel = dS_bulk - dS_bulk[0]
dS_singlet_rel = (dS_bulk - 1.5 * np.log(s5)) - (dS_bulk[0] - 1.5 * np.log(s5[0]))
kin_singlet_rel = -1.5 * np.log(s5 / s5[0])

s_line = np.geomspace(s5[0], s5[-1], 100)
sen_mixed_line = (106.0 / 45.0) * np.log(s_line / s5[0])
sen_singlet_line = (77.0 / 90.0) * np.log(s_line / s5[0])

ax3.axhline(0.0, color=BASELINE, linestyle=":", linewidth=0.7, zorder=0)
ax3.semilogx(s_line, sen_mixed_line, linestyle="--", color=PALETTE["orange"], label="Sen exact mixed ($+\\frac{106}{45}\\ln s$)")
ax3.semilogx(s5, dS_mixed_rel, marker="o", linestyle="none", color=PALETTE["red"], label="Path 2 discrete-harmonic mixed")
ax3.semilogx(s_line, sen_singlet_line, linestyle="-.", color=PALETTE["purple"], label="Sen exact singlet ($+\\frac{77}{90}\\ln s$)")
ax3.semilogx(s5, dS_singlet_rel, marker="s", linestyle="none", color=PALETTE["blue"], label="Path 2 glued $SU(2)$ singlet")
ax3.semilogx(s5, kin_singlet_rel, marker="^", linestyle=":", color=BASELINE, label="Kinematical singlet only ($-\\frac{3}{2}\\ln s$)")

ax3.set_xlim(220, 4.5e5)
ax3.set_ylim(-13, 21)
ax3.set_yticks([-10, -5, 0, 5, 10, 15, 20])
ax3.minorticks_off()
ax3.set_xlabel("Schwarzschild horizon area $s = \\pi (r_H/\\ell_P)^2$")
ax3.set_ylabel("Relative subleading $\\Delta S(s) - \\Delta S(s_0)$")
ax3.legend(loc="upper left", fontsize=5.6)

# =========================================================================
# Panel (d): Unified 4D Spin-2 Log Correction Scorecard Across Geometries
# =========================================================================
labels = [
    "dS $S^4$ total\n($-\\frac{571}{45}$)",
    "dS $S^4$ bulk\n($-\\frac{331}{45}$)",
    "Flat sphere\n($-\\frac{61}{45}$)",
    "Schw. singlet\n($+\\frac{77}{90}$)",
    "Schw. mixed\n($+\\frac{106}{45}$)",
]
exact_vals = np.array([-571.0 / 45.0, -331.0 / 45.0, -61.0 / 45.0, 77.0 / 90.0, 106.0 / 45.0])
# Path 2 measured/exact values:
# S^4: exact rational match (-571/45, -331/45); Flat sphere: -1.35523 (alpha=2.0);
# Schw: asymptotic r_H=160..320 slope (+0.85334 singlet, +2.35334 mixed) / fit (+0.84557, +2.34557)
p2_vals = np.array([-571.0 / 45.0, -331.0 / 45.0, -1.35523, 0.84557, 2.34557])

y_pos = np.arange(len(labels))
h = 0.34
ax4.axvline(0.0, color=BASELINE, linewidth=0.8)
ax4.barh(y_pos + h / 2, exact_vals, height=h, color=BASELINE, alpha=0.75, label="Exact 4D GR target")
ax4.barh(y_pos - h / 2, p2_vals, height=h, color=PALETTE["blue"], label="Path 2 fuzzy-sphere / glued")

ax4.set_yticks(y_pos)
ax4.set_yticklabels(labels, fontsize=6.5)
ax4.set_xlim(-15.5, 4.8)
ax4.set_xticks([-15, -10, -5, 0, 5])
ax4.set_xlabel("Logarithmic coefficient $c_{\\log}$ (exact convention)")
ax4.legend(loc="upper left", fontsize=6.0)

panel_labels([ax1, ax2, ax3, ax4])
out_stem = os.path.join(HERE, "path2_breakthrough_summary")
save(fig, out_stem)
print("Saved path2_breakthrough_summary.svg and .pdf")
