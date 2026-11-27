#!/usr/bin/env python3
"""
Generate publication-grade 2-panel diagnostic figure for the Path 2 Prototype:
  (a) Radial UV rescue: Un-normal-ordered clock floor E0_un (exploding as N^2) vs
      tapered normal-ordered Flanagan QEI floor |E_min^{NO}| and vacuum variance Var_0
      for a genuine C^infty_c((0,1)) bump function g(u) = exp(-1/(4u(1-u)) + 1) as N_rad -> 160.
  (b) Entropy & logarithmic correction diagnostic: Exact log Omega(s) - s for Path 2
      unconstrained embadons (c_log = 0) and SU(2)-singlet embadons (c_log = -3/2)
      vs Sen's 4D pure-gravity 1-loop prediction (+106/45 log s for mixed, +77/90 log s for singlet).
"""

import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from orx_figstyle import BASELINE, PALETTE, WIDE, figure_grid, panel_labels, save, use_style

use_style()

with open(os.path.join(HERE, "path2_prototype_results.json")) as f:
    p2 = json.load(f)
with open(os.path.join(HERE, "path2_cinfty_bump_results.json")) as f:
    cb = json.load(f)

fig, (ax1, ax2) = figure_grid(1, 2, width=WIDE, ratio=0.43)

# --- Panel (a): Radial UV rescue on C^infty_c bump ---
rows = cb["convergence_table"]
Ns = np.array([r["N_rad"] for r in rows], dtype=float)
E0_un = np.array([r["E0_un"] for r in rows], dtype=float)
abs_Emin_NO = np.array([abs(r["E_min_NO"]) for r in rows], dtype=float)
Var0_NO = np.array([r["Var0_NO"] for r in rows], dtype=float)

ax1.loglog(Ns, E0_un, marker="o", color=PALETTE["red"], label="Un-normal-ordered $E_{0,\\mathrm{un}}$ (diverges)")
ax1.loglog(Ns, abs_Emin_NO, marker="s", color=PALETTE["blue"], label="Tapered QEI floor $|E_{\\min}^{\\mathrm{NO}}|$")
ax1.loglog(Ns, Var0_NO, marker="^", linestyle="--", color=PALETTE["green"], label="Vacuum variance $\\mathrm{Var}_0(H_g^{\\mathrm{NO}})$")

ax1.set_xlim(8, 200)
ax1.set_ylim(1e-3, 1.5e4)
ax1.set_yticks([1e-3, 1e-1, 1e1, 1e3])
ax1.set_xticks([10, 40, 160])
ax1.get_xaxis().set_major_formatter(lambda x, _: f"{int(x)}")
ax1.minorticks_off()
ax1.set_xlabel("Radial mode cutoff $N_{\\mathrm{rad}}$ (modes)")
ax1.set_ylabel("Operator magnitude (dimensionless)")
ax1.legend(loc="center left", bbox_to_anchor=(0.02, 0.39), fontsize=6.6)

# --- Panel (b): Subleading entropy S(s) - s diagnostic vs Sen ---
s_B = np.array([r["s"] for r in p2["module1_corner"]["modelB_su2_singlet"]], dtype=float)
res_B = np.array([r["residual_minus_s"] for r in p2["module1_corner"]["modelB_su2_singlet"]], dtype=float)
s_A = np.array([r["s"] for r in p2["module1_corner"]["modelA_unconstrained"]], dtype=float)
res_A = np.array([r["diff_from_s"] for r in p2["module1_corner"]["modelA_unconstrained"]], dtype=float)

s_ref = np.geomspace(12.0, 5500.0, 200)
# Align Sen curves at s_0 = 20 for visual slope comparison (additive O(1) constant is scheme-dependent)
s0 = 20.0
sen_mixed = (106.0 / 45.0) * np.log(s_ref / s0)
sen_singlet = (77.0 / 90.0) * np.log(s_ref / s0)

ax2.axhline(0.0, color=BASELINE, linestyle=":", linewidth=0.7, zorder=0)
ax2.semilogx(s_ref, sen_mixed, color=PALETTE["orange"], linestyle="--", label="Sen 4D pure GR mixed ($+\\frac{106}{45}\\ln s$)")
ax2.semilogx(s_ref, sen_singlet, color=PALETTE["purple"], linestyle="-.", label="Sen 4D pure GR singlet ($+\\frac{77}{90}\\ln s$)")
ax2.semilogx(s_A, res_A, marker="o", color=PALETTE["green"], label="Path 2 embadons, all $J$ ($c_{\\log}=0$)")
ax2.semilogx(s_B, res_B, marker="s", color=PALETTE["blue"], label="Path 2 $SU(2)$ singlet ($c_{\\log}=-\\frac{3}{2}$)")

ax2.set_xlim(10, 6000)
ax2.set_ylim(-14, 16)
ax2.set_yticks([-10, -5, 0, 5, 10, 15])
ax2.minorticks_off()
ax2.set_xlabel("Dimensionless horizon area $s = \\mathscr{A}/(4 G_N \\hbar)$")
ax2.set_ylabel("Subleading entropy $S(s) - s$ (nats)")
ax2.legend(loc="upper left", fontsize=6.4)

panel_labels([ax1, ax2])
out_stem = os.path.join(HERE, "path2_prototype_results")
save(fig, out_stem)
print("Saved path2_prototype_results.svg and .pdf")
