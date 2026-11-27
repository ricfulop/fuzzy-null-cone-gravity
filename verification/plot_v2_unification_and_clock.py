#!/usr/bin/env python3
"""
Generate publication-grade 4-panel Fig. 5 for v2 of the APS PRD manuscript,
directly illustrating all four reviewer-response resolutions:
  (a) Explicit Matrix Cubic Vertex V_3^{(N)}[F, Z]: parity-symmetrized O(1/N^2) (p=2.05) vs.
      one-sided O(1/N) (p=1.06) convergence to the continuum spin-2 Gaunt integral.
  (b) Interacting Projected Dynamics H_phys = P_B H_FSG P_B built from V_3^{(N)}: O(1/N^2)
      convergence of both 1->2 graviton transition probability P_{1->2} (p=2.03) and 1-loop
      self-energy Delta E_{1L} (p=2.02) with zero sector leakage.
  (c) Unified B(s) = s Q_* Eq. (7) Reconciliation: strictly finite mode set preserves Re M ~ s^{-3/2}
      and Im M ~ ln(s)/s^2 with exact O(s^{-2}) tail bound |M_inf - M_B| vs. fixed-B=4 rational decay.
  (d) Exponential-Clock & Spectral-Projector Bounds inside B(s) = s Q_*: uniform 2nd-moment,
      4th-moment, and exact Balian-Brezin Fredholm Chernoff failure probability bounds vs. area s.
"""

import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from orx_figstyle import BASELINE, PALETTE, WIDE, figure_grid, panel_labels, save, use_style

use_style()

with open(os.path.join(HERE, "v2_test1_matrix_vertex_projected_results.json")) as f:
    v1 = json.load(f)
with open(os.path.join(HERE, "v2_test2_qei_exp_clock_results.json")) as f:
    v2 = json.load(f)

fig, axes = figure_grid(2, 2, width=WIDE, ratio=0.77)
ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

# =========================================================================
# Panel (a): Explicit Matrix Cubic Vertex V_3^{(N)}[F, Z] Convergence
# =========================================================================
sc_a = v1["part_a_explicit_matrix_vertex"]["N_scan"]
Ns_a = np.array([r["N"] for r in sc_a], dtype=float)
e_one = np.array([r["triples"]["(2,2,2;m=0)"]["err_onesided_vs_cont"] for r in sc_a], dtype=float)
e_sym = np.array([r["triples"]["(2,2,2;m=0)"]["err_sym_vs_cont"] for r in sc_a], dtype=float)
e_422 = np.array([r["triples"]["(4,2,2;m=0)"]["err_sym_vs_cont"] for r in sc_a if "(4,2,2;m=0)" in r["triples"]], dtype=float)
Ns_422 = np.array([r["N"] for r in sc_a if "(4,2,2;m=0)" in r["triples"]], dtype=float)

ax1.loglog(Ns_a, e_one, marker="^", markersize=3.8, linestyle="--", color=PALETTE["orange"], label=r"One-sided $\mathrm{Tr}_{N-2}(\iota(F) Z^\dagger Z)$ ($p=1.06$)")
ax1.loglog(Ns_a, e_sym, marker="o", markersize=3.8, color=PALETTE["blue"], linewidth=1.4, label=r"Sym. $\mathcal{V}_3^{(N)}$ $(2,2,2)$ ($p=2.05$)")
ax1.loglog(Ns_422, e_422, marker="s", markersize=3.4, color=PALETTE["green"], linewidth=1.3, label=r"Sym. $\mathcal{V}_3^{(N)}$ $(4,2,2)$ ($O(1/N^2)$)")

ax1.set_xlim(5.3, 23)
ax1.set_ylim(1.2e-3, 0.7)
ax1.set_xticks([6, 8, 12, 16, 20])
ax1.get_xaxis().set_major_formatter(lambda x, _: f"{int(x)}")
ax1.set_yticks([1e-2, 1e-1])
ax1.minorticks_off()
ax1.set_xlabel(r"Fuzzy-sphere matrix rank $N = J_{\max}+1$")
ax1.set_ylabel(r"Cubic vertex error $|C^{(N)} - C^{\mathrm{cont}}|$")
ax1.legend(loc="lower left", fontsize=5.1, labelspacing=0.25)

# =========================================================================
# Panel (b): Interacting Projected H_phys = P_B H_FSG P_B Convergence
# =========================================================================
sc_b = v1["part_b_projected_dynamics_and_eq7"]["interacting_N_scan"]
Ns_b = np.array([r["N"] for r in sc_b], dtype=float)
e_dE = np.array([r["dE_1loop_err_vs_cont"] for r in sc_b], dtype=float)
e_P = np.array([r["P_1to2_err_vs_cont"] for r in sc_b], dtype=float)

ax2.loglog(Ns_b, e_dE, marker="o", markersize=3.8, color=PALETTE["blue"], linewidth=1.4, label=r"1-loop self-energy $|\Delta E_{\mathrm{1L}}^{(N)} - \Delta E_{\mathrm{1L}}^{\mathrm{cont}}|$ ($p=2.02$)")
ax2.loglog(Ns_b, e_P, marker="s", markersize=3.5, color=PALETTE["green"], linewidth=1.4, label=r"$1\to 2$ prob. $|P_{1\to 2}^{(N)} - P_{1\to 2}^{\mathrm{cont}}|$ ($p=2.03$)")

ax2.set_xlim(5.2, 37)
ax2.set_ylim(1.5e-8, 8e-4)
ax2.set_xticks([6, 8, 12, 16, 24, 32])
ax2.get_xaxis().set_major_formatter(lambda x, _: f"{int(x)}")
ax2.set_yticks([1e-7, 1e-5, 1e-3])
ax2.minorticks_off()
ax2.set_xlabel(r"Fuzzy-sphere matrix rank $N = J_{\max}+1$")
ax2.set_ylabel(r"Projected $H_{\mathrm{phys}}$ observable error")
ax2.legend(loc="lower left", fontsize=5.0, labelspacing=0.28)

# =========================================================================
# Panel (c): Unified B(s) = s Q_* Eq. (7) Reconciliation & O(s^{-2}) Tail
# =========================================================================
sc_c = v1["part_b_projected_dynamics_and_eq7"]["unified_eq7_scan"]
s_c = np.array([r["s"] for r in sc_c], dtype=float)
re_B = np.array([r["Re_M_B"] for r in sc_c], dtype=float)
im_B = np.array([abs(r["Im_M_B"]) for r in sc_c], dtype=float)
diff_B = np.array([r["tail_diff_abs"] for r in sc_c], dtype=float)
fixed_B4 = 1.483012 / (s_c**2)

ax3.loglog(s_c, re_B, marker="o", markersize=3.8, color=PALETTE["blue"], linewidth=1.4, label=r"$\mathrm{Re}\,\mathcal{M}_{B(s)}$ ($B=sQ_*$, $\propto s^{-3/2}$)")
ax3.loglog(s_c, im_B, marker="s", markersize=3.5, color=PALETTE["orange"], linewidth=1.3, label=r"$|\mathrm{Im}\,\mathcal{M}_{B(s)}|$ ($\propto \ln s / s^2$)")
ax3.loglog(s_c, diff_B, marker="^", markersize=3.4, linestyle=":", color=PALETTE["green"], linewidth=1.3, label=r"Tail $|\mathcal{M}_\infty - \mathcal{M}_{B(s)}|$ ($0.33/s^2$)")
ax3.loglog(s_c, fixed_B4, linestyle="--", color=PALETTE["red"], linewidth=1.1, label=r"Fixed $B=4$ ($\mathrm{Re} = 1.483/s^2$)")

ax3.set_xlim(75, 4.2e4)
ax3.set_ylim(1.5e-10, 6e-3)
ax3.set_xticks([1e2, 1e3, 1e4])
ax3.set_yticks([1e-9, 1e-6, 1e-3])
ax3.minorticks_off()
ax3.set_xlabel(r"Corner area $s = \mathcal{A}/(4 G_N \hbar)$")
ax3.set_ylabel(r"Projected 4D amplitude magnitude")
ax3.legend(loc="lower left", fontsize=5.0, labelspacing=0.22)

# =========================================================================
# Panel (d): Exponential-Clock & Spectral-Projector Failure Bounds
# =========================================================================
sc_d = v2["unified_Bs_exponential_clock_scan"]
s_d = np.array([r["s"] for r in sc_d], dtype=float)
f2 = np.array([r["spectral_proj_fail_2nd_moment"] for r in sc_d], dtype=float)
f4 = np.array([r["spectral_proj_fail_4th_moment"] for r in sc_d], dtype=float)
fc = np.array([10.0 ** r["chernoff_log10_fail_prob"] for r in sc_d], dtype=float)

ax4.loglog(s_d, f2, marker="^", markersize=3.8, linestyle="--", color=PALETTE["orange"], label=r"2nd-moment bound ($O(1/s)$, uniform $N_{\mathrm{rad}}\to\infty$)")
ax4.loglog(s_d, fc, marker="s", markersize=3.5, linestyle="-.", color=PALETTE["green"], linewidth=1.3, label=r"Balian–Brezin Chernoff ($B=sQ_*$, $\lambda_{\mathrm{crit}}\gg 2$)")
ax4.loglog(s_d, f4, marker="o", markersize=3.8, color=PALETTE["blue"], linewidth=1.4, label=r"4th-moment bound ($O(1/s^2)$, uniform $N_{\mathrm{rad}}\to\infty$)")

ax4.set_xlim(16, 300)
ax4.set_ylim(2e-10, 2e-3)
ax4.set_xticks([20, 40, 80, 160, 240])
ax4.get_xaxis().set_major_formatter(lambda x, _: f"{int(x)}")
ax4.set_yticks([1e-9, 1e-7, 1e-5, 1e-3])
ax4.minorticks_off()
ax4.set_xlabel(r"Horizon area $s$ ($N_{\mathrm{rad}} = B(s) = 0.5\,s$)")
ax4.set_ylabel(r"Clock failure bound $\|(1-P_{\mathcal{D}})|\Omega_s\rangle\|^2$")
ax4.legend(loc="lower left", fontsize=5.0, labelspacing=0.24)

panel_labels([ax1, ax2, ax3, ax4])
out_stem = os.path.join(HERE, "fig5_v2_unification_and_clock")
save(fig, out_stem)
print("Saved fig5_v2_unification_and_clock.svg and .pdf")
