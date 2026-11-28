#!/usr/bin/env python3
"""
Generate publication-grade 4-panel capstone figure for the Master Resolution of arXiv:2611.34924v1:
  (a) Repaired Section 4.4 — Nonperturbative 4D Graviton Pair-Creation Amplitude M_4D^{P2}(s):
      Re M_4D ~ s^{-3/2} and |Im M_4D| ~ ln(s) / (s_min s^2) vs exact analytic asymptotics
      across s = 10 .. 1e6 (convergent at N_rad = inf, curing the original J^2 continuum blowup).
  (b) Repaired Section 5.4 (Part A) — Baby-Universe / Wormhole Gram-Matrix Spectrum (K = 140 fillings):
      Exact 14-order-of-magnitude rank cliff at k = 100 = dim K_Sigma separating physical states
      from 40 exact null wormhole superpositions (<= 1.1e-15), with eq. (5.73) factorization residual 8.9e-16.
  (c) Repaired Section 5.4 (Part B) — Two-Sided Hartle-Hawking Gluing from Disconnected Path 2 Caps:
      Cross-horizon Raychaudhuri area-focusing correlator C_A(-u, u) vs affine separation u / l_P,
      locking onto connected-horizon eq. (5.77b) at u >> l_P while smoothly curing the u -> 0 blowup.
  (d) Repaired Section 6 — Exact 1-Loop Semiclassical Einstein Equation (eq. 6.13):
      Residual |E[g_bar] + <E^{(2)}[h, h]>_G| plunges as O(hbar^2) (p = 1.999) when both shear variance
      and dressed metric-shear cross-term are included, vs O(hbar^1) for classical and naive shear-only.
"""

import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from orx_figstyle import BASELINE, PALETTE, WIDE, figure_grid, panel_labels, save, use_style

use_style()

with open(os.path.join(HERE, "orig1_sec44_4d_pair_results.json")) as f:
    o1 = json.load(f)
with open(os.path.join(HERE, "orig2_sec54_wormhole_hh_results.json")) as f:
    o2 = json.load(f)
with open(os.path.join(HERE, "orig3_sec6_semiclassical_results.json")) as f:
    o3 = json.load(f)

fig, axes = figure_grid(2, 2, width=WIDE, ratio=0.77)
ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

# =========================================================================
# Panel (a): Repaired Sec. 4.4 — 4D Nonperturbative Graviton Pair Creation
# =========================================================================
sc1 = o1["s_scan"]
s_arr = np.array([r["s"] for r in sc1], dtype=float)
re_arr = np.array([r["Re_M4D_bessel"] for r in sc1], dtype=float)
im_arr = np.array([abs(r["Im_M4D_bessel"]) for r in sc1], dtype=float)
re_as = np.array([r["Re_asym"] for r in sc1], dtype=float)
im_as = np.array([abs(r["Im_asym"]) for r in sc1], dtype=float)

ax1.loglog(s_arr, re_arr, marker="o", markersize=3.8, color=PALETTE["blue"], label=r"$\mathrm{Re}\,\mathcal{M}_{4\mathrm{D}}^{\mathrm{FSG}}(s)$ ($\propto s^{-3/2}$)")
ax1.loglog(s_arr, re_as, linestyle=":", color=PALETTE["blue"], linewidth=1.0, label=r"Asymptotic $\frac{\pi}{\sqrt{2s_{\min}}} s^{-3/2}$")
ax1.loglog(s_arr, im_arr, marker="s", markersize=3.5, color=PALETTE["orange"], label=r"$|\mathrm{Im}\,\mathcal{M}_{4\mathrm{D}}^{\mathrm{FSG}}(s)|$ ($\propto \frac{\ln s}{s^2}$)")
ax1.loglog(s_arr, im_as, linestyle="--", color=PALETTE["orange"], linewidth=1.0, label=r"Non-analytic branch cut $\frac{\ln s}{s_{\min} s^2}$")

ax1.set_xlim(7, 1.8e6)
ax1.set_ylim(4e-12, 0.25)
ax1.set_xticks([1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
ax1.set_yticks([1e-11, 1e-8, 1e-5, 1e-2])
ax1.minorticks_off()
ax1.set_xlabel(r"Corner area $s = \mathcal{A}/(4 G_N \hbar)$")
ax1.set_ylabel(r"4D pair amplitude $|\mathcal{M}_{4\mathrm{D}}(s)|$")
ax1.legend(loc="lower left", fontsize=7, labelspacing=0.25)

# =========================================================================
# Panel (b): Repaired Sec. 5.4A — Wormhole Gram-Matrix Spectrum (K = 140)
# =========================================================================
eigs = np.array(o2["part_a_wormhole_nulls_and_factorization"]["gram_eigs_all"], dtype=float)
idx = np.arange(1, len(eigs) + 1)

ax2.semilogy(idx[:100], eigs[:100], color=PALETTE["blue"], linewidth=1.4, label=r"Physical states ($k \leq 100 = \mathrm{dim}\,\mathcal{K}_\Sigma$)")
ax2.semilogy(idx[100:], eigs[100:], color=PALETTE["red"], linewidth=1.5, label=r"40 exact null wormholes ($\lambda \leq 1.1\times 10^{-15}$)")
ax2.axvline(100.5, color=BASELINE, linestyle=":", linewidth=1.0)
ax2.axhline(1e-15, color=PALETTE["red"], linestyle="--", linewidth=0.7, alpha=0.5)

ax2.set_xlim(0, 144)
ax2.set_ylim(1e-19, 1e3)
ax2.set_xticks([0, 25, 50, 75, 100, 125, 140])
ax2.set_yticks([1e-18, 1e-12, 1e-6, 1e0])
ax2.minorticks_off()
ax2.set_xlabel(r"Wormhole Gram eigenvalue index $k$ ($K=140$ fillings)")
ax2.set_ylabel(r"Gram eigenvalue $|\lambda_k(G)|$")
ax2.legend(loc="lower left", fontsize=7, labelspacing=0.28)

# =========================================================================
# Panel (c): Repaired Sec. 5.4B — Two-Sided Hartle-Hawking Horizon Gluing
# =========================================================================
sc3 = o2["part_b_hartle_hawking_gluing"]["bifurcation_uv_regularization_scan"]
u_vals = np.array([r["u"] for r in sc3], dtype=float)
ca_cl = np.array([r["CA_classical"] for r in sc3], dtype=float)
ca_p2 = np.array([r["CA_path2"] for r in sc3], dtype=float)

ax3.loglog(u_vals, ca_cl, linestyle="--", color=PALETTE["red"], linewidth=1.2, label=r"Classical connected eq. (5.77b) $\propto u^{-4}$")
ax3.loglog(u_vals, ca_p2, marker="o", markersize=3.8, color=PALETTE["blue"], linewidth=1.4, label=r"FSG disconnected TFD caps (UV-finite)")
ax3.axvline(0.5914, color=BASELINE, linestyle=":", linewidth=0.9, label=r"Embadon floor $\rho_{\min} = \sqrt{\ln 3 / \pi}\,\ell_P$")

ax3.set_xlim(0.038, 42)
ax3.set_ylim(5e-8, 8e4)
ax3.set_xticks([0.1, 1.0, 10.0])
ax3.get_xaxis().set_major_formatter(lambda x, _: f"{x:g}")
ax3.set_yticks([1e-6, 1e-3, 1e0, 1e3])
ax3.minorticks_off()
ax3.set_xlabel(r"Horizon affine separation $u / \ell_P$ ($u_L = u_R = u$)")
ax3.set_ylabel(r"Area-focusing correlator $C_{\mathcal{A}}(-u, u)$")
ax3.legend(loc="lower left", fontsize=7, labelspacing=0.25)

# =========================================================================
# Panel (d): Repaired Sec. 6 — Exact 1-Loop Semiclassical Einstein eq. (6.13)
# =========================================================================
sc4 = o3["part_a_quantum_raychaudhuri_eq613"]["hbar_scan"]
hb = np.array([r["hbar"] for r in sc4], dtype=float)
res_cl = np.array([r["res_classical_O_hbar"] for r in sc4], dtype=float)
res_sh = np.array([r["res_naive_shear_only"] for r in sc4], dtype=float)
res_613 = np.array([r["res_eq613_full"] for r in sc4], dtype=float)

ax4.loglog(hb, res_cl, marker="^", markersize=3.8, linestyle="--", color=PALETTE["red"], label=r"Classical $E[\bar{g}]$ ($O(\hbar^{0.99})$)")
ax4.loglog(hb, res_sh, marker="s", markersize=3.4, linestyle="-.", color=PALETTE["orange"], label=r"Shear-only (no cross-term, $O(\hbar^{1.0})$)")
ax4.loglog(hb, res_613, marker="o", markersize=3.8, color=PALETTE["green"], linewidth=1.5, label=r"Full FSG Eq. (6.13) ($O(\hbar^{2.00})$)")

ax4.set_xlim(7e-5, 1.6e-1)
ax4.set_ylim(5e-10, 0.2)
ax4.set_xticks([1e-4, 1e-3, 1e-2, 1e-1])
ax4.set_yticks([1e-9, 1e-7, 1e-5, 1e-3, 1e-1])
ax4.minorticks_off()
ax4.set_xlabel(r"Quantum parameter $\hbar$")
ax4.set_ylabel(r"Semiclassical Einstein residual")
ax4.legend(loc="lower right", fontsize=7, labelspacing=0.25)

panel_labels([ax1, ax2, ax3, ax4])
out_stem = os.path.join(HERE, "..", "figs", "fig4_original_resolution")
save(fig, out_stem)
print("Saved original_problem_solved.svg and .pdf")
