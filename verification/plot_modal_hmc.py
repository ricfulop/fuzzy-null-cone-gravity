#!/usr/bin/env python3
"""
Generate publication-grade 4-panel figure summarizing the 44-ensemble Modal NVIDIA B200
nonperturbative FAHMC simulation of 4D Fuzzy-Sphere Null-Cone Matrix Gravity (I_FSG^{(N)}):
  (a) Nonperturbative spin-2 graviton dispersion E_eff^2(j) / (2 t_j^2) across j=2..47 at N=48, g_P=3.3821
      and large-N suppression of quadrupole dressing 1 - Z(j=2) across N=8..48;
  (b) Quantum area matrix rho^2 = exp((g/sqrt(N)) Phi) eigenvalue PDF rho(lambda) at g_P=3.3821 across N=8..48,
      showing strictly positive caustic-free horizon concentrating as (1 - lambda_min) ~ N^{-0.44};
  (c) Nonperturbative monopole susceptibility chi_phi0(g) (>12-sigma stiffening) and specific heat C_V/dof
      matching exact analytic scaling 0.5 + (g^2/3) chi_phi0 across g=0..4.5;
  (d) DeWitt-Fujikawa invariant measure vs. flat-measure graviton-induced horizon collapse (180x collapse
      in <rho^2> and N^2 blowup of C_V/dof without the DeWitt measure across N=8..32 at g=2.5).
"""

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from orx_figstyle import WIDE, PALETTE, BASELINE, figure_grid, panel_labels, save, use_style

use_style()

t1 = json.load(open(os.path.join(HERE, "modal_hmc_tier1_results.json")))
t2 = json.load(open(os.path.join(HERE, "modal_hmc_tier2_results.json")))
all_runs = t1["runs"] + t2["runs"]
main_runs = [r for r in all_runs if r["params"]["mu2"] == 0.5 and "flatmeas" not in r["tag"]]
flat_runs = sorted([r for r in all_runs if "flatmeas" in r["tag"]], key=lambda r: r["params"]["N"])

fig, axes = figure_grid(2, 2, width=WIDE, ratio=0.77)
ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

# =========================================================================
# Panel (a): Full 46-Harmonic Spin-2 Dispersion at N=48, g_P=3.3821 & Quadrupole Scaling
# =========================================================================
r48_p = next(r for r in main_runs if r["params"]["N"] == 48 and abs(r["params"]["g"] - 3.3821) < 1e-3)
r48_0 = next(r for r in main_runs if r["params"]["N"] == 48 and abs(r["params"]["g"] - 0.0) < 1e-3)
r8_p = next(r for r in main_runs if r["params"]["N"] == 8 and abs(r["params"]["g"] - 3.3821) < 1e-3)

js_48 = np.array([m["j"] for m in r48_p["spin2_modes"]], dtype=float)
z0_48p = np.array([m["Z_ratio_0"] for m in r48_p["spin2_modes"]], dtype=float)
z0_48p_se = np.array([m["G_meas_0_se"] / m["G_free_0"] for m in r48_p["spin2_modes"]], dtype=float)
z0_480 = np.array([m["Z_ratio_0"] for m in r48_0["spin2_modes"]], dtype=float)
js_8 = np.array([m["j"] for m in r8_p["spin2_modes"]], dtype=float)
z0_8p = np.array([m["Z_ratio_0"] for m in r8_p["spin2_modes"]], dtype=float)

ax1.axhline(1.0, color=BASELINE, linestyle=":", linewidth=0.9)
ax1.plot(js_48, z0_480, color=BASELINE, linestyle="--", linewidth=1.0, label=r"$N=48,\ g=0$ (free benchmark)")
ax1.errorbar(js_48, z0_48p, yerr=z0_48p_se, fmt="o-", markersize=2.6, linewidth=1.1, elinewidth=0.7, color=PALETTE["blue"], label=r"$N=48,\ g_P=3.38$ (1.33M dof)")
ax1.plot(js_8, z0_8p, marker="s", markersize=3.5, linewidth=1.2, color=PALETTE["orange"], label=r"$N=8,\ g_P=3.38$ (microscopic)")

ax1.set_xlim(1, 49)
ax1.set_ylim(0.942, 1.022)
ax1.set_xticks([2, 10, 20, 30, 40, 47])
ax1.set_yticks([0.95, 0.97, 0.99, 1.00, 1.02])
ax1.set_xlabel(r"Spin-2 graviton angular momentum $j$ ($2t_j^2 = 4 \to 2254$)")
ax1.set_ylabel(r"Dressed propagator ratio $G_{\mathrm{meas}}(j)/G_{\mathrm{free}}(j)$")
ax1.legend(loc="lower right", fontsize=5.6, labelspacing=0.25)

# =========================================================================
# Panel (b): Quantum Area Eigenvalue Spectrum rho(lambda) at g_P = 3.3821
# =========================================================================
specs_b = [
    (8, PALETTE["red"], "--", r"$N=8$ ($s=70$, $\lambda \in [0.41, 2.53]$)"),
    (16, PALETTE["orange"], "-.", r"$N=16$ ($s=281$, $\lambda \in [0.53, 1.97]$)"),
    (32, PALETTE["green"], "-", r"$N=32$ ($s=1125$, $\lambda \in [0.66, 1.54]$)"),
    (48, PALETTE["blue"], "-", r"$N=48$ ($s=2531$, $\lambda \in [0.73, 1.39]$)"),
]
for N_v, col, ls, lab in specs_b:
    r = next(x for x in main_runs if x["params"]["N"] == N_v and abs(x["params"]["g"] - 3.3821) < 1e-3)
    bc = np.array(r["rho2_eig_histogram"]["bin_centers"], dtype=float)
    pdf = np.array(r["rho2_eig_histogram"]["pdf"], dtype=float)
    ax2.plot(bc, pdf, linestyle=ls, color=col, linewidth=1.3, label=lab)

ax2.axvline(1.0, color=BASELINE, linestyle=":", linewidth=0.9)
ax2.axvline(0.0, color=PALETTE["red"], linestyle=":", linewidth=0.8)
ax2.set_xlim(0.0, 2.75)
ax2.set_ylim(0.0, 3.35)
ax2.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
ax2.set_yticks([0.0, 1.0, 2.0, 3.0])
ax2.set_xlabel(r"Quantum area matrix eigenvalue $\lambda(\widehat{\rho}^2 / r_0^2)$")
ax2.set_ylabel(r"Spectral density $\rho(\lambda)$ at $g_P = 3.382$")
ax2.legend(loc="upper right", fontsize=5.3, labelspacing=0.24)

# =========================================================================
# Panel (c): Monopole Stiffening chi_phi0(g) & Specific Heat vs Analytic Scaling
# =========================================================================
g_vals = [0.0, 0.75, 1.5, 2.5, 3.3821, 4.5]
chi_w, chi_se = [], []
cv_w, cv_se = [], []
for g_v in g_vals:
    sub = [r for r in main_runs if abs(r["params"]["g"] - g_v) < 1e-3 and r["params"]["N"] <= 24]
    c_m = np.array([r["observables"]["chi_phi0"]["mean"] for r in sub])
    c_s = np.array([r["observables"]["chi_phi0"]["se"] for r in sub])
    w = 1.0 / np.maximum(c_s**2, 1e-12)
    chi_w.append(float(np.sum(c_m * w) / np.sum(w)))
    chi_se.append(float(1.0 / math.sqrt(np.sum(w))))

    v_m = np.array([r["observables"]["C_V_per_dof"]["mean"] for r in sub])
    v_s = np.array([r["observables"]["C_V_per_dof"]["se"] for r in sub])
    wv = 1.0 / np.maximum(v_s**2, 1e-12)
    cv_w.append(float(np.sum(v_m * wv) / np.sum(wv)))
    cv_se.append(float(1.0 / math.sqrt(np.sum(wv))))

g_arr = np.array(g_vals)
chi_arr = np.array(chi_w)
cv_arr = np.array(cv_w)
cv_pred = 0.5 + (g_arr**2 / 3.0) * chi_arr

ax3_r = ax3.twinx()
l1 = ax3.errorbar(g_arr, chi_arr, yerr=chi_se, fmt="o-", color=PALETTE["blue"], linewidth=1.35, markersize=4.2, label=r"Monopole $\chi_{\phi_0}(g)$ ($12\sigma$ shift)")
l2 = ax3_r.errorbar(g_arr, cv_arr, yerr=cv_se, fmt="s-", color=PALETTE["orange"], linewidth=1.3, markersize=3.8, label=r"Specific heat $C_V / N_{\mathrm{dof}}$")
l3 = ax3_r.plot(g_arr, cv_pred, linestyle="--", color=PALETTE["green"], linewidth=1.1, label=r"Analytic $\frac{1}{2} + \frac{g^2}{3}\chi_{\phi_0}$")

ax3.axvline(3.3821, color=BASELINE, linestyle=":", linewidth=0.9)
ax3.set_xlim(-0.2, 4.75)
ax3.set_ylim(0.95, 2.38)
ax3_r.set_ylim(0.0, 11.5)
ax3.set_xticks([0.0, 1.0, 2.0, 3.3821, 4.5])
ax3.set_xticklabels(["0", "1", "2", "$g_P$", "4.5"])
ax3.set_xlabel(r"Gravitational coupling $g$ ($g_P = \sqrt{4\pi/\ln 3}$)")
ax3.set_ylabel(r"Monopole susceptibility $\chi_{\phi_0}$", color=PALETTE["blue"])
ax3_r.set_ylabel(r"Specific heat $C_V / N_{\mathrm{dof}}$", color=PALETTE["orange"])
ax3_r.spines["top"].set_visible(False)
ax3_r.spines["right"].set_visible(True)
ax3_r.tick_params(axis="y", colors=PALETTE["orange"])
ax3.tick_params(axis="y", colors=PALETTE["blue"])
handles = [l1, l2, l3[0]]
labels_c = [h.get_label() for h in handles]
ax3.legend(handles, labels_c, loc="upper left", fontsize=5.1, labelspacing=0.22)

# =========================================================================
# Panel (d): DeWitt-Fujikawa Measure vs. Flat-Measure Horizon Collapse (g=2.5)
# =========================================================================
Ns_d = np.array([r["params"]["N"] for r in flat_runs], dtype=float)
rho_flat = np.array([r["observables"]["rho2_eig_mean"]["mean"] for r in flat_runs], dtype=float)
cv_flat = np.array([r["observables"]["C_V_per_dof"]["mean"] for r in flat_runs], dtype=float)
rho_dw = np.array(
    [next(x for x in main_runs if x["params"]["N"] == int(n) and abs(x["params"]["g"] - 2.5) < 1e-3)["observables"]["rho2_eig_mean"]["mean"] for n in Ns_d],
    dtype=float,
)
cv_dw = np.array(
    [next(x for x in main_runs if x["params"]["N"] == int(n) and abs(x["params"]["g"] - 2.5) < 1e-3)["observables"]["C_V_per_dof"]["mean"] for n in Ns_d],
    dtype=float,
)

ax4.semilogy(Ns_d, rho_dw, marker="o", color=PALETTE["blue"], linewidth=1.4, label=r"DeWitt measure $\langle\widehat{\rho}^2\rangle \approx 1.0$ (geometric)")
ax4.semilogy(Ns_d, rho_flat, marker="s", linestyle="--", color=PALETTE["red"], linewidth=1.3, label=r"Flat measure $\langle\widehat{\rho}^2\rangle \approx 0.006$ ($180\times$ collapse)")
ax4.semilogy(Ns_d, cv_flat, marker="^", linestyle="-.", color=PALETTE["orange"], linewidth=1.2, label=r"Flat measure $C_V/N_{\mathrm{dof}} \propto N^2$")
ax4.semilogy(Ns_d, cv_dw, marker="D", linestyle=":", color=PALETTE["green"], linewidth=1.2, label=r"DeWitt measure $C_V/N_{\mathrm{dof}} \approx 4.2$ (bounded)")

ax4.set_xlim(6, 34)
ax4.set_ylim(1.8e-3, 1.4e3)
ax4.set_xticks([8, 16, 24, 32])
ax4.set_yticks([1e-2, 1e0, 1e2])
ax4.minorticks_off()
ax4.set_xlabel(r"Fuzzy-sphere matrix rank $N = J_{\max} + 1$ ($g = 2.5$)")
ax4.set_ylabel(r"Horizon area $\langle\widehat{\rho}^2\rangle$ and $C_V/N_{\mathrm{dof}}$")
ax4.legend(loc="center right", bbox_to_anchor=(0.99, 0.30), fontsize=5.1, labelspacing=0.20, frameon=True, facecolor="white", framealpha=0.92, edgecolor="none")

panel_labels([ax1, ax2, ax3, ax4])
out_stem = os.path.join(HERE, "fig6_modal_b200_hmc")
save(fig, out_stem)
print("Saved fig6_modal_b200_hmc.svg and .pdf")
