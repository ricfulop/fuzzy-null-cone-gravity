#!/usr/bin/env python3
"""
Generate final publication-grade 4-panel Fig. 6 summarizing the converged 28-GPU Modal NVIDIA B200
dual IR + UV nonperturbative FAHMC simulation (all 28 runs R-hat <= 1.005, Schwinger-Dyson <= 0.16%):
  (a) Dressed spin-2 graviton propagator G_meas(j)/G_free(j) across j=2..47 at N=48 (mu^2=0.05, g_P=3.382)
      and strong-coupling fixed-g_N=2.0 at N=32, showing UV finiteness across 2t_j^2 = 4..2254;
  (b) Quantum area matrix rho^2 = exp((g/sqrt(N)) Phi) eigenvalue density rho(lambda) at mu^2=0.05, g_P=3.382
      across N=8..48, showing strictly positive caustic-free horizon concentration onto lambda = 1;
  (c) Fixed-g_N UV strong-coupling scan (N=8, 16, 24, 32): N-independent nonperturbative horizon-anisotropy
      suppression Q_2 / Q_{2,free} down to 0.522 (-47.8%) at g_N=2.0 with bounded intrinsic C_V^{int};
  (d) DeWitt-Fujikawa invariant measure vs. flat-measure graviton-induced 180x horizon collapse (g=2.5).
"""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from orx_figstyle import BASELINE, PALETTE, WIDE, figure_grid, panel_labels, save, use_style

use_style()

fin = json.load(open(os.path.join(HERE, "modal_hmc_final_results.json")))
t1 = json.load(open(os.path.join(HERE, "modal_hmc_tier1_results.json")))
t2 = json.load(open(os.path.join(HERE, "modal_hmc_tier2_results.json")))
fu_ana_path = os.path.join(HERE, "modal_hmc_followup_analysis.json")
fu_b2_path = os.path.join(HERE, "modal_hmc_followup_b2_both_results.json")
fu_ana = json.load(open(fu_ana_path)) if os.path.exists(fu_ana_path) else None
fu_b2 = json.load(open(fu_b2_path)) if os.path.exists(fu_b2_path) else None

fin_runs = fin["runs"]
old_all = t1["runs"] + t2["runs"]
old_main = [r for r in old_all if r["params"]["mu2"] == 0.5 and "flatmeas" not in r["tag"]]
flat_runs = sorted([r for r in old_all if "flatmeas" in r["tag"]], key=lambda r: r["params"]["N"])

fig, axes = figure_grid(2, 2, width=WIDE, ratio=0.77)
ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

# =========================================================================
# Panel (a): Full Spin-2 Graviton Propagator at N=48 (IR, mu2=0.05) & N=32 (UV, g_N=2.0)
# =========================================================================
r48_ir = next(r for r in fin_runs if r["tag"] == "final_IR_mu005_N48")
r32_uv = next(r for r in fin_runs if r["tag"] == "final_UV_gN2.00_N32")
r8_ir = next(r for r in fin_runs if r["tag"] == "final_IR_mu005_N8")

for r_cur, col, mk, lab in [
    (r48_ir, PALETTE["blue"], "o", r"IR $N=48,\ \mu^2=0.05,\ g_P=3.38$"),
    (r32_uv, PALETTE["red"], "^", r"UV $N=32,$ fixed $g_N=2.0$"),
    (r8_ir, PALETTE["orange"], "s", r"IR $N=8,\ \mu^2=0.05,\ g_P=3.38$"),
]:
    js = np.array([m["j"] for m in r_cur["spin2_modes"]], dtype=float)
    z0 = np.array([m["Z_ratio_0"] for m in r_cur["spin2_modes"]], dtype=float)
    zse = np.array([m["G_meas_0_se"] / m["G_free_0"] for m in r_cur["spin2_modes"]], dtype=float)
    ax1.errorbar(js, z0, yerr=zse, fmt=f"{mk}-", markersize=2.5, linewidth=1.1, elinewidth=0.65, color=col, label=lab)

ax1.axhline(1.0, color=BASELINE, linestyle=":", linewidth=0.9)
ax1.set_xlim(1, 49)
ax1.set_ylim(0.935, 1.055)
ax1.set_xticks([2, 10, 20, 30, 40, 47])
ax1.set_yticks([0.94, 0.97, 1.00, 1.03])
ax1.set_xlabel(r"Spin-2 graviton angular momentum $j$ ($2t_j^2 = 4 \to 2254$)")
ax1.set_ylabel(r"Dressed propagator $G_{\mathrm{meas}}(j)/G_{\mathrm{free}}(j)$")
ax1.legend(loc="lower right", fontsize=5.3, labelspacing=0.22)

# =========================================================================
# Panel (b): Quantum Area Eigenvalue Spectrum rho(lambda) at mu^2=0.05, g_P=3.382
# =========================================================================
specs_b = [
    (8, PALETTE["red"], "--", r"$N=8$ ($s=70$, $\lambda \in [0.43, 2.66]$)"),
    (16, PALETTE["orange"], "-.", r"$N=16$ ($s=281$, $\lambda \in [0.53, 1.98]$)"),
    (32, PALETTE["green"], "-", r"$N=32$ ($s=1125$, $\lambda \in [0.64, 1.52]$)"),
    (48, PALETTE["blue"], "-", r"$N=48$ ($s=2531$, $\lambda \in [0.73, 1.38]$)"),
]
for N_v, col, ls, lab in specs_b:
    if N_v == 48 and fu_b2 is not None:
        n48_ws = [x for x in fu_b2["runs"] if "N48" in x["tag"]]
        bc = np.array(n48_ws[0]["rho2_eig_histogram"]["bin_centers"], dtype=float)
        pdf = np.mean([np.array(x["rho2_eig_histogram"]["pdf"], dtype=float) for x in n48_ws], axis=0)
    else:
        r = next(x for x in fin_runs if x["tag"] == f"final_IR_mu005_N{N_v}")
        bc = np.array(r["rho2_eig_histogram"]["bin_centers"], dtype=float)
        pdf = np.array(r["rho2_eig_histogram"]["pdf"], dtype=float)
    ax2.plot(bc, pdf, linestyle=ls, color=col, linewidth=1.3, label=lab)

ax2.axvline(1.0, color=BASELINE, linestyle=":", linewidth=0.9)
ax2.set_xlim(0.0, 2.75)
ax2.set_ylim(0.0, 3.35)
ax2.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
ax2.set_yticks([0.0, 1.0, 2.0, 3.0])
ax2.set_xlabel(r"Quantum area eigenvalue $\lambda(\widehat{\rho}^2 / r_0^2)$")
ax2.set_ylabel(r"Spectral density $\rho(\lambda)$ ($\mu^2=0.05$)")
ax2.legend(loc="upper right", fontsize=5.2, labelspacing=0.22)

# =========================================================================
# Panel (c): Fixed-g_N UV Strong-Coupling Horizon-Anisotropy Suppression Q2/Q2_free
#            with MEASURED zero-coupling controls (mu^2=0.25, lam4=0.1) & 12x8 check
# =========================================================================
gN_grid = [0.0, 0.5, 1.0, 1.5, 2.0]
specs_c = [
    (8, PALETTE["red"], "o", "--", r"$N=8$ ($8\times 4$)"),
    (16, PALETTE["orange"], "s", "-.", r"$N=16$ ($8\times 4$)"),
    (24, PALETTE["green"], "^", ":", r"$N=24$ ($8\times 4$)"),
    (32, PALETTE["blue"], "D", "-", r"$N=32$ ($12\times 4$)"),
]
for N_v, col, mk, ls, lab in specs_c:
    ys, yses = [], []
    for gN_v in gN_grid:
        if gN_v == 0.0 and fu_ana is not None:
            row = next(r for r in fu_ana["primary_4rank_table"] if r["N"] == N_v)
            ys.append(row["Q2_g0_mean"] / row["Q2_free"])
            yses.append(1.96 * row["Q2_g0_se"] / row["Q2_free"])
        elif gN_v == 0.0:
            ys.append(1.0)
            yses.append(0.0)
        else:
            r = next(x for x in fin_runs if x["tag"] == f"final_UV_gN{gN_v:.2f}_N{N_v}")
            ys.append(r["observables"]["Q2_anisotropy"]["mean"] / r["observables"]["Q2_free"])
            yses.append(1.96 * r["observables"]["Q2_anisotropy"]["se"] / r["observables"]["Q2_free"])
    ax3.errorbar(gN_grid, ys, yerr=yses, fmt=f"{mk}{ls}", color=col, markersize=3.5, linewidth=1.2, capsize=1.8, label=lab)

if fu_ana is not None:
    lat = fu_ana["lattice_12x8_test"]
    q2f_lat = lat["lat_12x8_g0"]["params"]["Q2_free"]
    y_lat = [
        lat["lat_12x8_g0"]["diagnostics"]["Q2"]["mean"] / q2f_lat,
        lat["lat_12x8_g2"]["diagnostics"]["Q2"]["mean"] / q2f_lat,
    ]
    ye_lat = [
        1.96 * max(lat["lat_12x8_g0"]["diagnostics"]["Q2"]["se_chains"], lat["lat_12x8_g0"]["diagnostics"]["Q2"]["se_ess"]) / q2f_lat,
        1.96 * max(lat["lat_12x8_g2"]["diagnostics"]["Q2"]["se_chains"], lat["lat_12x8_g2"]["diagnostics"]["Q2"]["se_ess"]) / q2f_lat,
    ]
    ax3.errorbar([0.04, 2.04], y_lat, yerr=ye_lat, fmt="*", color=PALETTE["purple"], markersize=5.4, markerfacecolor="white", markeredgewidth=0.9, capsize=1.8, label=r"$N=24$ ($12\times 8$)")

ax3.axhline(1.0, color=BASELINE, linestyle=":", linewidth=0.8)
ax3.set_xlim(-0.10, 2.15)
ax3.set_ylim(0.45, 1.06)
ax3.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax3.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax3.set_xlabel(r"Fixed microscopic coupling $g_N$ ($\mu^2 = 0.25,\ \lambda_4 = 0.1$)")
ax3.set_ylabel(r"Roughness ratio $Q_2(g_N) / Q_{2,\mathrm{free}}$")
ax3.legend(loc="lower left", ncol=2, fontsize=5.1, labelspacing=0.20, columnspacing=0.6)

# =========================================================================
# Panel (d): DeWitt-Fujikawa Measure vs. Flat-Measure Horizon Collapse (g=2.5)
# =========================================================================
Ns_d = np.array([r["params"]["N"] for r in flat_runs], dtype=float)
rho_flat = np.array([r["observables"]["rho2_eig_mean"]["mean"] for r in flat_runs], dtype=float)
cv_flat = np.array([r["observables"]["C_V_per_dof"]["mean"] for r in flat_runs], dtype=float)
rho_dw = np.array(
    [next(x for x in old_main if x["params"]["N"] == int(n) and abs(x["params"]["g"] - 2.5) < 1e-3)["observables"]["rho2_eig_mean"]["mean"] for n in Ns_d],
    dtype=float,
)
cv_int_dw = np.array(
    [next(x for x in fin_runs if x["tag"] == f"final_IR_mu050_N{int(n)}")["observables"]["C_V_intrinsic_per_dof"]["mean"] for n in Ns_d],
    dtype=float,
)

ax4.semilogy(Ns_d, rho_dw, marker="o", color=PALETTE["blue"], linewidth=1.4, label=r"DeWitt measure $\langle\widehat{\rho}^2\rangle \approx 1.0$")
ax4.semilogy(Ns_d, rho_flat, marker="s", linestyle="--", color=PALETTE["red"], linewidth=1.3, label=r"Flat measure $\langle\widehat{\rho}^2\rangle \approx 0.006$ ($180\times$)")
ax4.semilogy(Ns_d, cv_flat, marker="^", linestyle="-.", color=PALETTE["orange"], linewidth=1.2, label=r"Flat measure $C_V/N_{\mathrm{dof}} \propto N^2$")
ax4.semilogy(Ns_d, cv_int_dw, marker="D", linestyle=":", color=PALETTE["green"], linewidth=1.2, label=r"DeWitt $C_V^{\mathrm{int}}/N_{\mathrm{dof}} \approx 0.46$")

ax4.set_xlim(6, 34)
ax4.set_ylim(1.8e-3, 1.4e3)
ax4.set_xticks([8, 16, 24, 32])
ax4.set_yticks([1e-2, 1e0, 1e2])
ax4.minorticks_off()
ax4.set_xlabel(r"Fuzzy-sphere matrix rank $N = J_{\max} + 1$")
ax4.set_ylabel(r"Horizon area $\langle\widehat{\rho}^2\rangle$ and $C_V/N_{\mathrm{dof}}$")
ax4.legend(loc="center right", bbox_to_anchor=(0.99, 0.25), fontsize=5.1, labelspacing=0.20, frameon=True, facecolor="white", framealpha=0.92, edgecolor="none")

panel_labels([ax1, ax2, ax3, ax4])
out_stem = os.path.join(HERE, "fig6_modal_b200_hmc")
save(fig, out_stem)
print("Saved final fig6_modal_b200_hmc.svg and .pdf")
