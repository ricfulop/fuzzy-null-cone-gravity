#!/usr/bin/env python3
"""
Generate publication-grade 4-panel figure for the Path 2 Frontier (Unified SU(N) Fuzzy-Sphere Null-Cone Gravity):
  (a) Nonlinear SU(N) matrix backreaction convergence: one-sided O(1/N) vs parity-symmetrized exact O(1/N^2)
      for slow-Kerr O(a^2) quadrupole (p=2.00) and nonlinear GW spin-2 shear quadrupole/hexadecapole (p=2.03).
  (b) Second variation delta^2 I_FSG^{(N)} around Schwarzschild (j=2): Regge-Wheeler (-6M/r^3 Weyl tidal),
      Zerilli Darboux partner, and minimal scalar without Weyl tidal coupling (+2M/r^3).
  (c) Exponential restoration of finite-rapidity SO(3,1) Lorentz boosts on spin-2 graviton wavepackets
      across N = 6..44 for rapidities eta = 0.5, 1.0, 1.5 (plunging 12-15 orders of magnitude to 1e-15).
  (d) Nonlinear Raychaudhuri caustic bounce: classical rho -> 0 singularity at v_c = 0.7405 vs Path 2
      Ermakov-Pinney embadon quantum bounce at s_min = ln 3 (rho_min = 0.5914).
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

with open(os.path.join(HERE, "p2_frontier1_sun_action_results.json")) as f:
    f1 = json.load(f)
with open(os.path.join(HERE, "p2_frontier3_so31_caustic_results.json")) as f:
    f3 = json.load(f)

fig, axes = figure_grid(2, 2, width=WIDE, ratio=0.77)
ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

# =========================================================================
# Panel (a): Nonlinear SU(N) Matrix Product Convergence (O(1/N) vs O(1/N^2))
# =========================================================================
rows1 = f1["nonlinear_matrix_product_convergence"]
tg1 = f1["continuum_gaunt_targets"]
Ns1 = np.array([r["N"] for r in rows1], dtype=float)
err_k_1side = np.array([abs(r["kerr_c2_onesided"] - tg1["kerr_c2"]) for r in rows1], dtype=float)
err_k_sym = np.array([r["kerr_c2_err"] for r in rows1], dtype=float)
err_gw_sym = np.array([r["gw_c2_err"] for r in rows1], dtype=float)
Ns_c4 = np.array([r["N"] for r in rows1 if r["gw_c4_err"] is not None], dtype=float)
err_c4_sym = np.array([r["gw_c4_err"] for r in rows1 if r["gw_c4_err"] is not None], dtype=float)

# Reference 1/N^2 guide line
N_ref = np.array([6.0, 24.0])
ax1.loglog(Ns1, err_k_1side, marker="o", linestyle="--", color=PALETTE["orange"], label="One-sided $T^\\dagger T$ ($O(1/N)$)")
ax1.loglog(Ns_c4, err_c4_sym, marker="^", color=PALETTE["green"], label="Sym. GW $j=4$ ($O(1/N^2)$)")
ax1.loglog(Ns1, err_gw_sym, marker="s", color=PALETTE["red"], label="Sym. GW $j=2$ ($p=2.03$)")
ax1.loglog(Ns1, err_k_sym, marker="D", color=PALETTE["blue"], label="Sym. Kerr $O(a^2)$ ($p=2.00$)")
ax1.loglog(N_ref, 0.7 * (N_ref / 6.0) ** (-2.0) * 0.0187 / 0.7, linestyle=":", color=BASELINE, linewidth=0.9)

ax1.set_xlim(5.3, 27.5)
ax1.set_ylim(7e-4, 4e-1)
ax1.set_xticks([6, 8, 12, 16, 24])
ax1.get_xaxis().set_major_formatter(lambda x, _: f"{int(x)}")
ax1.set_yticks([1e-3, 1e-2, 1e-1])
ax1.minorticks_off()
ax1.set_xlabel("Fuzzy-sphere matrix rank $N = J_{\\max}+1$")
ax1.set_ylabel(r"Nonlinear error $|c_j(N) - c_j^{\mathrm{GR}}|$")
ax1.legend(loc="lower left", fontsize=5.6, labelspacing=0.25)

# =========================================================================
# Panel (b): Second Variation delta^2 I_FSG -> RW, Zerilli, & Weyl Tidal Shift
# =========================================================================
M = 1.0
rs_p = np.linspace(-10.0, 22.0, 400)
r_p = np.where(rs_p > 4.0, rs_p, 2.0 + 2.0 * np.exp(0.5 * (rs_p - 2.0)))
r_p = np.maximum(r_p, 2.0 + 1e-10)
for _ in range(25):
    u = np.maximum(r_p / 2.0 - 1.0, 1e-14)
    r_p = np.maximum(r_p - (r_p + 2.0 * np.log(u) - rs_p) / (1.0 + 1.0 / u), 2.0 + 1e-12)
f_p = 1.0 - 2.0 * M / r_p
j = 2
lam = 0.5 * (j - 1) * (j + 2)  # = 2
V_RW_p = f_p * (j * (j + 1) / (r_p**2) - 6.0 * M / (r_p**3))
V_Z_p = f_p * (
    2 * lam**2 * (lam + 1) * r_p**3 + 6 * lam**2 * M * r_p**2 + 18 * lam * M**2 * r_p + 18 * M**3
) / ((r_p**3) * ((lam * r_p + 3.0 * M) ** 2))
V_sc_p = f_p * (j * (j + 1) / (r_p**2) + 2.0 * M / (r_p**3))

ax2.plot(rs_p, V_sc_p, linestyle=":", color=BASELINE, linewidth=1.2, label="No Weyl tidal ($+\\frac{2M}{r^3}$)")
ax2.plot(rs_p, V_RW_p, color=PALETTE["blue"], linewidth=1.3, label="Axial Regge-Wheeler ($-\\frac{6M}{r^3}$)")
ax2.plot(rs_p, V_Z_p, linestyle="--", color=PALETTE["orange"], linewidth=1.2, label="Polar Zerilli (Darboux partner)")

ax2.set_xlim(-10, 22)
ax2.set_ylim(-0.01, 0.285)
ax2.set_xticks([-10, 0, 10, 20])
ax2.set_yticks([0.0, 0.1, 0.2])
ax2.set_xlabel("Tortoise coordinate $r_* / M$ (horizon $r_* \\to -\\infty$)")
ax2.set_ylabel("Effective potential $M^2 V_{j=2}(r_*)$")
ax2.legend(loc="upper right", fontsize=5.7, labelspacing=0.25)

# =========================================================================
# Panel (c): Exponential Restoration of SO(3,1) Finite-Rapidity Boosts
# =========================================================================
p2_boost = f3["part2_wavepacket_boosts"]
specs_c = [
    ("eta_0.5", "o", PALETTE["blue"], "Rapidity $\\eta = 0.5$"),
    ("eta_1.0", "s", PALETTE["orange"], "Rapidity $\\eta = 1.0$"),
    ("eta_1.5", "^", PALETTE["red"], "Rapidity $\\eta = 1.5$"),
]
for key, mk, col, lab in specs_c:
    rws = p2_boost[key]
    Ns_b = np.array([r["N"] for r in rws], dtype=float)
    infs = np.maximum(np.array([r["infidelity"] for r in rws], dtype=float), 1e-16)
    ax3.semilogy(Ns_b, infs, marker=mk, color=col, label=lab)

ax3.set_xlim(4, 46)
ax3.set_ylim(5e-16, 2.0)
ax3.set_xticks([6, 14, 24, 32, 44])
ax3.set_yticks([1e-15, 1e-10, 1e-5, 1e0])
ax3.minorticks_off()
ax3.set_xlabel("Fuzzy-sphere matrix rank $N = J_{\\max}+1$")
ax3.set_ylabel("Boost infidelity $1 - |\\langle\\psi_N(\\eta)|\\psi_{\\mathrm{exact}}(\\eta)\\rangle|^2$")
ax3.legend(loc="lower left", fontsize=5.9, labelspacing=0.25)

# =========================================================================
# Panel (d): Nonlinear Raychaudhuri Caustic Bounce (Classical vs Path 2)
# =========================================================================
r0 = 10.0
sigma0 = 3.0
omega = sigma0 / math.sqrt(2.0)
rho_min = math.sqrt(math.log(3.0) / math.pi)
vs_d = np.linspace(0.0, 1.0, 250)
rho_cl_d = r0 * np.cos(omega * vs_d)
rho_Q_d = np.sqrt((r0 * np.cos(omega * vs_d)) ** 2 + (rho_min * np.sin(omega * vs_d)) ** 2)

ax4.axhline(0.0, color=BASELINE, linestyle=":", linewidth=0.7)
ax4.axhline(rho_min, color=PALETTE["green"], linestyle="-.", linewidth=0.9, label="Min embadon $\\rho_{\\min}=\\sqrt{\\ln 3/\\pi}$")
ax4.plot(vs_d, rho_cl_d, linestyle="--", color=PALETTE["red"], linewidth=1.2, label="Classical Raychaudhuri (caustic)")
ax4.plot(vs_d, rho_Q_d, color=PALETTE["blue"], linewidth=1.35, label="FSG quantum bounce (Ermakov)")

ax4.set_xlim(0.0, 1.0)
ax4.set_ylim(-6.5, 11.2)
ax4.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
ax4.set_yticks([-5, 0, 5, 10])
ax4.set_xlabel("Null affine parameter $v$")
ax4.set_ylabel("Areal radius $\\rho(v) / \\ell_P$")
ax4.legend(loc="lower left", fontsize=5.6, labelspacing=0.25)

panel_labels([ax1, ax2, ax3, ax4])
out_stem = os.path.join(HERE, "path2_frontier_unified")
save(fig, out_stem)
print("Saved path2_frontier_unified.svg and .pdf")
