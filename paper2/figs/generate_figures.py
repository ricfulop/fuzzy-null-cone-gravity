#!/usr/bin/env python3
"""Generate publication-quality vector figures for Paper II (v0.2) and the Paper I audit.

Reads only verified numerical outputs from verification/*.json.
Writes PDF + SVG pairs into figs/ using orx_figstyle (built at exact TEXT = 5.5 in width)
and copies the SVGs to the durable artifacts folder.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from orx_figstyle import BASELINE, MUTED, PALETTE, TEXT, figure_grid, panel_labels, save, use_style


def load_json(name: str) -> dict:
    with open(os.path.join(ROOT, "verification", name)) as fh:
        return json.load(fh)


def make_fig1_scattering_and_ablations(d_bs: dict):
    """Fig 1: (a) 4-graviton amplitude residual across channels & ablations;
    (b) Moyal noncommutative scaling of eps_4 vs theta_patch."""
    use_style()
    fig, (ax1, ax2) = figure_grid(nrows=1, ncols=2, width=TEXT, ratio=0.47)

    # Panel (a): Horizontal log bar / dot chart of eps_4 across channels & ablations
    c0 = d_bs["corrected_action_F0_0"]
    c1 = d_bs["corrected_action_F0_1"]
    ab = d_bs["ablations"]

    labels = [
        "Corrected BS ($--++$)",
        "Corrected BS ($-+-+$)",
        "Corrected BS ($++++$)",
        "Corrected BS ($-+++$)",
        "Corrected BS ($F_0=1$)",
        "No SL(2,R) quartic",
        "Paper I Eq. (A.1) ($U(1)$)",
        "No Prop. 1 constr.",
    ]
    vals = [
        c0["eps4_mmpp"],
        c0["eps4_mpmp"],
        c0["eps4_pppp"],
        c0["eps4_mppp"],
        c1["eps4_mmpp"],
        ab["drop_sl2_quartic_mmpp"],
        ab["paper1_eq_A1_u1_symmetry_mmpp"],
        ab["drop_prop1_constraint_mmpp"],
    ]
    colors = [
        PALETTE["blue"],
        PALETTE["blue"],
        PALETTE["green"],
        PALETTE["green"],
        PALETTE["cyan"],
        PALETTE["orange"],
        PALETTE["red"],
        PALETTE["red"],
    ]

    y = np.arange(len(labels))
    # Draw stems from 1e-15 floor to value (log scale)
    floor = 1e-15
    for yi, v, col in zip(y, vals, colors):
        ax1.plot([floor, v], [yi, yi], color=col, linewidth=1.6, zorder=2)
        ax1.scatter([v], [yi], color=col, s=22, zorder=3)

    ax1.axvline(1e-11, color=BASELINE, linestyle=":", linewidth=0.8, zorder=1)
    ax1.text(2.5e-11, 0.35, "FP64 match", fontsize=6.5, color=BASELINE, va="bottom", ha="left")
    ax1.set_xscale("log")
    ax1.set_xlim(floor, 3e4)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=6.5)
    ax1.invert_yaxis()
    ax1.set_xlabel(r"4-graviton residual $\epsilon_4$ (dimensionless)")
    ax1.grid(True, axis="x", color="#D9D9D9", linewidth=0.4)
    ax1.grid(False, axis="y")

    # Panel (b): Moyal noncommutative deformation scaling vs theta_patch
    scan = [r for r in d_bs["moyal_noncommutative_scan"] if r["theta_patch"] > 0]
    th = np.array([r["theta_patch"] for r in scan])
    e_sym_mhv = np.array([r["eps4_sym_mmpp"] for r in scan])
    e_sym_pppp = np.array([r["eps4_sym_pppp"] for r in scan])
    e_one_mhv = np.array([r["eps4_onesided_mmpp"] for r in scan])

    ax2.loglog(th, e_one_mhv, "s-", color=PALETTE["red"], markersize=3.5, label=r"One-sided ($--++$)")
    ax2.loglog(th, e_sym_mhv, "o-", color=PALETTE["blue"], markersize=3.5, label=r"Symmetrized ($--++$)")
    ax2.loglog(th, e_sym_pppp, "^--", color=PALETTE["orange"], markersize=3.5, label=r"Symmetrized ($++++$)")

    # Reference slopes O(theta) and O(theta^2)
    th_ref = np.array([2e-4, 3e-2])
    ax2.loglog(th_ref, 43.0 * th_ref, ":", color=BASELINE, linewidth=0.8)
    ax2.text(1.5e-4, 8.5e-2, r"$\propto \theta_{\mathrm{patch}}$", fontsize=6.5, color=BASELINE)
    ax2.loglog(th_ref, 21.0 * (th_ref**2), "--", color=BASELINE, linewidth=0.8)
    ax2.text(1.5e-4, 5.5e-5, r"$\propto \theta_{\mathrm{patch}}^2$", fontsize=6.5, color=BASELINE)

    ax2.set_xlabel(r"Noncommutativity $\theta_{\mathrm{patch}} E_{\mathrm{ref}}^2$")
    ax2.set_ylabel(r"Amplitude residual $\epsilon_4$")
    ax2.set_xlim(6e-5, 2.0)
    ax2.set_ylim(5e-8, 8e1)
    ax2.legend(loc="lower right", fontsize=6.5)

    panel_labels((ax1, ax2))
    save(fig, os.path.join(HERE, "fig1_scattering_and_ablations"))


def make_fig2_fuzzy_scaling_resolutions(d_sc: dict):
    """Fig 2: (a) Exact Wigner-6j vertex error across four asymptotic sequences;
    (b) Band-limited commutative Gaunt product associativity defect vs fuzzy star error."""
    use_style()
    fig, (ax1, ax2) = figure_grid(nrows=1, ncols=2, width=TEXT, ratio=0.56)

    seq = d_sc["scaling_sequences"]

    # Panel (a): Vertex relative error vs matrix rank N
    fix = seq["fixed_j2_sequence"]
    ax1.loglog(
        [r["N"] for r in fix], [r["err_sym"] for r in fix],
        "o-", color=PALETTE["blue"], markersize=3.2,
        label=rf"Fixed $j=2$ sym ($p={seq['fixed_j2_exponents']['p_sym']:.2f}$)",
    )
    ax1.loglog(
        [r["N"] for r in fix], [r["err_onesided"] for r in fix],
        "v:", color=PALETTE["cyan"], markersize=3.0,
        label=rf"Fixed $j=2$ 1-sided ($p={seq['fixed_j2_exponents']['p_onesided']:.2f}$)",
    )

    meso_eq = seq["mesoscopic_sqrtN_sequence"]["equilateral_c1"]
    meso_col = seq["mesoscopic_sqrtN_sequence"]["collinear_c1"]
    ax1.loglog(
        [r["N"] for r in meso_eq], [r["err_sym"] for r in meso_eq],
        "s-", color=PALETTE["orange"], markersize=3.2,
        label=r"Mesoscopic $j=\sqrt{N}$ equil.",
    )
    ax1.loglog(
        [r["N"] for r in meso_col], [r["err_sym"] for r in meso_col],
        "d--", color=PALETTE["green"], markersize=3.0,
        label=r"Mesoscopic $j=\sqrt{N}$ collin.",
    )

    sat10 = seq["fixed_physical_momentum_saturated_cap"]["alpha_0.10"]
    ax1.loglog(
        [r["N"] for r in sat10], [r["err_sym"] for r in sat10],
        "^-", color=PALETTE["red"], markersize=3.4,
        label=r"Fixed $p_\perp$ ($j=0.10N$)",
    )

    dec = seq["decoupled_regulator_j_N13"]
    ax1.loglog(
        [r["N"] for r in dec], [r["err_sym"] for r in dec],
        "P--", color=PALETTE["purple"], markersize=3.4,
        label=r"Decoupled $j \propto N^{1/3}$",
    )

    ax1.set_xlabel(r"Matrix rank $N$")
    ax1.set_ylabel(r"Vertex error $|V_3^{(N)}/V_3^{\mathrm{cont}} - 1|$")
    ax1.set_xlim(7, 620)
    ax1.set_ylim(1.5e-5, 5.0)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.48, -0.24), fontsize=6.0, ncol=2, columnspacing=0.8)

    # Panel (b): Route C band-limited Gaunt product associativity defect vs j/N (N=24)
    gaunt = d_sc["resolution_comparison"]["route_C_gaunt_associativity_scan"]
    jn_g = np.array([r["j_over_N"] for r in gaunt])
    def_g = np.array([max(r["associativity_defect"], 1e-14) for r in gaunt])

    sat05 = seq["fixed_physical_momentum_saturated_cap"]["alpha_0.05"]
    sat15 = seq["fixed_physical_momentum_saturated_cap"]["alpha_0.15"]
    jn_fuz = [2.0 / 192.0]
    err_fuz = [next(r["err_sym"] for r in fix if r["N"] == 192)]
    for arr in (sat05, sat10, sat15):
        row = next(r for r in arr if r["N"] == 200)
        jn_fuz.append(row["alpha_actual"])
        err_fuz.append(row["err_sym"])

    ax2.semilogy(jn_g, def_g, "o-", color=PALETTE["green"], markersize=3.5, label=r"Band-limited Gaunt ($N=24$)")
    ax2.semilogy(jn_fuz, err_fuz, "s--", color=PALETTE["red"], markersize=3.5, label=r"Fuzzy star ($N=200$, equil.)")
    ax2.axvline(0.5, color=BASELINE, linestyle=":", linewidth=0.9)
    ax2.text(0.53, 3e-11, r"$j = N/2$" "\n(Planck scale)", fontsize=6.5, color=BASELINE, va="bottom")

    ax2.set_xlabel(r"Angular ratio $j / N$")
    ax2.set_ylabel(r"Algebraic / vertex defect")
    ax2.set_xlim(0.0, 0.98)
    ax2.set_ylim(2e-15, 6.0)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.50, -0.24), fontsize=6.0, ncol=1)

    panel_labels((ax1, ax2))
    save(fig, os.path.join(HERE, "fig2_fuzzy_scaling_resolutions"))


def make_fig3_projection_defects(d_pr: dict):
    """Fig 3: (a) Prop 2 resolvent self-energy & missed energy shift vs budget B;
    (b) Prop 3 spin-2 Lorentz boost compression defect (1-body boundary vs 2-body leakage)."""
    use_style()
    fig, (ax1, ax2) = figure_grid(nrows=1, ncols=2, width=TEXT, ratio=0.47)

    # Panel (a): Prop 2 vs budget B
    bscan = d_pr["prop2_occupation_budget_resolvent"]["budget_scan"]
    Bs = np.array([r["B"] for r in bscan])
    sig = np.array([r["sigma_norm"] for r in bscan])
    dE = np.array([r["missed_energy_shift"] for r in bscan])
    bnd = np.array([r["eq16_upper_bound"] for r in bscan])

    ax1.semilogy(Bs, bnd, "--", color=BASELINE, linewidth=1.0, label=r"Prop. 2 bound $\|PHQ\|^2/d$")
    ax1.semilogy(Bs, sig, "o-", color=PALETTE["blue"], markersize=3.5, label=r"Self-energy $\|\Sigma_P(E_0)\|$")
    ax1.semilogy(Bs, dE, "^:", color=PALETTE["orange"], markersize=3.2, label=r"Shift $|E_{\mathrm{exact}} - E_{PHP}|$")
    ax1.set_xlabel(r"Occupation budget $B$ (dimensionless)")
    ax1.set_ylabel(r"Virtual-state shift / bound")
    ax1.set_xlim(5.2, 24.8)
    ax1.set_ylim(1.2e-2, 1.5)
    ax1.legend(loc="lower left", fontsize=6.3)

    # Panel (b): Prop 3 spin-2 boost compression defect
    sh = d_pr["prop3_lorentz_compression_defect"]["single_body_shell_scan"]
    tb = d_pr["prop3_lorentz_compression_defect"]["two_body_total_budget_leakage"]

    j_cuts = np.array([r["j_cut"] for r in sh])
    b_def = np.array([r["boundary_shell_defect"] for r in sh])
    ax2.plot(j_cuts, b_def, "o-", color=PALETTE["red"], markersize=3.8, label=r"1-body boundary $j=j_{\mathrm{cut}}$")

    # Overlay 2-body total-budget leakage onto lowest mode j1=2 (plotted on top axis or inset?
    # Better: plot 1-body boundary defect as main curve, and mark 2-body (j1=2, j2=3) defect
    # at B=9, 12, 13 cleanly via a crisp horizontal reference / scatter using a shared y-axis!
    # Wait: even cleaner — use shared y-axis defect norm and plot a grouped comparison or twin x-axis?
    # Avoid twin axes (visually confusing); instead plot 2-body leakage vs B directly in (b)
    # and annotate the 1-body boundary defect values! Let's check: BOTH 1-body j_cut in [3..8]
    # and 2-body B in [9, 12, 13, 16] fit cleanly on ONE integer/dimensionless axis? No, different x variables!
    # Cleanest solution: plot 2-body total-budget leakage vs B on ax2, OR plot 1-body boundary defect vs j_cut
    # with an explicit horizontal dashed marker showing 2-body (j1=2, j2=3) leakage at B=9 (7.14) and B=12 (5.14)!
    ax2.axhline(tb[0]["defect_on_j1_2_j2_3"], color=PALETTE["purple"], linestyle="--", linewidth=1.0,
                label=r"2-body $(j_1=2,j_2=3)$ at $B=9$")
    ax2.axhline(tb[1]["defect_on_j1_2_j2_3"], color=PALETTE["blue"], linestyle=":", linewidth=1.0,
                label=r"2-body $(j_1=2,j_2=3)$ at $B=12$")
    ax2.scatter(j_cuts, [0.0] * len(j_cuts), marker="s", color=PALETTE["green"], s=18,
                label=r"1-body interior ($j < j_{\mathrm{cut}}$)")

    ax2.set_xlabel(r"Angular truncation shell $j_{\mathrm{cut}}$")
    ax2.set_ylabel(r"Compression defect $\|D_{xy}\|_2$")
    ax2.set_xlim(2.6, 8.4)
    ax2.set_ylim(-2.5, 40.0)
    ax2.legend(loc="upper left", fontsize=6.0)

    panel_labels((ax1, ax2))
    save(fig, os.path.join(HERE, "fig3_projection_defects"))


def main():
    d_bs = load_json("bondi_sachs_4pt_results.json")
    d_sc = load_json("fuzzy_scaling_results.json")
    d_pr = load_json("projection_defects_results.json")

    make_fig1_scattering_and_ablations(d_bs)
    make_fig2_fuzzy_scaling_resolutions(d_sc)
    make_fig3_projection_defects(d_pr)

    art_dir = os.environ.get("PAPER2_ARTIFACTS_DIR")
    if art_dir:
        os.makedirs(art_dir, exist_ok=True)
        for stem in (
            "fig1_scattering_and_ablations",
            "fig2_fuzzy_scaling_resolutions",
            "fig3_projection_defects",
        ):
            shutil.copy2(os.path.join(HERE, f"{stem}.svg"), os.path.join(art_dir, f"{stem}.svg"))
            shutil.copy2(os.path.join(HERE, f"{stem}.pdf"), os.path.join(art_dir, f"{stem}.pdf"))
        print(f"Copied vector figures to {art_dir}")


if __name__ == "__main__":
    main()
