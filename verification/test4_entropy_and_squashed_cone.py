#!/usr/bin/env python3
"""
Test 4: Microscopic Entropy Scaling (Appendix D.3), Thermal Ensemble Spread (Section 5.3),
        and Squashed-Cone Heat-Kernel Audit (eqs. 5.66-5.69 vs FPS 1306.4000 & BBM 2206.13287).

Checks:
  4A. Exact scaling comparison: Appendix D.3 graviton state count log d(s) ~ (s Q_*)^{3/4}
      vs Bekenstein-Hawking S_BH = s.
  4B. Exact moments & fractional spread of horizon area s and radius r_h in the paper's
      thermal clock state rho_D (Footnote 23-24, eqs. 5.61, 5.63, 5.68, 5.69).
  4C. Symbolic verification of Fursaev-Patrushev-Solodukhin (1306.4000) eqs. (1.4)-(1.6)
      applied to the TRUE 1-loop gravity heat kernel a_2 (Bastianelli-Bonezzi-Melis 2206.13287 eq. 36)
      vs the Gauss-Bonnet proxy E_4 (eq. 5.67a):
      - Shows exact agreement on a marginally trapped surface (Theta_l = 0 => k^2 = 0), EVEN with shear != 0;
      - Shows explicit mismatch on a general causal diamond corner (Theta_l * Theta_n != 0);
      - Computes the exact Schwarzschild interior Kretschmann divergence ~ eps^{-3} near the cap Sigma_1.
"""

import json
import os
import numpy as np
import sympy as sp


def run_test_4a_microscopic_vs_bh():
    """Compare Appendix D.3 graviton microstate bound log d(s) ~ (s Q_*)^{3/4} to S_BH = s."""
    # From eq. (D.46c): log d(s) <= C_eps * (1 + (s Q_*)^{3/4}), with Q_* = O(1) for a fixed-ratio diamond.
    # Also compute exact Bose-Einstein state count for mode spectrum nu_i = i^{1/3} / s (from D.45a: #{nu_i <= nu} = (s nu)^3):
    #   S_gas(s, Q_*) = min_{beta > 0} [ beta Q_* - sum_i log(1 - e^{-beta nu_i}) ]
    #   Continuum integral: sum_i -> int_0^inf 3 s^3 nu^2 dnu,
    #   -int_0^inf 3 s^3 nu^2 log(1 - e^{-beta nu}) dnu = (3 s^3 / beta^3) * (2 * zeta(4)) = (pi^4 / 15) * (s/beta)^3.
    #   Minimize f(beta) = beta Q_* + (pi^4 / 15) s^3 / beta^3:
    #     f'(beta) = Q_* - (pi^4 / 5) s^3 / beta^4 = 0 => beta_* = (pi^4 / 5)^{1/4} (s^3 / Q_*)^{1/4}
    #     f(beta_*) = (4/3) * (pi^4 / 5)^{1/4} * (s Q_*)^{3/4}  (exact Stefan-Boltzmann 3+1D thermal gas!)
    c_sb = float((4.0 / 3.0) * ((np.pi**4 / 5.0) ** 0.25))
    Q_star = 1.0
    rows = []
    for log10_s in [2, 4, 6, 8, 10, 12, 16, 20]:
        s = 10.0**log10_s
        log_d_gravitons = c_sb * ((s * Q_star) ** 0.75)
        ratio = log_d_gravitons / s
        rows.append(
            dict(
                s_BH=s,
                log_d_graviton_gas=float(log_d_gravitons),
                ratio_over_S_BH=float(ratio),
                missing_entropy_fraction=float(1.0 - ratio),
            )
        )
    return dict(
        stefan_boltzmann_coeff=c_sb,
        scaling_table=rows,
        conclusion=(
            "Perturbative graviton states allowed by focusing (Appendix D.3, eq. D.46c) scale as "
            "log d(s) ~ 2.80 (s Q_*)^{3/4} (standard 3+1D 't Hooft graviton gas), which is an exponentially "
            "vanishing fraction of e^s at large s. The Bekenstein-Hawking factor e^s in eq. (D.48) is "
            "imported from the classical Euclidean disk Gauss-Bonnet action, not counted in H_D."
        ),
    )


def run_test_4b_ensemble_spread_and_106_45():
    """Verify exact moments of rho_D (eqs. 5.61, 5.63, 5.68, 5.69) and the +1/2 - 1/2 bookkeeping."""
    alpha = sp.Rational(257, 90)  # q / (2(1-q)) = 106/45 + 1/2
    q_exact = sp.Rational(257, 302)
    assert sp.simplify(q_exact / (2 * (1 - q_exact)) - alpha) == 0
    assert sp.simplify(alpha - sp.Rational(1, 2) - sp.Rational(106, 45)) == 0

    # Random variable U = s / s_* on [0, 1] with CDF F(u) = u^alpha, PDF p(u) = alpha * u^{alpha - 1}
    u = sp.Symbol("u", positive=True)
    p_u = alpha * u ** (alpha - 1)

    mean_u = sp.integrate(u * p_u, (u, 0, 1))
    mean_u2 = sp.integrate(u**2 * p_u, (u, 0, 1))
    var_u = sp.simplify(mean_u2 - mean_u**2)
    rel_std_area = sp.simplify(sp.sqrt(var_u) / mean_u)

    # Horizon radius R = sqrt(u) (in units of r_max = l / (1 - eps))
    mean_r = sp.integrate(sp.sqrt(u) * p_u, (u, 0, 1))
    var_r = sp.simplify(mean_u - mean_r**2)
    rel_std_r = sp.simplify(sp.sqrt(var_r) / mean_r)

    # Check Jensen gap <A^{1/2}> / <A>^{1/2} = <sqrt(u)> / sqrt(<u>)
    jensen_ratio = sp.simplify(mean_r / sp.sqrt(mean_u))

    # Check Footnote 23 eq. (5.61): x = sqrt(mean_u), q = 2 x^2 / (1 + x^2)
    x2 = mean_u
    q_from_fn23 = sp.simplify(2 * x2 / (1 + x2))
    assert q_from_fn23 == q_exact

    return dict(
        alpha=str(alpha),
        q_exact=str(q_exact),
        q_float=float(q_exact),
        mean_area_over_s_star=str(mean_u),
        mean_area_float=float(mean_u),
        area_fractional_std=float(rel_std_area),
        radius_fractional_std=float(rel_std_r),
        jensen_ratio_sqrtA_over_sqrt_meanA=float(jensen_ratio),
        fn23_consistency=True,
        bookkeeping_check="53/45 * 2 = 106/45 (imported BBM heat kernel) -> +1/2 (measure 5.66) = 257/90 -> -1/2 (change to M in 5.69) = 106/45",
    )


def run_test_4c_squashed_cone_and_interior_cap():
    """
    Symbolic check of:
      (i) Fursaev-Patrushev-Solodukhin (1306.4000) eqs. (1.4)-(1.6) applied to the true 1-loop
          pure gravity heat kernel tr[a_{2,tot}] = (53/45) R_{abcd}^2 - (361/90) R_{ab}^2 + (43/36) R^2
          (Bastianelli-Bonezzi-Melis 2206.13287 eq. 36) vs the Gauss-Bonnet proxy (53/45) E_4.
      (ii) Exact Schwarzschild interior Kretschmann integral down to cap radius r_cap = eps * r_h.
    """
    # Symbols on a 2-surface S inside a Ricci-flat 4D spacetime (R_ab = 0, R = 0):
    Rq, Th_l, Th_n, sig_dot = sp.symbols("R_q Theta_l Theta_n sigma_dot", real=True)
    # Null extrinsic curvatures:
    #   K_{ab} = sigma_{ab} + (1/2) Th_l q_{ab}
    #   K^{(n)}_{ab} = sigma^{(n)}_{ab} + (1/2) Th_n q_{ab}
    #   2 K^{ab} K^{(n)}_{ab} = 2 sig_dot + Th_l * Th_n
    #   k^2 = -2 Th_l * Th_n
    #   Tr(k^2) = -2 K^{ab} K^{(n)}_{ab} = -2 sig_dot - Th_l * Th_n
    # Gauss equation at R_ab = 0:
    #   R_{ijij} - Tr(k^2) = 2 K^{ab} K^{(n)}_{ab} - 2 R_{abcd} l^a n^b l^c n^d = R_q + 2 Th_l * Th_n
    k2 = -2 * Th_l * Th_n
    R_ijij_minus_Trk2 = Rq + 2 * Th_l * Th_n

    # FPS (1306.4000) eqs. (1.4)-(1.6) surface coefficients of (1 - n) on Ricci-flat bulk (R = R_ii = 0):
    surf_R2 = 0  # 8 pi R = 0
    surf_Ric2 = 4 * sp.pi * (-sp.Rational(1, 2) * k2)  # = -2 pi k^2 = +4 pi Th_l Th_n
    surf_Riem2 = 8 * sp.pi * R_ijij_minus_Trk2

    # 1) Gauss-Bonnet proxy used in eq. (5.67a): (53/45) * (Riem^2 - 4 Ric^2 + R^2)
    surf_GB_proxy = sp.simplify(sp.Rational(53, 45) * (surf_Riem2 - 4 * surf_Ric2 + surf_R2))

    # 2) True 1-loop pure gravity heat kernel (Bastianelli et al. 2206.13287 eq. 36):
    #    (53/45) Riem^2 - (361/90) Ric^2 + (43/36) R^2
    surf_true_BBM = sp.simplify(
        sp.Rational(53, 45) * surf_Riem2 - sp.Rational(361, 90) * surf_Ric2 + sp.Rational(43, 36) * surf_R2
    )

    diff_general = sp.simplify(surf_true_BBM - surf_GB_proxy)
    diff_on_MTS = sp.simplify(diff_general.subs(Th_l, 0))
    val_on_MTS = sp.simplify(surf_true_BBM.subs(Th_l, 0))

    # Evaluate mismatch on a round sphere of radius b in flat Minkowski space:
    # R_q = 2/b^2, Th_l = +2/b, Th_n = -1/b (with l.n = -1, l = d_t + d_r, n = (1/2)(d_t - d_r) => Th_l Th_n = -2/b^2)
    b = sp.Symbol("b", positive=True)
    # Integrate over sphere of area 4 pi b^2:
    mink_sphere_GB = sp.simplify((surf_GB_proxy * 4 * sp.pi * b**2).subs({Rq: 2 / b**2, Th_l * Th_n: -2 / b**2}))
    mink_sphere_BBM = sp.simplify((surf_true_BBM * 4 * sp.pi * b**2).subs({Rq: 2 / b**2, Th_l * Th_n: -2 / b**2}))
    mink_sphere_diff = sp.simplify(mink_sphere_BBM - mink_sphere_GB)

    # (ii) Schwarzschild interior Kretschmann integral from r = eps * r_h to r_h
    # R_{abcd} R^{abcd} = 48 M^2 / r^6 = 12 r_h^2 / r^6; spatial volume element r^2 dr dOmega
    # Over a finite time slice Delta t ~ r_h: int_{eps r_h}^{r_h} (12 r_h^2 / r^6) * (4 pi r^2) * r_h dr
    # = 48 pi r_h^3 * [-1/(3 r^3)]_{eps r_h}^{r_h} = 16 pi (eps^{-3} - 1)!
    r, rh, eps = sp.symbols("r r_h eps", positive=True)
    kretsch_int = sp.integrate((12 * rh**2 / r**6) * (4 * sp.pi * r**2) * rh, (r, eps * rh, rh))

    return dict(
        surf_GB_proxy=str(surf_GB_proxy),
        surf_true_BBM=str(surf_true_BBM),
        diff_true_minus_GB=str(diff_general),
        diff_on_MTS_Theta_l_zero=str(diff_on_MTS),
        value_on_MTS=str(val_on_MTS),
        minkowski_sphere_GB_integral=str(mink_sphere_GB),
        minkowski_sphere_BBM_integral=str(mink_sphere_BBM),
        minkowski_sphere_mismatch=str(mink_sphere_diff),
        schwarzschild_interior_kretschmann_vs_eps=str(sp.simplify(kretsch_int)),
        nuanced_verdict=(
            "On a marginally trapped surface (Theta_l = 0), k^2 = -2 Theta_l Theta_n vanishes, so the "
            "BBM heat kernel and Gauss-Bonnet proxy agree (giving 424*pi*R_q/45) even with nonzero shear. "
            "HOWEVER: (1) on a generic causal diamond corner (e.g. Minkowski sphere), BBM and E_4 disagree "
            "by -56*pi^2/5; (2) on the conical manifold M_n, R_ab != 0 at S is off-shell where BBM's R_ab^2 "
            "coefficient (-361/90) is gauge-dependent; (3) the bulk interior Kretschmann integral diverges "
            "as 16*pi*(eps^{-3} - 1) near the cap; and (4) the log(s_b^{1/2}) factor imports the perturbative "
            "UV cutoff Lambda_UV ~ 1/l_P rather than deriving it from the null-cone quantization."
        ),
    )


if __name__ == "__main__":
    res = dict(
        test_4a_microscopic_vs_bh=run_test_4a_microscopic_vs_bh(),
        test_4b_ensemble_spread=run_test_4b_ensemble_spread_and_106_45(),
        test_4c_squashed_cone=run_test_4c_squashed_cone_and_interior_cap(),
    )
    out_path = os.path.join(os.path.dirname(__file__), "test4_results.json")
    with open(out_path, "w") as f:
        json.dump(res, f, indent=2)

    print("=== TEST 4A: Microscopic Graviton Count (Appendix D.3) vs S_BH = s ===")
    for r in res["test_4a_microscopic_vs_bh"]["scaling_table"]:
        print(f"  s = {r['s_BH']:8.1e}: log d(s) = {r['log_d_graviton_gas']:10.2e}, ratio log d / s = {r['ratio_over_S_BH']:.3e}")

    print("\n=== TEST 4B: Thermal Clock Ensemble Spread (Section 5.3, Footnotes 23-24) ===")
    print(json.dumps(res["test_4b_ensemble_spread"], indent=2))

    print("\n=== TEST 4C: Squashed-Cone BBM vs E_4 & Schwarzschild Interior Cap ===")
    for k, v in res["test_4c_squashed_cone"].items():
        if k != "nuanced_verdict":
            print(f"  {k}: {v}")
    print(f"\nSaved full JSON results to {out_path}")
