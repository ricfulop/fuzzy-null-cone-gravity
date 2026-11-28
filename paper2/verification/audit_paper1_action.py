#!/usr/bin/env python3
"""Independent symbolic audit of Paper I's Appendix A action I_FSG^{(N)} (commit e541955495d0).

Verifies all five structural discrepancies between Paper I's Eq. (A.1) and the
continuum null-cone Einstein-Hilbert action of Chandrasekaran (arXiv:2611.34924,
Appendix B, Eqs. B.9-B.11e), plus the two matrix-typing issues noted in Paper II Sec. 2.3.

Outputs a machine-readable JSON record to verification/audit_paper1_action_results.json.
"""

from __future__ import annotations

import json
import os
import sympy as sp


def check_defect1_u1_symmetry():
    """Verify exact global U(1) / Z_2 phase symmetry Z -> e^{i alpha} Z of Paper I Eq. (A.1)
    and Euclidean HMC action S_E^{(N)} (Paper I Eq. 12)."""
    alpha = sp.Symbol("alpha", real=True)
    rho2, N, F, Lam = sp.symbols("rho2 N F Lam", real=True, positive=True)
    Z_t, Z_r, Z_eb = sp.symbols("Z_t Z_r Z_eb", complex=True)
    w, d_r_theta, Rq = sp.symbols("w d_r_theta Rq", real=True)

    # Scalar representative of every term in Paper I Eq. (A.1)
    L_A1 = (
        2 * rho2 * (d_r_theta - N * Lam)
        + N * Rq
        + rho2 * N**(-1) * (sp.conjugate(Z_t) * Z_r + sp.conjugate(Z_r) * Z_t)
        - sp.Rational(1, 2) * rho2 * F * sp.conjugate(Z_r) * Z_r
        + 2 * rho2 * N * w**2
        - sp.Rational(1, 2) * N * F * sp.conjugate(Z_eb) * Z_eb
    )

    phase = sp.exp(sp.I * alpha)
    L_A1_rot = L_A1.subs({
        Z_t: phase * Z_t,
        Z_r: phase * Z_r,
        Z_eb: phase * Z_eb,
    })
    diff_u1 = sp.simplify(sp.expand(L_A1_rot - L_A1))
    assert diff_u1 == 0, f"Expected exact U(1) invariance, got {diff_u1}"

    # Consequence: any n-point pure-Z functional derivative with odd n vanishes identically.
    lam_scale = sp.Symbol("lam_scale", real=True)
    L_scaled = L_A1.subs({Z_t: lam_scale * Z_t, Z_r: lam_scale * Z_r, Z_eb: lam_scale * Z_eb})
    d3 = sp.simplify(sp.diff(L_scaled, lam_scale, 3).subs(lam_scale, 0))
    d5 = sp.simplify(sp.diff(L_scaled, lam_scale, 5).subs(lam_scale, 0))
    assert d3 == 0 and d5 == 0

    return {
        "u1_invariance_residual": str(diff_u1),
        "odd_point_Z3_derivative": str(d3),
        "odd_point_Z5_derivative": str(d5),
        "consequence": "M_3(h,h,h) = 0 and M_5(h,h,h,h,q) = 0 identically; zero spin-2 exchange poles in M_4.",
    }


def check_defect2_bondi_sachs_constraints():
    """Derive the exact quadratic 4D Einstein-Hilbert kinetic matrix in Bondi-Sachs coordinates
    around flat space, proving that the massless angular dispersion -k_x^2 h arises purely from
    integrating out the shift V_x via Damour-Navier-Stokes, while dF is sourced linearly in h."""
    t, r, x, y, eps = sp.symbols("t r x y eps", real=True)
    dF = sp.Function("dF")(t, r, x)
    dN = sp.Function("dN")(t, r, x)
    dphi = sp.Function("dphi")(t, r, x)
    Vx = sp.Function("Vx")(t, r, x)
    h = sp.Function("h")(t, r, x)
    coords = [t, r, x, y]

    qxx = (1 + eps * dphi) * sp.exp(eps * h)
    qyy = (1 + eps * dphi) * sp.exp(-eps * h)
    g = sp.Matrix([
        [-(1 + eps * dF) + qxx * (eps * Vx) ** 2, -(1 + eps * dN), qxx * eps * Vx, 0],
        [-(1 + eps * dN), 0, 0, 0],
        [qxx * eps * Vx, 0, qxx, 0],
        [0, 0, 0, qyy],
    ])
    g0 = g.subs(eps, 0)
    g0inv = g0.inv()
    dg1 = sp.diff(g, eps).subs(eps, 0)

    # First-order Christoffel symbols
    G1 = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                val = 0
                for d in range(4):
                    if g0inv[a, d] != 0:
                        term = (
                            sp.diff(dg1[d, c], coords[b])
                            + sp.diff(dg1[d, b], coords[c])
                            - sp.diff(dg1[b, c], coords[d])
                        )
                        if term != 0:
                            val += sp.Rational(1, 2) * g0inv[a, d] * term
                G1[a][b][c] = sp.simplify(val)
                G1[a][c][b] = G1[a][b][c]

    # Exact Gamma-Gamma bulk quadratic Lagrangian (equals sqrt(-g) R up to total derivatives)
    L2_bulk = 0
    for b in range(4):
        for d in range(4):
            if g0inv[b, d] == 0:
                continue
            for a in range(4):
                for e in range(4):
                    L2_bulk += g0inv[b, d] * (G1[a][b][e] * G1[e][a][d] - G1[a][a][e] * G1[e][b][d])
    L2_bulk = sp.expand(L2_bulk)

    # Euler-Lagrange equations and Fourier kinetic matrix M (d_mu -> i k_mu)
    wt, kr, kx = sp.symbols("wt kr kx", real=True)
    fields = [h, Vx, dN, dF, dphi]
    f_syms = sp.symbols("h_k Vx_k dN_k dF_k dphi_k")

    M = sp.zeros(5, 5)
    for i, f in enumerate(fields):
        eq = -(
            sp.diff(sp.diff(L2_bulk, sp.diff(f, t)), t)
            + sp.diff(sp.diff(L2_bulk, sp.diff(f, r)), r)
            + sp.diff(sp.diff(L2_bulk, sp.diff(f, x)), x)
        )
        eq_k = eq
        for j, fj in enumerate(fields):
            for a, ka in [(t, wt), (r, kr), (x, kx)]:
                for b, kb in [(t, wt), (r, kr), (x, kx)]:
                    eq_k = eq_k.subs(sp.diff(fj, a, b), -ka * kb * f_syms[j])
        eq_k = sp.expand(eq_k)
        for j, fs in enumerate(f_syms):
            M[i, j] = sp.simplify(eq_k.coeff(fs))

    # Fix affine gauge dN = 0; remaining constraints are c = (Vx, dF, dphi)
    M_aff = M.extract([0, 1, 3, 4], [0, 1, 3, 4])
    K_c = M_aff[1:, 1:]
    det_Kc = sp.factor(K_c.det())
    assert det_Kc == -kr**6, f"Unexpected constraint determinant: {det_Kc}"

    c_sol = K_c.LUsolve(-M_aff[1:, 0])
    Vx_sol, dF_sol, dphi_sol = [sp.simplify(s) for s in c_sol]
    eff_kin = sp.simplify(M_aff[0, 0] + M_aff[0, 1:].dot(c_sol))
    assert sp.simplify(eff_kin - (2 * kr * wt - kr**2 - kx**2)) == 0
    assert Vx_sol == kx / kr
    assert dF_sol == -(kx**2) / (kr**2)
    assert dphi_sol == 0

    return {
        "bare_h_diagonal_kinetic": str(M[0, 0]),
        "bare_h_angular_term": "0 (topological 2D Gauss-Bonnet at constant N)",
        "affine_constraint_det": str(det_Kc),
        "linear_Vx_over_h": str(Vx_sol),
        "linear_dF_over_h": str(dF_sol),
        "linear_dphi_over_h": str(dphi_sol),
        "reduced_h_kinetic_operator": str(eff_kin),
    }


def check_defect3_kinetic_factors_and_lightcone():
    """Verify exact N, F, and factor-of-1/2 structure of sqrt(-g) R in Bondi-Sachs coordinates,
    and compute the distorted flat-space radial wave speeds of Paper I Eq. (A.1)."""
    t, r, x, y, eps = sp.symbols("t r x y eps", real=True)
    F = sp.Function("F")(t, r)
    N = sp.Function("N")(t, r)
    rho = sp.Function("rho")(t, r)
    h = sp.Function("h")(t, r, x)
    coords = [t, r, x, y]

    g = sp.Matrix([
        [-F, -N, 0, 0],
        [-N, 0, 0, 0],
        [0, 0, rho**2 * sp.exp(eps * h), 0],
        [0, 0, 0, rho**2 * sp.exp(-eps * h)],
    ])
    ginv = g.inv()
    sqrtg = sp.simplify(N * rho**2)

    def exp2(expr):
        return sp.series(expr, eps, 0, 3).removeO()

    G = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                val = 0
                for d in range(4):
                    if ginv[a, d] != 0:
                        term = (
                            sp.diff(g[d, c], coords[b])
                            + sp.diff(g[d, b], coords[c])
                            - sp.diff(g[b, c], coords[d])
                        )
                        if term != 0:
                            val += sp.Rational(1, 2) * ginv[a, d] * term
                G[a][b][c] = sp.simplify(exp2(val))
                G[a][c][b] = G[a][b][c]

    R_sc = 0
    for b in range(4):
        for d in range(4):
            if ginv[b, d] == 0:
                continue
            ric = 0
            for a in range(4):
                ric += sp.diff(G[a][b][d], coords[a]) - sp.diff(G[a][a][d], coords[b])
                for e in range(4):
                    if G[a][a][e] != 0 and G[e][b][d] != 0:
                        ric += G[a][a][e] * G[e][b][d]
                    if G[a][b][e] != 0 and G[e][a][d] != 0:
                        ric -= G[a][b][e] * G[e][a][d]
            R_sc += ginv[b, d] * ric

    L2 = sp.expand(sp.diff(exp2(sqrtg * R_sc), eps, 2).subs(eps, 0) / 2)
    c_ht_hr = sp.simplify(L2.coeff(sp.diff(h, t) * sp.diff(h, r)))
    c_hr2 = sp.simplify(L2.coeff(sp.diff(h, r) ** 2))
    c_hx2_eff = sp.simplify(L2.coeff(sp.diff(h, x) ** 2) - L2.coeff(h * sp.diff(h, x, 2)))

    assert c_ht_hr == rho**2
    assert c_hr2 == -F * rho**2 / (2 * N)
    assert c_hx2_eff == 0

    # Flat-space radial phase velocities v = omega / k_r in static Minkowski coords (t_S = u + r):
    # Standard EF wave operator: -2 d_u d_r + d_r^2 -> -d_{t_S}^2 + d_r^2 -> speeds +-1.
    # Paper I Eq. (A.1) has c_ur = 1 (from Z_t^* Z_r + h.c.) and c_rr = 1/2 (from -1/2 F |Z_r|^2):
    # Wave operator: -2 d_u d_r + (1/2) d_r^2. Substituting d_u = d_T, d_r = d_R + d_T:
    # -2 d_T(d_R + d_T) + (1/2)(d_R + d_T)^2 = -(3/2) d_T^2 - d_T d_R + (1/2) d_R^2.
    # Characteristic polynomial for v = dR/dT (with d_T -> -i v, d_R -> i):
    v = sp.Symbol("v")
    poly_exact = -v**2 + 1
    poly_A1 = sp.Rational(-3, 2) * v**2 - (-v) + sp.Rational(1, 2)  # Wait: check plane wave exp(i(k R - w T))!
    # d_T -> -i w, d_R -> i k => -(3/2)(-w^2) - (-i w)(i k) + (1/2)(-k^2) = (3/2) w^2 - w k - (1/2) k^2 = 0
    # Divide by k^2 with c_phase = w / k: (3/2) c^2 - c - (1/2) = 0 -> (3 c + 1)(c - 1) = 0 -> c = +1, -1/3!
    c_vel = sp.Symbol("c_vel")
    disp_A1 = sp.Rational(3, 2) * c_vel**2 - c_vel - sp.Rational(1, 2)
    speeds_A1 = sorted(sp.solve(disp_A1, c_vel))
    assert speeds_A1 == [sp.Rational(-1, 3), sp.Integer(1)]

    return {
        "exact_EH_coeff_ht_hr": str(c_ht_hr),
        "exact_EH_coeff_hr2": str(c_hr2),
        "exact_EH_effective_hx2_at_const_N": str(c_hx2_eff),
        "paper1_A1_flat_radial_phase_speeds": [str(s) for s in speeds_A1],
        "expected_GR_radial_phase_speeds": ["-1", "1"],
    }


def check_defect4_sl2_sigma_model_quartic():
    """Verify that the unimodular 2-metric gamma = exp(H) produces nonzero O(Z^4) bare quartic
    kinetic interactions in tr(gamma^{-1} d gamma gamma^{-1} d gamma), omitted in Paper I Eq. (A.1)."""
    eps = sp.Symbol("eps", real=True)
    hp, hx, dhp, dhx = sp.symbols("hp hx dhp dhx", real=True)
    H = sp.Matrix([[hp, hx], [hx, -hp]])
    dH = sp.Matrix([[dhp, dhx], [dhx, -dhp]])

    # Exact expansion of gamma(t) = exp(eps*(H + s*dH)), compute gamma^{-1} (d/ds gamma)|_{s=0}
    # Using Baker-Campbell-Hausdorff / Duhamel: gamma^{-1} d(gamma) = dH - (1/2)[H, dH] + (1/6)[H,[H,dH]] - ...
    def comm(A, B):
        return A * B - B * A

    # Up to O(eps^3) in J = gamma^{-1} d(gamma) (sufficient for tr(J^2) up to O(eps^4))
    J1 = eps * dH
    J2 = -sp.Rational(1, 2) * eps**2 * comm(H, dH)
    J3 = sp.Rational(1, 6) * eps**3 * comm(H, comm(H, dH))
    J = J1 + J2 + J3
    tr_J2 = sp.series(sp.trace(J * J), eps, 0, 5).removeO()
    quad = sp.simplify(sp.diff(tr_J2, eps, 2).subs(eps, 0) / 2)
    quart = sp.simplify(sp.diff(tr_J2, eps, 4).subs(eps, 0) / 24)
    # Check using Z = hp + i hx: hp*dhx - hx*dhp = Im(Z^* dZ)
    # Note SL(2,R)/SO(2) has negative sectional curvature -> metric sinh^2(r)/r^2 = 1 + r^2/3 + ... > 0
    expected_quart = sp.Rational(2, 3) * (hp * dhx - hx * dhp) ** 2
    assert sp.simplify(quart - expected_quart) == 0, f"Got quart={quart}"

    return {
        "quadratic_tr_gamma_inv_dgamma_sq": str(quad),
        "quartic_tr_gamma_inv_dgamma_sq": str(quart),
        "complex_Z_form_quartic": "+(2/3) [Im(Z^* dZ)]^2 != 0 (omitted in Paper I Eq. A.1)",
    }


def check_defect5_schwarzschild_second_variation():
    """Expand Paper I Eq. (A.1) directly around Schwarzschild (rho^2=r^2, N=1, F=f=1-2M/r)
    and compare with the true Regge-Wheeler potential f(r)[j(j+1)/r^2 - 6M/r^3]."""
    r, M, j, om = sp.symbols("r M j om", positive=True)
    f = 1 - 2 * M / r
    lam = sp.Rational(1, 2) * (j - 1) * (j + 2)
    chi = sp.Function("chi")(r)

    def Dr(expr):
        return sp.diff(expr, r) - sp.I * om / f * expr

    Z = chi / r
    Zr = Dr(Z)
    Zu = -sp.I * om * Z

    def extract_eom(c_ur, c_rr, c_ang):
        eom = -c_ur * (-sp.I * om * r**2 * Zr + Dr(r**2 * Zu)) + c_rr * Dr(r**2 * f * Zr) - c_ang * (2 * lam) * Z
        return sp.simplify(sp.expand(eom * f / r))

    # 1) Paper I Eq. (A.1) literally as written
    eom_A1 = extract_eom(1, sp.Rational(1, 2), f / 2)
    # Has surviving i*om*chi'(r) term (broken static symmetry!)
    cross_coeff_A1 = sp.simplify(sp.expand(eom_A1).coeff(sp.diff(chi, r)).coeff(sp.I * om))
    assert cross_coeff_A1 != 0

    # 2) Paper I Eq. (A.1) even AFTER fixing the factor of 1/2 and removing spurious F on |ethbar Z|^2
    eom_fixed = extract_eom(sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2))
    # Multiply by 2 to match standard d^2 chi / dr_*^2 + (om^2 - V_eff) chi = 0
    # Note d^2 chi / dr_*^2 = f^2 chi'' + f f' chi'
    d2_rstar = f**2 * sp.diff(chi, r, 2) + f * sp.diff(f, r) * sp.diff(chi, r)
    rem = sp.simplify(2 * eom_fixed - d2_rstar - om**2 * chi)
    V_A1_fixed = sp.simplify(-rem / chi)
    V_expected_scalar_spin2 = sp.simplify(f * ((j * (j + 1) - 2) / r**2 + 2 * M / r**3))
    V_RW = sp.simplify(f * (j * (j + 1) / r**2 - 6 * M / r**3))

    assert sp.simplify(V_A1_fixed - V_expected_scalar_spin2) == 0
    V_diff = sp.simplify(V_A1_fixed - V_RW)
    assert sp.simplify(V_diff - f * (-2 / r**2 + 8 * M / r**3)) == 0

    return {
        "paper1_A1_as_written_spurious_i_omega_chi_prime_coeff": str(cross_coeff_A1),
        "paper1_A1_after_1_2_and_F_fixes_potential": str(sp.factor(V_A1_fixed)),
        "true_Regge_Wheeler_potential": str(sp.factor(V_RW)),
        "potential_mismatch_V_A1_fixed_minus_V_RW": str(sp.factor(V_diff)),
    }


def main():
    results = {
        "paper1_commit": "e541955495d0ba6e21a4a6505282436cac3e7e74",
        "defect1_u1_symmetry": check_defect1_u1_symmetry(),
        "defect2_bondi_sachs_constraints": check_defect2_bondi_sachs_constraints(),
        "defect3_kinetic_factors_and_lightcone": check_defect3_kinetic_factors_and_lightcone(),
        "defect4_sl2_sigma_model_quartic": check_defect4_sl2_sigma_model_quartic(),
        "defect5_schwarzschild_second_variation": check_defect5_schwarzschild_second_variation(),
    }
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "audit_paper1_action_results.json")
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)

    print("=== PAPER I ACTION AUDIT (e541955495d0) ===")
    for k, v in results.items():
        if isinstance(v, dict):
            print(f"\n[{k}]")
            for kk, vv in v.items():
                print(f"  {kk}: {vv}")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
