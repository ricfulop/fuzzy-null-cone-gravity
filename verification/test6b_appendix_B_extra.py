#!/usr/bin/env python3
"""
Test 6B (Part 2): Audit of the remaining Appendix B angular & propagation equations:
  - Eq. (B.4a): Damour-Navier-Stokes angular constraint d_v Pi_A on general metric (B.1) with y-dependence
  - Eq. (B.5c): Propagation of transverse connection d_v mu_c
  - Eq. (B.7): Divergence identity reducing bulk Hamiltonian (3.47) to pure surface charge (3.54)
  - Eq. (B.8f): Inaffinity kappa_m = m^a nabla_a log N - d_v(s_a s^a)/(2 N) on displaced null boundary!
"""

import json
import os
import sympy as sp


def audit_b7_divergence_identity():
    """
    Check eq. (B.7) algebraically from eq. (B.5a):
      (B.5a): d_v(rho^2 vartheta) = -N rho^2 [ (1/2) R[q] - div(varpi) - varpi_A varpi^A - Lambda ]
      varpi_A = w_A - (1/2) nabla_A f,   f = -log N  =>  nabla_A f = -(nabla_A N)/N
      Check whether 2 d_v(rho^2 vartheta) - 2 d_A(N rho^2 varpi^A) EQUALS
        N rho^2 [ -R[q] - (1/2) nabla_A f nabla^A f + 2 w_A w^A + 2 Lambda ]!
    """
    y = sp.Symbol("y", real=True)
    N = sp.Function("N")(y)
    A = sp.Function("A")(y)      # rho^2
    w = sp.Function("w")(y)      # w_y (1D angular test with q_{yy} = 1 for clean divergence test, or general q(y)!)
    q = sp.Function("q")(y)      # q_{yy} -> metric factor; sqrt(q) = A in 2D, wait:
    # In 2D angular coordinates with sqrt(det q) = rho^2 = A:
    # For a 2D vector field X^A, the covariant divergence is:
    #   nabla_A X^A = (1/A) d_A(A X^A)  =>  A nabla_A(varpi^A) = d_A(A varpi^A)!
    # Let us compute:
    #   2 d_v(A vartheta) - 2 d_A(N A varpi^A)
    #   Substitute d_v(A vartheta) from (B.5a):
    #     = -N A R[q] + 2 N A nabla_A(varpi^A) + 2 N A varpi_A varpi^A + 2 N A Lambda
    #       - 2 (d_A N) A varpi^A - 2 N d_A(A varpi^A)
    #   Since d_A(A varpi^A) = A nabla_A(varpi^A), the +2 N A nabla_A varpi^A and -2 N d_A(A varpi^A) CANCEL!
    #   Now substitute d_A N = -N nabla_A f (since f = -log N => nabla_A f = -(d_A N)/N):
    #     Remaining terms divided by (N A):
    #     = -R[q] + 2 varpi_A varpi^A + 2 Lambda + 2 (nabla_A f) varpi^A!
    #   Now substitute varpi_A = w_A - (1/2) nabla_A f:
    df, w_val, Rq, Lam = sp.symbols("df w Rq Lam", real=True)
    varpi = w_val - sp.Rational(1, 2) * df
    remaining = sp.expand(-Rq + 2 * varpi**2 + 2 * Lam + 2 * df * varpi)
    paper_B7 = sp.expand(-Rq - sp.Rational(1, 2) * df**2 + 2 * w_val**2 + 2 * Lam)
    res_B7 = sp.simplify(remaining - paper_B7)
    assert res_B7 == 0, f"B.7 mismatch: {res_B7}"
    return dict(status="PASS", eq="B.7", remaining=str(remaining), paper_B7=str(paper_B7), residual=str(res_B7))


def audit_b8f_displaced_inaffinity():
    """
    Check eq. (B.8f) directly from 4D Christoffel symbols of metric (B.1) with N(t,r), F(t,r), b(t)
    (and also angular dependence b(t,y) at N(t,r,y), V=0, q=rho(t,r,y)^2 I_2!):
      s_a = N nabla_a(v - b), m^a = -g^{ab} s_b
      m^b nabla_b m^a = kappa_m m^a
      Check eq. (B.8f): kappa_m == m^a nabla_a log N - d_v(s_a s^a) / (2 N)  on s_a s^a = 0!
    """
    t, r, y, zc = sp.symbols("t r y zc", real=True)
    coords = [t, r, y, zc]
    F = sp.Function("F")(t, r, y)
    N = sp.Function("N")(t, r, y)
    rho = sp.Function("rho")(t, r, y)
    b = sp.Function("b")(t, y)

    q22 = rho**2
    g = sp.Matrix([
        [-F, -N, 0, 0],
        [-N, 0, 0, 0],
        [0, 0, q22, 0],
        [0, 0, 0, q22],
    ])
    ginv = sp.simplify(g.inv())
    Gamma = [
        [
            [
                sp.simplify(
                    sp.Rational(1, 2)
                    * sum(
                        ginv[a, d]
                        * (
                            sp.diff(g[d, c], coords[b])
                            + sp.diff(g[d, b], coords[c])
                            - sp.diff(g[b, c], coords[d])
                        )
                        for d in range(4)
                    )
                )
                for c in range(4)
            ]
            for b in range(4)
        ]
        for a in range(4)
    ]

    # s_a = N nabla_a(r - b(t, y)) = (-N b_t, N, -N b_y, 0)
    s_cov = sp.Matrix([-N * sp.diff(b, t), N, -N * sp.diff(b, y), 0])
    m = sp.simplify(-ginv * s_cov)
    s_sq = sp.simplify((s_cov.T * ginv * s_cov)[0])  # = F + 2 N b_t + N^2 b_y^2 / rho^2

    # Compute (m^b nabla_b m^0) / m^0 (since m^0 = 1, kappa_m = m^b nabla_b m^0 = m^b Gamma^0_{bc} m^c!)
    acc_0 = sp.simplify(
        sum(m[b_i] * (sp.diff(m[0], coords[b_i]) + sum(Gamma[0][b_i][c_i] * m[c_i] for c_i in range(4))) for b_i in range(4))
    )
    # Paper eq. (B.8f): kappa_m = m^a d_a log N - d_r(s_a s^a) / (2 N)  (with d_r at fixed b, b_t, b_y)
    m_dlogN = sp.simplify(sum(m[a] * sp.diff(N, coords[a]) / N for a in range(4)))
    dr_ssq = sp.diff(s_sq, r)
    kappa_paper = sp.simplify(m_dlogN - dr_ssq / (2 * N))

    # Compare ON the null surface s_a s^a = 0 (eliminate b_t = -(F + N^2 b_y^2 / rho^2)/(2 N)):
    sub_null = {sp.diff(b, t): -(F + N**2 * sp.diff(b, y) ** 2 / rho**2) / (2 * N)}
    res_B8f = sp.simplify((acc_0 - kappa_paper).subs(sub_null))

    # Also check component 2 (y-component): m^b nabla_b m^y - kappa_m m^y on the null surface!
    # Wait: s_a s^a = 0 holds ALONG the hypersurface r = b(t, y), so its tangential derivatives vanish:
    # (d_t + b_t d_r)(s_a s^a) = 0 and (d_y + b_y d_r)(s_a s^a) = 0!
    # We check acc_0 directly because m^0 = 1 has NO derivatives of b_t or b_y!
    return dict(
        status="PASS" if res_B8f == 0 else "FAIL",
        eq="B.8f",
        residual_on_null_surface=str(res_B8f),
    )


def audit_b4a_and_b5c():
    """
    Check eq. (B.4a) (Damour-Navier-Stokes R_{ry} = 0) and eq. (B.5c) (propagation of mu_c)
    on a metric (B.1) with N(t,r,y), V^y(t,r), rho(r) = r * rhoS(t), F(t,r)!
    """
    t, r, y, zc = sp.symbols("t r y zc", real=True)
    coords = [t, r, y, zc]
    # Part 1: check (B.4a) with y-dependent N(r,y), sh(r,y), V(r) (fixed t, F=0, rho = r*rhoS(y)? Take rho=r for simplicity, so Theta_l = 2/r, R[q]=0)
    N = sp.Function("N")(r, y)
    V = sp.Function("V")(r, y)
    sh = sp.Function("sh")(r, y)
    q22 = r**2 * sp.exp(sh)
    q33 = r**2 * sp.exp(-sh)
    g = sp.Matrix([
        [q22 * V**2, -N, q22 * V, 0],
        [-N, 0, 0, 0],
        [q22 * V, 0, q22, 0],
        [0, 0, 0, q33],
    ])
    ginv = sp.simplify(g.inv())
    Gamma = [
        [
            [
                sp.simplify(
                    sp.Rational(1, 2)
                    * sum(
                        ginv[a, d]
                        * (
                            sp.diff(g[d, c], coords[b])
                            + sp.diff(g[d, b], coords[c])
                            - sp.diff(g[b, c], coords[d])
                        )
                        for d in range(4)
                    )
                )
                for c in range(4)
            ]
            for b in range(4)
        ]
        for a in range(4)
    ]

    def Ric(b, d):
        val = 0
        for a in range(4):
            val += sp.diff(Gamma[a][b][d], coords[a]) - sp.diff(Gamma[a][a][d], coords[b])
            for e in range(4):
                val += Gamma[a][a][e] * Gamma[e][b][d] - Gamma[a][b][e] * Gamma[e][a][d]
        return sp.simplify(val)

    R_ry = Ric(1, 2)
    # Paper eq. (B.4a):
    #   Pi_y = -2 rho^2 w_y = (rho^2 / N) q_{yy} d_r V
    #   d_r Pi_y == rho^2 [ 2 nabla_B sigma^B_y + nabla_y(d_r f) - nabla_y Theta_l - Theta_l nabla_y f ]
    # Compute 2D Christoffel of q = diag(q22, q33) in (y, zc):
    #   Gamma2^y_{yy} = (1/2) sh_y, Gamma2^z_{yz} = -(1/2) sh_y, Gamma2^y_{zz} = (1/2) e^{-2 sh} sh_y
    #   sigma^y_y = (1/2) sh_r, sigma^z_z = -(1/2) sh_r (diagonal!)
    #   nabla_B sigma^B_y = d_y sigma^y_y + Gamma2^y_{yy} sigma^y_y + Gamma2^z_{yz} sigma^y_y - Gamma2^y_{yy} sigma^y_y - Gamma2^z_{zy} sigma^z_z
    #                     = (1/2) sh_{ry} + (-(1/2) sh_y)(1/2 sh_r) - (-(1/2) sh_y)(-1/2 sh_r) = (1/2) sh_{ry} - (1/2) sh_y sh_r?
    #   Wait: let us compute nabla_B sigma^B_y = (1/sqrt(q)) d_B(sqrt(q) sigma^B_y) - (1/2)(d_y q_{CD}) sigma^{CD}!
    #   Since sqrt(q) = r^2 (independent of y!), (1/sqrt(q)) d_y(r^2 sigma^y_y) = (1/2) sh_{ry}.
    #   And -(1/2)(d_y q22 sigma^{22} + d_y q33 sigma^{33}) = -(1/2)(q22 sh_y * (1/2) sh_r / q22 + q33 (-sh_y) * (-1/2) sh_r / q33)
    #     = -(1/2)(sh_y sh_r / 2 + sh_y sh_r / 2) = -(1/2) sh_y sh_r!
    #   Wait: check the sign! For diagonal q = diag(q22, q33), nabla_B T^B_A = (1/sqrt(q)) d_B(sqrt(q) T^B_A) - (1/2)(d_A q_{CD}) T^{CD}:
    #     Here T^2_2 = +1/2 sh_r, T^{22} = +(1/2) sh_r / q22, d_y q22 = +q22 sh_y -> prod = +(1/2) sh_r sh_y.
    #     T^3_3 = -1/2 sh_r, T^{33} = -(1/2) sh_r / q33, d_y q33 = -q33 sh_y -> prod = +(1/2) sh_r sh_y.
    #     Sum = sh_r sh_y => -(1/2) sh_r sh_y.
    #   So 2 nabla_B sigma^B_y = sh_{ry} - sh_y sh_r!
    #   Wait: let's double check with explicit 2D covariant derivative formula!
    q2d = sp.diag(q22, q33)
    q2inv = q2d.inv()
    c2 = [y, zc]
    G2 = [
        [
            [
                sp.simplify(
                    sp.Rational(1, 2)
                    * sum(q2inv[a, d] * (sp.diff(q2d[d, c], c2[b]) + sp.diff(q2d[d, b], c2[c]) - sp.diff(q2d[b, c], c2[d])) for d in range(2))
                )
                for c in range(2)
            ]
            for b in range(2)
        ]
        for a in range(2)
    ]
    sig_mixed = sp.diag(sp.Rational(1, 2) * sp.diff(sh, r), -sp.Rational(1, 2) * sp.diff(sh, r))
    div_sig_y = 0
    for B in range(2):
        div_sig_y += sp.diff(sig_mixed[B, 0], c2[B])
        for C in range(2):
            div_sig_y += G2[B][B][C] * sig_mixed[C, 0] - G2[C][B][0] * sig_mixed[B, C]
    div_sig_y = sp.simplify(div_sig_y)

    f = -sp.log(N)
    Theta_l = 2 / r
    Pi_y = (r**2 / N) * q22 * sp.diff(V, r)
    lhs_B4a = sp.diff(Pi_y, r)
    rhs_B4a = r**2 * (2 * div_sig_y + sp.diff(f, r, y) - sp.diff(Theta_l, y) - Theta_l * sp.diff(f, y))
    # Compare (lhs_B4a - rhs_B4a) / r^2 with 2 * R_{ry} (or -2 * R_{ry})!
    res_B4a_plus = sp.simplify((lhs_B4a - rhs_B4a) / r**2 + 2 * R_ry)
    res_B4a_minus = sp.simplify((lhs_B4a - rhs_B4a) / r**2 - 2 * R_ry)

    # Part 2: Check (B.5c) on (t, r)-dependent metric (N(t,r), F(t,r), V(t,r), rho(t,r), sh(t,r), no y-dep so R[q]=0, df_A=0)
    F_tr = sp.Function("F")(t, r)
    N_tr = sp.Function("N")(t, r)
    V_tr = sp.Function("V")(t, r)
    rho_tr = sp.Function("rho")(t, r)
    sh_tr = sp.Function("sh")(t, r)
    q22_tr = rho_tr**2 * sp.exp(sh_tr)
    q33_tr = rho_tr**2 * sp.exp(-sh_tr)
    g_tr = sp.Matrix([
        [-F_tr + q22_tr * V_tr**2, -N_tr, q22_tr * V_tr, 0],
        [-N_tr, 0, 0, 0],
        [q22_tr * V_tr, 0, q22_tr, 0],
        [0, 0, 0, q33_tr],
    ])
    g_tr_inv = sp.simplify(g_tr.inv())
    G_tr = [
        [
            [
                sp.simplify(
                    sp.Rational(1, 2)
                    * sum(
                        g_tr_inv[a, d]
                        * (sp.diff(g_tr[d, c], coords[b]) + sp.diff(g_tr[d, b], coords[c]) - sp.diff(g_tr[b, c], coords[d]))
                        for d in range(4)
                    )
                )
                for c in range(4)
            ]
            for b in range(4)
        ]
        for a in range(4)
    ]

    def Ric_tr(b, d):
        val = 0
        for a in range(4):
            val += sp.diff(G_tr[a][b][d], coords[a]) - sp.diff(G_tr[a][a][d], coords[b])
            for e in range(4):
                val += G_tr[a][a][e] * G_tr[e][b][d] - G_tr[a][b][e] * G_tr[e][a][d]
        return sp.simplify(val)

    # Compute R_{ab} l^a n^b on this metric!
    # l = (0, 1, 0, 0), n = (1/N, -F/(2 N^2), -V/N, 0)
    l_v = [0, 1, 0, 0]
    n_v = [1 / N_tr, -F_tr / (2 * N_tr**2), -V_tr / N_tr, 0]
    R_ln = sp.simplify(sum(Ric_tr(a, b) * l_v[a] * n_v[b] for a in range(4) for b in range(4)))
    # Also compute (1/2) q^{AB} R_{AB}
    half_qR = sp.simplify(sp.Rational(1, 2) * (Ric_tr(2, 2) / q22_tr + Ric_tr(3, 3) / q33_tr))
    # Recall on Rab = Lambda g_ab: R_{ab} l^a n^b = -Lambda AND (1/2) q^{AB} R_{AB} = +Lambda,
    # so R_{ln} + (1/2) q^{AB} R_{AB} = 0 INDEPENDENTLY OF LAMBDA!
    # Check eq. (B.5c):
    #   d_r mu_c == sigma^{AB} eta_{AB} - (1/2) Theta_l vartheta + 3 N w_A w^A   (when R[q]=0, nabla_A f = 0)
    # Build quantities from (B.2a-g):
    f_tr = -sp.log(N_tr)
    mu_c = sp.diff(f_tr, t) + sp.diff(F_tr, r) / (2 * N_tr)
    Th_l_tr = 2 * sp.diff(rho_tr, r) / rho_tr
    Th_t_tr = 2 * sp.diff(rho_tr, t) / rho_tr
    vartheta = Th_t_tr - (F_tr / (2 * N_tr)) * Th_l_tr
    # sigma_{22} = (1/2) q22 sh_r, sigma_{33} = -(1/2) q33 sh_r
    # eta_{22} = (1/2) q22 (sh_t - (F/(2N)) sh_r), eta_{33} = -(1/2) q33 (sh_t - (F/(2N)) sh_r)
    sig_eta = sp.Rational(1, 2) * sp.diff(sh_tr, r) * (sp.diff(sh_tr, t) - (F_tr / (2 * N_tr)) * sp.diff(sh_tr, r))
    w_sq_tr = (q22_tr * sp.diff(V_tr, r) ** 2) / (4 * N_tr**2)

    lhs_B5c = sp.diff(mu_c, r)
    rhs_B5c = sig_eta - sp.Rational(1, 2) * Th_l_tr * vartheta + 3 * N_tr * w_sq_tr
    # Check relation between (lhs_B5c - rhs_B5c) and N_tr * (R_ln + half_qR)!
    # Wait: let us test BOTH on-shell (R_ln = -Lam, half_qR = +Lam => R_ln + half_qR = 0)
    diff_B5c_vs_Ric = sp.simplify((lhs_B5c - rhs_B5c) - N_tr * (R_ln + half_qR))
    diff_B5c_vs_negRic = sp.simplify((lhs_B5c - rhs_B5c) + N_tr * (R_ln + half_qR))
    # Also check with R_rr term in case rho_{rr} != 0 (remember in (B.5c) rho = v * rho_S(t), so rho_{rr} = 0!):
    rhoS_t = sp.Function("rhoS")(t)
    diff_B5c_area_plus = sp.simplify(diff_B5c_vs_Ric.subs(rho_tr, r * rhoS_t).doit())
    diff_B5c_area_minus = sp.simplify(diff_B5c_vs_negRic.subs(rho_tr, r * rhoS_t).doit())

    return dict(
        eq_B4a=dict(
            status="PASS" if (res_B4a_plus == 0 or res_B4a_minus == 0) else "CHECK",
            lhs_minus_rhs_over_r2_plus_2Rry=str(res_B4a_plus),
            lhs_minus_rhs_over_r2_minus_2Rry=str(res_B4a_minus),
        ),
        eq_B5c=dict(
            status="PASS" if (diff_B5c_area_plus == 0 or diff_B5c_area_minus == 0) else "CHECK",
            diff_area_coord_minus_N_Rln_halfqR=str(diff_B5c_area_plus),
            diff_area_coord_plus_N_Rln_halfqR=str(diff_B5c_area_minus),
        ),
    )


if __name__ == "__main__":
    res = dict(
        b7_divergence=audit_b7_divergence_identity(),
        b8f_inaffinity=audit_b8f_displaced_inaffinity(),
        b4a_and_b5c=audit_b4a_and_b5c(),
    )
    out_path = os.path.join(os.path.dirname(__file__), "test6b_appendix_B_extra_results.json")
    with open(out_path, "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))
    print(f"Saved results to {out_path}")
