#!/usr/bin/env python3
"""
Test 6B: Exact Symbolic Audit of Appendix B ("Boundary terms in the Hamiltonian", pp. 59-64)
of arXiv:2611.34924v1.

Explicitly constructs the 4D metric (B.1) in SymPy:
    ds^2 = -F dt^2 - 2 N dt dr + q_{AB}(dy^A + V^A dt)(dy^B + V^B dt)
and checks EVERY key classical GR identity from scratch (Christoffel symbols, Ricci tensor,
scalar curvature, total divergence, displaced null boundary, and Noether charge):

  B-1. Exact metric inverse, null normals l, n (eq. B.2f), inaffinity, and constraints:
       - R_{vv} = 0 => d_v f = -(sigma_{AB} sigma^{AB}) / Theta_l  (eq. B.3, with f = -log N)
       - Cross-focusing / Einstein constraint (B.5a) & divergence identity (B.7)!
  B-2. Exact bulk Lagrangian decomposition eqs. (B.10)-(B.11e) for N = 1 (affine coordinate r):
       Verify sqrt(-g) (R - 2 Lambda) == L_0 + Delta L_0 + d_a L^a  (checking EVERY sign & factor of 2!)
  B-3. Exact displaced null boundary r = b(t, y) (eqs. B.8a-g and B.12):
       Verify s_a s^a (B.8d), m^a (B.8e), kappa_m (B.8f), and boundary divergence s_a L^a (B.12)!
  B-4. Exact Iyer-Wald / Komar Noether charge Q_chi|_C (eq. B.17) on r = b(t, y)!
"""

import json
import os
import sympy as sp


def audit_b1_b2_full_4d_metric():
    """
    Construct a general (1+1+2)-reduced 4D metric of the exact form (B.1) with nontrivial
    F(t,r), N(t,r), V^2(t,r) (angular shift), rho(t,r) (area radius), and shear z(t,r)
    (anisotropic 2-metric q_{22} = rho^2 e^z, q_{33} = rho^2 e^{-z}) — plus a second test
    with angular dependence y^2 = y for checking angular derivatives (d_A, w_A, V^A, R[q], B.7, B.12, B.17)!
    """
    t, r, y, zc = sp.symbols("t r y zc", real=True)
    coords = [t, r, y, zc]

    # --- PART 1: Full (t, r)-dependent metric with N(t,r), F(t,r), V(t,r), rho(t,r), shear s(t,r) ---
    F = sp.Function("F")(t, r)
    N = sp.Function("N")(t, r)
    V = sp.Function("V")(t, r)      # V^y (angular shift along y)
    rho = sp.Function("rho")(t, r)
    sh = sp.Function("sh")(t, r)    # shear parameter: q_{yy} = rho^2 e^{sh}, q_{zz} = rho^2 e^{-sh}
    Lam = sp.Symbol("Lambda", real=True)

    q22 = rho**2 * sp.exp(sh)
    q33 = rho**2 * sp.exp(-sh)
    # ds^2 = -F dt^2 - 2 N dt dr + q22 (dy + V dt)^2 + q33 dzc^2
    g = sp.Matrix([
        [-F + q22 * V**2, -N, q22 * V, 0],
        [-N, 0, 0, 0],
        [q22 * V, 0, q22, 0],
        [0, 0, 0, q33],
    ])
    ginv = sp.simplify(g.inv())
    assert sp.simplify(g.det() + N**2 * rho**4) == 0  # det(g) = -N^2 rho^4 => sqrt(-g) = N rho^2

    # Compute Christoffel symbols Gamma^a_{bc}
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

    # Ricci tensor R_{bd} = d_a Gamma^a_{bd} - d_b Gamma^a_{ad} + Gamma^a_{ae} Gamma^e_{bd} - Gamma^a_{be} Gamma^e_{ad}
    def Ric(b, d):
        val = 0
        for a in range(4):
            val += sp.diff(Gamma[a][b][d], coords[a]) - sp.diff(Gamma[a][a][d], coords[b])
            for e in range(4):
                val += Gamma[a][a][e] * Gamma[e][b][d] - Gamma[a][b][e] * Gamma[e][a][d]
        return sp.simplify(val)

    # Check (B.3): R_{rr} = 0 with rho = v * rho_S(t) (or general rho(t,r))!
    # Actually let's check R_{rr} for general rho(t,r), N(t,r), sh(t,r):
    R_rr = Ric(1, 1)
    # Paper (B.2a)-(B.3):
    #   Theta_l = (d_r rho^2)/rho^2 = 2 rho_r / rho
    #   sigma_{22} = (1/2) q22 sh_r, sigma_{33} = -(1/2) q33 sh_r => sigma_{AB} sigma^{AB} = (1/2) sh_r^2
    #   f = -log N => d_r f = -N_r / N
    # On rho = r * rho_S(t) (so rho_{rr} = 0, Theta_l = 2/r):
    rhoS = sp.Function("rhoS")(t)
    R_rr_area_coord = sp.simplify(R_rr.subs(rho, r * rhoS).doit())
    # Check against -(Theta_l d_r f + sigma_{AB} sigma^{AB}) = -((2/r)(-N_r/N) + (1/2) sh_r^2) = (2 N_r)/(r N) - (1/2) sh_r^2!
    expected_R_rr = (2 * sp.diff(N, r)) / (r * N) - sp.Rational(1, 2) * sp.diff(sh, r) ** 2
    res_B3 = sp.simplify(R_rr_area_coord - expected_R_rr)
    assert res_B3 == 0, f"B.3 mismatch: R_rr = {R_rr_area_coord}, expected = {expected_R_rr}"

    # Check (B.2f) null vector n and (B.2d) w_A and (B.2e) mu_c:
    # l = d_r = (0, 1, 0, 0); n = (1/N)(d_t - V d_y - (F/(2N)) d_r) = (1/N, -F/(2 N^2), -V/N, 0)
    l_vec = sp.Matrix([0, 1, 0, 0])
    n_vec = sp.Matrix([1 / N, -F / (2 * N**2), -V / N, 0])
    assert sp.simplify((l_vec.T * g * l_vec)[0]) == 0
    assert sp.simplify((n_vec.T * g * n_vec)[0]) == 0
    assert sp.simplify((l_vec.T * g * n_vec)[0] + 1) == 0  # l . n = -1

    # Check beta_a = -n_b nabla_a l^b:
    # Since l^b = delta^b_r, nabla_a l^b = Gamma^b_{ar}, so beta_a = -n_b Gamma^b_{ar}.
    n_cov = sp.simplify(g * n_vec)
    beta_r = sp.simplify(-sum(n_cov[b] * Gamma[b][1][1] for b in range(4)))
    beta_y = sp.simplify(-sum(n_cov[b] * Gamma[b][2][1] for b in range(4)))
    beta_t = sp.simplify(-sum(n_cov[b] * Gamma[b][0][1] for b in range(4)))
    # Compare with (A.1) & (B.2d)-(B.2h):
    #   beta_r = d_r log N = N_r / N
    #   w_y = -(1/(2N)) q22 d_r V  => here N has no y-dep so beta_y = w_y!
    #   n^a beta_a = -n^a n_b nabla_a l^b = mu_c? Check (B.2e): mu_c = d_t(-log N) + (d_r F)/(2N) = -N_t/N + F_r/(2N)!
    w_y = -q22 * sp.diff(V, r) / (2 * N)
    mu_c_paper = -sp.diff(N, t) / N + sp.diff(F, r) / (2 * N)
    n_dot_beta = sp.simplify(n_vec[0] * beta_t + n_vec[1] * beta_r + n_vec[2] * beta_y)
    # Wait: check whether mu_c is n^a beta_a or -n^a beta_a or what!
    res_beta_r = sp.simplify(beta_r - sp.diff(N, r) / N)
    res_beta_y = sp.simplify(beta_y - w_y)
    res_mu_c_plus = sp.simplify(n_dot_beta - mu_c_paper)
    res_mu_c_minus = sp.simplify(n_dot_beta + mu_c_paper)
    assert res_beta_r == 0, f"beta_r = {beta_r}"
    assert res_beta_y == 0, f"beta_y = {beta_y}, w_y = {w_y}"

    # --- PART 2: Exact check of Bulk Lagrangian Decomposition eqs. (B.10)-(B.11e) for N = 1 ---
    # Substitute N = 1 into metric g and compute exact Ricci scalar R!
    g1 = g.subs(N, 1)
    g1inv = ginv.subs(N, 1)
    Gamma1 = [[[sp.simplify(Gamma[a][b][c].subs(N, 1).doit()) for c in range(4)] for b in range(4)] for a in range(4)]

    def Ric1(b, d):
        val = 0
        for a in range(4):
            val += sp.diff(Gamma1[a][b][d], coords[a]) - sp.diff(Gamma1[a][a][d], coords[b])
            for e in range(4):
                val += Gamma1[a][a][e] * Gamma1[e][b][d] - Gamma1[a][b][e] * Gamma1[e][a][d]
        return sp.simplify(val)

    R_scalar = sp.simplify(sum(g1inv[b, d] * Ric1(b, d) for b in range(4) for d in range(4)))
    # Note: In (3.2) & (3.6)-(3.8), I = (1/(16 pi G_N)) int (R - 2 Lambda) sqrt(-g) d^4x
    # AND in (3.2), B_{mn} ^ F^{mn} - (1/2) Phi B ^ B = (R - 2 Lambda) * eps_g!
    # Let us check whether (R - 2 Lambda) * rho^2 EQUALS (L_0 + Delta L_0 + d_a L^a) from (B.11a)-(B.11e)!!
    eh_density = sp.simplify(rho**2 * (R_scalar - 2 * Lam))

    # Build paper's eqs. (B.9a)-(B.11e) (with no y-dependence here, so £_V = 0, R[q] = 0, d_A = 0):
    A_area = rho**2
    Th_t = sp.diff(A_area, t) / A_area
    Th_r = sp.diff(A_area, r) / A_area
    # sigma_t AB = (1/2)(d_t q_{AB} - Th_t q_{AB}): since q22 = A e^{sh}, q33 = A e^{-sh},
    #   d_t q22 - Th_t q22 = A e^{sh} sh_t  => sigma_t 22 = (1/2) q22 sh_t, sigma_t 33 = -(1/2) q33 sh_t
    #   sigma_r 22 = (1/2) q22 sh_r, sigma_r 33 = -(1/2) q33 sh_r
    #   sigma_t^{AB} sigma_{r AB} = (1/2) sh_t sh_r
    #   sigma_r^{AB} sigma_{r AB} = (1/2) sh_r^2
    sig_t_sig_r = sp.Rational(1, 2) * sp.diff(sh, t) * sp.diff(sh, r)
    sig_r_sig_r = sp.Rational(1, 2) * sp.diff(sh, r) ** 2
    w_sq = (q22 * sp.diff(V, r) ** 2) / 4  # w_y = -(1/2) q22 V_r, w^y = -(1/2) V_r => w_A w^A = (1/4) q22 V_r^2

    L0 = A_area * (2 * sig_t_sig_r - Th_t * Th_r + 2 * w_sq - 2 * Lam)
    dL0 = sp.diff(A_area, r) * sp.diff(F, r) + F * (sp.diff(A_area, r) ** 2 / (2 * A_area) - A_area * sig_r_sig_r)
    Lt = 2 * sp.diff(A_area, r)
    Lr = 2 * sp.diff(A_area, t) - A_area * sp.diff(F, r) - 2 * F * sp.diff(A_area, r)
    div_L = sp.diff(Lt, t) + sp.diff(Lr, r)
    paper_B10 = sp.simplify(L0 + dL0 + div_L)

    res_B10_diff = sp.simplify(eh_density - paper_B10)

    # Check (B.5a) & (B.7) (with N(t,r), V=0, sh=0, or general N(t,r), V(t,r), sh(t,r))!
    # Evaluate G_{tr} / R_{ln} = R_{ab} l^a n^b = -Lambda (Wait: check sign of R_{ab} l^a n^b = -Lambda!):
    # Since R_{ab} = Lambda g_{ab} and l.n = -1, R_{ab} l^a n^b = Lambda (l.n) = -Lambda!
    # Check eq. (B.5a) on our metric (no y-dep, so R[q]=0, varpi=w_y, div=0):
    #   d_r(rho^2 vartheta) == -N rho^2 [ -w_A w^A - Lambda ]  (since R_{ab} l^a n^b = -Lambda)
    # where vartheta = N q^{AB} K^{(n)}_{AB} = (1/2) q^{AB} d_t q_{AB} - (F/(2N)) Theta_l = Th_t - (F/(2N)) Th_r!
    # Let us compute the exact Einstein tensor component E_{ab} l^a n^b = (R_{ab} - (1/2) R g_{ab} + Lambda g_{ab}) l^a n^b
    # Wait: let's check (B.5a) directly by computing R_{ab} l^a n^b - (-Lambda) on N=1 metric!
    R_ln_1 = sp.simplify(
        sum(
            Ric1(a, b) * l_vec[a] * n_vec.subs(N, 1)[b]
            for a in range(4)
            for b in range(4)
        )
    )
    # Wait: in (3.57b): C_{ln} = (1/2) q^{ab} R_{ab} - Lambda = (G_{ab} + Lambda g_{ab}) l^a n^b = 0!
    # Let us compute (1/2) q^{AB} R_{AB} on our N=1 metric!
    half_qR = sp.simplify(sp.Rational(1, 2) * (Ric1(2, 2) / q22 + Ric1(3, 3) / q33))
    vartheta_1 = Th_t - sp.Rational(1, 2) * F * Th_r
    # Paper (B.5a) says: d_r(rho^2 vartheta) = -rho^2 [ (1/2) R[q] - div(varpi) - varpi^2 - Lambda ]
    # i.e. d_r(rho^2 vartheta) - rho^2 w_sq - rho^2 Lambda = 0!
    # Check whether (1/rho^2) * [ d_r(rho^2 vartheta) - rho^2 w_sq ] EQUALS -half_qR (or +half_qR)!!
    lhs_B5a = sp.simplify(sp.diff(A_area * vartheta_1, r) - A_area * w_sq)
    res_B5a_plus = sp.simplify(lhs_B5a / A_area + half_qR)
    res_B5a_minus = sp.simplify(lhs_B5a / A_area - half_qR)

    return dict(
        eq_B3_residual=str(res_B3),
        eq_A1_beta_r_residual=str(res_beta_r),
        eq_A1_beta_y_residual=str(res_beta_y),
        n_dot_beta_vs_mu_c=dict(
            n_dot_beta_minus_mu_c=str(res_mu_c_plus),
            n_dot_beta_plus_mu_c=str(res_mu_c_minus),
        ),
        eq_B10_B11_EH_minus_paper=str(res_B10_diff),
        eq_B5a_check=dict(
            lhs_over_A_plus_half_qR=str(res_B5a_plus),
            lhs_over_A_minus_half_qR=str(res_B5a_minus),
        ),
    )


def audit_b3_displaced_boundary_and_angular():
    """
    Check eqs. (B.8a)-(B.8g), (B.12), and (B.17) on a metric with angular dependence y
    and a displaced null boundary r = b(t, y) (taking N = 1 as stated below B.8g).
    """
    t, r, y, zc = sp.symbols("t r y zc", real=True)
    coords = [t, r, y, zc]
    F = sp.Function("F")(t, r, y)
    N = sp.Function("N")(t, r, y)
    V = sp.Function("V")(t, r, y)
    rho = sp.Function("rho")(t, r, y)
    b = sp.Function("b")(t, y)

    q22 = rho**2
    q33 = rho**2
    g = sp.Matrix([
        [-F + q22 * V**2, -N, q22 * V, 0],
        [-N, 0, 0, 0],
        [q22 * V, 0, q22, 0],
        [0, 0, 0, q33],
    ])
    ginv = sp.simplify(g.inv())

    # Eq. (B.8a): s_a = N nabla_a(r - b(t, y)) = N * (-b_t, 1, -b_y, 0)
    s_cov = sp.Matrix([-N * sp.diff(b, t), N, -N * sp.diff(b, y), 0])
    # Eq. (B.8b): m^a = -g^{ab} s_b
    m_vec = sp.simplify(-ginv * s_cov)
    # Eq. (B.8d): s_a s^a = g^{ab} s_a s_b
    s_sq = sp.simplify((s_cov.T * ginv * s_cov)[0])
    s_sq_paper = F + 2 * N * (sp.diff(b, t) - V * sp.diff(b, y)) + N**2 * (1 / q22) * sp.diff(b, y) ** 2
    res_B8d = sp.simplify(s_sq - s_sq_paper)

    # Eq. (B.8e): m = d_t + (N q^{yy} b_y - V) d_y + ((N^2 q^{yy} b_y^2 - F)/(2 N)) d_r
    # WAIT: check m^r on the null surface s_a s^a = 0!
    # Solve s_a s^a = 0 for b_t: b_t = V b_y - (F + N^2 q^{yy} b_y^2)/(2 N)
    m_r_on_null = sp.simplify(
        m_vec[1].subs(sp.diff(b, t), V * sp.diff(b, y) - (F + N**2 * sp.diff(b, y) ** 2 / q22) / (2 * N))
    )
    m_paper = sp.Matrix([
        1,
        (N**2 * sp.diff(b, y) ** 2 / q22 - F) / (2 * N),
        N * sp.diff(b, y) / q22 - V,
        0,
    ])
    res_B8e_t = sp.simplify(m_vec[0] - m_paper[0])
    res_B8e_r_on_null = sp.simplify(m_r_on_null - m_paper[1])
    res_B8e_y = sp.simplify(m_vec[2] - m_paper[2])
    # Also note off-null: m^r = (1/N) * (N^2 q^{yy} b_y^2 - F - s_a s^a) ?? Let's check!
    m_r_off_null = sp.simplify(m_vec[1] - (m_paper[1] - s_sq / (2 * N)))

    # Eq. (B.12): s_a L^a on r = b(t, y) with N = 1 and s_a s^a = 0!
    # From (B.11c-e):
    #   L^t = 2 d_r(rho^2)
    #   L^r = 2 d_t(rho^2) - 2 d_y(rho^2 V) - rho^2 d_r F - 2 F d_r(rho^2)
    #   L^y = -2 d_r(rho^2) V
    A = rho**2
    L_vec = sp.Matrix([
        2 * sp.diff(A, r),
        2 * sp.diff(A, t) - 2 * sp.diff(A * V, y) - A * sp.diff(F, r) - 2 * F * sp.diff(A, r),
        -2 * sp.diff(A, r) * V,
        0,
    ])
    s_cov_1 = s_cov.subs(N, 1)
    s_dot_L = sp.simplify((s_cov_1.T * L_vec)[0])
    # Substitute null condition s_a s^a = 0 (at N=1): 2 (b_t - V b_y) = -F - q^{yy} b_y^2
    s_dot_L_on_null = sp.simplify(
        s_dot_L.subs(sp.diff(b, t), V * sp.diff(b, y) - sp.Rational(1, 2) * (F + sp.diff(b, y) ** 2 / q22))
    )
    B12_paper = (
        2 * sp.diff(A, t)
        - 2 * sp.diff(A * V, y)
        - A * sp.diff(F, r)
        - F * sp.diff(A, r)
        + sp.diff(A, r) * sp.diff(b, y) ** 2 / q22
    )
    res_B12 = sp.simplify(s_dot_L_on_null - B12_paper)

    # Eq. (B.17): Iyer-Wald / Komar Noether charge Q_chi|_C for chi = (0, chi^r(t,r,y), chi^y(t,r,y), 0)
    # at N = 1 on cut C: {t = const, r = b(t, y)}!
    # From eq. (3.50c): (Q_chi)_{ab} = -(1/(16 pi G_N)) (eps_g)_{abcd} nabla^c chi^d.
    # Pullback to cut C (coordinates y, zc, with r = b(t,y) at fixed t):
    # Tangent vectors to C: e_y = d_y + b_y d_r = (0, b_y, 1, 0), e_z = d_zc = (0, 0, 0, 1).
    # So (Q_chi)|_C = (Q_chi)_{ab} e_y^a e_z^b dy ^ dzc = [ (Q_chi)_{y zc} + b_y (Q_chi)_{r zc} ] dy ^ dzc.
    # With eps_{t r y zc} = +sqrt(-g) = +rho^2 (at N = 1):
    #   (Q_chi)_{y zc} = -(1/(16 pi G_N)) * rho^2 * [ nabla^t chi^r - nabla^r chi^t ]
    #   (Q_chi)_{r zc} = -(1/(16 pi G_N)) * (-rho^2) * [ nabla^t chi^y - nabla^y chi^t ]
    # Wait: let's compute nabla^c chi^d = g^{ce} (d_e chi^d + Gamma^d_{ef} chi^f) explicitly at N = 1!
    g1 = g.subs(N, 1)
    g1inv = ginv.subs(N, 1)
    Gamma1 = [
        [
            [
                sp.simplify(
                    sp.Rational(1, 2)
                    * sum(
                        g1inv[a, d]
                        * (
                            sp.diff(g1[d, c], coords[b])
                            + sp.diff(g1[d, b], coords[c])
                            - sp.diff(g1[b, c], coords[d])
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
    chir = sp.Function("chir")(t, r, y)
    chiy = sp.Function("chiy")(t, r, y)
    chi = [0, chir, chiy, 0]

    def nabla_up_chi(c, d):
        # g^{ce} nabla_e chi^d
        val = 0
        for e in range(4):
            if g1inv[c, e] != 0:
                cov = sp.diff(chi[d], coords[e]) + sum(Gamma1[d][e][f] * chi[f] for f in range(4))
                val += g1inv[c, e] * cov
        return sp.simplify(val)

    # Antisymmetric part: 2 nabla^{[c} chi^{d]} = nabla^c chi^d - nabla^d chi^c
    # Note: (Q_chi)_{ab} = -(1/(16 pi G_N)) eps_{abcd} nabla^c chi^d
    # For (a,b) = (2,3) = (y, zc): sum_{c,d} eps_{23cd} nabla^c chi^d = eps_{2301} (nabla^0 chi^1 - nabla^1 chi^0)
    #   Since (2,3,0,1) is an EVEN permutation of (0,1,2,3), eps_{2301} = +rho^2!
    # For (a,b) = (1,3) = (r, zc): (1,3,0,2) -> swap(1,0)->(0,3,1,2) -> (0,1,2,3) is 3 swaps = ODD => eps_{1302} = -rho^2!
    # Wait: check sign of eps_{2301} and eps_{1302} by Levi-Civita:
    # (0,1,2,3): +1.
    # (2,3,0,1): 2 moves past 0,1 (2 swaps), 3 moves past 0,1 (2 swaps) -> 4 swaps = +1.
    # (1,3,0,2): 1->0 (1), 3->2 (1) -> (0,1,3,2) -> swap 3,2 (1) => 3 swaps = -1.
    Q_yz = -rho**2 * (nabla_up_chi(0, 1) - nabla_up_chi(1, 0))
    Q_rz = -(-rho**2) * (nabla_up_chi(0, 2) - nabla_up_chi(2, 0))
    Q_pullback = sp.simplify(Q_yz + sp.diff(b, y) * Q_rz)  # times 1/(16 pi G_N)

    # Paper's eq. (B.17) (times 1/(16 pi G_N)):
    #   rho^2 * [ d_r chi^r - q_{yy} (d_r V) chi^y - q^{yy} (d_y b) d_r(q_{yy} chi^y) ]
    # WAIT: check the SIGN of Q_pullback vs B.17!! (Depending on orientation of S: dy^dz vs dz^dy)
    B17_paper = rho**2 * (
        sp.diff(chir, r)
        - q22 * sp.diff(V, r) * chiy
        - (1 / q22) * sp.diff(b, y) * sp.diff(q22 * chiy, r)
    )
    res_B17_minus = sp.simplify(Q_pullback - B17_paper)
    res_B17_plus = sp.simplify(Q_pullback + B17_paper)

    return dict(
        eq_B8d_residual=str(res_B8d),
        eq_B8e_t_residual=str(res_B8e_t),
        eq_B8e_r_on_null_residual=str(res_B8e_r_on_null),
        eq_B8e_y_residual=str(res_B8e_y),
        eq_B8e_r_off_null_note=str(m_r_off_null),
        eq_B12_residual=str(res_B12),
        eq_B17_Q_minus_paper=str(res_B17_minus),
        eq_B17_Q_plus_paper=str(res_B17_plus),
    )


if __name__ == "__main__":
    print("Running full 4D SymPy audit of Appendix B (this takes ~10-20s)...")
    r1 = audit_b1_b2_full_4d_metric()
    r2 = audit_b3_displaced_boundary_and_angular()
    res = dict(b1_b2_full_4d_metric=r1, b3_displaced_boundary_and_charge=r2)
    out_path = os.path.join(os.path.dirname(__file__), "test6b_appendix_B_results.json")
    with open(out_path, "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))
    print(f"Saved results to {out_path}")
