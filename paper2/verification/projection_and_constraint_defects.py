#!/usr/bin/env python3
"""Verification of Paper II Section 4 diagnostics (Propositions 2 and 3) on BOTH:
  (A) The exact finite-matrix examples from Paper II (Eqs. 17-20, 35) for full self-contained reproducibility;
  (B) Actual gravitational operators from the finite-state causal-diamond model:
      1. Proposition 2 (Resolvent / virtual-state self-energy Sigma_P(z)):
         Interacting 1 <-> 2 graviton Hamiltonian with Paper I's exact cost operator
         Q_hat = sum omega_{k,j} n_{k,j} (omega_{k,j} = k + (j-1)(j+2)/(2k)), exact Wigner-6j
         cubic vertices, and budget projector P_B = 1_{[0,B]}(Q_hat) across B = 6..24!
         Quantifies the discarded virtual-state self-energy ||Sigma_P(z)|| and energy shift
         |E_exact - E_PHP| as a function of B, verifying the resolvent bound ||PHQ|| ||QHP|| / d.
      2. Proposition 3 (Constraint & Lorentz-boost compression defect):
         Evaluates [P K_i P, P K_j P] - P [K_i, K_j] P = -P K_i Q K_j P + P K_j Q K_i P
         directly on Paper I's spin-2 Gel'fand-Naimark Lorentz boost generators (L_i, K_i)
         under BOTH:
           (i) single-graviton angular-shell projection j <= j_max(B) (showing exact boundary-shell
               localization on j = j_max and exact zero in the interior j <= j_max - 1), AND
          (ii) two-graviton occupation-budget projection omega_1 + omega_2 <= B (demonstrating
               that total-budget coupling leaks the boost-compression defect into LOW-j single-particle
               modes whenever a spectator graviton saturates the remaining budget!).

Outputs: verification/projection_defects_results.json
"""

from __future__ import annotations

import json
import math
import os
import sys
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fuzzy_scaling_resolution import exact_wigner_3j_float, exact_wigner_6j_float


# =============================================================================
# 1. Exact Symbolic Verification of Paper II Draft Diagnostics (Eqs. 13, 17-20, 35)
# =============================================================================
def verify_draft_exact_diagnostics() -> dict:
    # Prop 1 scalar example: S = (1/2) Delta c^2 + kappa alpha c q^2 -> -kappa^2 alpha^2 q^4 / (2 Delta)
    kap, alpha, Delta = sp.symbols("kap alpha Delta", positive=True)
    q, c = sp.symbols("q c", real=True)
    S_aux = sp.Rational(1, 2) * Delta * c**2 + kap * alpha * c * q**2
    c_sol = sp.solve(sp.diff(S_aux, c), c)[0]
    S_red = sp.simplify(S_aux.subs(c, c_sol))
    assert sp.simplify(S_red - (-kap**2 * alpha**2 * q**4 / (2 * Delta))) == 0

    # Prop 2 2x2 example (Eqs. 17-18)
    g, D, z = sp.symbols("g D z", positive=True)
    H2 = sp.Matrix([[0, g], [g, D]])
    ev = sorted(H2.eigenvals().keys(), key=lambda e: float(e.subs({g: 0.1, D: 1.0})))
    E_minus = sp.simplify(ev[0])
    E_ser = sp.series(E_minus, g, 0, 7).removeO()
    assert sp.simplify(E_ser - (-g**2 / D + g**4 / D**3 - 2 * g**6 / D**5)) == 0
    # Numerical resolvent identity check at random Hermitian 6x6 matrix
    rng = np.random.default_rng(11)
    A = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
    H6 = 0.5 * (A + A.conj().T)
    P6 = np.diag([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    Q6 = np.eye(6) - P6
    z_num = 1.7 + 0.4j
    lhs_res = (P6 @ np.linalg.inv(z_num * np.eye(6) - H6) @ P6)[:3, :3]
    PHP = H6[:3, :3]
    PHQ = H6[:3, 3:]
    QHP = H6[3:, :3]
    QHQ = H6[3:, 3:]
    Sigma = PHQ @ np.linalg.inv(z_num * np.eye(3) - QHQ) @ QHP
    rhs_res = np.linalg.inv(z_num * np.eye(3) - PHP - Sigma)
    res_id_err = float(np.linalg.norm(lhs_res - rhs_res))
    assert res_id_err < 1e-14

    # Prop 3 3x3 example (Eqs. 19-20, 35)
    C1 = sp.Matrix([[0, 1, 1], [1, 0, 2], [1, 2, 0]])
    C2 = C1**2
    P3 = sp.diag(1, 1, 0)
    Q3 = sp.eye(3) - P3
    assert C1 * C2 - C2 * C1 == sp.zeros(3, 3)
    PC1P = (P3 * C1 * P3)[:2, :2]
    PC2P = (P3 * C2 * P3)[:2, :2]
    comm_comp = PC1P * PC2P - PC2P * PC1P
    assert comm_comp == sp.Matrix([[0, 3], [-3, 0]])
    # Check Eq. (19) identity explicitly
    lhs_19 = (P3 * C1 * P3) * (P3 * C2 * P3) - (P3 * C2 * P3) * (P3 * C1 * P3) - P3 * (C1 * C2 - C2 * C1) * P3
    rhs_19 = -P3 * C1 * Q3 * C2 * P3 + P3 * C2 * Q3 * C1 * P3
    assert lhs_19 == rhs_19

    return {
        "prop1_scalar_S_red": str(S_red),
        "prop2_E_minus_series": str(E_ser),
        "prop2_6x6_resolvent_identity_err": res_id_err,
        "prop3_commutator_full_zero": True,
        "prop3_compressed_commutator_2x2": [[int(comm_comp[i, j]) for j in range(2)] for i in range(2)],
        "prop3_eq19_identity_exact": True,
    }


# =============================================================================
# 2. Proposition 2 on an Interacting Bondi-Sachs / Wigner-6j Occupation-Budget Sector
# =============================================================================
def evaluate_prop2_occupation_budget_resolvent() -> dict:
    """Build an interacting 1 <-> 2 graviton Fock Hamiltonian H = H_0 + g * V_3^{(N)}
    with Paper I's exact cost operator omega_{k,j} = k + (j-1)(j+2)/(2k) (spin-2) and
    omega_{k,j}^{(0)} = k + j(j+1)/(2k) (scalar constraint/focusing mode), exact Wigner-6j
    cubic couplings at N = 32, and extended Fock cutoff B_ext = 26.
    Sweep the physical occupation budget B in {6, 8, 10, 12, 14, 16, 18, 20, 22, 24} and measure:
      - Exact lowest physical eigenvalue E_exact(H) vs compressed eigenvalue E_PHP(P_B H P_B)
      - Exact resolvent self-energy ||Sigma_P(z)|| at z = E_0 = omega_{1,2} = 3.0 (in the 1-graviton sector)
      - Paper II Eq. (16) upper bound ||PHQ|| ||QHP|| / dist(z, spec QHQ)
      - Exact verification of Prop 2 identity P(z-H)^{-1}P = (z - PHP - Sigma_P(z))^{-1}."""
    N = 32
    B_ext = 26.0
    g_coup = 0.08  # weak gravitational coupling (kappa / R ~ 1/sqrt(s))

    def om_gr(j, k):
        return k + 0.5 * (j - 1) * (j + 2) / k

    def om_sc(j, k):
        return k + 0.5 * j * (j + 1) / k

    # Precompute exact Wigner-6j cubic vertex C^{(N)}(j1, 0; 2, 0; j2, 0) for scalar j1 + graviton j2 <-> graviton (2,0)
    two_jL = N + 1
    two_jR = N - 3
    c6_cache = {}
    for j1 in range(2, 7):
        for j2 in range(2, 7):
            if abs(j1 - j2) <= 2 <= j1 + j2 and ((j1 + j2 + 2) % 2 == 0):
                w3 = exact_wigner_3j_float(j1, 2, j2, 0, 0, 0)
                if abs(w3) > 1e-12:
                    s6m = math.sqrt(N - 2) * exact_wigner_6j_float(2 * j1, 4, 2 * j2, two_jL, two_jR, two_jR)
                    s6p = math.sqrt(N + 2) * exact_wigner_6j_float(2 * j1, 4, 2 * j2, two_jR, two_jL, two_jL)
                    pref = math.sqrt((2 * j1 + 1) * 5 * (2 * j2 + 1)) * w3
                    val = 0.5 * pref * (s6m + s6p)  # parity even since j1+j2+2 even
                    if abs(val) > 1e-12:
                        c6_cache[(j1, j2)] = abs(val)

    # Enumerate Fock states: 1-graviton probe |Z(j=2, k=1)> (cost 3.0) + all coupled 2-particle states
    # |F(j1, k1), Z(j2, k2)> with total cost <= B_ext!
    states = [("in", 2, 1, None, None, om_gr(2, 1))]
    for (j1, j2), cval in sorted(c6_cache.items()):
        for k1 in range(1, 20):
            e1 = om_sc(j1, k1)
            if e1 >= B_ext:
                continue
            for k2 in range(1, 20):
                e2 = om_gr(j2, k2)
                etot = e1 + e2
                if etot <= B_ext + 1e-12:
                    states.append(("pair", j1, k1, j2, k2, etot))

    n_ext = len(states)
    H_ext = np.zeros((n_ext, n_ext), dtype=float)
    for i, st in enumerate(states):
        H_ext[i, i] = st[5]
    # Bondi-Sachs 2-derivative cubic vertex scales as g * C_6j * sqrt(k1 * k2) / (omega_1 * omega_2)^{1/2}
    for i in range(1, n_ext):
        _, j1, k1, j2, k2, _ = states[i]
        e1 = om_sc(j1, k1)
        e2 = om_gr(j2, k2)
        v_elem = g_coup * c6_cache[(j1, j2)] * math.sqrt(float(k1 * k2)) / math.sqrt(e1 * e2)
        H_ext[i, 0] = v_elem
        H_ext[0, i] = v_elem

    # Exact unprojected ground-state energy and full resolvent computed ONCE!
    evals_full = np.linalg.eigvalsh(H_ext)
    E_exact_1g = float(evals_full[0])
    z_probe = 3.0
    z_c = 3.0 + 0.2j
    G_full_c = np.linalg.inv(z_c * np.eye(n_ext) - H_ext)

    # Sweep physical budget B
    budget_scan = []
    for B in (6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0):
        in_B = np.array([st[5] <= B + 1e-12 for st in states], dtype=bool)
        idx_P = np.where(in_B)[0]
        idx_Q = np.where(~in_B)[0]
        nP, nQ = len(idx_P), len(idx_Q)

        H_PP = H_ext[np.ix_(idx_P, idx_P)]
        H_PQ = H_ext[np.ix_(idx_P, idx_Q)]
        H_QP = H_ext[np.ix_(idx_Q, idx_P)]
        diag_Q = np.diag(H_ext)[idx_Q]  # H_QQ is strictly diagonal (1 <-> 2 coupling)!

        E_php = float(np.linalg.eigvalsh(H_PP)[0])
        energy_shift_missed = abs(E_exact_1g - E_php)

        gap_d = float(np.min(np.abs(diag_Q - z_probe)))
        Sigma_P = H_PQ @ ((1.0 / (z_probe - diag_Q))[:, None] * H_QP)
        sigma_norm = float(np.linalg.norm(Sigma_P, 2))
        sigma_00 = float(abs(Sigma_P[0, 0]))
        norm_PHQ = float(np.linalg.norm(H_PQ, 2))
        eq16_bound = (norm_PHQ**2) / gap_d

        # Verify exact resolvent identity (Eq. 14) at z_c = 3.0 + 0.2j
        lhs = G_full_c[np.ix_(idx_P, idx_P)]
        Sig_c = H_PQ @ ((1.0 / (z_c - diag_Q))[:, None] * H_QP)
        rhs = np.linalg.inv(z_c * np.eye(nP) - H_PP - Sig_c)
        res_err = float(np.linalg.norm(lhs - rhs))
        assert res_err < 1e-12
        assert sigma_norm <= eq16_bound + 1e-12

        budget_scan.append({
            "B": B,
            "n_P": nP,
            "n_Q": nQ,
            "gap_d": gap_d,
            "norm_PHQ": norm_PHQ,
            "sigma_norm": sigma_norm,
            "sigma_00": sigma_00,
            "eq16_upper_bound": eq16_bound,
            "E_PHP": E_php,
            "E_exact_ext": E_exact_1g,
            "missed_energy_shift": energy_shift_missed,
            "resolvent_identity_err": res_err,
        })

    return {
        "N": N,
        "B_ext": B_ext,
        "n_ext_states": n_ext,
        "g_coupling": g_coup,
        "bare_1g_energy": z_probe,
        "exact_ext_1g_energy": E_exact_1g,
        "budget_scan": budget_scan,
    }


# =============================================================================
# 3. Proposition 3 on Paper I's Spin-2 Lorentz Boost Generators (1-body & 2-body Budget!)
# =============================================================================
def build_spin2_gelfand_naimark_generators(j_max: int):
    """Exact unitary Gel'fand-Naimark principal-series spin-2 (helicity s=2, c=0)
    SO(3,1) generators (L_x, L_y, L_z, K_x, K_y, K_z) truncated to angular shells
    2 <= j <= j_max (matching Paper I Sec. III.D!).

    Recall for unitary massless/principal-series spin s=2 (Casimirs C1 = L^2 - K^2 = s^2 - 1 + ...;
    for Naimark's (k_0=2, c=0) representation on S^2 spin-2 harmonics):
      L_z |j,m> = m |j,m>,  L_+- |j,m> = sqrt((j -+ m)(j +- m + 1)) |j, m +- 1>
      K_z |j,m> = C_j * sqrt(j^2 - m^2) |j-1, m> + A_j * m |j, m> - C_{j+1} * sqrt((j+1)^2 - m^2) |j+1, m>
    with k_0 = 2, c = 0 (unitary principal series with s=2):
      A_j = i * k_0 * c / (j(j+1)) = 0
      C_j = (i / j) * sqrt((j^2 - k_0^2)(j^2 - c^2) / (4 j^2 - 1))  [wait: Hermitian K_z => real or pure imag!]
    Check Hermiticity of K_z: <j-1, m | K_z | j, m> = C_j sqrt(j^2 - m^2),
    <j, m | K_z | j-1, m> = -C_j sqrt(j^2 - m^2) if C_j is purely imaginary!
    With C_j = +i * (1/j) * sqrt((j^2 - 4) * j^2 / (4 j^2 - 1)) = i * sqrt((j^2 - 4) / (4 j^2 - 1)),
    then (-C_j) = +C_j^* => K_z is EXACTLY HERMITIAN!
    Check [K_z, K_x] = -i L_y using K_x = i [L_y, K_z] and K_y = -i [L_x, K_z]!
    """
    basis = [(j, m) for j in range(2, j_max + 1) for m in range(-j, j + 1)]
    idx_of = {jm: i for i, jm in enumerate(basis)}
    dim = len(basis)

    Lx = np.zeros((dim, dim), dtype=complex)
    Ly = np.zeros((dim, dim), dtype=complex)
    Lz = np.zeros((dim, dim), dtype=complex)
    Kz = np.zeros((dim, dim), dtype=complex)

    for (j, m), col in idx_of.items():
        Lz[col, col] = m
        if m + 1 <= j:
            row = idx_of[(j, m + 1)]
            lp = math.sqrt((j - m) * (j + m + 1))
            Lx[row, col] += 0.5 * lp
            Ly[row, col] += -0.5j * lp
        if m - 1 >= -j:
            row = idx_of[(j, m - 1)]
            lm = math.sqrt((j + m) * (j - m + 1))
            Lx[row, col] += 0.5 * lm
            Ly[row, col] += 0.5j * lm

        # K_z coupling j <-> j-1 and j <-> j+1 (for k_0 = 2, c = 0):
        # C_j = i * sqrt((j^2 - 4) / (4 j^2 - 1))  (vanishes at j=2 -> no j=1 leakage!)
        if j - 1 >= 2 and abs(m) <= j - 1:
            row = idx_of[(j - 1, m)]
            Cj = 1j * math.sqrt((j * j - 4.0) / (4.0 * j * j - 1.0))
            Kz[row, col] += Cj * math.sqrt(j * j - m * m)
        if j + 1 <= j_max and abs(m) <= j + 1:
            row = idx_of[(j + 1, m)]
            Cjp1 = 1j * math.sqrt(((j + 1) ** 2 - 4.0) / (4.0 * (j + 1) ** 2 - 1.0))
            Kz[row, col] += -Cjp1 * math.sqrt((j + 1) ** 2 - m * m)

    # Exact SU(2) vector construction: [L_a, K_b] = i eps_{abc} K_c => K_x = i [L_y, K_z], K_y = -i [L_x, K_z]
    Kx = 1j * (Ly @ Kz - Kz @ Ly)
    Ky = -1j * (Lx @ Kz - Kz @ Lx)

    # Check Hermiticity of all 6 generators
    for G in (Lx, Ly, Lz, Kx, Ky, Kz):
        assert np.linalg.norm(G - G.conj().T) < 1e-13

    return basis, idx_of, (Lx, Ly, Lz), (Kx, Ky, Kz)


def evaluate_prop3_lorentz_compression_defect() -> dict:
    """Evaluate Paper II Proposition 3 (Eq. 19):
      D_{xy} := [P K_x P, P K_y P] - P [K_x, K_y] P = -P K_x Q K_y P + P K_y Q K_x P
    on BOTH:
      (1) Single-particle spin-2 space with angular truncation P_{j_cut} (j <= j_cut inside parent j_max = 10):
          shows exact boundary-shell localization at j = j_cut (interior j <= j_cut - 1 has ZERO defect!).
      (2) Two-particle Fock states |(k1, j1, m1), (k2, j2, m2)> under Paper I's TOTAL occupation budget
          omega(k1,j1) + omega(k2,j2) <= B (Eq. 3-4 of Paper II!):
          proves that total-budget coupling leaks a NONZERO O(1) boost-algebra compression defect
          directly onto LOW-ANGULAR-MOMENTUM j1 = 2 gravitons whenever a spectator graviton (k2, j2)
          places the pair within one boost step of the budget ceiling B!"""
    # --- Part (1): Single-particle shell-resolved compression defect ---
    j_parent = 10
    basis1, idx1, (Lx, Ly, Lz), (Kx, Ky, Kz) = build_spin2_gelfand_naimark_generators(j_parent)
    dim1 = len(basis1)

    # First verify parent [K_x, K_y] + i L_z == 0 on all interior shells j <= j_parent - 1!
    full_boost_err = (Kx @ Ky - Ky @ Kx) + 1j * Lz
    int_mask = np.array([j <= j_parent - 1 for (j, m) in basis1])
    parent_interior_err = float(np.linalg.norm(full_boost_err[:, int_mask]))
    assert parent_interior_err < 1e-11, f"Parent boost closure err={parent_interior_err}"

    shell_scan = []
    for j_cut in (3, 4, 5, 6, 7, 8):
        in_P = np.array([j <= j_cut for (j, m) in basis1], dtype=float)
        P = np.diag(in_P)
        Q = np.eye(dim1) - P

        PKxP = P @ Kx @ P
        PKyP = P @ Ky @ P
        D_lhs = (PKxP @ PKyP - PKyP @ PKxP) - P @ (Kx @ Ky - Ky @ Kx) @ P
        D_rhs = -P @ Kx @ Q @ Ky @ P + P @ Ky @ Q @ Kx @ P
        id_err = float(np.linalg.norm(D_lhs - D_rhs))
        assert id_err < 1e-11

        # Resolve norm per angular shell j = 2..j_cut!
        per_shell = {}
        for j_val in range(2, j_cut + 1):
            m_idx = [i for i, (j, m) in enumerate(basis1) if j == j_val]
            per_shell[str(j_val)] = float(np.linalg.norm(D_lhs[np.ix_(m_idx, m_idx)], 2))
        shell_scan.append({
            "j_cut": j_cut,
            "total_defect_op_norm": float(np.linalg.norm(D_lhs, 2)),
            "interior_max_defect": max(per_shell[str(jv)] for jv in range(2, j_cut)),
            "boundary_shell_defect": per_shell[str(j_cut)],
            "per_shell_norm": per_shell,
            "eq19_identity_err": id_err,
        })

    # --- Part (2): Two-particle total-budget leakage into LOW-j1 (j1 = 2) modes! ---
    # Consider distinguishable two-mode sector (probe graviton 1 at k1=1, j1 in 2..5; spectator graviton 2 at k2=1, j2 in 2..5)
    # where total cost is omega(1, j1) + omega(1, j2) <= B.
    # For k = 1, omega(1, j) = 1 + (j-1)(j+2)/2:
    #   j=2 -> 3.0,  j=3 -> 6.0,  j=4 -> 10.0,  j=5 -> 15.0!
    # Total boost generator on the tensor product is K_a^{tot} = K_a (x) I + I (x) K_a.
    j_sub = 5
    b_sub, idx_sub, _, (Kx_s, Ky_s, Kz_s) = build_spin2_gelfand_naimark_generators(j_sub)
    d_sub = len(b_sub)
    I_sub = np.eye(d_sub, dtype=complex)
    Kx_2 = np.kron(Kx_s, I_sub) + np.kron(I_sub, Kx_s)
    Ky_2 = np.kron(Ky_s, I_sub) + np.kron(I_sub, Ky_s)

    def om_k1(j):
        return 1.0 + 0.5 * (j - 1) * (j + 2)

    costs_2 = np.array([om_k1(j1) + om_k1(j2) for (j1, _) in b_sub for (j2, _) in b_sub])
    j1_arr = np.array([j1 for (j1, _) in b_sub for (j2, _) in b_sub])
    j2_arr = np.array([j2 for (j1, _) in b_sub for (j2, _) in b_sub])

    # Evaluate defect on LOW-j probe states (j1 = 2, j2 = 3 -> cost 3 + 6 = 9.0) vs budget B!
    # Notice: at B = 9.0 or B = 10.0, state |(j1=2), (j2=3)> IS IN P_B (cost 9 <= B),
    # YET boosting particle 1 from j1=2 -> j1=3 raises cost to 6 + 6 = 12 > B (OUTSIDE P_B!),
    # AND boosting particle 2 from j2=3 -> j2=4 raises cost to 3 + 10 = 13 > B (OUTSIDE P_B!)!
    # Therefore EVEN THE LOWEST GRAVITON MODE j1 = 2 suffers an O(1) boost-algebra compression defect!
    two_body_scan = []
    for B_2 in (9.0, 12.0, 13.0, 16.0):
        in_P2 = (costs_2 <= B_2 + 1e-12).astype(float)
        P2 = np.diag(in_P2)
        Q2 = np.eye(len(costs_2)) - P2
        D2_lhs = (P2 @ Kx_2 @ P2) @ (P2 @ Ky_2 @ P2) - (P2 @ Ky_2 @ P2) @ (P2 @ Kx_2 @ P2) - P2 @ (Kx_2 @ Ky_2 - Ky_2 @ Kx_2) @ P2
        D2_rhs = -P2 @ Kx_2 @ Q2 @ Ky_2 @ P2 + P2 @ Ky_2 @ Q2 @ Kx_2 @ P2
        assert np.linalg.norm(D2_lhs - D2_rhs) < 1e-12

        # Defect restricted specifically to the (j1=2, j2=3) subspace (lowest physical graviton j1=2!)
        low_probe_idx = np.where((j1_arr == 2) & (j2_arr == 3))[0]
        defect_j1_2_j2_3 = float(np.linalg.norm(D2_lhs[np.ix_(low_probe_idx, low_probe_idx)], 2))
        # Also (j1=2, j2=2) subspace (cost 6.0)
        low_22_idx = np.where((j1_arr == 2) & (j2_arr == 2))[0]
        defect_j1_2_j2_2 = float(np.linalg.norm(D2_lhs[np.ix_(low_22_idx, low_22_idx)], 2))
        two_body_scan.append({
            "B": B_2,
            "n_P2": int(np.sum(in_P2)),
            "defect_on_j1_2_j2_2": defect_j1_2_j2_2,
            "defect_on_j1_2_j2_3": defect_j1_2_j2_3,
            "max_defect_all_P2": float(np.linalg.norm(D2_lhs, 2)),
        })

    return {
        "parent_j_max": j_parent,
        "parent_interior_boost_closure_err": parent_interior_err,
        "single_body_shell_scan": shell_scan,
        "two_body_total_budget_leakage": two_body_scan,
    }


def main():
    diag_res = verify_draft_exact_diagnostics()
    prop2_res = evaluate_prop2_occupation_budget_resolvent()
    prop3_res = evaluate_prop3_lorentz_compression_defect()

    results = {
        "draft_exact_diagnostics": diag_res,
        "prop2_occupation_budget_resolvent": prop2_res,
        "prop3_lorentz_compression_defect": prop3_res,
    }
    out_path = os.path.join(HERE, "projection_defects_results.json")
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)

    print("=== PROJECTION & CONSTRAINT-COMPRESSION DIAGNOSTICS ===")
    print(f"[1] Draft exact diagnostics (Eqs. 13, 17-20, 35): ALL VERIFIED "
          f"(6x6 resolvent err = {diag_res['prop2_6x6_resolvent_identity_err']:.2e})")

    print(f"\n[2] Prop 2 on interacting Wigner-6j occupation-budget sector (N={prop2_res['N']}, n_ext={prop2_res['n_ext_states']}):")
    for r in prop2_res["budget_scan"]:
        if r["B"] in (6.0, 10.0, 14.0, 18.0, 22.0, 24.0):
            print(
                f"    B={r['B']:4.1f} (nP={r['n_P']:3d}, d={r['gap_d']:4.1f}): "
                f"||Sigma_P||={r['sigma_norm']:.3e} <= bound={r['eq16_upper_bound']:.3e} | "
                f"|E_exact - E_PHP|={r['missed_energy_shift']:.3e}"
            )

    print(f"\n[3] Prop 3 on spin-2 Lorentz boost generators (parent interior closure err={prop3_res['parent_interior_boost_closure_err']:.2e}):")
    for r in prop3_res["single_body_shell_scan"]:
        print(f"    1-body j_cut={r['j_cut']}: interior max defect={r['interior_max_defect']:.2e}, "
              f"boundary shell j={r['j_cut']} defect={r['boundary_shell_defect']:.4f}")
    print("    2-body total-budget leakage onto LOWEST graviton mode j1=2:")
    for r in prop3_res["two_body_total_budget_leakage"]:
        print(f"    B={r['B']:4.1f} (nP={r['n_P2']:3d}): defect(j1=2,j2=2)={r['defect_on_j1_2_j2_2']:.4f}, "
              f"defect(j1=2,j2=3)={r['defect_on_j1_2_j2_3']:.4f}")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
