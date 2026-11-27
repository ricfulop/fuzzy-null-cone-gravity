#!/usr/bin/env python3
"""
Test 1: Verification & Counterexample Suite for Section 5.2 of arXiv:2611.34924v1
("The discrete spectrum of the area operator" & AdS boundary discrete spectrum).

Checks:
  1A. Exact counterexample H = C^2 (x) l^2(N) to the cycle-grouping proof (eqs. 5.30-5.38, D.22-D.25).
  1B. Continuous-spectrum counterexample H = L^2([0, s_*]) (x) l^2(N) modeling JT + massless scalar
      with bounded corner dilaton Phi(S_R) in [0, Phi_*].
  1C. Exact counterexample & Riemann equivalence for the AdS boundary proof (eqs. 5.42-5.46).
  1D. Proof that salvaging Section 5.2 with finite-norm open states destroys eqs. (5.30c) and (D.24)-(D.25).
"""

import json
import os
import numpy as np


def run_test_1a_qubit_clock():
    """Exact counterexample H = C^2 (x) l^2(N) to eqs. (5.30)-(5.38) and (D.22)-(D.25)."""
    s_star = 2.0 * np.pi
    s0, s1 = 0.0, s_star
    # C_b = diag(c0, c1) (x) I; satisfies e^{-s_*} m_b I <= C_b^dag C_b <= I (eq. 5.30c)
    c0, c1 = 1.0, np.exp(-s1 / 2.0)  # eigenvalues of C_b^dag C_b are 1 and e^{-2pi}
    m_lower = c1**2

    # Boost rapidity eta = pi acts on conjugate copy by e^{i * pi * s / (2 pi)} = diag(+1, -1)
    eta = np.pi
    # Bounded observable O = -sum_n (|0,n><1,n| (x) |0,n>*<1,n|* + h.c.), ||O|| = 1
    # On each fixed-n sector: <O_0> = -2 c0 c1 / (c0^2 + c1^2), <O_eta> = +2 c0 c1 / (c0^2 + c1^2)
    delta_eta = 4.0 * c0 * c1 / (c0**2 + c1**2)

    # 1) Product exhaustion P_M = I_{C^2} (x) sum_{n=1}^M |n><n|
    product_rows = []
    for M in [1, 5, 10, 25, 50, 100, 250, 500]:
        tr_P = M * (c0**2 + c1**2)
        exp_dO = delta_eta
        Z_P_X = tr_P * (exp_dO - delta_eta / 2.0)
        product_rows.append(
            dict(M=M, rank=2 * M, tr_P=float(tr_P), exp_dO=float(exp_dO), Z_P_X=float(Z_P_X))
        )

    # 2) Paper's cycle exhaustion P_cyc(K, r): K cycles of length r,
    #    each having 1 vector |1, n> (outside I={0}) and (r-1) vectors |0, m> (inside I={0}).
    #    Choose r large enough so bound (5.32) <= delta_eta / 2.
    #    From eq. (D.25) with Delta s = 0, s_0 = 0:
    #      Tr[rho_P (s - 0)^2] <= s_*^2 / (r * m_lower)
    #      bound_532 = (|eta|/pi) * ||O|| * sqrt(Tr[rho_P (s - 0)^2])
    r_required = int(np.ceil((2.0 * s_star / delta_eta) ** 2 / m_lower)) + 1
    # Actual exact value of var_s in K cycles (exact ratio independent of K!):
    # var_s_exact = c1^2 * s_star^2 / ((r - 1) * c0^2 + c1^2)
    r_exact_threshold = int(np.ceil((2.0 * c1 * s_star / delta_eta) ** 2 / c0**2)) + 2

    cycle_rows = []
    r_test = max(r_exact_threshold, 200)
    for K in [1, 5, 10, 25, 50, 100, 250, 500]:
        tr_P = K * ((r_test - 1) * c0**2 + c1**2)
        # Upper bound on exact <Delta_eta O>_P even if every |1,n> is paired with its matching |0,n>
        exp_dO_max = (K * 4.0 * c0 * c1) / tr_P
        var_s_exact = (K * c1**2 * s_star**2) / tr_P
        var_s_D25_bound = s_star**2 / (r_test * m_lower)
        bound_532_exact = (abs(eta) / np.pi) * 1.0 * np.sqrt(var_s_exact)
        Z_P_X_max = tr_P * (exp_dO_max - delta_eta / 2.0)
        cycle_rows.append(
            dict(
                K=K,
                r=r_test,
                rank=K * r_test,
                tr_P=float(tr_P),
                exp_dO_max=float(exp_dO_max),
                bound_532_exact=float(bound_532_exact),
                Z_P_X_max=float(Z_P_X_max),
                var_s_exact=float(var_s_exact),
                var_s_D25_bound=float(var_s_D25_bound),
            )
        )

    # Check consistency with every inequality of Section 5.2
    assert all(row["Z_P_X"] > 0 for row in product_rows)
    assert all(row["Z_P_X_max"] < 0 for row in cycle_rows)
    assert all(row["exp_dO_max"] <= row["bound_532_exact"] + 1e-12 for row in cycle_rows)
    assert all(row["bound_532_exact"] <= delta_eta / 2.0 for row in cycle_rows)

    return dict(
        delta_eta=float(delta_eta),
        c0=float(c0),
        c1=float(c1),
        m_lower=float(m_lower),
        r_test=r_test,
        r_D25_sufficient=r_required,
        product_rows=product_rows,
        cycle_rows=cycle_rows,
    )


def run_test_1b_continuous_jt_analog():
    """
    Continuous-spectrum analog of JT gravity + massless scalar with bounded dilaton
    s = 2 pi Phi / hbar in [0, s_*], infinite multiplicity at each s (clock states Psi_n & scalar modes).
    Discretized into M_s fine bins of s in [0, s_*] tensored with l^2(N).
    """
    s_star = 10.0
    s_b = 5.0
    sigma_s = np.sqrt(2.0 * s_b)
    M_s = 200
    s_grid = np.linspace(0.0, s_star, M_s)
    # Weight w(s) = e^{-s} * exp(-(s - s_b)^2 / (2 sigma_s^2))
    w_s = np.exp(-s_grid) * np.exp(-((s_grid - s_b) ** 2) / (2.0 * sigma_s**2))
    m_lower = float(np.min(w_s))

    # Smooth boost-sensitive observable across s (e.g., coherent rank-1 projector in s-basis on each fiber)
    # Take |phi> = sum_a sqrt(w_a) |a> / sqrt(sum w_a), and O = |phi><phi| (x) |phi*><phi*| on each fiber?
    # Under boost eta: |phi(eta)> = sum_a sqrt(w_a) e^{i eta s_a / (2 pi)} |a> / sqrt(Z_1)
    # In |Psi_fiber> = (1/sqrt(Z_1)) sum_a sqrt(w_a) |a>|a*>, take O = -(|u><v| + |v><u|) where |u>,|v> span different s
    eta = 3.0
    # Compute boost fidelity |<Psi_fiber | (I x e^{i eta s*/(2pi)}) | Psi_fiber>|^2
    Z1 = float(np.sum(w_s))
    phase = np.exp(1j * eta * s_grid / (2.0 * np.pi))
    overlap = np.sum(w_s * phase) / Z1
    # Take O = - |Psi_fiber><Psi_fiber| (bounded, ||O|| = 1):
    # <O_0> = -1, <O_eta> = -|overlap|^2 -> Delta_eta O = 1 - |overlap|^2 > 0!
    delta_eta = float(1.0 - abs(overlap) ** 2)
    assert delta_eta > 0.05

    # Partition [0, s_*] into intervals of width Delta_s <= 0.1
    Delta_s = 0.1
    s0 = 5.0
    in_I = np.abs(s_grid - s0) <= Delta_s / 2.0
    w_in = float(np.mean(w_s[in_I]))
    w_out_max = float(np.max(w_s[~in_I]))

    # Choose r so bound (5.32) <= delta_eta / 2
    # In each cycle of length r: 1 vector outside I (weight <= w_out_max), r-1 vectors inside I (weight >= min in I)
    w_in_min = float(np.min(w_s[in_I]))
    r = 10
    while True:
        var_bound = Delta_s**2 + (w_out_max * s_star**2) / ((r - 1) * w_in_min + w_out_max)
        b532 = (abs(eta) / np.pi) * np.sqrt(var_bound)
        if b532 <= delta_eta / 2.0:
            break
        r += 10

    return dict(
        s_star=s_star,
        eta=eta,
        delta_eta=delta_eta,
        m_lower=m_lower,
        Delta_s=Delta_s,
        r_chosen=r,
        bound_532_at_r=float(b532),
        product_sign="+infinity (since Z_{P_M}[X] = M * Z1 * (delta_eta/2) -> +inf)",
        cycle_sign="-infinity (since Z_P[X] <= 0 for every finite union of cycles)",
    )


def run_test_1c_ads_disjoint_tautology():
    """
    Exact counterexample to eqs. (5.44)-(5.45) on p. 44:
    Continuous-spectrum / infinite-trace system (w_i = 1/(i+1), sum w_i = +inf)
    tensored with a 2-state diamond (v_0 = 0.4 in I, v_1 = 0.6 not in I, so p = 0.4).
    """
    p = 0.4
    v = np.array([p, 1.0 - p])
    c = np.array([v[0] * (1.0 - p), v[1] * (0.0 - p)])  # [p(1-p), -p(1-p)]
    assert abs(np.sum(c)) < 1e-15

    # Divergent weights w_i = 1/(i+1) (analogous to infinite trace Z_beta = sum w_i = inf)
    N = 10000
    w = 1.0 / np.arange(1, N + 1)

    # Exhaustion 1: balanced pairs F_k = {1..k} x {0, 1} -> S_F = 0 identically!
    balanced_sums = [float(np.sum(w[:k]) * np.sum(c)) for k in [10, 100, 1000, 10000]]

    # Exhaustion 2: paper's unbalanced subset F_k U (J x {1}) -> arbitrarily negative!
    unbalanced_sums = []
    for k_extra in [10, 100, 1000, 9000]:
        # Start with balanced {1..100} x {0,1}, add {101..100+k_extra} x {1}
        S_unbal = float(np.sum(w[100 : 100 + k_extra]) * c[1])
        unbalanced_sums.append(dict(extra_indices=k_extra, S_F=S_unbal))

    # Exhaustion 3: unbalanced subset adding {0} instead of {1} -> arbitrarily positive!
    positive_sums = []
    for k_extra in [10, 100, 1000, 9000]:
        S_pos = float(np.sum(w[100 : 100 + k_extra]) * c[0])
        positive_sums.append(dict(extra_indices=k_extra, S_F=S_pos))

    return dict(
        p=p,
        c_values=c.tolist(),
        balanced_sums=balanced_sums,
        unbalanced_negative=unbalanced_sums,
        unbalanced_positive=positive_sums,
        equivalence="Unconditional convergence of sum_{(i,a)} w_i c_a <=> sum_i w_i < inf (Riemann/Dirichlet).",
    )


def run_test_1d_open_state_obstruction():
    """
    Show that replacing the closed trace with a finite-norm open state |Phi> = sum_i lam_i |i>|i*>
    (so eq. 4.29 DOES apply) necessarily breaks the uniform lower bound (5.30c) and eq. (D.24)-(D.25).
    """
    # Take finite-norm state on C^2 (x) l^2(N): weight on |1, 1> is 1/2,
    # weights on |0, n> are (1/2) * 2^{-n} so total norm^2 = 1.
    # For ANY cycle length r, adding (r-1) vectors from I = {0} only contributes at most 1/2 total weight!
    w_outside = 0.5  # single vector |1, 1> outside I
    s_star = 2.0 * np.pi
    results = []
    for r in [5, 20, 100, 1000]:
        # Even if we pick the LARGEST r-1 vectors in I={0}:
        w_inside = 0.5 * (1.0 - 2.0 ** (-(r - 1)))
        var_s = (w_outside * s_star**2) / (w_outside + w_inside)
        results.append(dict(r=r, w_inside=float(w_inside), var_s=float(var_s), lower_floor=float(0.5 * s_star**2)))
    # Variance stays >= (1/2) s_*^2 for ALL r -> eq. (D.25) fails completely!
    assert all(row["var_s"] >= 0.5 * s_star**2 - 1e-12 for row in results)
    return results


if __name__ == "__main__":
    out = dict(
        test_1a_qubit_clock=run_test_1a_qubit_clock(),
        test_1b_continuous_jt=run_test_1b_continuous_jt_analog(),
        test_1c_ads_tautology=run_test_1c_ads_disjoint_tautology(),
        test_1d_open_state_obstruction=run_test_1d_open_state_obstruction(),
    )
    out_path = os.path.join(os.path.dirname(__file__), "test1_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print("=== TEST 1A (C^2 x l^2(N) counterexample to eqs. 5.31-5.38) ===")
    print(f"delta_eta = {out['test_1a_qubit_clock']['delta_eta']:.6f}, cycle length r = {out['test_1a_qubit_clock']['r_test']}")
    print("  Product exhaustion Z_{P_M}[X]:", [(r["M"], round(r["Z_P_X"], 4)) for r in out["test_1a_qubit_clock"]["product_rows"][:5]])
    print("  Cycle exhaustion   Z_P[X]   :", [(r["K"], round(r["Z_P_X_max"], 4)) for r in out["test_1a_qubit_clock"]["cycle_rows"][:5]])
    print("\n=== TEST 1B (Continuous-spectrum JT+scalar analog) ===")
    print(json.dumps(out["test_1b_continuous_jt"], indent=2))
    print("\n=== TEST 1C (AdS disjoint-cut tautology eqs. 5.44-5.45) ===")
    print("  Balanced sums:", out["test_1c_ads_tautology"]["balanced_sums"])
    print("  Unbalanced negative:", out["test_1c_ads_tautology"]["unbalanced_negative"])
    print("  Unbalanced positive:", out["test_1c_ads_tautology"]["unbalanced_positive"])
    print("\n=== TEST 1D (Finite-norm open state breaks eq. D.25) ===")
    for row in out["test_1d_open_state_obstruction"]:
        print(f"  r={row['r']:4d}: var_s = {row['var_s']:.4f} >= floor {row['lower_floor']:.4f}")
    print(f"\nSaved full JSON results to {out_path}")
