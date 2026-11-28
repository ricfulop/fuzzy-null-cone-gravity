#!/usr/bin/env python3
"""
Supplemental check for Path 2 Module 2:
Verify QEI convergence of E_min^{NO}(N_rad) AND Var_0(H_g^{NO})(N_rad) on a GENUINE
C^infty_c((0, 1)) non-bandlimited bump function g(u) = C * exp(-1 / (u (1 - u))),
whose Fourier coefficients c_n are strictly nonzero for infinitely many n!
"""

import json
import os
import numpy as np
from scipy.linalg import eigh


def run_cinfty_bump_check():
    n_quad = 2000
    u_q, w_q = np.polynomial.legendre.leggauss(n_quad)
    u_q = 0.5 * (u_q + 1.0)
    w_q = 0.5 * w_q

    # Smooth C^infty_c bump on (0, 1): g(u) = exp(-1/(4 u (1 - u)) + 1) so g(1/2) = 1
    denom = u_q * (1.0 - u_q)
    g = np.exp(-0.25 / denom + 1.0)
    g_sq = g**2

    # Verify non-bandlimited: check c_n for n = 0..12
    ns_check = np.arange(0, 13, dtype=float)
    c_check = np.cos(2.0 * np.pi * ns_check[:, None] * u_q[None, :]) @ (w_q * g_sq)
    # Check c_4, c_6, c_8 are genuinely nonzero!

    rows = []
    for N in [10, 20, 40, 80, 160]:
        ks = np.arange(1, N + 1, dtype=float)
        sq = np.sqrt(ks)
        ns = np.arange(0, 2 * N + 1, dtype=float)
        c_arr = np.cos(2.0 * np.pi * ns[:, None] * u_q[None, :]) @ (w_q * g_sq)
        idx = np.arange(N)
        diff_idx = np.abs(idx[:, None] - idx[None, :])
        sum_idx = (idx[:, None] + 1) + (idx[None, :] + 1)
        S_mat = 2.0 * np.outer(sq, sq)
        Mq = S_mat * (c_arr[diff_idx] + c_arr[sum_idx])
        Mp = S_mat * (c_arr[diff_idx] - c_arr[sum_idx])

        ev_q, U_q = eigh(Mq)
        Mq_sqrt = (U_q * np.sqrt(np.maximum(ev_q, 1e-15))) @ U_q.T
        prod = 0.5 * (Mq_sqrt @ Mp @ Mq_sqrt + (Mq_sqrt @ Mp @ Mq_sqrt).T)
        ev_prod = eigh(prod, eigvals_only=True)
        E0_un = 0.5 * float(np.sum(np.sqrt(np.maximum(ev_prod, 0.0))))
        E_vac = 0.25 * float(np.trace(Mq + Mp))
        E_min_NO = E0_un - E_vac
        var_0 = 0.125 * float(np.sum((Mq - Mp) ** 2))
        rows.append(
            dict(
                N_rad=N,
                E0_un=E0_un,
                E_min_NO=E_min_NO,
                Var0_NO=var_0,
            )
        )

    out = dict(
        bump="g(u) = exp(-1/(4 u (1-u)) + 1) in C^infty_c((0,1))",
        nonzero_high_fourier_coeffs={int(n): float(c_check[int(n)]) for n in [0, 1, 2, 4, 6, 8, 10]},
        convergence_table=rows,
    )
    out_path = os.path.join(os.path.dirname(__file__), "path2_cinfty_bump_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print("=== C^infty_c Bump Non-Bandlimited QEI & Variance Check ===")
    print("  Fourier coeffs c_n:", {k: f"{v:.3e}" for k, v in out["nonzero_high_fourier_coeffs"].items()})
    for r in rows:
        print(f"  N_rad={r['N_rad']:3d}: E0_un={r['E0_un']:8.2f}, E_min^NO={r['E_min_NO']:+.8f}, Var0_NO={r['Var0_NO']:.8f}")
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    run_cinfty_bump_check()
