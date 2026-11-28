#!/usr/bin/env python3
"""
Test 3: Interacting C*-Algebra, Pure-State Extension Non-Uniqueness,
        Reduced Hamiltonian H_c Singularities, and the Section 4.4 "Check" (eqs. 4.32-4.34).

Checks:
  3A. Exact asymptotic expansion of the Section 4.4 free 1D sum Li_2(e^{-2i/s}) / (2 s^2)
      and demonstration that including 4D angular harmonics (j, m) makes it quadratically
      UV-divergent even at finite G_N (and that it never touches any vertex of H_c).
  3B. Exact divergence of the lapse factor 1/N(eps, y) = exp[Xi_y(1) - Xi_y(eps)] appearing
      in the reduced Hamiltonian H_c (eqs. 3.54-3.55), showing H_c violates the bounded-generator
      hypothesis of Appendix C.2 eq. (C.7).
  3C. Concrete finite-dimensional C*-algebra witness showing non-uniqueness of Kadison-Singer /
      Hahn-Banach pure-state extensions omega from A_in to A: two pure extensions agree identically
      on A_in = M_A yet give completely different transition amplitudes <a|U_M|b> = omega(a^dag U_M b).
"""

import json
import os
import numpy as np
import sympy as sp


def run_test_3a_pair_creation_check():
    """Verify eqs. (4.32)-(4.34) and test what happens when restoring 4D angular harmonics."""
    # Exact numerical evaluation of S_N(s) = (1/(2 s^2)) sum_{k=1}^N exp(-2 i k / s) / k^2
    # Compare with exact Li_2(e^{-2i/s}) / (2 s^2) and its small-(1/s) asymptotic expansion:
    #   Li_2(e^{-2i/s}) / (2 s^2) = pi^2/(12 s^2) - pi/(2 s^3) - i log(s)/s^3 + i (log(2)-1)/s^3 + O(s^{-4})
    s_vals = [10.0, 100.0, 1000.0, 10000.0]
    rows_1d = []
    for s in s_vals:
        # High-precision sum up to N = 500_000 (tail error ~ 1/(2 s^2 N) < 1e-8)
        k = np.arange(1, 500001, dtype=np.float64)
        exact_num = np.sum(np.exp(-2j * k / s) / (k**2)) / (2.0 * s**2)
        asymp = (
            (np.pi**2) / (12.0 * s**2)
            - np.pi / (2.0 * s**3)
            + 1j * (np.log(2.0) - 1.0 - np.log(s)) / (s**3)
        )
        # Note: at s = infinity (G_N = 0), normalized by 2 s^2, the sum ALREADY converges to pi^2/6!
        s2_times_exact = (2.0 * s**2) * exact_num
        rows_1d.append(
            dict(
                s=s,
                exact_real=float(exact_num.real),
                exact_imag=float(exact_num.imag),
                asymp_real=float(asymp.real),
                asymp_imag=float(asymp.imag),
                rel_err_imag=float(abs(exact_num.imag - asymp.imag) / abs(exact_num.imag)),
                two_s2_times_sum_real=float(s2_times_exact.real),
                pi2_over_6=float(np.pi**2 / 6.0),
            )
        )

    # Now restore 4D angular harmonics j = 2..J (multiplicity 2j+1 per j, 2 polarizations):
    # Each angular harmonic (j,m,p) has the same radial variance (Delta q_k)^2 = 1/(2 s k)
    # -> sum over all angular modes multiplies the 1D result by 2 * ((J+1)^2 - 4) -> +inf as J -> inf!
    s_fixed = 100.0
    base_val = abs(rows_1d[1]["exact_real"] + 1j * rows_1d[1]["exact_imag"])
    angular_divergence = []
    for J in [2, 10, 100, 1000, 10000]:
        mult = 2 * ((J + 1) ** 2 - 4)
        angular_divergence.append(
            dict(J=J, angular_modes=mult, amplitude_mag_at_s100=float(mult * base_val))
        )

    return dict(
        asymptotic_1d_check=rows_1d,
        angular_mode_divergence_at_finite_GN=angular_divergence,
        diagnosis=(
            "Eq. (4.33) is a 1D free-oscillator sum with a manually injected phase shift 1/s. "
            "Normalized by G_N^2 (i.e. 1/s^2), it already converges at G_N=0 to pi^2/6; the log N "
            "divergence in eq. (4.34) is merely d/dx [sum e^{ikx}/k^2]|_{x=0} = i sum 1/k. "
            "Restoring the 2*((J+1)^2-4) angular graviton modes on S^2 makes it quadratically UV divergent "
            "even at finite G_N, and no nonlinear interaction vertex of H_c enters."
        ),
    )


def run_test_3b_Hc_lapse_singularity():
    """
    Use Test 2's exact Williamson lower bound on Xi_y(1) - Xi_y(eps) to evaluate the
    operator lower bound on 1/N(eps, y) = N(1, y)^{-1} exp[Xi_y(1) - Xi_y(eps)] entering
    H_Sigma in eqs. (3.54)-(3.55).
    """
    t2_path = os.path.join(os.path.dirname(__file__), "test2_results.json")
    with open(t2_path) as f:
        t2 = json.load(f)
    c_rad = t2["asymptotic_c_rad"]
    s_macro = t2["s_macro"]

    rows = []
    for (N, J) in [(10, 10), (100, 100), (500, 500), (1000, 1000)]:
        n_ang = (J + 1) ** 2 - 4
        E0_Xi = (n_ang * c_rad * (N**2)) / s_macro
        # Since Xi_y(eps) = 0, eq. (3.55) gives N(eps, y)/N(1, y) = exp(-Xi_y(1))
        # and the inverse lapse factor in eq. (3.54) is 1/N(eps, y) = N(1, y)^{-1} exp(+Xi_y(1))!
        inv_lapse_min = float(np.exp(min(E0_Xi, 700.0)))
        lapse_max = float(np.exp(-min(E0_Xi, 700.0)))
        rows.append(
            dict(
                N=N,
                J=J,
                E0_Xi=float(E0_Xi),
                lapse_ratio_upper_bound=lapse_max,
                inv_lapse_factor_lower_bound=inv_lapse_min,
            )
        )
    return dict(
        s_macro=s_macro,
        lapse_collapse_table=rows,
        conclusion=(
            "In eq. (3.54)-(3.55), H_c contains (d_v F)/(2 N(v,y)) with 1/N(eps,y) = N(1,y)^{-1} exp(Xi_y(1)). "
            "As the UV cutoff N, J -> inf, N(eps,y) -> 0 and 1/N(eps,y) -> +inf on every state in K_A, "
            "violating the bounded, self-adjoint generator hypothesis of Appendix C.2 eq. (C.7)."
        ),
    )


def run_test_3c_pure_extension_nonuniqueness():
    """
    Exact operator-algebraic demonstration of Kadison-Singer / Hahn-Banach pure-state
    extension non-uniqueness (Section 4.2, p. 29 & Footnote 20, p. 38).

    Let A = M_4(C) (representing the full algebra with two area sectors A_1, A_2, each dim 2),
    and let A_in = M_2(C) (+) M_2(C) be the area-preserving block-diagonal subalgebra.
    Fix a pure state |Omega_A> = (1, 0)^T in sector A_1 on A_in:
      omega_in(diag(X_1, X_2)) = <Omega_A | X_1 | Omega_A> = (X_1)_{00}.

    Wait: what if instead A_in = M_2(C) (x) I_2 (diagonally embedded, representing area-preserving
    operators acting identically across sectors, or a unital subalgebra where U_M maps A_1 -> A_2)?
    Then ANY unit vector |v> = (alpha, beta)^T in C^2 defines a pure state |Omega_A> (x) |v> on A
    whose restriction to A_in = M_2(C) (x) I_2 is IDENTICAL to <Omega_A | . | Omega_A>,
    yet for an area-changing unitary U_M = I_2 (x) sigma_x (mapping A_1 <-> A_2):
      omega_v(U_M) = 2 Re(alpha^* beta) ranges over the ENTIRE interval [-1, +1]!
    """
    # Explicit matrices on C^2 (fixed-area Fock) (x) C^2 (area sectors)
    I2 = np.eye(2, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    U_M = np.kron(I2, sx)  # area-changing unitary in A \ A_in
    # Random test observable c in A_in = M_2(C) (x) I_2
    rng = np.random.default_rng(42)
    X_in = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    c_in = np.kron(X_in, I2)

    omega_A = np.array([1.0, 0.0], dtype=complex)  # |Omega_A>
    expected_in = complex(omega_A.conj() @ X_in @ omega_A)

    extensions = []
    for theta in [0.0, np.pi / 6, np.pi / 4, np.pi / 2, 3 * np.pi / 4]:
        v = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
        psi = np.kron(omega_A, v)
        val_in = complex(psi.conj() @ c_in @ psi)
        val_UM = complex(psi.conj() @ U_M @ psi)
        assert abs(val_in - expected_in) < 1e-14
        extensions.append(
            dict(
                theta_rad=float(theta),
                omega_on_A_in_diff=float(abs(val_in - expected_in)),
                transition_amp_omega_UM=float(val_UM.real),
            )
        )

    return dict(
        extensions=extensions,
        conclusion=(
            "All pure extensions omega_theta agree identically on A_in (diff = 0), "
            "generate irreducible Type I factors B(H_omega), and yet assign arbitrary transition "
            "amplitudes omega(U_M) in [-1, 1] to area-changing unitaries. Thus pure-state extension "
            "alone cannot fix gravitational dynamics or equate omega(a^dag U_M b) with the path integral (4.5)."
        ),
    )


if __name__ == "__main__":
    res = dict(
        test_3a_pair_creation=run_test_3a_pair_creation_check(),
        test_3b_Hc_lapse=run_test_3b_Hc_lapse_singularity(),
        test_3c_pure_extension=run_test_3c_pure_extension_nonuniqueness(),
    )
    out_path = os.path.join(os.path.dirname(__file__), "test3_results.json")
    with open(out_path, "w") as f:
        json.dump(res, f, indent=2)

    print("=== TEST 3A: Section 4.4 Li_2(e^{-2i/s})/(2s^2) Check & 4D Angular Divergence ===")
    for r in res["test_3a_pair_creation"]["asymptotic_1d_check"]:
        print(
            f"  s={r['s']:7.1f}: 2s^2*Re(S)={r['two_s2_times_sum_real']:.6f} (pi^2/6={r['pi2_over_6']:.6f}), "
            f"Im(S)={r['exact_imag']:.4e} vs asymp={r['asymp_imag']:.4e} (rel_err={r['rel_err_imag']:.2e})"
        )
    print("  Restoring 4D angular harmonics at fixed s=100:")
    for r in res["test_3a_pair_creation"]["angular_mode_divergence_at_finite_GN"]:
        print(f"    J={r['J']:5d} ({r['angular_modes']:10d} modes): |amp| = {r['amplitude_mag_at_s100']:.4e}")

    print("\n=== TEST 3B: Reduced Hamiltonian H_c Inverse-Lapse Singularity (eqs. 3.54-3.55) ===")
    for r in res["test_3b_Hc_lapse"]["lapse_collapse_table"]:
        print(
            f"  N=J={r['N']:4d}: E0_Xi={r['E0_Xi']:10.3e}, N(eps)/N(1) <= {r['lapse_ratio_upper_bound']:.3e}, "
            f"1/N(eps) >= {r['inv_lapse_factor_lower_bound']:.3e}"
        )

    print("\n=== TEST 3C: Pure-State Extension Non-Uniqueness (Section 4.2, p. 29) ===")
    for r in res["test_3c_pure_extension"]["extensions"]:
        print(
            f"  theta={r['theta_rad']:.3f}: |omega(c_in) - <Omega|c_in|Omega>| = {r['omega_on_A_in_diff']:.1e}, "
            f"omega(U_M) = {r['transition_amp_omega_UM']:+.6f}"
        )
    print(f"\nSaved full JSON results to {out_path}")
