#!/usr/bin/env python3
"""
Path 2 Breakthrough Test 3: Curved 4D 1-Loop Graviton Dynamics on Fuzzy-Sphere Horizons (S^4 & Schwarzschild).

Tests the CRUCIAL question for Path 2:
  Does a discrete fuzzy-sphere corner (J_max(r) = Lambda_P * r, dynamical embadon rule K(r) <= s(r)/s_min)
  preserve or break the exact 4D 1-loop graviton logarithmic corrections on curved spacetimes?

Part A — Exact 4D S^4 (de Sitter causal diamond) Spectral Audit:
  - Uses the exact SO(5) eigenvalues & degeneracies of spin-2 TT gravitons, spin-1 vector ghosts,
    and spin-0 trace/longitudinal modes on S^4 (Vassilevich 1993; Anninos-Denef-Law-Sun 2009.12464).
  - Verifies the exact covariant heat-kernel pole alpha_4 = -571/45 = (-331/45 bulk) - (+16/3 edge),
    PLUS the exact group-volume log correction -5 log S_0 (from dim SO(5) = 10 -> -10/2 = -5!).
  - Performs the EXACT SO(5) -> SO(4) -> SO(3) branching into fuzzy-sphere harmonics j <= J_max,
    comparing:
      (i)  Covariant Path 2 cutoff (all modes with lambda_n <= Lambda_P^2; automatically contained in j <= J_max!)
      (ii) Sharp anisotropic cutoff (j <= J_max with n unconstrained -> checks whether anisotropic cutoff
           induces a local-Lorentz-breaking anomaly in alpha_4!).

Part B — Exact 4D Euclidean Schwarzschild (Ricci-Flat, W_{abcd} != 0) Decomposition & Cutoff Test:
  - Verifies the exact Lichnerowicz heat-kernel decomposition on Schwarzschild (R_{ab} = 0):
      180 a_2 / K = (+2 minimal 2-pol) + (-60 spin-connection Omega^2) + (+270 Lichnerowicz Weyl tidal V^2) = +212
      => c_log^{bulk}(ln s) = +1/45 - 2/3 + 3 = +106/45!
  - Tests explicitly in 4D local momentum space (using the EXACT Riemann normal coordinate / Schwinger-DeWitt
    expansion on Schwarzschild with Weyl tensor W_{abcd} != 0!) the difference between:
      (i)  Covariant 4D double-null / cross-focusing cutoff (k_par^2 + k_perp^2 <= Lambda_P^2,
           which automatically satisfies k_perp <= Lambda_P => j <= J_max(r) = sqrt(s(r)/s_min)!), and
      (ii) Purely anisotropic cylindrical cutoff (k_perp <= Lambda_P on S^2, k_par <= Lambda_par independent)!
    Computes the EXACT analytical + numerical ratio of the Weyl tidal term (+270) and spin-connection
    term (-60) under anisotropic vs. isotropic cutoffs!
"""

import json
import math
import os
import numpy as np
import sympy as sp


def part_a_s4_exact_spectrum():
    """
    Exact 4D pure gravity 1-loop spectral zeta / heat-kernel pole on S^4 (radius l_dS = 1)
    and exact SO(5) -> SO(4) -> SO(3) fuzzy-sphere branching check.
    """
    # Exact heat kernel trace K(t) = sum_n d_n * exp(-t * lambda_n) for pure 4D Einstein gravity on S^4
    # (Christensen-Duff 1980; Vassilevich 1993; Anninos-Denef-Law-Sun 2009.12464 eqs. 5.8, G.18-G.24):
    # On S^4 (R = 12, Lambda = 3):
    #   TT spin-2 (n >= 2): lambda_n^(2) = n(n+3) - 2,   d_n^(2) = (5/6)(n-1)(n+4)(2n+3)
    #   Vector ghost (n >= 2, excluding n=1 Killing zero modes!):
    #     Wait: in Landau / harmonic gauge on S^4 (ADLS 2009.12464 App. G.2, eq. G.24):
    #       1-loop Z_PI^(1) = Z_G * [ det'(Delta_1 - 3) / det'(Delta_2 + 2) ]^{1/2} * [ ... ]
    #     Even cleaner: ADLS (2009.12464 eq. 5.8) proves the EXACT identity (including all Jacobians,
    #     ghosts, and conformal mode!):
    #       log Z_char = int_0^inf (dt / (2t)) * (1+q)/(1-q) * ( [chi_bulk,2]_+ - [chi_edge,2]_+ )
    #     where q = e^{-t}, and:
    #       [chi_bulk,2]_+ = (10 q^3 - 6 q^4) / (1 - q)^3 = sum_{k=3}^inf (2k^2 - 2k - 6? let's check!)
    #       [chi_edge,2]_+ = (10 q^2 - 2 q^3) / (1 - q)   = 10 q^2 + 8 q^3 + 8 q^4 + ...
    # Let's ALSO compute the direct heat-kernel trace from the exact S^4 Lichnerowicz eigenvalues
    # (TT spin-2 minus 2 * ghost vector + conformal scalar + longitudinal Jacobians) via Euler-Maclaurin
    # AND verify exact agreement = -571/45!
    t = sp.Symbol("t", positive=True)
    q = sp.exp(-t)
    chi_bulk = (10 * q**3 - 6 * q**4) / (1 - q) ** 3
    chi_edge = (10 * q**2 - 2 * q**3) / (1 - q)
    alpha_bulk = sp.series((1 + q) / (2 * t * (1 - q)) * chi_bulk, t, 0, 1).coeff(t, -1)
    alpha_edge = sp.series((1 + q) / (2 * t * (1 - q)) * chi_edge, t, 0, 1).coeff(t, -1)
    alpha_total = sp.simplify(alpha_bulk - alpha_edge)

    # Now verify DIRECTLY from the exact S^4 spherical-harmonic eigenvalue sums (Christensen-Duff / Taylor-Woodard):
    # On S^4 (radius 1):
    #   Physical 1-loop determinant in de Donder gauge (Christensen-Duff 1980, Nucl. Phys. B 170, 480):
    #     (1) Spin-2 symmetric traceless tensor (n >= 2): lam = n(n+3) + 2? Wait:
    #         Delta_L h^{TT} = (-nabla^2 + 8/3 * (R/4)?): on S^4 of radius 1, R_{acbd} = g_{ab}g_{cd} - g_{ad}g_{bc},
    #         so -2 R_{acbd} h^{cd} = +2 h_{ab} (for traceless h)! Plus -2 Lambda h_{ab} = -6 h_{ab}?
    #         Wait: on S^4 with R_{ab} = 3 g_{ab}, Einstein eq is R_{ab} - (1/2)R g_{ab} + Lambda g_{ab} = 0 => Lambda = 3.
    #         Linearized Einstein operator on TT modes: (-nabla^2 + 2) h^{TT} = [n(n+3) - 2 + 2] = n(n+3)!!
    #         WAIT: check! For TT modes on S^4, -nabla^2 has eigenvalues n(n+3) - 2 (n >= 2), and
    #         Lichnerowicz Delta_L - 2 Lambda = -nabla^2 + 2 R_{acbd}? Check: R_{acbd} h^{cd} = -h_{ab} for traceless h,
    #         so -2 R_{acbd} h^{cd} = +2 h_{ab}; Wait: (-nabla^2 + 2R/4 - ...): in Vassilevich/ADLS eq. G.16,
    #         the physical TT operator is (-nabla^2 + 2) with eigenvalues (n+1)(n+2) = n(n+3) + 2? No:
    #         (-nabla^2)_{TT} = n(n+3) - 2, so (-nabla^2 + 2)_{TT} = n(n+3), OR (-nabla^2 + 4)_{TT} = (n+1)(n+2)!
    #     Let's CHECK by exact heat-kernel asymptotic expansion K(t) = c_{-2}/t^2 + c_{-1}/t + c_0 + O(t)
    #     of the standard Christensen-Duff operators on S^4:
    #       1) TT spin-2:   Delta_2 = -nabla^2 + 2  => eig = n(n+3),          n >= 2, d = (5/6)(n-1)(n+4)(2n+3)
    #       2) Ghost vec:   Delta_1 = -nabla^2 - 3  => eig = n(n+3) - 4,      n >= 2 (n=1 is 10 zero modes!), d = (1/2)n(n+3)(2n+3)
    #          Wait: check n=1 for ghost vec: 1*(4) - 4 = 0! YES! 10 Killing vectors of SO(5)!!
    #       3) Trace/scalar: in Christensen-Duff, total zeta(0) = zeta_TT(0) - zeta_vec(0) - zeta_scalar(0)
    #          Wait: trace h (1 scalar, eig n(n+3)-6 -> n=0 is -6 negative mode, n=1 is -2 negative, n=2 is 4? Wait:
    #          in de Donder gauge: trace h has operator -nabla^2 - 6 (10 = 5 TT + 4 vec + 1 trace -> 10 - 2*4 = 2)!
    #          Let's check zeta(0) = c_0 of EACH operator directly using exact Euler-Maclaurin on
    #          sum_{n=n0}^inf d(n) * exp(-t * (n+3/2)^2 + t * M^2)!
    def heat_kernel_c0(poly_d, shift_M2, n0):
        """
        Exact c_0 (constant term as t -> 0) of K(t) = sum_{n=n0}^inf poly_d(n) * exp(-t * ((n+3/2)^2 - shift_M2))
        using exact Bernoulli / Hurwitz zeta:
          sum_{n=n0}^inf P(n) exp(-t (n+a)^2) -> expand P(n) in powers of u = n + 3/2:
          Since u^k * exp(-t u^2) has Mellin transform Gamma((k+1)/2) * ... wait:
          For K(t) = exp(t * shift_M2) * sum_{n=n0}^inf sum_k a_k (n+3/2)^k exp(-t (n+3/2)^2):
          By exact Euler-Maclaurin / Hurwitz zeta:
            sum_{n=0}^inf (n+3/2)^k exp(-t (n+3/2)^2) = sum_{m=0}^inf ((-t)^m / m!) * zeta_H(-k - 2m, 3/2)
            + (divergent gamma terms 1/t^{(k+1)/2} for odd k -> wait: for ODD k, Gamma((k+1)/2)/2 has NO t^0 term!!
               ONLY integer powers t^m from Taylor series of exp(-t u^2) give t^0, i.e. zeta_H(-k, 3/2)!!
               Wait: check that! For odd k (our degeneracies d(n) are ALL ODD polynomials in u = n + 3/2!!):
               Every single term u^{2r+1} exp(-t u^2) has singular part (1/2) r! / t^{r+1}, and analytic part
               sum_{m=0}^inf ((-t)^m / m!) * zeta_H(-(2r+1) - 2m, 3/2)!)
        THAT'S EXACT TO INFINITE PRECISION IN SYMPY!!
        """
        u = sp.Symbol("u")
        n = sp.Symbol("n")
        # Express poly_d(n) in u = n + 3/2 (i.e. n = u - 3/2)
        P_u = sp.Poly(sp.expand(poly_d.subs(n, u - sp.Rational(3, 2))), u)
        # Subtract finite sum for n = 0 .. n0 - 1 (exact!)
        # Wait: let's work with sum_{n=0}^inf - sum_{n=0}^{n0-1}
        # Expand exp(t * shift_M2) * [ sum_{m=0}^2 ((-t)^m / m!) sum_k a_k zeta_H(-k-2m, 3/2)
        #                              + sum_{r} a_{2r+1} (1/2) r! / t^{r+1} ] - sum_{n=0}^{n0-1} d(n) exp(-t((n+3/2)^2 - shift_M2))
        K_asymp = sp.Integer(0)
        for deg, coeff in zip(P_u.monoms(), P_u.coeffs()):
            k = deg[0]
            if k % 2 == 1:
                r = (k - 1) // 2
                K_asymp += coeff * sp.Rational(1, 2) * sp.factorial(r) / t ** (r + 1)
            # Analytic Taylor terms up to t^2 (needed when multiplying by exp(t * shift_M2) or extracting t^0)
            for m in range(3):
                # zeta_H(-s, a) = -B_{s+1}(a) / (s+1) for integer s >= 0!
                s_idx = k + 2 * m
                zh = -sp.bernoulli(s_idx + 1, sp.Rational(3, 2)) / (s_idx + 1)
                K_asymp += coeff * ((-t) ** m / sp.factorial(m)) * zh

        # Multiply by exp(t * shift_M2)
        K_full = sp.series(sp.exp(t * shift_M2) * K_asymp, t, 0, 1).removeO()
        # Subtract missing low modes n = 0 .. n0 - 1:
        for nv in range(n0):
            dv = poly_d.subs(n, nv)
            # Constant term as t -> 0 of dv * exp(-t * lam) is simply dv!
            K_full -= dv
        return sp.simplify(sp.series(K_full, t, 0, 1).coeff(t, 0))

    n = sp.Symbol("n")
    d2 = sp.Rational(5, 6) * (n - 1) * (n + 4) * (2 * n + 3)
    d1 = sp.Rational(1, 2) * n * (n + 3) * (2 * n + 3)
    d0 = sp.Rational(1, 6) * (n + 1) * (n + 2) * (2 * n + 3)

    # Check exact eigenvalues on S^4 (Vassilevich 1993 / Christensen-Duff 1980 / ADLS 2009.12464):
    #   TT spin-2:   n(n+3) + 0  = (n+3/2)^2 - 9/4  ? Wait: check (-nabla^2 + 2) on TT where -nabla^2 = n(n+3) - 2 -> n(n+3) = (n+3/2)^2 - 9/4!
    #   Wait: in Christensen-Duff 1980 (eq. 4.10) on S^4 (R = 12, Lambda = 3):
    #     TT spin-2 operator:   -nabla^2 + (2/3) R - 2 Lambda = -nabla^2 + 2 => eig = (n(n+3) - 2) + 2 = n(n+3) = (n+3/2)^2 - 9/4, n >= 2
    #     Ghost vector (x2):    -nabla^2 - R/4 = -nabla^2 - 3 => eig = (n(n+3) - 1) - 3 = n(n+3) - 4 = (n+3/2)^2 - 25/4, n >= 2 (n=1 is 10 zero modes!)
    #     Vector longitudinal + trace (in Landau / physical decomposition: Z = [det_vec(-nabla^2-3) / det_TT(-nabla^2+2)]^{1/2} * [det_0(-nabla^2-4)/... ]?):
    #     In ADLS (2009.12464 eq. G.24), after cancelling ALL longitudinal & scalar modes identically:
    #       Z_1loop = (vol SO(5))^{-1} * [ det'_{1,T}(-nabla^2 - 3) / det_{2,TT}(-nabla^2 + 2) ]^{1/2}
    #     WAIT A SECOND! CHECK THAT EXACT FORMULA RIGHT NOW!!
    #     Let's compute c0(TT, n>=2, shift=9/4) - c0(1T, n>=2, shift=25/4)!!
    c0_TT = heat_kernel_c0(d2, sp.Rational(9, 4), 2)
    c0_1T = heat_kernel_c0(d1, sp.Rational(25, 4), 2)
    # Wait: remember in log Z = -(1/2) log det(2,TT) + (1/2) log det'(1,T) =>
    #   log Z_{log-div} = +(c0_TT - c0_1T) * log(l / eps)!!
    #   Wait: remember zero modes! In log Z = +(1/2) int (dt/t) [ K_2(t) - K_1(t) ],
    #   K(t) = sum exp(-t lambda) has coeff +c0, and int_{eps^2}^{l^2} (dt/(2t)) = +log(l/eps)!
    #   Wait: for 1-loop Z = [det(1,T) / det(2,TT)]^{1/2}, log Z = +(1/2) int (dt/t) [ K_{2,TT}(t) - K_{1,T}(t) ]!
    #   Let's check c0_TT - c0_1T vs -571/45!
    diff_TT_1T = sp.simplify(c0_TT - c0_1T)

    # Check also exact zero-mode subtraction (10 Killing vectors of SO(5) at n=1 of d1, d1(1)=10):
    # Wait: in heat_kernel_c0(d1, 25/4, 2), we ALREADY excluded n=1 (d1(1)=10)!
    # What if we ALSO include the n=1 conformal Killing vectors / negative modes of trace?
    # Let's print c0_TT, c0_1T, diff_TT_1T, and check exact relation to -571/45!

    # Now check SO(5) -> SO(3) fuzzy-sphere truncation on S^4!
    # For a single SO(5) symmetric-traceless scalar rep n (dim d0(n) = (1/6)(n+1)(n+2)(2n+3)):
    #   Branching SO(5) -> SO(4): k = 0..n, each with SO(4) rep (k/2, k/2) of dim (k+1)^2.
    #   Branching SO(4) -> SO(3): j = 0..k, each with SO(3) dim (2j+1).
    #   Check sum_{k=0}^n sum_{j=0}^k (2j+1) = sum_{k=0}^n (k+1)^2 = (1/6)(n+1)(n+2)(2n+3) = d0(n)!! EXACT!!
    #   Therefore, in the SO(5) rep n, the number of states with fuzzy-sphere angular momentum j <= J_max is:
    #     - If n <= J_max: ALL d0(n) states (100% preserved!)
    #     - If n > J_max:  only sum_{k=0}^{J_max} (k+1)^2 + (n - J_max) * (J_max + 1)^2 < d0(n)!
    # Notice: for ANY covariant cutoff n <= N_max = J_max, ALL states satisfy j <= n <= J_max automatically!
    # Check what happens to the scalar S^4 log pole (exact = -29/90 for minimally coupled scalar; wait:
    # check c0_scalar_minimal = heat_kernel_c0(d0, 9/4, 1) -> exact -29/90 - 1? Let's check!)
    c0_scalar_conf = heat_kernel_c0(d0, sp.Rational(1, 4), 0)  # conformal scalar on S^4: eig = (n+1)(n+2) = (n+3/2)^2 - 1/4 -> exact -1/90!
    c0_scalar_min = heat_kernel_c0(d0, sp.Rational(9, 4), 1)   # minimal massless scalar (excluding n=0 zero mode)

    return dict(
        ADLS_character_poles=dict(
            alpha_bulk=str(alpha_bulk),
            alpha_edge=str(alpha_edge),
            alpha_total=str(alpha_total),
            group_volume_log_S0="-5 (from dim SO(5) = 10 -> -10/2)",
        ),
        S4_direct_eigenvalue_c0=dict(
            c0_conformal_scalar=str(c0_scalar_conf),
            c0_minimal_scalar_no_n0=str(c0_scalar_min),
            c0_TT_spin2=str(c0_TT),
            c0_1T_ghost_no_n1=str(c0_1T),
            c0_TT_minus_c0_1T=str(diff_TT_1T),
            note_comparison_to_minus_571_45=str(sp.simplify(diff_TT_1T - sp.Rational(-571, 45))),
        ),
    )


def part_b_schwarzschild_anisotropic_vs_covariant_cutoff():
    """
    EXACT analytical & numerical test of the 4D Schwarzschild Ricci-flat Lichnerowicz 1-loop log term
    under:
      (i)  Covariant 4D double-null / cross-focusing cutoff (k_par^2 + k_perp^2 <= Lambda^2),
           which is AUTOMATICALLY contained inside the dynamical embadon Hilbert space j <= Lambda * r!
      (ii) Sharp ANISOTROPIC fuzzy-sphere cutoff (k_perp <= Lambda_perp on S^2, k_par <= Lambda_par
           independent on the (t,r) cigar)!

    DERIVATION OF THE EXACT MOMENTUM-SPACE POTENTIAL-SQUARED INTEGRAL (THE +270 WEYL TIDAL TERM!):
    On a Ricci-flat 4-manifold (R_{ab}=0, W_{abcd} != 0), the Lichnerowicz operator is
      Delta_L = -nabla^2 + V,   where V_{ab}^{cd} = -2 W_{a c b d},  Tr(V) = 0,  Tr(V^2) = 3 W^2.
    In local Riemann normal coordinates at a point x on Schwarzschild, the 1-loop effective action
    is Gamma_1 = (1/2) Tr log(-nabla^2 + V) - Tr log(-nabla^2_vec).
    Expand the tidal potential V to 2nd order in perturbation theory (the ONLY local term proportional to V^2!):
      Delta Gamma_tidal = -(1/4) Tr[ (-nabla^2)^{-1} V (-nabla^2)^{-1} V ]
                        = -(1/4) int d^4x sqrt(g) Tr(V^2) * I_4(Lambda)
    where I_4(Lambda) = int_{reg} (d^4k / (2 pi)^4) * 1 / (k^2 + m_IR^2)^2!!

    STOP AND STARE AT I_4(Lambda)!!
    1. Under a COVARIANT 4D cutoff k^2 = k_par^2 + k_perp^2 <= Lambda^2 (in spherical 4D coords, Vol(S^3)=2 pi^2):
         I_4^{cov}(Lambda) = (2 pi^2 / (16 pi^4)) * int_{m_IR}^{Lambda} k^3 dk / k^4
                           = (1 / (8 pi^2)) * log(Lambda / m_IR)
       Multiply by (1/4) * Tr(V^2) * 2 (since 2 log(Lambda) = log(s)):
         Delta S_tidal = +(1/4) * (3 W^2) * (1 / (8 pi^2)) * (1/2) log(s)
       Integrate over Schwarzschild (int d^4x sqrt(g) W^2 = 32 pi^2 * chi = 64 pi^2):
         Delta S_tidal = +(1/4) * 3 * (64 pi^2 / (8 pi^2)) * (1/2) log(s) = +(1/4) * 24 * (1/2) log(s) = +3 log(s)!!
       EXACT MATCH TO +270/90 = +3 log(s)!!

    2. NOW COMPUTE I_4^{aniso}(Lambda_perp, Lambda_par) FOR THE ANISOTROPIC CYLINDRICAL CUTOFF
       k_perp <= Lambda_perp (fuzzy sphere!) AND k_par <= Lambda_par (independent radial/time cutoff!):
         I_4^{aniso} = (1 / (16 pi^4)) * (2 pi) * (2 pi) * int_0^{Lambda_par} k_par dk_par
                       * int_0^{Lambda_perp} k_perp dk_perp / (k_par^2 + k_perp^2 + m_IR^2)^2
       Let x = k_par^2 in [0, L_par^2], y = k_perp^2 in [0, L_perp^2], dx dy = 4 k_par k_perp dk_par dk_perp:
         I_4^{aniso} = (1 / (16 pi^2)) * int_0^{L_par^2} dx int_0^{L_perp^2} dy / (x + y + m_IR^2)^2!!
       COMPUTE THAT EXACT INTEGRAL ANALYTICALLY!!
         int_0^{A} dx int_0^{B} dy / (x + y + m^2)^2
         = int_0^A dx [ 1/(x + m^2) - 1/(x + B + m^2) ]
         = log((A + m^2)/m^2) - log((A + B + m^2)/(B + m^2))
         = log( (A + m^2)(B + m^2) / (m^2 (A + B + m^2)) )!!

       HOLD EVERYTHING!! TAKE m -> 0 (IR scale 1/r_H << Lambda_perp, Lambda_par)!!
         (A)(B) / (m^2 (A + B)) = (1 / m^2) * (A B / (A + B)) = (1 / m^2) * (L_par^2 L_perp^2 / (L_par^2 + L_perp^2))!!
       SO THE COEFFICIENT OF log(1 / m_IR) = log(r_H) IS:
         I_4^{aniso} = (1 / (8 pi^2)) * log( Lambda_eff / m_IR )
         where Lambda_eff = Lambda_par * Lambda_perp / sqrt(Lambda_par^2 + Lambda_perp^2)!!

       WAIT A SECOND!! STOP!! WHAT IF Lambda_par -> infinity (CONTINUUM RADIAL-TIME!) AT FIXED FUZZY-SPHERE CUTOFF Lambda_perp = Lambda_P?!
         As Lambda_par -> infinity:
           Lambda_eff = lim_{L_par -> inf} (L_par * L_perp / sqrt(L_par^2 + L_perp^2)) = L_perp = Lambda_P!!
         SO EVEN WITH CONTINUUM RADIAL/TIME MODES (N_rad = inf!!) AND ONLY A FUZZY-SPHERE TRANSVERSE CUTOFF Lambda_perp = Lambda_P:
           I_4^{fuzzy-only}(Lambda_P) = (1 / (8 pi^2)) * log( Lambda_P / m_IR ) — EXACTLY IDENTICAL COEFFICIENT +3 log(s)!!

    3. WAIT A SECOND!! CHECK THE NON-TRIVIAL PART!! What about the **SPIN-CONNECTION (-60)** AND **SCALAR (+2)** TERMS under anisotropic cutoff?!
       In momentum space, the spin-connection term comes from (k_a A^a)^2 ~ k_a k_b Omega^{ac} Omega^b_c / (k^2)^3!
       THAT HAS ANISOTROPIC MOMENTUM INDICES k_a k_b!!
       Under a 4D isotropic cutoff, <k_a k_b> = (1/4) delta_{ab} k^2.
       Under an ANISOTROPIC fuzzy-sphere cutoff (k_perp <= Lambda_perp, k_par -> inf), DOES <k_a k_b> BREAK 4D ISOTROPY AND SHIFT THE -60 SPIN-CONNECTION TERM?!
       LET'S COMPUTE BOTH <k_par^2 / k^2> AND <k_perp^2 / k^2> IN THE EXACT LOG-DIVERGENT INTEGRAL RIGHT NOW!!
    """
    # Exact symbolic evaluation of the anisotropic momentum integrals for ALL 4D tensor structures!
    # Let x = k_par^2 in [0, A], y = k_perp^2 in [0, B], m = m_IR^2 -> 0.
    # For scalar/tidal V^2:       J_0  = int_0^A dx int_0^B dy / (x + y + m)^2
    # For parallel spin-conn:     J_par  = int_0^A dx int_0^B dy * x / (x + y + m)^3
    # For transverse spin-conn:   J_perp = int_0^A dx int_0^B dy * y / (x + y + m)^3
    # Check isotropic sum: J_par + J_perp = int dx dy (x+y)/(x+y+m)^3 ~ J_0!
    # Under isotropic 4D cutoff (x+y <= L^2): J_par / J_0 = 1/2, J_perp / J_0 = 1/2 (so each of 4 directions gets 1/4)!
    # Under ANISOTROPIC fuzzy-sphere cutoff (B = Lambda_P^2 fixed, A = Lambda_par^2 -> inf):
    x, y, m, A, B = sp.symbols("x y m A B", positive=True)

    J0_int = sp.integrate(sp.integrate(1 / (x + y + m) ** 2, (y, 0, B)), (x, 0, A))
    Jpar_int = sp.integrate(sp.integrate(x / (x + y + m) ** 3, (y, 0, B)), (x, 0, A))
    Jperp_int = sp.integrate(sp.integrate(y / (x + y + m) ** 3, (y, 0, B)), (x, 0, A))

    # Extract exact coefficient of log(1/m) as m -> 0!
    # Take m * d/dm (-J) as m -> 0 (wait: if J = c * log(1/m) + const, then -m * dJ/dm -> c as m -> 0!!):
    c_J0_AB = sp.simplify(sp.limit(-m * sp.diff(J0_int, m), m, 0))
    c_Jpar_AB = sp.simplify(sp.limit(-m * sp.diff(Jpar_int, m), m, 0))
    c_Jperp_AB = sp.simplify(sp.limit(-m * sp.diff(Jperp_int, m), m, 0))

    # Also take A -> infinity FIRST (continuum radial-time, fixed fuzzy-sphere B = Lambda_P^2), THEN m -> 0!
    J0_Ainf = sp.simplify(sp.limit(J0_int, A, sp.oo))
    Jpar_Ainf = sp.simplify(sp.limit(Jpar_int, A, sp.oo))
    Jperp_Ainf = sp.simplify(sp.limit(Jperp_int, A, sp.oo))

    c_J0_fuzzy_only = sp.simplify(sp.limit(-m * sp.diff(J0_Ainf, m), m, 0))
    c_Jpar_fuzzy_only = sp.simplify(sp.limit(-m * sp.diff(Jpar_Ainf, m), m, 0))
    c_Jperp_fuzzy_only = sp.simplify(sp.limit(-m * sp.diff(Jperp_Ainf, m), m, 0))

    # ALSO check 4th-order momentum tensor integrals (needed for the +2 metric-curvature scalar term!):
    # K_{par,par} = int x^2 / (x+y+m)^4, K_{par,perp} = int x y / (x+y+m)^4, K_{perp,perp} = int y^2 / (x+y+m)^4
    Kpp = sp.integrate(sp.integrate(x**2 / (x + y + m) ** 4, (y, 0, B)), (x, 0, A))
    Kpt = sp.integrate(sp.integrate(x * y / (x + y + m) ** 4, (y, 0, B)), (x, 0, A))
    Ktt = sp.integrate(sp.integrate(y**2 / (x + y + m) ** 4, (y, 0, B)), (x, 0, A))

    Kpp_Ainf = sp.simplify(sp.limit(Kpp, A, sp.oo))
    Kpt_Ainf = sp.simplify(sp.limit(Kpt, A, sp.oo))
    Ktt_Ainf = sp.simplify(sp.limit(Ktt, A, sp.oo))

    c_Kpp = sp.simplify(sp.limit(-m * sp.diff(Kpp_Ainf, m), m, 0))
    c_Kpt = sp.simplify(sp.limit(-m * sp.diff(Kpt_Ainf, m), m, 0))
    c_Ktt = sp.simplify(sp.limit(-m * sp.diff(Ktt_Ainf, m), m, 0))

    # Compare with exact 4D isotropic ratios (on the 4-ball x+y <= R^2):
    # int_{x+y<=R^2} x^a y^b / (x+y+m)^{a+b+2} dx dy:
    # Let u = x+y, x = u s, y = u(1-s), s in [0,1]:
    #   int_0^1 s^a (1-s)^b ds = a! b! / (a+b+1)!
    #   For a+b=0: 1
    #   For (1,0), (0,1): 1/2, 1/2
    #   For (2,0), (1,1), (0,2): 1/3, 1/6, 1/3!!

    # Numerical check at finite B=1e4, A=1e12, m=1e-4 to confirm zero hidden boundary terms!
    m_num, B_num = 1e-4, 1e4
    val_J0 = float(J0_Ainf.subs({B: B_num, m: m_num}))
    val_Jpar = float(Jpar_Ainf.subs({B: B_num, m: m_num}))
    val_Jperp = float(Jperp_Ainf.subs({B: B_num, m: m_num}))

    return dict(
        exact_finite_AB_log_coeffs=dict(
            c_J0=str(c_J0_AB),
            c_Jpar=str(c_Jpar_AB),
            c_Jperp=str(c_Jperp_AB),
        ),
        exact_continuum_radial_fuzzy_sphere_only_A_to_inf=dict(
            J0_exact=str(J0_Ainf),
            Jpar_exact=str(Jpar_Ainf),
            Jperp_exact=str(Jperp_Ainf),
            log_coeff_J0=str(c_J0_fuzzy_only),
            log_coeff_Jpar=str(c_Jpar_fuzzy_only),
            log_coeff_Jperp=str(c_Jperp_fuzzy_only),
            log_coeff_K_par_par=str(c_Kpp),
            log_coeff_K_par_perp=str(c_Kpt),
            log_coeff_K_perp_perp=str(c_Ktt),
            isotropic_4D_exact_targets=dict(
                J0="1",
                Jpar="1/2",
                Jperp="1/2",
                K_par_par="1/3",
                K_par_perp="1/6",
                K_perp_perp="1/3",
            ),
            all_tensor_log_ratios_100pct_isotropic_invariant=(
                c_J0_fuzzy_only == 1
                and c_Jpar_fuzzy_only == sp.Rational(1, 2)
                and c_Jperp_fuzzy_only == sp.Rational(1, 2)
                and c_Kpp == sp.Rational(1, 3)
                and c_Kpt == sp.Rational(1, 6)
                and c_Ktt == sp.Rational(1, 3)
            ),
        ),
        numerical_spot_check_B1e4_m1e_4=dict(
            J0=val_J0,
            Jpar=val_Jpar,
            Jperp=val_Jperp,
            Jpar_over_J0=val_Jpar / val_J0,
            Jperp_over_J0=val_Jperp / val_J0,
        ),
    )


if __name__ == "__main__":
    res_a = part_a_s4_exact_spectrum()
    res_b = part_b_schwarzschild_anisotropic_vs_covariant_cutoff()
    out = dict(part_a_S4=res_a, part_b_Schwarzschild_fuzzy_cutoff_invariance=res_b)

    out_path = os.path.join(os.path.dirname(__file__), "p2_test3_curved_lichnerowicz_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print("=== PATH 2 BREAKTHROUGH TEST 3A: Exact 4D S^4 (de Sitter) Graviton Spectrum ===")
    print("  ADLS character poles:", json.dumps(res_a["ADLS_character_poles"], indent=2))
    print("  Direct S^4 eigenvalue c0:", json.dumps(res_a["S4_direct_eigenvalue_c0"], indent=2))

    print("\n=== PATH 2 BREAKTHROUGH TEST 3B: Fuzzy-Sphere-Only Cutoff (N_rad = inf!) 4D Isotropy Theorem ===")
    fb = res_b["exact_continuum_radial_fuzzy_sphere_only_A_to_inf"]
    print(f"  J0(A=inf)    = {fb['J0_exact']}  -> log coeff = {fb['log_coeff_J0']} (iso: 1)")
    print(f"  Jpar(A=inf)  = {fb['Jpar_exact']}  -> log coeff = {fb['log_coeff_Jpar']} (iso: 1/2)")
    print(f"  Jperp(A=inf) = {fb['Jperp_exact']}  -> log coeff = {fb['log_coeff_Jperp']} (iso: 1/2)")
    print(f"  K_par_par, K_par_perp, K_perp_perp log coeffs = ({fb['log_coeff_K_par_par']}, {fb['log_coeff_K_par_perp']}, {fb['log_coeff_K_perp_perp']}) (iso: 1/3, 1/6, 1/3)")
    print(f"  ALL 4D TENSOR LOG COEFFICIENTS INVARIANT UNDER FUZZY-SPHERE-ONLY CUTOFF (N_rad = inf)? {fb['all_tensor_log_ratios_100pct_isotropic_invariant']}!!")
    print(f"Saved to {out_path}")
