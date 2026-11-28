#!/usr/bin/env python3
"""Second angular-shift preservation at Paper I's vacuum witness.

Verifies the geometry and exact contraction in a differentiated continuum
Noether identity. The analytic identity and its hypotheses are documented in
the accompanying report; this is not a generic Dirac classification.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import subprocess
import time

import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[2]
PREVIOUS = ROOT / 'verification/dirac_stabilization_20261127/stabilize_witness.py'
spec = importlib.util.spec_from_file_location('first_preservation', PREVIOUS)
first = importlib.util.module_from_spec(spec)
spec.loader.exec_module(first)
x = first.x
OUT = ROOT / 'reports/higher-order-consistency-2026-11-27'


def lie_metric(metric, vector, coords):
    return s.Matrix(3, 3, lambda i, j: s.factor(sum(
        vector[k] * s.diff(metric[i, j], coords[k])
        + metric[k, j] * s.diff(vector[k], coords[i])
        + metric[i, k] * s.diff(vector[k], coords[j])
        for k in range(3))))


def angular_defect():
    """Check full reconstruction and its changing tangent with arbitrary w(r)."""
    r, t, phi, kappa = s.symbols('r t phi kappa', real=True)
    w = s.Function('w')(r)
    d = 1 - x*x
    xi = 3*w*x*d
    lam = 1 + 2*kappa*t*x*x
    gxx = r*r/d + 2*kappa*t*r*r
    gpp = r*r*d
    V = 2*kappa*t*x*d/r
    sigma = s.Matrix([[lam + gxx*V*V, gxx*V, 0],
                      [gxx*V, gxx, 0], [0, 0, gpp]])
    delta_lam = xi*s.diff(lam, x)
    delta_gxx = xi*s.diff(gxx, x) + 2*gxx*s.diff(xi, x)
    delta_gpp = xi*s.diff(gpp, x)
    delta_V = s.diff(xi, r) + xi*s.diff(V, x) - V*s.diff(xi, x)

    def tangent(dl, dg, dp, dv):
        return s.Matrix([[dl + dg*V*V + 2*gxx*V*dv, dg*V + gxx*dv, 0],
                         [dg*V + gxx*dv, dg, 0], [0, 0, dp]])

    full_variation = tangent(delta_lam, delta_gxx, delta_gpp, delta_V)
    geometric_variation = lie_metric(sigma, [0, xi, 0], [r, x, phi])
    assert all(s.factor(z) == 0 for z in full_variation-geometric_variation)
    zero_tests = [z.subs(t, 0) for z in [delta_lam, delta_gxx, delta_gpp, delta_V]]
    tangent_dot = tangent(*zero_tests).diff(t).subs(t, 0)
    sigma_dot = sigma.diff(t).subs(t, 0)
    U = (lie_metric(sigma_dot, [0, xi, 0], [r, x, phi]) - tangent_dot).applyfunc(s.factor)
    coefficient_route = tangent(*[s.diff(z, t).subs(t, 0)
                                  for z in [delta_lam, delta_gxx, delta_gpp, delta_V]]).subs(t, 0)
    assert all(s.factor(z) == 0 for z in U-coefficient_route)
    assert all(not z.has(s.diff(w, r)) for z in U)
    # Convention: f, g, v_theta/sqrt(1-x^2), c/(1-x^2).
    components = [U[0, 0],
                  (d*U[1, 1]/r**2 + U[2, 2]/(r**2*d))/2,
                  -U[0, 1]/r,
                  (d*U[1, 1]/r**2 - U[2, 2]/(r**2*d))/(2*d)]
    unit = [s.factor(z/(kappa*w)) for z in components]
    expected = [12*x*x*d, 6*d*(1-3*x*x), 0, 6*(1-3*x*x)]
    assert all(s.expand(a-b) == 0 for a, b in zip(unit, expected))
    assert s.factor(s.diff(delta_V, t)) == 0
    return unit, {
        'descriptor': 'xi^x = 3*w(r)*x*(1-x^2); xi^r=xi^phi=0',
        'coordinate_tangent_defect': [[str(z) for z in row] for row in U.tolist()],
        'unit_component_defect': [str(z) for z in unit],
        'coefficient_transform_matches_full_Lie_derivative': True,
        'two_derivations_of_defect_agree': True,
        'radial_profile_derivatives_cancel': True,
        'coefficient_V_variation': str(s.factor(delta_V)),
    }


def adm_chain_rule():
    """Check the initial zero coefficients in the differentiated EL chain rule."""
    a, b, lam, V, W = s.symbols('a b lambda V W', real=True)
    alpha = s.sqrt(-a+b*b/lam)
    beta = [b/lam, W-b*V/lam]
    witness = {a: -1, b: 0, W: 0}
    spatial_jacobians = [s.diff(z, q).subs(witness)
                         for z in [alpha, *beta] for q in [lam, V]]
    assert all(z == 0 for z in spatial_jacobians)
    assert alpha.subs(witness) == 1
    assert all(z.subs(witness) == 0 for z in beta)
    return {'alpha_at_fixed_temporal_witness': '1 for every spatial coefficient',
            'beta_at_fixed_temporal_witness': '0 for every spatial coefficient',
            'spatial_derivatives_of_alpha_beta_vanish': True,
            'other_chain_rule_terms': 'multiply initially zero full continuum constraints'}


def verify():
    start = time.monotonic()
    U, geometry = angular_defect()
    chain = adm_chain_rule()
    all_m = first.prior.all_m_hessian()
    rows = []
    for L in [2, 3]:
        mu = s.Rational(32, 9)*s.legendre(2, x)
        R = first.Q([b+mu*k for b, k in zip(first.B, first.K)], L)
        # The first preservation equations really hold at the point tested.
        for ell in range(L+1):
            T = first.Q([s.legendre(ell, x)*k for k in first.K], L)
            assert first.pair(T, R) == 0
        scalar = next(z for z in all_m if z['L'] == L)
        kernel = [label for label, norm in zip(scalar['labels'], scalar['discarded_norm_diagonal'])
                  if s.sympify(norm) == 0]
        assert len(kernel) == (L-1)**2
        assert all(ell <= L-2 for ell, m in kernel)
        assert len(scalar['labels'])-len(kernel) == 4*L
        QU = first.Q(U, L)
        contraction = first.pair(R, U)
        assert first.pair(R, QU) == contraction
        assert contraction == -s.Rational(1024, 1575)
        # E_sigma = -G_DW(R, .)/2; dot C_W = -E_sigma[Lie_xi sigma]
        # at alpha=1, beta=0. Differentiating restricted spatial E=0 gives
        # ddot C_W = G_DW(R,U)/2; see report for discarded-constraint terms.
        obstruction = contraction/2
        assert obstruction == -s.Rational(512, 1575)
        integrand = s.factor(first.pair_integrand(R, U))
        integrand_fun = s.lambdify(x, integrand, 'numpy')
        quadrature = []
        for n in [24, 48, 96]:
            points, weights = np.polynomial.legendre.leggauss(n)
            measured = float(weights@integrand_fun(points)/4)
            error = abs(measured-float(obstruction))
            assert error < 1e-12
            quadrature.append({'nodes': n, 'C_W_second_derivative_coefficient': measured,
                               'absolute_error': error})
        rows.append({'N': L+1, 'L': L,
                     'first_secondary_preservation_solvable': True,
                     'lapse_rate_kernel_all_m': kernel,
                     'normal_acceleration_residue': [str(z) for z in R],
                     'discarded_changing_shift_tangent': [str(z) for z in QU],
                     'unit_DeWitt_integrand': str(integrand),
                     'unit_DeWitt_contraction': str(contraction),
                     'second_preservation_coefficient': str(obstruction),
                     'radial_and_curvature_factor': 'kappa^3 * integral r^2*w(r) dr',
                     'depends_on_free_first_temporal_rates': False,
                     'depends_on_second_temporal_rates': False,
                     'smooth_projected_solution_through_fixed_witness': False,
                     'quadrature_checks': quadrature})
    result = {
        'normalization': '(16piG)/(4pi) times Einstein-Hilbert action; angular average dOmega/(4pi)',
        'geometry': geometry, 'ADM_chain_rule': chain, 'checks': rows,
        'meaning': 'No sufficiently smooth retained Euler-Lagrange trajectory through this fixed witness at N=3,4.',
        'regularity': 'enough time and radial derivatives for the twice differentiated secondary and differentiated spatial equations',
        'final_constraint_surface_classified': False,
        'generic_physical_DOF_count_established': False,
        'all_initial_data_excluded': False,
        'full_nonlinear_goal_complete': False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'checks.json').write_text(json.dumps(result, indent=2)+'\n')
    sources = [Path(__file__), PREVIOUS, first.PRIOR, ROOT/'sections/nonlinear_action.tex']
    manifest = {
        'base_commit_at_execution': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'python': platform.python_version(), 'platform': platform.platform(),
        'dependencies': {k: importlib.metadata.version(k) for k in ['numpy', 'sympy']},
        'seconds': round(time.monotonic()-start, 3), 'GPU_jobs': 0, 'experiment_launches': 0,
        'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    }
    (OUT/'reproducibility.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print('PASS: exact angular diffeomorphism of full coefficient reconstruction; arbitrary radial profile retained.')
    print('PASS: geometric and coefficient derivatives give the same changing tangent; radial-profile derivatives cancel.')
    print('PASS: first-preservation equations and all-m homogeneous lapse-rate kernel checked at N=3,4.')
    for row in rows:
        print(f"N={row['N']}: C_W second derivative = -512/1575 * kappa^3 * integral r^2*w(r) dr.")
    print('This is nonzero for kappa != 0 and a nonnegative, nonzero compact radial bump.')
    print('Maximum independent angular quadrature error:', max(q['absolute_error'] for row in rows for q in row['quadrature_checks']))
    print('CONCLUSION: the specified witness is excluded by second preservation; no smooth retained trajectory through it.')
    print('Other initial data, final Dirac classification, physical DOF, and the full nonlinear objective remain open.')


if __name__ == '__main__':
    verify()
