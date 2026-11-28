#!/usr/bin/env python3
"""Check the replacement Einstein-Hilbert angular Galerkin prescription.

The finite matrices store harmonic coefficients; nonlinear operations act on
reconstructed functions. This is not the legacy fuzzy matrix-product action.
The general curvature routine accepts an arbitrary symbolic metric. The test
family is a finite-amplitude, angularly anisotropic conformal Schwarzschild
metric with band-limited metric coefficients, not just a quadratic expansion.
"""
from pathlib import Path
import json
import math
import sys
import numpy as np
import sympy as sp
from scipy.special import sph_harm_y

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
from p2_frontier1_sun_fuzzy_action_and_saddles import build_fuzzy_monopole_basis, su2_generators

OUT = HERE.parents[2] / 'reports/repair-audit-2026-11-27'


def reduce(expr):
    expr = sp.expand_trig(sp.simplify(expr))
    expr = expr.replace(lambda z: z.func == sp.tan,
                        lambda z: sp.sin(z.args[0])/sp.cos(z.args[0]))
    expr = expr.replace(lambda z: z.func == sp.cot,
                        lambda z: sp.cos(z.args[0])/sp.sin(z.args[0]))
    return sp.factor(sp.trigsimp(sp.cancel(expr)))


def ricci_scalar(metric, coords):
    """Direct Levi-Civita calculation; no field equations or target R inserted."""
    n = len(coords)
    inv = metric.inv()
    gamma = [[[reduce(sum(inv[a, d] * (
        sp.diff(metric[d, c], coords[b]) + sp.diff(metric[d, b], coords[c])
        - sp.diff(metric[b, c], coords[d])) for d in range(n)) / 2)
        for c in range(n)] for b in range(n)] for a in range(n)]
    ric = sp.Matrix(n, n, lambda b, d: reduce(sum(
        sp.diff(gamma[a][b][d], coords[a]) - sp.diff(gamma[a][a][d], coords[b])
        + sum(gamma[a][a][e] * gamma[e][b][d]
              - gamma[a][b][e] * gamma[e][a][d] for e in range(n))
        for a in range(n))))
    return reduce(sum(inv[a, b] * ric[a, b] for a in range(n) for b in range(n)))


def scalar_symbol(matrix, basis, theta, phi):
    """Isometric reconstruction in harmonics normalized by int dOmega/(4 pi)."""
    N = matrix.shape[0]
    return sum((np.trace(T.conj().T @ matrix) / N)
               * math.sqrt(4 * math.pi) * sph_harm_y(j, m, theta, phi)
               for (j, m), T in basis.items())


def main():
    t, r, theta, phi, mass = sp.symbols('t r theta phi M', positive=True)
    A = sp.Function('A')(theta)
    f = 1 - 2 * mass / r
    metric = sp.diag(-A * f, A / f, A * r**2, A * r**2 * sp.sin(theta)**2)
    R = ricci_scalar(metric, (t, r, theta, phi))
    lapA = sp.diff(A, theta, 2) + sp.cot(theta) * sp.diff(A, theta)
    expected = (-3 * lapA / A**2 + sp.Rational(3, 2) * sp.diff(A, theta)**2 / A**3) / r**2
    assert reduce(R - expected) == 0
    assert reduce(R.subs(A, 1).doit()) == 0
    # sqrt(-g) R / (sin(theta)) = -3 Lap(A) + 3 |grad A|^2/(2 A).
    density = reduce(A**2 * r**2 * R)
    assert reduce(density + 3 * lapA - sp.Rational(3, 2) * sp.diff(A, theta)**2 / A) == 0
    # A=1+epsilon Y20 lies in every scalar band N>=3. Products/inverses do not.
    x = sp.symbols('x', real=True)
    Y = sp.sqrt(5) * (3 * x*x - 1) / 2
    gradY2 = (1 - x*x) * sp.diff(Y, x)**2
    quadratic = sp.integrate(sp.Rational(3, 4) * gradY2, (x, -1, 1))
    cubic = sp.integrate(-sp.Rational(3, 4) * Y * gradY2, (x, -1, 1))
    assert quadratic == 9
    assert sp.simplify(cubic + 9 * sp.sqrt(5) / 7) == 0
    rows = []
    xq, wq = np.polynomial.legendre.leggauss(160)
    y = np.sqrt(5) * (3*xq*xq - 1) / 2
    grad = 45 * xq*xq * (1-xq*xq)
    th = np.arccos(xq)
    for eps in [-0.2, 0.1, 0.3]:
        aq = 1 + eps*y
        assert np.min(aq) > 0
        # Direct Ricci density and its integration-by-parts reduction.
        direct = np.dot(wq/2, 18*eps*y + 1.5*eps**2*grad/aq)
        reduced = np.dot(wq/2, 1.5*eps**2*grad/aq)
        assert abs(direct-reduced) < 1e-12
        # First variation from the Einstein trace, -G^mu_mu=R (Lambda=0).
        einstein_variation = np.dot(wq/2, y * (18*eps*y/aq + 1.5*eps**2*grad/aq**2))
        differentiated_action = np.dot(wq/2, 1.5*grad*(2*eps/aq-eps**2*y/aq**2))
        assert abs(einstein_variation-differentiated_action) < 1e-11
        representations = []
        for N in [3, 4, 6, 8]:
            basis = build_fuzzy_monopole_basis(N, 0)['basis']
            mat = np.eye(N) + eps*basis[(2, 0)]
            symbol = scalar_symbol(mat, basis, th, np.zeros_like(th))
            err = float(np.max(abs(symbol-aq)))
            assert err < 2e-11
            # A deliberately different fuzzy-product regulator, reported only
            # to expose that nonlinear symbol reconstruction is not trace substitution.
            fuzzy_grad = np.zeros((N, N), complex)
            for J in su2_generators(N-1)[:3]:
                comm = J @ mat - mat @ J
                fuzzy_grad += comm.conj().T @ comm
            fuzzy_value = float((1.5*np.trace(np.linalg.solve(mat, fuzzy_grad))/N).real)
            representations.append({'N':N, 'symbol_error':err,
                                    'distinct_fuzzy_product_density':fuzzy_value})
        rows.append({'epsilon':eps, 'exact_angular_density':float(reduced),
                     'direct_Ricci_vs_integrated_density_error':float(abs(direct-reduced)),
                     'Einstein_first_variation_error':float(abs(einstein_variation-differentiated_action)),
                     'representations':representations})
    result = {'curvature_symbolic_residual':'0', 'Schwarzschild_R_residual':'0',
              'conformal_Ricci_scalar':str(R), 'quadratic_coefficient':str(quadratic),
              'cubic_coefficient':str(cubic), 'nonlinear_checks':rows,
              'scope':'Full angular Galerkin EH action is defined by pullback; these tests exercise one nonlinear family. No closed finite-N constraint algebra, null-boundary gluing, quantization, or HMC equivalence is proved.'}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'nonlinear_galerkin_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
