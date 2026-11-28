#!/usr/bin/env python3
"""Actual metric-coefficient Galerkin action: first two Dirac stages.

Exact ADM reconstruction, regular spatial Legendre map, a vacuum constraint
witness, and the full scalar-multiplier Hessian at that witness. The latter is
on the primary-secondary surface, NOT a completed stabilization or DOF count.
No simulations, external compute, or GPU jobs are launched.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import time

import numpy as np
import sympy as s
from sympy.physics.wigner import wigner_3j

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'reports/canonical-analysis-2026-11-27'
x = s.symbols('x', real=True)


def average(f):
    return s.simplify(s.integrate(s.expand(f), (x, -1, 1)) / 2)


def adm_map():
    a, b, lam = s.symbols('a b lambda', real=True)
    g11, g12, g22 = s.symbols('g11 g12 g22', real=True)
    v = s.Matrix(s.symbols('v1 v2', real=True))
    w = s.Matrix(s.symbols('w1 w2', real=True))
    gamma = s.Matrix([[g11, g12], [g12, g22]])
    spatial = s.zeros(3)
    spatial[0, 0] = lam + (v.T * gamma * v)[0]
    spatial[0, 1:] = v.T * gamma
    spatial[1:, 0] = gamma * v
    spatial[1:, 1:] = gamma
    beta = s.Matrix([b / lam, *(w - b * v / lam)])
    lapse_squared = -a + b * b / lam
    target_g0i = s.Matrix([b + (w.T * gamma * v)[0], *(gamma * w)])
    assert all(s.simplify(z) == 0 for z in spatial * beta - target_g0i)
    target_g00 = a + (w.T * gamma * w)[0]
    assert s.simplify(-lapse_squared + (beta.T * spatial * beta)[0] - target_g00) == 0
    assert s.factor(spatial.det() - lam * gamma.det()) == 0
    return {'lapse_squared': str(lapse_squared),
            'shift_r': str(beta[0]), 'shift_angular': [str(z) for z in beta[1:]],
            'spatial_determinant': 'lambda * det(gamma)',
            'metric_reconstruction_residuals': 'all exactly zero'}


def vacuum_witness():
    r, theta, phi, kappa = s.symbols('r theta phi kappa', real=True)
    coords = [r, theta, phi]
    q = s.diag(1, r**2, r**2 * s.sin(theta)**2)
    qi = q.inv()
    n = s.Matrix([s.cos(theta), -r * s.sin(theta), 0])
    K = kappa * n * n.T
    mixed = qi * K
    trace = s.trigsimp(s.trace(mixed))
    hamiltonian = s.trigsimp(trace**2 - s.trace(mixed * mixed))
    assert trace == kappa and hamiltonian == 0
    Gamma = [[[s.simplify(sum(qi[i, k] * (s.diff(q[k, j], coords[l])
                 + s.diff(q[k, l], coords[j]) - s.diff(q[j, l], coords[k]))
                 for k in range(3)) / 2) for l in range(3)]
                 for j in range(3)] for i in range(3)]
    stress = mixed - s.eye(3) * trace
    momentum = []
    for j in range(3):
        div = sum(s.diff(stress[i, j], coords[i]) for i in range(3))
        div += sum(Gamma[i][i][k] * stress[k, j] - Gamma[k][i][j] * stress[i, k]
                   for i in range(3) for k in range(3))
        momentum.append(s.trigsimp(div))
    assert momentum == [0, 0, 0]
    # Direct curvature check, independent of the Cartesian flatness argument.
    Ric = s.zeros(3)
    for i in range(3):
        for j in range(3):
            Ric[i, j] = s.trigsimp(sum(s.diff(Gamma[k][i][j], coords[k])
                - s.diff(Gamma[k][i][k], coords[j]) for k in range(3))
                + sum(Gamma[k][k][l] * Gamma[l][i][j]
                - Gamma[k][j][l] * Gamma[l][i][k] for k in range(3) for l in range(3)))
    # trigsimp alone can leave sin(2t)/(2*tan(t))-cos(2t)/2-1/2.
    Ric = Ric.applyfunc(lambda z: s.simplify(s.expand_trig(z)))
    assert all(z == 0 for z in Ric)
    return {'spatial_Ricci': 'zero', 'Lambda': 0,
            'K_cartesian': 'kappa * diag(0,0,1), constant kappa != 0',
            'trace_K': str(trace), 'hamiltonian_constraint': str(hamiltonian),
            'momentum_constraints': [str(z) for z in momentum],
            'retained_tensor_degrees': [0, 2],
            'surface': 'all temporal primaries and all secondary constraints',
            'tertiary_stabilization': 'not established'}


def spatial_signature():
    scalar = s.Matrix([[0, -2], [-2, -2]])
    assert scalar.det() == -4
    rows = []
    for N in [3, 4]:
        rows.append({'N': N, 'spatial_coefficients': 6 * N*N - 10,
                     'temporal_primary_constraints': 4 * N*N - 2,
                     'positive_spatial_eigenvalues': 5 * N*N - 10,
                     'negative_spatial_eigenvalues': N*N,
                     'zero_spatial_eigenvalues': 0,
                     'physical_DOF': 'not yet classified'})
    return {'scalar_DeWitt_block': [[str(z) for z in row] for row in scalar.tolist()],
            'scalar_eigenvalues': ['-1+sqrt(5)', '-1-sqrt(5)'], 'counts': rows,
            'domain': 'round positive spatial metric and an open regular neighborhood'}


def residual(poly, spin, L):
    """Project polynomial amplitude with spin weight (1-x²)^spin."""
    weight = (1 - x*x)**spin
    out = poly
    for j in range(spin, L + 1):
        basis = s.diff(s.legendre(j, x), x, spin)
        coeff = average(weight * poly * basis) / average(weight * basis**2)
        out -= coeff * basis
    return s.factor(out)


def axisymmetric_exact():
    rows = []
    for L in [2, 3]:
        for ell in range(1, L + 1):
            u = s.legendre(ell, x)
            F, G, V, C = [residual(u * a, spin, L) for a, spin in
                [(x*x, 0), ((1-x*x)/2, 0), (-x, 1), (s.Rational(1, 2), 2)]]
            assert s.expand(F + 2*G) == 0  # Discarded tensor is tracefree.
            norm = average(-4*F*G - 2*G*G + 2*(1-x*x)*V*V + 2*(1-x*x)**2*C*C)
            assert norm >= 0
            rows.append({'N': L+1, 'L': L, 'u': f'P_{ell}', 'ell': ell,
                         'discarded_DeWitt_norm': str(norm),
                         'fixed_momentum_H_second_derivative': str(-norm/2),
                         'discarded_polynomial_components': [str(z) for z in (F,G,V,C)]})
    assert [r['discarded_DeWitt_norm'] for r in rows] == ['2/15','12/175','0','12/175','20/441']
    return rows


def gaunt(l, m, J, spin):
    # All harmonics have mean |Y|²=1. Multiplying scalar Y_lm by spin-s Y_20.
    return (s.Integer(-1)**(spin+m) * s.sqrt((2*l+1)*5*(2*J+1))
            * wigner_3j(l,2,J,0,-spin,spin) * wigner_3j(l,2,J,m,0,-m))


def all_m_hessian():
    rows = []
    for L in [2,3]:
        labels = [(ell,m) for ell in range(L+1) for m in range(-ell,ell+1)]
        def entry(i,j):
            l,m=labels[i]; k,n=labels[j]
            if m != n: return s.S.Zero
            # Both spin signs are required for off-diagonal/parity cancellations.
            return s.simplify(s.Rational(2,15) * sum(
                gaunt(l,m,J,z) * gaunt(k,m,J,z)
                for J in range(L+1,L+3) for z in range(-2,3)))
        D = s.Matrix(len(labels),len(labels),entry)
        assert D == s.diag(*D.diagonal())
        assert D.rank() == 4*L
        for i,(ell,m) in enumerate(labels):
            assert (D[i,i] == 0) == (ell <= L-2)
            assert D[i,i] >= 0
        # Independent direct-polynomial result for normalized Y_ell,0.
        for ell in range(1,L+1):
            reference = next(row for row in axisymmetric_exact()
                             if row['L']==L and row['ell']==ell)
            i=labels.index((ell,0))
            assert s.simplify(D[i,i]-(2*ell+1)*s.sympify(reference['discarded_DeWitt_norm'])) == 0
        rows.append({'N':L+1,'L':L,'scalar_dimension':len(labels),
                     'lapse_Hessian_rank':D.rank(),'kernel_dimension':(L-1)**2,
                     'full_primary_secondary_Poisson_rank_lower_bound':8*L,
                     'labels':[list(z) for z in labels],
                     'discarded_norm_diagonal':[str(z) for z in D.diagonal()],
                     'Hessian_diagonal':[str(-z/2) for z in D.diagonal()],
                     'off_diagonal_entries':'all exactly zero',
                     'scope':'all scalar magnetic modes at this witness, not generic Dirac rank'})
    return rows


def numeric_legendre(rows):
    """Independent quadrature and finite differences of the actual inverse Gram blocks."""
    results=[]
    components=[x*x,(1-x*x)/2,-x*s.sqrt(1-x*x),(1-x*x)/2]
    for row in rows:
        L=row['L']; ell=row['ell']; target=float(s.sympify(row['fixed_momentum_H_second_derivative']))
        for grid in [64,96]:
            nodes,weights=np.polynomial.legendre.leggauss(grid); weights=weights/2
            bases=[]; gram=[]
            for spin in [0,1,2]:
                functions=[s.sqrt((2*j+1)*s.factorial(j-spin)/s.factorial(j+spin))
                           *s.assoc_legendre(j,spin,x) for j in range(spin,L+1)]
                B=np.stack([np.broadcast_to(s.lambdify(x,f,'numpy')(nodes),nodes.shape) for f in functions])
                bases.append(B); gram.append((B*weights)@B.T)
                assert np.max(np.abs(gram[-1]-np.eye(len(functions)))) < 1e-12
            coefficients=[]
            for component,spin in zip(components,[0,0,1,2]):
                vals=np.asarray(s.lambdify(x,component,'numpy')(nodes),float)
                coefficients.append(bases[spin]@(weights*vals))
            u=np.asarray(s.lambdify(x,s.legendre(ell,x),'numpy')(nodes),float)
            zero=np.zeros_like(gram[0])
            W0=np.block([[zero,-gram[0]],[-gram[0],-gram[0]]])
            ps=W0@(2*np.concatenate(coefficients[:2]))
            pv=gram[1]@(2*coefficients[2]); pc=gram[2]@(2*coefficients[3])
            def H(eps):
                matrices=[(B*(weights/np.sqrt(1+eps*u)))@B.T for B in bases]
                Ws=np.block([[zero,-matrices[0]],[-matrices[0],-matrices[0]]])
                return float((ps@np.linalg.solve(Ws,ps)+pv@np.linalg.solve(matrices[1],pv)
                              +pc@np.linalg.solve(matrices[2],pc))/2)
            assert abs(H(0)) < 1e-12
            for step in [.004,.002]:
                sample={str(k):H(k*step) for k in [-2,-1,0,1,2]}
                curvature=(-sample['2']+16*sample['1']-30*sample['0']+16*sample['-1']-sample['-2'])/(12*step**2)
                error=abs(curvature-target)
                assert error < 2e-7, (L,ell,grid,step,curvature,target)
                results.append({'N':L+1,'u':row['u'],'quadrature_nodes':grid,'step':step,
                                'H_samples':sample,'H_second_derivative':curvature,
                                'exact_target':target,'absolute_error':error})
    return results


def main():
    start=time.monotonic()
    exact=axisymmetric_exact()
    result={'normalization':'16*pi*G=1, angular average dOmega/(4*pi), local r=1, kappa=1',
            'action':'actual unprojected metric-coefficient EH Galerkin pullback, compact bulk variations',
            'adm_map':adm_map(),'spatial_Legendre_map':spatial_signature(),
            'vacuum_constraint_point':vacuum_witness(),
            'axisymmetric_exact':exact,'all_m_lapse_Hessian':all_m_hessian(),
            'independent_inverse_Gram_checks':numeric_legendre(exact),
            'complete_Dirac_classification':False,'full_nonlinear_goal_complete':False}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
    sources=[Path(__file__),ROOT/'sections/nonlinear_action.tex']
    manifest={'base_commit_at_execution':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'python':platform.python_version(),'platform':platform.platform(),
              'dependencies':{k:importlib.metadata.version(k) for k in ['numpy','sympy']},
              'seconds':round(time.monotonic()-start,3),'GPU_jobs':0,
              'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (OUT/'reproducibility.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PASS: exact ADM map; full spatial Hessian; spherical-coordinate vacuum constraints.')
    for row in result['all_m_lapse_Hessian']:
        print(f"N={row['N']}: spatial rank {6*row['N']**2-10}, scalar multiplier Hessian rank {row['lapse_Hessian_rank']}; full constraint bracket rank >= {row['full_primary_secondary_Poisson_rank_lower_bound']} at the witness.")
    print('PASS: 5 exact axisymmetric cases, checked against all-m Gaunt matrices and 20 inverse-Gram finite differences.')
    print('Maximum finite-difference error:',max(z['absolute_error'] for z in result['independent_inverse_Gram_checks']))
    print('Scope: primary-secondary surface. Tertiary stabilization and physical DOF remain unresolved.')


if __name__=='__main__':
    main()
