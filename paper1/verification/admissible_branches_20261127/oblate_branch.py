#!/usr/bin/env python3
"""An exact admissible Paper I family with a nonzero shift bracket.

Checks finite-band oblate Minkowski geometry, the leading exact fixed-momentum
shift Hessian, and finite-deformation tangent projections. Does not classify
the full final constraint surface or prove a regular surrounding stratum.
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
PRIOR = ROOT/'verification/dirac_stabilization_20261127/stabilize_witness.py'
spec = importlib.util.spec_from_file_location('first_checkpoint', PRIOR)
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
x = prior.x
OUT = ROOT/'reports/admissible-branches-2026-11-27'


def geometry():
    r, eps, phi = s.symbols('r epsilon phi', real=True)
    d = 1-x*x
    rho, Sigma = r*r+eps, r*r+eps*x*x
    q = s.diag(Sigma/rho, Sigma/d, rho*d)
    # Cylindrical radius squared F and Cartesian axial coordinate Z provide
    # an independent Euclidean pullback without trigonometric simplification.
    F, Z = rho*d, r*x
    coords = [r, x, phi]
    pullback = s.zeros(3)
    for i in range(2):
        for j in range(2):
            pullback[i,j] = s.factor(s.diff(F,coords[i])*s.diff(F,coords[j])/(4*F)
                                         +s.diff(Z,coords[i])*s.diff(Z,coords[j]))
    pullback[2,2] = F
    assert all(s.factor(z)==0 for z in pullback-q)
    assert s.factor(q.det()-Sigma**2)==0
    P2 = s.legendre(2,x)
    lam = (r*r+eps/3)/rho+2*eps*P2/(3*rho)
    area = r*r+2*eps/3+eps*P2/3
    Cxx, Cpp = -eps*s.diff(P2,x,2)/6, eps*d*d*s.diff(P2,x,2)/6
    assert s.factor(lam-q[0,0])==0
    assert s.factor(area/d+Cxx-q[1,1])==0
    assert s.factor(area*d+Cpp-q[2,2])==0
    assert s.factor(d*Cxx+Cpp/d)==0
    # Direct Ricci calculation independent of the explicit Euclidean map.
    qi=q.inv()
    Gamma=[[[s.factor(sum(qi[i,k]*(s.diff(q[k,j],coords[l])
        +s.diff(q[k,l],coords[j])-s.diff(q[j,l],coords[k]))
        for k in range(3))/2) for l in range(3)] for j in range(3)] for i in range(3)]
    Ric=s.Matrix(3,3,lambda i,j:s.factor(
        sum(s.diff(Gamma[k][i][j],coords[k])-s.diff(Gamma[k][i][k],coords[j]) for k in range(3))
        +sum(Gamma[k][k][l]*Gamma[l][i][j]-Gamma[k][j][l]*Gamma[l][i][k]
             for k in range(3) for l in range(3))))
    assert all(z==0 for z in Ric)
    return q, r, eps, {
        'metric_diagonal':[str(z) for z in q.diagonal()],
        'lambda_harmonics':str(lam), 'area_harmonics':str(area),
        'C_xx':str(Cxx), 'C_phiphi':str(Cpp),
        'maximum_angular_degree':2, 'minimum_N':3,
        'Euclidean_coordinate_pullback':'exact', 'spatial_Ricci':'all components exactly zero',
        'spatial_volume_density':str(Sigma),
        'temporal_coefficients':'a=-1; b=W=0; spatial V=0',
        'K_and_retained_momentum':'identically zero',
        'trajectory':'stationary, exact continuum vacuum and retained solution for all time',
        'domain':'compact shell with r>0, r^2+epsilon>0, r^2+epsilon*x^2>0',
    }


def leading_defect(q,r,eps):
    d=1-x*x
    w=s.Function('w')(r)
    xi=3*w*x*d
    U=s.Matrix([[xi*s.diff(q[0,0],x),q[1,1]*s.diff(xi,r),0],
                [q[1,1]*s.diff(xi,r),xi*s.diff(q[1,1],x)+2*q[1,1]*s.diff(xi,x),0],
                [0,0,xi*s.diff(q[2,2],x)]])
    # Mixed terms are exactly in the actual tangent via delta V^x=partial_r xi.
    Vtangent=s.Matrix([[0,q[1,1]*s.diff(xi,r),0],
                       [q[1,1]*s.diff(xi,r),0,0],[0,0,0]])
    diagonal=U-Vtangent
    assert all(not z.has(s.diff(w,r)) for z in diagonal)
    round_source=diagonal.subs(eps,0)
    linear=diagonal.diff(eps).subs(eps,0)
    def components(tensor):
        return [tensor[0,0],(d*tensor[1,1]/r**2+tensor[2,2]/(r**2*d))/2,
                -tensor[0,1]/r,
                (d*tensor[1,1]/r**2-tensor[2,2]/(r**2*d))/(2*d)]
    unit=[s.factor(z*r*r/w) for z in components(linear)]
    expected=[6*x*x*d,3*x*x*(1-3*x*x),0,9*x*x]
    assert all(s.expand(a-b)==0 for a,b in zip(unit,expected))
    for L in [2,3]:
        assert prior.Q(components(round_source),L)==[0,0,0,0]
    return unit, {'unit_linear_defect':[str(z) for z in unit],
                  'radial_derivative_terms':'exactly tangent for arbitrary w(r), all epsilon',
                  'round_shift_source':'retained for both cutoffs'}


def finite_projection(L,epsilon,nodes):
    """Full diagonal m=0 tangent Gram solve at r=1 and unit angular smearing.

    Axisymmetry decouples m!=0. All off-diagonal tensor tangent fields are
    pointwise DeWitt orthogonal to this diagonal source. Thus this diagonal
    block suffices for the full normal projection of the specified source.
    """
    points,weights=np.polynomial.legendre.leggauss(nodes)
    d=1-points*points; rho=1+epsilon; Sigma=1+epsilon*points*points
    metric=np.stack([Sigma/rho,Sigma/d,rho*d],axis=-1)
    xi=3*points*d; xi_prime=3-9*points*points
    source=np.stack([xi*2*epsilon*points/rho,
                     xi*2*points*rho/(d*d)+2*Sigma*xi_prime/d,
                     -2*rho*points*xi],axis=-1)
    basis=[]; labels=[]
    def evaluate(expr):
        v=s.lambdify(x,expr,'numpy')(points)
        return np.broadcast_to(v,points.shape)
    for ell in range(L+1):
        P=evaluate(s.legendre(ell,x))
        basis.append(np.stack([P,np.zeros_like(P),np.zeros_like(P)],axis=-1));labels.append(['lambda',ell])
        basis.append(np.stack([np.zeros_like(P),P/d,d*P],axis=-1));labels.append(['X',ell])
    for ell in range(2,L+1):
        D2=evaluate(s.diff(s.legendre(ell,x),x,2))
        basis.append(np.stack([np.zeros_like(D2),D2,-d*d*D2],axis=-1));labels.append(['C_polar',ell])
    basis=np.asarray(basis)
    raised=basis/metric[None,:,:]; source_raised=source/metric
    tr=raised.sum(axis=-1); source_trace=source_raised.sum(axis=-1)
    measure=weights*Sigma/2
    gram=np.einsum('n,ani,bni->ab',measure,raised,raised)-np.einsum('n,an,bn->ab',measure,tr,tr)
    rhs=np.einsum('n,ani,ni->a',measure,raised,source_raised)-np.einsum('n,an,n->a',measure,tr,source_trace)
    coefficients=np.linalg.solve(gram,rhs)
    normal=source-np.einsum('a,ani->ni',coefficients,basis)
    normal_raised=normal/metric
    normal_trace=normal_raised.sum(axis=-1)
    norm=float(measure@(np.sum(normal_raised**2,axis=-1)-normal_trace**2))
    orthogonality=np.einsum('n,ani,ni->a',measure,raised,normal_raised)-np.einsum('n,an,n->a',measure,tr,normal_trace)
    assert np.max(np.abs(orthogonality))<1e-10
    eig=np.linalg.eigvalsh(gram)
    assert np.count_nonzero(eig<0)==L+1
    assert np.count_nonzero(eig>0)==2*L
    raw_norm=float(measure@(np.sum(source_raised**2,axis=-1)-source_trace**2))
    schur_hessian=float((rhs@coefficients-raw_norm)/2)
    direct_hessian=-norm/2
    assert abs(schur_hessian-direct_hessian)<1e-10
    return {'L':L,'epsilon_over_r_squared':epsilon,'quadrature_nodes':nodes,
            'diagonal_tangent_dimension':len(labels),
            'diagonal_Gram_positive':int(np.count_nonzero(eig>0)),
            'diagonal_Gram_negative':int(np.count_nonzero(eig<0)),
            'Gram_condition_number':float(np.linalg.cond(gram)),
            'normal_orthogonality_error':float(np.max(np.abs(orthogonality))),
            'shift_Hessian':direct_hessian,'shift_Hessian_over_epsilon_squared':direct_hessian/epsilon**2,
            'Schur_complement_Hessian':schur_hessian,
            'Schur_vs_direct_error':abs(schur_hessian-direct_hessian)}


def verify():
    start=time.monotonic()
    q,r,eps,geom=geometry()
    unit,defect=leading_defect(q,r,eps)
    exact=[]
    for L in [2,3]:
        normal=prior.Q(unit,L)
        norm=prior.pair(normal,normal)
        assert norm==s.Rational(192,1225)
        assert s.factor(normal[0]+2*normal[1])!=0
        exact.append({'N':L+1,'L':L,'normal_linear_components':[str(z) for z in normal],
                      'trace_of_normal_linear_term':str(s.factor(normal[0]+2*normal[1])),
                      'normal_DeWitt_norm':str(norm),
                      'leading_shift_Hessian_coefficient':str(-norm/2),
                      'smearing_and_radial_factor':'epsilon^2 * integral w(r)^2/r^2 dr',
                      'bracket_rank_lower_bound_on_exact_trajectory':2})
    numerical=[]
    for L in [2,3]:
        for amplitude in [0.0001,0.001,0.01,0.05,0.1,0.2]:
            rows=[finite_projection(L,amplitude,n) for n in [64,128]]
            assert abs(rows[0]['shift_Hessian']-rows[1]['shift_Hessian'])<1e-12
            if amplitude<=0.001:
                assert rows[1]['shift_Hessian']<0
                assert abs(rows[1]['shift_Hessian_over_epsilon_squared']+96/1225)<0.001
            numerical.extend(rows)
    result={'normalization':'(16piG)/(4pi) I; dOmega/(4pi) angular average',
            'geometry':geom,'leading_defect':defect,'exact_checks':exact,
            'finite_deformation_checks':numerical,
            'exact_dynamically_admissible_family':True,
            'nonzero_primary_secondary_bracket_on_exact_trajectory':True,
            'regular_full_constraint_stratum_proved':False,
            'generic_full_bracket_rank_proved':False,
            'physical_DOF_count_established':False,
            'approaches_excluded_fixed_nonzero_kappa_witness':False,
            'full_nonlinear_goal_complete':False}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
    sources=[Path(__file__),PRIOR,prior.PRIOR,ROOT/'sections/nonlinear_action.tex']
    manifest={'base_commit_at_execution':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'python':platform.python_version(),'platform':platform.platform(),
              'dependencies':{k:importlib.metadata.version(k) for k in ['numpy','sympy']},
              'seconds':round(time.monotonic()-start,3),'GPU_jobs':0,'experiment_launches':0,
              'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (OUT/'reproducibility.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PASS: exact Euclidean pullback and independent zero Ricci tensor; coefficient band j<=2.')
    print('PASS: stationary vacuum trajectory at N=3,4 with a=-1,b=V=W=0 and K=p=0.')
    print('PASS: arbitrary radial shift-profile derivatives lie exactly in the retained tangent.')
    print('Exact leading discarded DeWitt norm: 192/1225 (trace contribution retained).')
    print('H_WW = -96/1225 * epsilon^2 * integral w(r)^2/r^2 dr + O(epsilon^3).')
    for L in [2,3]:
        small=next(z for z in numerical if z['L']==L and z['epsilon_over_r_squared']==0.0001 and z['quadrature_nodes']==128)
        print(f"N={L+1}: finite Gram solve H_WW/epsilon^2 at epsilon/r^2=1e-4: {small['shift_Hessian_over_epsilon_squared']:.15g}")
    print('Maximum tangent orthogonality error:',max(z['normal_orthogonality_error'] for z in numerical))
    print('CONCLUSION: primary-secondary Poisson rank at least 2 at an exact admissible point, for sufficiently small nonzero epsilon.')
    print('Full constant-rank stratum, generic physical DOF, and nonlinear completion remain open.')


if __name__=='__main__':
    verify()
