#!/usr/bin/env python3
"""First secondary-preservation step for the actual Paper I vacuum witness.

Keeps the nonlinear coefficient embedding and radial derivatives. This checks
instantaneous solvability of Cdot=0, not final Dirac closure or an evolution.
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
PRIOR = ROOT/'verification/canonical_analysis_20261127/actual_adm.py'
spec = importlib.util.spec_from_file_location('canonical_checkpoint', PRIOR)
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
x = prior.x
OUT = ROOT/'reports/dirac-stabilization-2026-11-27'
# Entries represent f,g,v_theta/sqrt(1-x²),c/(1-x²), in an orthonormal frame.
K = [x*x, (1-x*x)/2, -x, s.Rational(1,2)]
B = [8*x*x*(1-x*x), s.S.Zero, -8*x*(1-x*x), s.S.Zero]
SPINS = [0,0,1,2]


def Q(fields, L):
    return [prior.residual(z, spin, L) for z,spin in zip(fields,SPINS)]


def pair_integrand(U,V):
    f,g,v,c=U; F,G,Vv,C=V
    return s.expand(-2*f*G-2*g*F-2*g*G+2*(1-x*x)*v*Vv+2*(1-x*x)**2*c*C)


def pair(U,V):
    return prior.average(pair_integrand(U,V))


def covariant(fields,r):
    f,g,v,c=fields; d=1-x*x
    return s.Matrix([[f,-r*v,0],[-r*v,r*r*(g+d*c)/d,0],
                     [0,0,r*r*d*(g-d*c)]])


def embedding_check():
    r,t=s.symbols('r t',positive=True)
    d=1-x*x
    # Coordinate x=cos(theta) avoids coordinate trigonometric simplification.
    lam=1+2*x*x*t
    gammaxx=r*r/d+2*r*r*t
    gammaphi=r*r*d
    Vx=2*x*d*t/r
    sigma=s.Matrix([[lam+gammaxx*Vx**2,gammaxx*Vx,0],
                    [gammaxx*Vx,gammaxx,0],[0,0,gammaphi]])
    velocity=sigma.diff(t).subs(t,0)
    acceleration=sigma.diff(t,2).subs(t,0)
    assert all(s.simplify(z)==0 for z in velocity-2*covariant(K,r))
    assert all(s.simplify(z)==0 for z in acceleration-covariant(B,r))
    return {'embedding_acceleration_components':[str(z) for z in B],
            'linear_coefficient_path_velocity':'2K exactly',
            'linear_coefficient_path_acceleration':'B exactly'}


def coordinate_constraint_derivative(R):
    """Differentiate continuum constraint densities on the actual Galerkin jet.

    qdot=2K is covariantly constant; all connection time derivatives vanish.
    Kdot=K0+R/2 at kappa=1 after subtracting alphadot. No radial freeze.
    """
    r,phi=s.symbols('r phi',positive=True)
    coords=[r,x,phi]; d=1-x*x
    q=s.diag(1,r*r/d,r*r*d); qi=q.inv()
    Kcov=covariant(K,r); qdot=2*Kcov; qidot=-qi*qdot*qi
    Rcov=covariant(R,r); Kdot=Kcov+Rcov/2
    Gamma=[[[s.factor(sum(qi[i,k]*(s.diff(q[k,j],coords[l])
                 +s.diff(q[k,l],coords[j])-s.diff(q[j,l],coords[k]))
                 for k in range(3))/2) for l in range(3)] for j in range(3)] for i in range(3)]
    Gdot=[[[s.factor(sum(qidot[i,k]*(s.diff(q[k,j],coords[l])
                 +s.diff(q[k,l],coords[j])-s.diff(q[j,l],coords[k]))
                 +qi[i,k]*(s.diff(qdot[k,j],coords[l])+s.diff(qdot[k,l],coords[j])
                 -s.diff(qdot[j,l],coords[k])) for k in range(3))/2)
                 for l in range(3)] for j in range(3)] for i in range(3)]
    assert all(z==0 for plane in Gdot for row in plane for z in row)
    mixed=qi*Kcov; mixed_dot=qidot*Kcov+qi*Kdot
    Hdot=s.factor(2*s.trace(mixed*mixed_dot)-2*s.trace(mixed)*s.trace(mixed_dot))
    assert s.simplify(Hdot-pair_integrand(K,R))==0
    stress=mixed-s.eye(3)*s.trace(mixed)
    stress_dot=mixed_dot-s.eye(3)*s.trace(mixed_dot)
    Mdot=[]
    for j in range(3):
        z=sum(s.diff(stress_dot[i,j],coords[i]) for i in range(3))
        z+=sum(Gamma[i][i][k]*stress_dot[k,j]-Gamma[k][i][j]*stress_dot[i,k]
               +Gdot[i][i][k]*stress[k,j]-Gdot[k][i][j]*stress[i,k]
               for i in range(3) for k in range(3))
        Mdot.append(s.factor(-2*z))
    assert Mdot[2]==0
    assert s.diff(s.factor(r*Mdot[0]),r)==0
    assert s.diff(Mdot[1],r)==0
    return {'Hdot':Hdot,'Mdot':Mdot,'r':r,
            'spatial_connection_time_derivative':'all components zero',
            'spatial_curvature_time_derivative':'zero (qdot parallel and spatial Ricci zero)'}


def verify():
    start=time.monotonic(); embedding=embedding_check()
    prior_matrices=prior.all_m_hessian()
    rows=[]
    phi=s.symbols('phi',real=True)
    for L in [2,3]:
        QB=Q(B,L)
        mu=s.Rational(32,9)*s.legendre(2,x)
        R=Q([b+mu*k for b,k in zip(B,K)],L)
        trace_R=s.factor(R[0]+2*R[1])
        assert s.expand(trace_R+s.Rational(64,35)*s.legendre(4,x))==0
        # Actual coefficient acceleration is retained despite the nonzero normal residue.
        qacc=[s.factor(2*k-mu*k-b+z) for k,b,z in zip(K,B,R)]
        assert Q(qacc,L)==[0,0,0,0]
        assert any(z!=0 for z in R)
        sources=[]
        for ell in range(L+1):
            T=Q([s.legendre(ell,x)*k for k in K],L)
            assert s.expand(T[0]+2*T[1])==0
            force=pair(T,QB); norm=pair(T,T)
            residual=pair(T,R)
            assert residual==0
            sources.append({'ell':ell,'forcing':str(force),'Gram_diagonal':str(norm),
                            'stabilization_residual':str(-residual/2)})
        assert sources[2]['forcing']=='-128/525'
        assert sources[2]['Gram_diagonal']=='12/175'
        full=next(z for z in prior_matrices if z['L']==L)
        labels=[tuple(z) for z in full['labels']]
        diagonal=[s.sympify(z) for z in full['discarded_norm_diagonal']]
        forcing=[]
        for ell,m in labels:
            azimuth=s.simplify(s.integrate(s.exp(s.I*m*phi),(phi,0,2*s.pi))/(2*s.pi))
            assert azimuth==(1 if m==0 else 0)
            forcing.append(s.sqrt(2*ell+1)*s.sympify(sources[ell]['forcing']) if m==0 else s.S.Zero)
        rates=[s.simplify(-f/D) if D!=0 else s.S.Zero for f,D in zip(forcing,diagonal)]
        assert all(s.simplify(f+D*rate)==0 for f,D,rate in zip(forcing,diagonal,rates))
        assert rates[labels.index((2,0))]==32/(9*s.sqrt(5))
        assert sum(z!=0 for z in rates)==1
        direct=coordinate_constraint_derivative(R)
        r=direct['r']; Hdot=direct['Hdot']; Mr,Mx,Mphi=direct['Mdot']
        assert s.expand(Hdot-s.Rational(256,385)*s.legendre(4,x)
                        -s.Rational(128,231)*s.legendre(6,x))==0
        assert s.simplify(Mr+s.Rational(64,21)*s.legendre(4,x)/r)==0
        assert s.expand(Mx-s.Rational(16,105)*s.diff(s.legendre(4,x),x))==0
        coordinate_tests=[]
        for ell in range(L+1):
            Y=s.legendre(ell,x)
            h=prior.average(Y*Hdot)
            radial=prior.average(Y*r*Mr)
            angular=prior.average((1-x*x)*s.diff(Y,x)*Mx)
            assert h==radial==angular==0
            coordinate_tests.append({'ell':ell,'lapse':str(h),'radial_shift':str(radial),
                                     'polar_angular_shift':str(angular)})
        # Independent numerical integral of the source and Gram denominator.
        numerical=[]
        T=Q([s.legendre(2,x)*k for k in K],L)
        source_fun=s.lambdify(x,pair_integrand(T,QB),'numpy')
        gram_fun=s.lambdify(x,pair_integrand(T,T),'numpy')
        for order in [24,48,96]:
            points,weights=np.polynomial.legendre.leggauss(order)
            forcing_num=float(weights@source_fun(points)/2)
            gram_num=float(weights@gram_fun(points)/2)
            rate=-forcing_num/gram_num
            assert abs(rate-32/9)<1e-12
            numerical.append({'quadrature_order':order,'source':forcing_num,'Gram':gram_num,
                              'mu_P2_coefficient':rate,'absolute_error':abs(rate-32/9)})
        rows.append({'N':L+1,'L':L,'axisymmetric_sources':sources,
                     'particular_mu_a':'(32/9) * kappa * P_2(cos(theta))',
                     'particular_alpha_dot':'-(16/9) * kappa * P_2(cos(theta))',
                     'homogeneous_lapse_rate_space':f'V_{L-2}',
                     'temporal_coefficients':4*(L+1)**2-2,
                     'full_temporal_Hessian_rank_at_witness':4*L,
                     'free_radial_rate_functions_this_step':4*(L+1)**2-2-4*L,
                     'all_m_labels':[list(z) for z in labels],
                     'all_m_forcing':[str(z) for z in forcing],
                     'all_m_particular_mu_a':[str(z) for z in rates],
                     'normal_acceleration_residue':[str(z) for z in R],
                     'multiplier_independent_trace_residue':str(trace_R),
                     'exact_Einstein_second_order_time_jet_in_ansatz':False,
                     'retained_coefficient_acceleration_components':[str(z) for z in qacc],
                     'direct_continuum_constraint_derivatives':{
                         'Hdot':str(Hdot),'Mdot_r':str(Mr),'Mdot_x':str(Mx),'Mdot_phi':str(Mphi),
                         'connection_dot':direct['spatial_connection_time_derivative']},
                     'retained_coordinate_projections':coordinate_tests,
                     'quadrature_checks':numerical})
    result={'normalization':'normalized action (16piG)/(4pi) I; local components kappa=1; radial derivatives retained',
            'embedding':embedding,'checks':rows,
            'first_secondary_preservation_solvable':True,
            'tertiary_obstruction_at_this_instant':False,
            'generic_tertiary_constraints_classified':False,
            'continued_solution_established':False,'physical_DOF_count_established':False,
            'full_nonlinear_goal_complete':False}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
    sources=[Path(__file__),PRIOR,ROOT/'sections/nonlinear_action.tex']
    manifest={'base_commit_at_execution':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'python':platform.python_version(),'platform':platform.platform(),
              'dependencies':{k:importlib.metadata.version(k) for k in ['numpy','sympy']},
              'seconds':round(time.monotonic()-start,3),'GPU_jobs':0,'experiment_launches':0,
              'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (OUT/'reproducibility.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PASS: nonlinear coefficient-embedding acceleration, all-m lapse preservation, and direct coordinate constraint derivatives.')
    for row in rows:
        print(f"N={row['N']}: temporal Hessian rank {row['full_temporal_Hessian_rank_at_witness']}; {row['free_radial_rate_functions_this_step']} free rate functions at this step; all retained Cdot=0.")
    print('Particular rate: adot=(32/9)*kappa*P2; alphadot=-(16/9)*kappa*P2.')
    print('Exact Einstein two-jet obstruction: tr(R)=-(64/35)*kappa^2*P4 for any retained temporal rates.')
    print('Maximum independent quadrature error:',max(z['absolute_error'] for row in rows for z in row['quadrature_checks']))
    print('Nonzero discarded acceleration remains. No claim of continued solution, final constraint classification, or physical DOF.')


if __name__=='__main__':
    verify()
