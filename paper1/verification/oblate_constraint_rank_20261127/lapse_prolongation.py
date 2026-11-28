#!/usr/bin/env python3
"""Exact leading normal-Hessian pairing after eliminating shift image.

The entries retain separate radial profile and first-derivative slots. They
are not a final Dirac count; radial integration and constraint independence
are distinct analytic steps.
"""
from pathlib import Path
import importlib.util
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import time
import sympy as s

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('all_modes',HERE/'all_modes.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
certificate=a.certificate
x,d=a.x,a.d
OUT=a.ROOT/'reports/oblate-lapse-consistency-2026-11-27'


def leading_hessians(L,m):
    fields=[];labels=[]
    for ell in range(m,L+1):
        Y=a.harmonic(ell,m)[0];DY=s.diff(Y,x)
        H0=s.Matrix([[x*d*DY,2*x*x*DY,2*s.I*m*Y],
                      [2*x*x*DY,-x*DY,0],[2*s.I*m*Y,0,-x*d*d*DY]])
        H1=s.Matrix([[-d*Y,-x*Y-x*x*DY,-s.I*m*Y],
                      [-x*Y-x*x*DY,Y,0],[-s.I*m*Y,0,d*d*Y]])
        fields.extend([H0,H1]);labels.extend([[ell,'f/r^2'],[ell,'f_prime/r']])
    return fields,labels


def verify_hessian_derivative(L,m):
    e=a.eps;rho=1+e;Sigma=1+e*x*x
    for ell in range(m,L+1):
        Y,P,_,_,_=a.harmonic(ell,m);DY=s.diff(Y,x)
        H0=s.Matrix([[e*x*d*DY/(rho*Sigma),-DY/Sigma,-s.I*m*Y/rho],
                      [-DY/Sigma,s.diff(Y,x,2)-x*rho*DY/(Sigma*d),s.I*m*(DY+x*Y/d)],
                      [-s.I*m*Y/rho,s.I*m*(DY+x*Y/d),-m*m*Y-rho*x*d*DY/Sigma]])
        H1=s.Matrix([[-e*d*Y/(rho*Sigma),DY-e*x*Y/Sigma,s.I*m*Y],
                      [DY-e*x*Y/Sigma,rho*Y/(Sigma*d),0],
                      [s.I*m*Y,0,rho*d*Y/Sigma]])
        T1=s.zeros(3)
        T1[0,1]=T1[1,0]=x*x*P[0]/d
        T1[0,2]=T1[2,0]=d*P[1]
        D0=H0.diff(e).subs(e,0)+T1
        D1=H1.diff(e).subs(e,0)-T1
        expected=leading_hessians(L,m)[0][2*(ell-m):2*(ell-m)+2]
        assert all(s.factor(z)==0 for z in D0-expected[0])
        assert all(s.factor(z)==0 for z in D1-expected[1])


def calculate():
    start=time.monotonic()
    OUT.mkdir(parents=True,exist_ok=True)
    certificate_path=a.OUT/'exact-leading.json'
    certificate_sources=certificate.source_fingerprints(HERE)
    prior,certificate_provenance=certificate.load_certificate(certificate_path,certificate_sources)
    rows=[]
    for L in [2,3]:
        for m in range(L+1):
            verify_hessian_derivative(L,m)
            targets,labels=leading_hessians(L,m)
            coeffs=[];grams=[]
            for ell in [L+1,L+2]:
                T=a.round_degree_tangents(ell,m)
                norm=s.factorial(ell+m)/(s.factorial(ell-m)*(2*ell+1))
                G=s.diag(1,1,2*ell*(ell+1),2*ell*(ell+1),
                         s.Rational(1,2)*(ell-1)*ell*(ell+1)*(ell+2),
                         s.Rational(1,2)*(ell-1)*ell*(ell+1)*(ell+2))*norm
                G[:2,:2]=s.Matrix([[0,-2],[-2,-2]])*norm
                coupling=s.Matrix([[a.round_pair(t,h) for h in targets] for t in T])
                coeffs.extend((G.inv()*coupling).tolist());grams.append(G)
            H=s.Matrix(coeffs);G=s.diag(*grams)
            old=next(z for z in prior if z['L']==L and z['m']==m)
            C=s.Matrix([[s.sympify(z) for z in row] for row in old['normal_coefficients']])
            selected=old['nonzero_principal_minor_indices']
            C=C[:,selected]
            image_gram=C.conjugate().T*G*C
            Z=(H-C*image_gram.inv()*C.conjugate().T*G*H).applyfunc(s.factor)
            Zgram=(Z.conjugate().T*G*Z).applyfunc(s.factor)
            assert C.conjugate().T*G*Z==s.zeros(len(selected),len(labels))
            radial_forms=[]
            for ell in range(m,L+1):
                indices=[2*(ell-m),2*(ell-m)+1]
                block=Zgram.extract(indices,indices)
                others=[i for i in range(len(labels)) if i not in indices]
                assert Zgram.extract(indices,others)==s.zeros(2,len(others))
                if ell<L-1:
                    assert block==s.zeros(2)
                    continue
                direction=s.Matrix([1,-1] if (L,ell)==(2,1) else
                                   [0,1] if ell==2 else [1,1])
                coefficient=block[0,0] if direction[0] else block[1,1]
                assert coefficient!=0
                assert block==coefficient*direction*direction.T
                radial_forms.append({'ell':ell,'m':m,'DeWitt_coefficient':str(coefficient),
                                     'profile_operator':[int(z) for z in direction],
                                     'radial_form':'epsilon^2 * integral r^-6 (a*f+b*r*f_prime)^*(a*g+b*r*g_prime) dr',
                                     'constraint_bracket_factor':'-1/2 times DeWitt form'})
            rows.append({'L':L,'m':m,'labels':labels,
                         'residual_normal_coefficients':[[str(z) for z in row] for row in Z.tolist()],
                         'residual_DeWitt_Gram':[[str(z) for z in row] for row in Zgram.tolist()],
                         'residual_map_rank':Z.rank(),'residual_Gram_rank':Zgram.rank(),
                         'nondegenerate_radial_forms':radial_forms})
            nonzero=[(labels[i],labels[j],str(Zgram[i,j])) for i in range(len(labels)) for j in range(i,len(labels)) if Zgram[i,j]!=0]
            print('LAPSE_BLOCK',L,m,'map_rank',Z.rank(),'Gram_rank',Zgram.rank(),'entries',nonzero,flush=True)
    (OUT/'lapse-leading.json').write_text(json.dumps(rows,indent=2)+'\n')
    kernels=gradient_kernel()
    totals=[]
    for L in [2,3]:
        active=sum((1 if z['m']==0 else 2)*len(z['nondegenerate_radial_forms']) for z in rows if z['L']==L)
        assert active==4*L
        totals.append({'N':L+1,'L':L,'leading_active_scalar_multiplets_real':active,
                       'exact_redundant_lapse_profile_functions':(L-1)**2,
                       'nonlinear_final_Dirac_classification':False})
    (OUT/'lapse-summary.json').write_text(json.dumps({'exact_gradient_kernel':kernels,
        'counts':totals,'analytic_step_required':'weighted first-derivative forms on compact radial profiles and small-parameter operator perturbation',
        'physical_DOF_count_established':False},indent=2)+'\n')
    if certificate.source_fingerprints(HERE)!=certificate_sources:
        raise RuntimeError('certificate producer sources changed during lapse calculation')
    sources=[Path(__file__),HERE/'all_modes.py',HERE/'leading_certificate.py']
    manifest={'base_commit_at_execution':subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.ROOT,text=True).strip(),
              'python':platform.python_version(),'dependencies':{k:importlib.metadata.version(k) for k in ['numpy','sympy']},
              'seconds':round(time.monotonic()-start,3),'GPU_jobs':0,'experiment_launches':0,
              'exact_leading_certificate':certificate_provenance,
              'sha256':{str(p.relative_to(a.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    manifest['sha256'][str(certificate_path.relative_to(a.ROOT))]=certificate_provenance['payload_sha256']
    manifest['sha256'][str(certificate.manifest_path(certificate_path).relative_to(a.ROOT))]=certificate_provenance['manifest_sha256']
    (OUT/'lapse-reproducibility.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PASS: exact leading lapse forms retain the radial profile and its first derivative.',flush=True)
    print('PASS: active scalar profile counts are 8/12; exact redundant gradient families have 1/4 profile functions.',flush=True)
    print('These are linearized prolongation results, not a nonlinear physical degree-of-freedom count.',flush=True)


def gradient_kernel():
    r,e=s.symbols('r e',real=True);rho=r*r+e;Sigma=r*r+e*x*x
    f0=r*r*x+e*x**3/3
    f1=s.sqrt(d)*(r*r+2*e/3+e*x*x/3)
    W0=s.Matrix([d*s.diff(f0,x)/Sigma,0])
    assert all(s.factor(z)==0 for z in W0-s.Matrix([d,0]))
    W1=s.Matrix([d*s.diff(f1,x)/Sigma,s.I*f1/(rho*d)])
    P1=a.harmonic(1,1)[1];A2=a.harmonic(2,1)[2]
    c1=-(r*r+5*e/6)/rho;c2=-s.I*e/(18*rho)
    # In this m=1 chart the components carry sqrt(d) and 1/sqrt(d).
    # Remove those nonzero factors for -1<x<1, then cancel rationally.
    # factor() alone can return unevaluated 0*I, whose structural ==0 is false.
    residual=W1-c1*P1-c2*A2
    rational_residual=[s.cancel(residual[0]/s.sqrt(d)),
                       s.cancel(residual[1]*s.sqrt(d))]
    assert rational_residual==[s.S.Zero,s.S.Zero]
    return {'N3':'arbitrary angle-independent F0(r)',
            'N4_additional_templates':[str(f0),str(f1)],
            'N4_templates_meaning':'axisymmetric template plus real and imaginary parts of m=1 template times exp(i phi); each multiplied by an arbitrary smooth radial function',
            'N4_gradient_m0':'polar degree 1',
            'N4_gradient_m1_coefficients':[str(c1),str(c2)],
            'N4_gradient_m1_rationalized_residual':[str(z) for z in rational_residual],
            'gradient_identity_domain':'-1<x<1; equality of smooth retained fields extends across the coordinate poles',
            'N4_gradient_m1_basis':['polar (ell,m)=(1,1)','axial (ell,m)=(2,1)'],
            'radial_actual_shift_coefficient':'b=partial_r f, retained through L'}


if __name__=='__main__':calculate()
