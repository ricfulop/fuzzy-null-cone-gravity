#!/usr/bin/env python3
"""All magnetic sectors of the oblate spatial Gram and temporal shift Hessian.

Complex Fourier blocks m>=0 represent real m=0 modes once and m>0 modes twice.
Numerical ranks are diagnostics; exact leading matrices are computed separately.
"""
from pathlib import Path
from functools import lru_cache
import importlib.util
import sys
import json
import numpy as np
import sympy as s


def _load_leading_certificate():
    path=Path(__file__).resolve().with_name('leading_certificate.py')
    spec=importlib.util.spec_from_file_location('oblate_leading_certificate',path)
    if spec is None or spec.loader is None:
        raise ImportError(f'Cannot load sibling certificate helper: {path}')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


certificate=_load_leading_certificate()

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/oblate-constraint-rank-2026-11-27'
x,eps=s.symbols('x epsilon',real=True)
d=1-x*x


@lru_cache(None)
def harmonic(ell,m):
    y=s.assoc_legendre(ell,m,x)
    polar=s.Matrix([d*s.diff(y,x),s.I*m*y/d])
    axial=s.Matrix([s.I*m*y,-s.diff(y,x)])
    E=s.zeros(2)
    E[0,0]=s.diff(y,x,2)-x*s.diff(y,x)/d+ell*(ell+1)*y/(2*d)
    E[0,1]=E[1,0]=s.I*m*(s.diff(y,x)+x*y/d)
    E[1,1]=-d*d*E[0,0]
    B=s.Matrix([[E[0,1]/d,-d*E[0,0]],[-d*E[0,0],-d*E[0,1]]])
    return y,polar,axial,E.applyfunc(s.factor),B.applyfunc(s.factor)


def fields(L,m):
    rho=1+eps; Sigma=1+eps*x*x
    lam=Sigma/rho; gx=Sigma/d; gp=rho*d
    tangent=[]; tlabels=[]; shifts=[]; slabels=[]; leading=[]
    for ell in range(m,L+1):
        y,P,A,E,B=harmonic(ell,m)
        tangent.extend([s.diag(y,0,0),s.diag(0,y/d,d*y)])
        tlabels.extend([['lambda',ell],['X',ell]])
        # b' term is 2*b'*Y in rr and is exactly a retained lambda variation.
        R=s.zeros(3)
        R[0,0]=-2*eps*d*y/(rho*Sigma)
        R[0,1]=R[1,0]=s.diff(y,x)-2*eps*x*y/Sigma
        R[0,2]=R[2,0]=s.I*m*y
        R[1,1]=2*rho*y/(Sigma*d);R[2,2]=2*rho*d*y/Sigma
        shifts.append(R);slabels.append(['b',ell])
        # Subtract derivative of the round tangential delta V=grad Y.
        D=R.diff(eps).subs(eps,0)
        if ell>=1:
            D[0,1]-=x*x*P[0]/d;D[1,0]=D[0,1]
            D[0,2]-=d*P[1];D[2,0]=D[0,2]
        leading.append(D.applyfunc(s.factor))
        if ell>=1:
            for name,v in [('polar',P),('axial',A)]:
                T=s.zeros(3)
                T[0,1]=T[1,0]=gx*v[0];T[0,2]=T[2,0]=gp*v[1]
                tangent.append(T);tlabels.append(['V_'+name,ell])
                U=s.zeros(3)
                U[0,0]=v[0]*s.diff(lam,x)
                U[1,1]=v[0]*s.diff(gx,x)+2*gx*s.diff(v[0],x)
                U[2,2]=v[0]*s.diff(gp,x)+2*gp*s.I*m*v[1]
                U[1,2]=U[2,1]=gx*s.I*m*v[0]+gp*s.diff(v[1],x)
                shifts.append(U);slabels.append(['W_'+name,ell])
                leading.append(U.diff(eps).subs(eps,0).applyfunc(s.factor))
        if ell>=2:
            for name,C in [('polar',E),('axial',B)]:
                T=s.zeros(3);T[1:,1:]=C
                tangent.append(T);tlabels.append(['C_'+name,ell])
    return tangent,tlabels,shifts,slabels,leading


@lru_cache(None)
def evaluators(L,m):
    T,tl,U,ul,_=fields(L,m)
    def compile_matrices(matrices):
        return [[s.lambdify((x,eps),entry,'numpy') for entry in mat] for mat in matrices]
    return compile_matrices(T),tl,compile_matrices(U),ul


def evaluate(functions,points,amplitude):
    return np.asarray([[np.broadcast_to(f(points,amplitude),points.shape)
                        for f in row] for row in functions],dtype=complex).transpose(0,2,1).reshape(len(functions),len(points),3,3)


def gram(A,B,metric,measure,dewitt=True):
    denominator=np.sqrt(metric[:,:,None]*metric[:,None,:])
    a=A/denominator[None,:,:,:];b=B/denominator[None,:,:,:]
    g=np.einsum('n,anij,bnij->ab',measure,a.conj(),b)
    if dewitt:
        ta=np.trace(a,axis1=-2,axis2=-1);tb=np.trace(b,axis1=-2,axis2=-1)
        g-=np.einsum('n,an,bn->ab',measure,ta.conj(),tb)
    return (g+g.conj().T)/2 if A is B else g


def numerical_block(L,m,amplitude,nodes):
    points,weights=np.polynomial.legendre.leggauss(nodes)
    D=1-points*points;rho=1+amplitude;Sigma=1+amplitude*points*points
    metric=np.stack([Sigma/rho,Sigma/D,rho*D],axis=-1)
    measure=weights*Sigma/2
    tf,tl,uf,ul=evaluators(L,m)
    T=evaluate(tf,points,amplitude);U=evaluate(uf,points,amplitude)
    G=gram(T,T,metric,measure)
    # Scale coefficient basis to improve conditioning without changing its span.
    positive=gram(T,T,metric,measure,dewitt=False)
    scales=np.sqrt(np.real(np.diag(positive)))
    T=T/scales[:,None,None,None]
    G=gram(T,T,metric,measure);rhs=gram(T,U,metric,measure)
    coeff=np.linalg.solve(G,rhs)
    normal=U-np.einsum('ab,anij->bnij',coeff,T)
    Dgram=gram(normal,normal,metric,measure)
    Egram=gram(normal,normal,metric,measure,dewitt=False)
    # Normalize shift columns by their round vector/scalar coefficient norm.
    norms=[]
    for name,ell in ul:
        n=float(s.factorial(ell+m)/(s.factorial(ell-m)*(2*ell+1)))
        if name!='b':n*=ell*(ell+1)
        norms.append(np.sqrt(n))
    norms=np.asarray(norms)
    Dgram=Dgram/norms[:,None]/norms[None,:]
    Egram=Egram/norms[:,None]/norms[None,:]
    de=np.linalg.eigvalsh(Dgram);ee=np.linalg.eigvalsh(Egram);ge=np.linalg.eigvalsh(G)
    orth=float(np.max(np.abs(gram(T,normal,metric,measure))))
    assert orth<1e-9
    # Round scalar indefiniteness persists in each sampled full m-block.
    assert np.count_nonzero(ge<0)==L-m+1
    assert len(ge)==len(tl)
    ranks={str(tol):{'DeWitt':int(np.count_nonzero(np.abs(de)>tol)),
                     'normal_map':int(np.count_nonzero(ee>tol))}
           for tol in [1e-8,1e-10,1e-12,1e-14]}
    return {'L':L,'m':m,'epsilon_over_r_squared':amplitude,'nodes':nodes,
            'shift_column_normalization':'round coefficient L2 norm; angular vectors include ell*(ell+1)',
            'spatial_labels':tl,'shift_labels':ul,'spatial_dimension_complex':len(tl),
            'shift_dimension_complex':len(ul),'spatial_Gram_positive':int(np.count_nonzero(ge>0)),
            'spatial_Gram_negative':int(np.count_nonzero(ge<0)),
            'spatial_Gram_condition':float(np.linalg.cond(G)),
            'normal_orthogonality_error':orth,
            'normal_DeWitt_eigenvalues':de.tolist(),'normal_positive_eigenvalues':ee.tolist(),
            'rank_at_absolute_threshold':ranks,
            'Dgram_real':Dgram.real.tolist(),'Dgram_imag':Dgram.imag.tolist()}


def exploratory():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    for L in [2,3]:
        for amplitude in [0.01,0.05,0.1]:
            totals={str(t):[0,0] for t in [1e-8,1e-10,1e-12,1e-14]}
            for m in range(L+1):
                z=numerical_block(L,m,amplitude,96);rows.append(z)
                for tol,r in z['rank_at_absolute_threshold'].items():
                    totals[tol][0]+=(1 if m==0 else 2)*r['DeWitt']
                    totals[tol][1]+=(1 if m==0 else 2)*r['normal_map']
                print('BLOCK',L,amplitude,m,z['rank_at_absolute_threshold'],flush=True)
            print('TOTAL',L,amplitude,totals,flush=True)
    (OUT/'exploratory-all-modes.json').write_text(json.dumps(rows,indent=2)+'\n')


def round_degree_tangents(ell,m):
    y,P,A,E,B=harmonic(ell,m)
    T=[s.diag(y,0,0),s.diag(0,y/d,d*y)]
    for v in [P,A]:
        z=s.zeros(3);z[0,1]=z[1,0]=v[0]/d;z[0,2]=z[2,0]=d*v[1]
        T.append(z)
    for C in [E,B]:
        z=s.zeros(3);z[1:,1:]=C;T.append(z)
    return T


@lru_cache(None)
def average_polynomial(expr):
    expr=s.factor(expr)
    # In each same-m product the associated-Legendre square roots cancel.
    expr=s.cancel(s.powsimp(expr,force=True))
    poly=s.Poly(expr,x)
    value=sum(coeff/(power[0]+1) for power,coeff in poly.terms() if power[0]%2==0)
    return s.factor(value)


def round_pair(A,B):
    inverse=[1,d,1/d]
    # Amplitudes are real functions of x except explicit powers of I.
    Ac=A.applyfunc(lambda z:z.subs(s.I,-s.I))
    expr=sum(Ac[i,j]*B[i,j]*inverse[i]*inverse[j] for i in range(3) for j in range(3))
    expr-=sum(Ac[i,i]*inverse[i] for i in range(3))*sum(B[j,j]*inverse[j] for j in range(3))
    return average_polynomial(expr)


def exact_leading():
    source_hashes=certificate.source_fingerprints(Path(__file__).resolve().parent)
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    for L in [2,3]:
        total=0
        for m in range(L+1):
            _,_,_,labels,U=fields(L,m)
            D=s.zeros(len(U));normal_coefficients=[]
            for ell in [L+1,L+2]:
                T=round_degree_tangents(ell,m)
                norm=s.factorial(ell+m)/(s.factorial(ell-m)*(2*ell+1))
                G=s.diag(1,1,2*ell*(ell+1),2*ell*(ell+1),
                         s.Rational(1,2)*(ell-1)*ell*(ell+1)*(ell+2),
                         s.Rational(1,2)*(ell-1)*ell*(ell+1)*(ell+2))*norm
                G[:2,:2]=s.Matrix([[0,-2],[-2,-2]])*norm
                coupling=s.Matrix([[round_pair(t,u) for u in U] for t in T])
                coefficients=G.inv()*coupling
                D+=coupling.conjugate().T*coefficients
                normal_coefficients.extend(coefficients.tolist())
            D=D.applyfunc(s.factor)
            rank=D.rank();normal_rank=s.Matrix(normal_coefficients).rank()
            expected=3 if m==L else 4
            assert rank==normal_rank==expected,(L,m,rank,normal_rank,expected)
            total+=(1 if m==0 else 2)*rank
            # Principal minor: independent columns plus their corresponding rows.
            independent=list(D.rref()[1])
            minor=D.extract(independent,independent)
            assert minor.det()!=0
            sylvester=[s.factor(minor[:k,:k].det()) for k in range(1,len(independent)+1)]
            assert all(z>0 for z in sylvester)
            rows.append({'L':L,'N':L+1,'m':m,'labels':labels,
                         'normal_row_labels':certificate.normal_labels(L),
                         'leading_DeWitt_Gram':[[str(z) for z in row] for row in D.tolist()],
                         'exact_leading_rank':rank,'exact_leading_normal_map_rank':normal_rank,
                         'nonzero_principal_minor_indices':independent,
                         'nonzero_principal_minor_determinant':str(s.factor(minor.det())),
                         'positive_Sylvester_minors':[str(z) for z in sylvester],
                         'normal_coefficients':[[str(z) for z in row] for row in normal_coefficients]})
            print('EXACT_LEADING_BLOCK',L,m,rank,'minor',s.factor(minor.det()),flush=True)
        assert total==8*L+2,(L,total)
        print('EXACT_LEADING_TOTAL',L,total,flush=True)
    if certificate.source_fingerprints(Path(__file__).resolve().parent)!=source_hashes:
        raise RuntimeError('producer sources changed during exact-leading calculation')
    certificate.write_certificate(OUT/'exact-leading.json',rows,source_hashes)


def exact_kernel():
    rows=[]
    for L in [2,3]:
        total=0
        for m in range(L+1):
            labels=[];wx=[]
            for ell in range(max(1,m),L+1):
                _,P,A,_,_=harmonic(ell,m)
                labels.extend([['W_polar',ell],['W_axial',ell]])
                wx.extend([P[0],A[0]])
            matrix=[]
            for j in [L,L+1]:
                y=harmonic(j,m)[0]
                norm=s.factorial(j+m)/(s.factorial(j-m)*(2*j+1))
                matrix.append([average_polynomial(y*z)/norm for z in wx])
            F=s.Matrix(matrix)
            assert F.rank()==2
            null=F.nullspace()
            radial=max(0,L-1-m)
            total+=(1 if m==0 else 2)*(len(null)+radial)
            rows.append({'L':L,'m':m,'angular_labels':labels,
                         'forbidden_scalar_degrees':[L,L+1],
                         'W_x_high_degree_map':[[str(z) for z in row] for row in F.tolist()],
                         'angular_kernel_dimension_complex':len(null),
                         'angular_kernel_basis':[[str(z) for z in v] for v in null],
                         'radial_kernel_dimension_complex':radial})
        assert total==3*L*L-2*L-1
        print('EXACT_KERNEL_TOTAL',L,total,flush=True)
    # Exact tensor identity used for the radial null fields, arbitrary f.
    r,e=s.symbols('r e',real=True);fx,fp=s.symbols('f_x f_phi')
    rho=r*r+e;Sigma=r*r+e*x*x;lam=Sigma/rho
    actual=s.Matrix([lam*d*fx/Sigma,lam*fp/(rho*d)])
    claimed=s.Matrix([d*fx/rho,fp/(rho*d)-e*fp/rho**2])
    assert all(s.factor(z)==0 for z in actual-claimed)
    (OUT/'exact-kernel.json').write_text(json.dumps(rows,indent=2)+'\n')


def verify_all():
    import hashlib,importlib.metadata,platform,subprocess,time
    start=time.monotonic()
    exact_leading();exact_kernel()
    numerical=[];summaries=[]
    for L in [2,3]:
        for amplitude in [0.01,0.05,0.1]:
            rank=normal_rank=positive=negative=dimension=0
            max_difference=0.
            for m in range(L+1):
                a=numerical_block(L,m,amplitude,64)
                b=numerical_block(L,m,amplitude,128)
                difference=float(np.max(np.abs(np.asarray(a['Dgram_real'])-np.asarray(b['Dgram_real']))))
                difference=max(difference,float(np.max(np.abs(np.asarray(a['Dgram_imag'])-np.asarray(b['Dgram_imag'])))))
                assert difference<1e-10
                max_difference=max(max_difference,difference)
                expected=3 if m==L else 4
                for z in [a,b]:
                    for tol in ['1e-08','1e-10','1e-12','1e-14']:
                        assert z['rank_at_absolute_threshold'][tol]=={'DeWitt':expected,'normal_map':expected}
                multiplier=1 if m==0 else 2
                rank+=multiplier*expected;normal_rank+=multiplier*expected
                positive+=multiplier*b['spatial_Gram_positive']
                negative+=multiplier*b['spatial_Gram_negative']
                dimension+=multiplier*b['spatial_dimension_complex']
                numerical.extend([a,b])
            N=L+1
            assert rank==8*L+2
            assert dimension==6*N*N-10 and positive==5*N*N-10 and negative==N*N
            summaries.append({'N':N,'L':L,'epsilon_over_r_squared':amplitude,
                              'shift_Hessian_rank':rank,'normal_map_rank':normal_rank,
                              'spatial_dimension':dimension,'spatial_positive':positive,'spatial_negative':negative,
                              'quadrature_matrix_max_difference':max_difference})
            print('ALL_MODES',N,amplitude,'shift_rank',rank,'spatial_signature',positive,negative,flush=True)
    result={'normalization':'normalized ADM action; r=1 numerical blocks; arbitrary radial derivatives removed only as exact tangent terms',
            'exact_small_nonzero_deformation_results':[
                {'N':L+1,'L':L,'temporal_coefficients':4*(L+1)**2-2,
                 'shift_coefficients':3*(L+1)**2-2,'shift_kernel_dimension':3*L*L-2*L-1,
                 'temporal_Hessian_rank':8*L+2,'primary_plus_secondary_Poisson_rank':16*L+4,
                 'final_Dirac_rank_established':False} for L in [2,3]],
            'sampled_full_spatial_and_shift_matrices':summaries,
            'numerical_blocks':numerical,
            'kernel_description':'b=lambda*f, f in V_(L-2), plus angular W with W(x) in V_(L-1)',
            'secondary_bracket_rank_uses_analytic_flat_Bianchi_identity':True,
            'regular_surrounding_constraint_stratum_proved':False,
            'physical_DOF_count_established':False,'full_nonlinear_goal_complete':False}
    (OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
    sources=[Path(__file__),Path(__file__).with_name('leading_certificate.py'),
             ROOT/'sections/nonlinear_action.tex']
    manifest={'base_commit_at_execution':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'python':platform.python_version(),'platform':platform.platform(),
              'dependencies':{k:importlib.metadata.version(k) for k in ['numpy','sympy']},
              'seconds':round(time.monotonic()-start,3),'GPU_jobs':0,'experiment_launches':0,
              'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (OUT/'reproducibility.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PASS: exact leading positive minors and matching exact nullspace dimensions certify shift ranks18/26.',flush=True)
    print('PASS: independent all-m full spatial signatures and finite-deformation shift ranks agree on two quadratures.',flush=True)
    print('Primary-plus-secondary Poisson ranks36/52 use the reported analytic mixed-bracket reduction.',flush=True)
    print('Further preservation, final Dirac rank, physical DOF, and constant-rank neighborhood remain open.',flush=True)


if __name__=='__main__':
    if '--verify' in sys.argv:verify_all()
    elif '--exact-leading' in sys.argv:exact_leading()
    else:exploratory()
