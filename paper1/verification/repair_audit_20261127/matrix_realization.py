#!/usr/bin/env python3
"""Type/reality checks and exact quadratic EH harmonic realization at finite N.

This constructs the bilinear matrix representation after the EH constraint
reduction. It does not supply a nonlinear noncommutative diffeomorphism algebra.
"""
import json,math,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve();sys.path.insert(0,str(HERE.parents[1]))
from p2_frontier1_sun_fuzzy_action_and_saddles import build_fuzzy_monopole_basis
from v2_test1_matrix_vertex_and_projected_dynamics import iota_map
OUT=HERE.parents[2]/'reports/repair-audit-2026-11-27'

def inner(X,Y,s,N):return np.trace(X.conj().T@Y)/math.sqrt(N*N-s*s)
def Smap(X,Y,s,N,b0,bm,bp):
    return ((N-s)*iota_map(X.conj().T@Y,bm,b0)+(N+s)*iota_map(Y@X.conj().T,bp,b0))/(2*math.sqrt(N*N-s*s))
def scalar_action(F,X,s,N,b0,bm,bp):
    return (iota_map(F,b0,bp)@X+X@iota_map(F,b0,bm))/2
def hermitian(rng,N):
    x=rng.normal(size=(N,N))+1j*rng.normal(size=(N,N))
    return (x+x.conj().T)/2

def run():
    rng=np.random.default_rng(104729);results=[]
    for N in [4,6,8]:
        b0=build_fuzzy_monopole_basis(N,0)
        for s in [1,2]:
            b=build_fuzzy_monopole_basis(N,s)
            bm=build_fuzzy_monopole_basis(N-s,0);bp=build_fuzzy_monopole_basis(N+s,0)
            X=rng.normal(size=(N+s,N-s))+1j*rng.normal(size=(N+s,N-s))
            Y=rng.normal(size=X.shape)+1j*rng.normal(size=X.shape)
            F=hermitian(rng,N)
            selfadjoint_error=abs(inner(X,scalar_action(F,Y,s,N,b0,bm,bp),s,N)-inner(scalar_action(F,X,s,N,b0,bm,bp),Y,s,N))
            contract=Smap(X,X,s,N,b0,bm,bp)
            contract_error=np.linalg.norm(contract-contract.conj().T)
            identity_error=np.linalg.norm(scalar_action(np.eye(N),X,s,N,b0,bm,bp)-X)
            assert max(selfadjoint_error,contract_error,identity_error)<1e-10
            results.append({'N':N,'spin':s,'rectangular_shape':list(X.shape),
                            'scalar_action_self_adjoint_error':selfadjoint_error,
                            'quadratic_contraction_hermiticity_error':float(contract_error),
                            'scalar_identity_action_error':float(identity_error)})
    # The constraint-reduced EH master fields can be represented on B_2 by
    # a unitary harmonic-coefficient map. Test bilinear, radial/time derivative,
    # and both potential terms with independently randomized coefficients.
    quadratic=[]
    for N in [4,6,8]:
        b=build_fuzzy_monopole_basis(N,2)['basis'];labels=list(b)
        coeff=rng.normal(size=len(labels))+1j*rng.normal(size=len(labels))
        ct=rng.normal(size=len(labels))+1j*rng.normal(size=len(labels))
        cr=rng.normal(size=len(labels))+1j*rng.normal(size=len(labels))
        def lift(c):return sum((v*b[key] for key,v in zip(labels,c)),np.zeros((N+2,N-2),complex))
        Z,Zt,Zr=lift(coeff),lift(ct),lift(cr)
        mass=1.;rad=5.;f=1-2*mass/rad
        rw=[];ze=[]
        for j,m in labels:
            lam=(j-1)*(j+2)/2
            rw.append(f*(j*(j+1)/rad**2-6*mass/rad**3))
            ze.append(2*f*(lam**2*(lam+1)*rad**3+3*lam**2*mass*rad**2+9*lam*mass**2*rad+9*mass**3)/(rad**3*(lam*rad+3*mass)**2))
        errs={}
        for name,V in [('axial',np.array(rw)),('polar',np.array(ze))]:
            matrix=.5*(inner(Zt,Zt,2,N)-inner(Zr,Zr,2,N)-inner(Z,lift(V*coeff),2,N))
            spectral=.5*(np.vdot(ct,ct)-np.vdot(cr,cr)-np.vdot(coeff,V*coeff))
            errs[name]=float(abs(matrix-spectral));assert errs[name]<1e-10
        quadratic.append({'N':N,'harmonics':len(labels),'action_difference':errs})
    return {'rectangular_maps':results,'quadratic_EH_harmonic_realization':quadratic,
            'scope':'Exact bilinear realization of the derived EH perturbation sectors. Nonlinear completion and boundary/constraint algebra remain unresolved.'}

if __name__=='__main__':
    result=run();OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'matrix_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
