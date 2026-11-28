#!/usr/bin/env python3
"""Finite boundary-algebra obstructions and exact compression defects.

The canonical-commutator and noncompact-unitarity no-go results are trace
arguments, given in the report. These computations verify the constructive
compression identity and a low-excitation diagnostic, not gravity gluing.
"""
from pathlib import Path
import json
import numpy as np
import sympy as s

OUT=Path(__file__).resolve().parents[2]/'reports/nonlinear-consistency-2026-11-27'


def comm(a,b):
    return a@b-b@a


def main():
    oscillator=[]
    for d in [2,4,8,16]:
        a=s.zeros(d)
        for n in range(1,d):a[n-1,n]=s.sqrt(n)
        edge=s.zeros(d);edge[d-1,d-1]=1
        assert comm(a,a.T)==s.eye(d)-d*edge
        q=(a+a.T)/s.sqrt(2)
        p=(a-a.T)/(s.I*s.sqrt(2))
        defect=comm(q,p)-s.I*s.eye(d)
        assert defect==-s.I*d*edge
        oscillator.append({'dimension':d,'trace_commutator':str(s.trace(comm(q,p))),
                           'CCR_defect_operator_norm':d,
                           'CCR_defect_on_vacuum':str(defect[0,0]),
                           'CCR_defect_on_top_state':str(defect[d-1,d-1])})
    rng=np.random.default_rng(27112026)
    checks=[]
    for ambient, retained in [(8,3),(12,5),(20,8)]:
        a=rng.normal(size=(ambient,ambient))+1j*rng.normal(size=(ambient,ambient))
        b=rng.normal(size=(ambient,ambient))+1j*rng.normal(size=(ambient,ambient))
        a=(a+a.conj().T)/2;b=(b+b.conj().T)/2
        P=np.diag([1.]*retained+[0.]*(ambient-retained));Q=np.eye(ambient)-P
        lhs=comm(P@a@P,P@b@P)
        rhs=P@comm(a,b)@P-P@a@Q@b@P+P@b@Q@a@P
        error=float(np.linalg.norm(lhs-rhs))
        assert error<1e-12
        psi=rng.normal(size=ambient)+1j*rng.normal(size=ambient)
        psi=P@psi;psi/=np.linalg.norm(psi)
        defect=lhs-P@comm(a,b)@P
        expectation=float(abs(psi.conj()@defect@psi))
        bound=float(2*np.linalg.norm(Q@a@psi)*np.linalg.norm(Q@b@psi))
        assert expectation<=bound+1e-12
        checks.append({'ambient_dimension':ambient,'retained_dimension':retained,
                       'compression_identity_residual':error,
                       'state_defect_magnitude':expectation,'leakage_bound':bound})
    result={'finite_CCR':'Impossible for d>0 and nonzero hbar: Tr[A,B]=0, Tr(i*hbar*I)=i*hbar*d.',
            'finite_unitary_sl2R':'For Hermitian J,K1,K2 with [J,K1]=iK2, [J,K2]=-iK1, [K1,K2]=-iJ, trace cyclicity forces Tr(K2^2)=-Tr(J^2); all generators vanish.',
            'projected_commutator':'[PAP,PBP]=P[A,B]P-PAQBP+PBQAP; Q=I-P.',
            'expectation_defect_bound':'2 ||QA psi|| ||QB psi|| for normalized psi in ran(P).',
            'oscillator_checks':oscillator,'general_compression_checks':checks,
            'scope':'A finite admissible subspace can carry compressed observables with quantified defects; it cannot preserve the exact canonical pair or a nontrivial unitary sl(2,R) algebra. The physical corner reduction, inner product, and gluing remain to be derived.'}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'boundary_compression.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
