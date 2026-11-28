#!/usr/bin/env python3
"""Exact spin-two kinetic test: angular truncation and Legendre transform.

For L=1/2<dot(q), alpha^(-1) dot(q)> in a retained harmonic band,
H=1/2 p^T (P alpha^(-1) P)^(-1) p. This differs from restricting
the continuum Hamiltonian 1/2<p,alpha p>. The check isolates the positive
traceless angular kinetic block, not the full ADM constraint surface.
"""
from pathlib import Path
import json
import numpy as np
import sympy as s

OUT=Path(__file__).resolve().parents[2]/'reports/nonlinear-consistency-2026-11-27'
x=s.symbols('x',real=True)
u=s.legendre(2,x)


def inner(a,b):return s.simplify(s.integrate(s.expand(a*b),(x,-1,1))/2)


def main():
    rows=[]
    nodes,weights=np.polynomial.legendre.leggauss(100)
    for L in [2,3,4,5]:
        basis=[s.sqrt((2*j+1)*s.factorial(j-2)/s.factorial(j+2))*s.assoc_legendre(j,2,x)
               for j in range(2,L+1)]
        I=s.Matrix([[inner(a,b) for b in basis] for a in basis])
        assert I==s.eye(len(basis))
        M=s.Matrix([[inner(a,u*b) for b in basis] for a in basis])
        T=s.Matrix([[inner(a,u*u*b) for b in basis] for a in basis])
        leakage=s.simplify(T-M*M)
        top=s.simplify(leakage[-1,-1]);assert top>0
        if L==2:assert top==s.Rational(3,49)
        B=np.stack([np.asarray(s.lambdify(x,b,'numpy')(nodes),float) for b in basis])
        uq=s.lambdify(x,u,'numpy')(nodes)
        p=np.zeros(len(basis));p[-1]=1
        def ham(e):
            K=(B*(weights/(2*(1+e*uq))))@B.T
            return float(p@np.linalg.solve(K,p)/2)
        step=0.0002
        curvature=(ham(step)-2*ham(0)+ham(-step))/(step*step)
        assert abs(curvature+float(top))<2e-7
        exact_M=np.array(M,dtype=float)
        finite=[]
        for e in [-0.2,0.1,0.3]:
            lag=ham(e)
            can=float(p@(np.eye(len(basis))+e*exact_M)@p/2)
            assert lag<can
            finite.append({'epsilon':e,'truncate_L_then_Legendre':lag,
                           'Legendre_then_truncate':can,'difference':lag-can})
        rows.append({'L':L,'matrix_rank_N':L+1,'top_mode_leakage_norm_squared':str(top),
                     'exact_lapse_Hamiltonian_second_derivative':str(-top),
                     'finite_difference_second_derivative':curvature,
                     'numerical_error':abs(curvature+float(top)),'finite_amplitude_checks':finite})
    result={'checks':rows,
            'series_identity':'(P alpha^-1 P)^-1=I+epsilon M+epsilon^2(M^2-T)+O(epsilon^3), M=P u P, T=P u^2 P.',
            'exact_identity':'(P alpha^-1 P)^-1=P alpha P-P alpha Q(Q alpha Q)^-1 Q alpha P, on ran(P), for positive invertible alpha.',
            'conclusion':'Velocity-space truncation generates nonlinear dependence on nonconstant lapse harmonics, even though the continuum kinetic Hamiltonian is linear in lapse.',
            'scope':'Positive traceless angular kinetic block at fixed spatial metric, with a lapse-harmonic ansatz. This does not classify the full constraints of the current metric-coefficient Galerkin action, whose ADM lapse is also field-dependent. It proves the two truncation procedures cannot be identified without analysis.'}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'legendre_truncation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
