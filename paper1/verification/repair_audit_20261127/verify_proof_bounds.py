#!/usr/bin/env python3
"""Independent numerical/algebraic checks of lemmas used in the clock proof."""
from pathlib import Path
import json,sys
import numpy as np
import sympy as sp
from scipy.linalg import eigh
HERE=Path(__file__).resolve();sys.path.insert(0,str(HERE.parent))
from clock_derivation import quadrature,build_massive_clock_matrices
OUT=HERE.parents[2]/'reports/repair-audit-2026-11-27'

def root(A):
    e,u=eigh(A)
    return (u*np.sqrt(np.maximum(e,0)))@u.T

rng=np.random.default_rng(65431)
worst=0.
for _ in range(30):
    a,b=rng.normal(size=(6,6)),rng.normal(size=(6,6))
    P=a@a.T+.1*np.eye(6);Q=b@b.T+.1*np.eye(6)
    E=.5*np.linalg.svd(root(P)@root(Q),compute_uv=False).sum()-.25*np.trace(P+Q)
    lower=-.25*np.linalg.svd(P-Q,compute_uv=False).sum()
    assert E>=lower-1e-12
    worst=max(worst,float(lower-E))

quad=quadrature(1000)
Kbound=json.loads((OUT/'clock_results.json').read_text())['independent_lower_bounds']['uniform_mass_and_radial_bound']['conservative_K_g_numerical']
rows=[]
for mu in [0.,.001,.1,1.,10.,100.,10000.]:
    for nr in [5,40,120]:
        Q,P,c=build_massive_clock_matrices(nr,mu,*quad)
        K=(Q-P)/2
        ks=np.arange(1,nr+1,dtype=float)
        k,l=ks[:,None],ks[None,:]
        D=np.sqrt((k*k+mu*mu)*(l*l+mu*mu))
        stable=-mu*mu*(k-l)**2*c[np.abs(k-l).astype(int)]/(np.sqrt(D)*(D+k*l+mu*mu))
        stable-=(D+k*l-mu*mu)*c[(k+l).astype(int)]/np.sqrt(D)
        err=float(np.max(np.abs(K-stable)))
        assert err<1e-8
        entry=float(np.sum(abs(stable)))
        assert entry<=2*Kbound+1e-9
        rows.append({'mu':mu,'radial_modes':nr,'entrywise_kernel_bound_value':entry,'kernel_identity_error':err})

J=sp.symbols('J',positive=True)
assert sp.factor(sp.Rational(8,3)-2*((J+1)**2-4)/J**2-sp.Rational(2,3)*(J-3)**2/J**2)==0
result={'positive_matrix_trace_bound_trials':30,'kernel_checks':rows,
        'area_rank_bound_symbolic_residual':'0',
        'scope':'Sanity checks supporting the written proof; finite samples are not substituted for the analytic argument.'}
(OUT/'proof_bound_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('30 noncommuting positive-matrix tests; 21 mass/cutoff tests; exact angular-rank identity: passed.')
