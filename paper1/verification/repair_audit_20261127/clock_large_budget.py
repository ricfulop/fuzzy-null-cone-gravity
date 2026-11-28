#!/usr/bin/env python3
"""Exact finite-budget second moments without constructing the many-body space.

Xi|0> has only two-particle components. Their norm is an explicit masked sum,
with the true pair cost omega_k+omega_l<=B. FFT quadrature supplies smooth
Fourier coefficients; low-mode Gauss and doubled-grid checks are recorded.
"""
import math,json,sys,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve();sys.path.insert(0,str(HERE.parent))
from clock_derivation import quadrature
OUT=HERE.parents[2]/'reports/repair-audit-2026-11-27'

def fourier(grid):
    u=np.arange(grid,dtype=float)/grid;f=np.zeros(grid)
    f[1:]=np.exp(2-1/(2*u[1:]*(1-u[1:])))
    return np.fft.rfft(f).real/grid

def moment(s,B,J,c):
    accum=0.;nr=math.floor(B);ks=np.arange(1,nr+1,dtype=float)
    for j in range(2,J+1):
        t2=(j-1)*(j+2)/2
        if 4*math.sqrt(t2)>B:continue
        cost=ks+t2/ks;mu=math.sqrt(t2)/(2*math.pi)
        lk=ks[None,:];lc=cost[None,:];oml=np.sqrt(lk*lk+mu*mu)
        subtotal=0.
        for start in range(0,nr,96):
            kk=ks[start:start+96,None];kc=cost[start:start+96,None]
            omk=np.sqrt(kk*kk+mu*mu);D=omk*oml;root=np.sqrt(D)
            diff=np.abs(kk-lk).astype(int);sm=(kk+lk).astype(int)
            K=-mu*mu*(kk-lk)**2*c[diff]/(root*(D+kk*lk+mu*mu))
            K-=(D+kk*lk-mu*mu)*c[sm]/root
            K=np.where(kc+lc<=B+1e-12,K,0.)
            subtotal+=float(np.sum(K*K))
        accum+=(2*j+1)*subtotal # degeneracy/2 = (2j+1)
    return accum/(s*s)

def main():
    c=fourier(32768);c2=fourier(65536)
    u,w,f=quadrature(1600);ns=np.arange(65)
    cg=np.cos(2*math.pi*ns[:,None]*u)@(w*f)
    gauss_error=float(np.max(abs(c[:65]-cg)))
    assert gauss_error<1e-12
    rows=[];tic=time.time()
    for s in [8.,12.,16.,18.,20.,40.,80.,160.,240.,5000.]:
        B=s/2;J=math.floor(math.sqrt(s/math.log(3)))
        v=moment(s,B,J,c)
        # At large budget a second full sum checks FFT coefficient sensitivity.
        v2=moment(s,B,J,c2)
        assert abs(v-v2)<1e-12*max(v,1e-10)
        row={'s':s,'B':B,'J_max':J,'Xi_variance':v,
             'vacuum_failure_probability_upper_bound_at_log2':v/math.log(2)**2,
             'vacuum_clock_acceptance_lower_bound_at_log2':max(0.,1-v/math.log(2)**2),
             'FFT_grid_doubling_absolute_change':abs(v-v2)}
        rows.append(row);print(json.dumps(row),flush=True)
    existing=json.loads((OUT/'clock_results.json').read_text())
    for row in rows[:5]:
        match=next(x for x in existing['occupation_projected_clock'] if x['s']==row['s'])
        assert abs(row['Xi_variance']-match['Xi_second_moment'])<1e-12
    result={'rows':rows,'low_mode_FFT_vs_Gauss_max_error':gauss_error,
            'runtime_seconds':time.time()-tic,
            'scope':'Second moment and Chebyshev bound for the compressed quadratic clock, not exponential moments or a nonlinear gravitational clock.'}
    (OUT/'clock_large_budget_results.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
