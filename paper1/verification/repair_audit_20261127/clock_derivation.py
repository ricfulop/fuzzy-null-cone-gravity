#!/usr/bin/env python3
"""Occupation-budget clock and independent quadratic-operator lower bounds.

Every projected calculation uses sum(n_i*(k_i+t_j^2/k_i)) <= B.
For B<12, vacuum-connected states have zero or two quanta, all in one
angular/polarization channel. Identical channels are combined into their exact
permutation-symmetric bright state, with sqrt(degeneracy) vacuum coupling.
This is an exact reduction of the vacuum spectral measure, not a Fock cutoff.
"""
import json, math, sys
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import simpson
from scipy.optimize import minimize_scalar
from scipy.special import zeta
import sympy as sp

HERE=Path(__file__).resolve()
sys.path.insert(0,str(HERE.parents[1]))
from v2_test2_fewster_eveson_and_exponential_clock import build_massive_clock_matrices
OUT=HERE.parents[2]/'reports/repair-audit-2026-11-27'

def quadrature(n=1600):
    x,w=np.polynomial.legendre.leggauss(n)
    u=(x+1)/2; w=w/2
    g=np.exp(-1/(4*u*(1-u))+1)
    return u,w,g*g

def pair_basis(cost,B):
    return [(k,l) for k in range(len(cost)) for l in range(k,len(cost)) if cost[k]+cost[l]<=B]

def pair_number_matrix(A,pairs):
    index={p:i for i,p in enumerate(pairs)}
    H=np.zeros((len(pairs),len(pairs)))
    for col,(k,l) in enumerate(pairs):
        occupation=np.zeros(A.shape[0],dtype=int); occupation[k]+=1; occupation[l]+=1
        for b in np.flatnonzero(occupation):
            temp=occupation.copy(); amp=math.sqrt(temp[b]);temp[b]-=1
            for a in range(A.shape[0]):
                final=temp.copy(); factor=amp*math.sqrt(final[a]+1);final[a]+=1
                state=tuple(np.repeat(np.arange(len(final)),final))
                row=index.get(state)
                if row is not None:H[row,col]+=A[a,b]*factor
    assert np.max(np.abs(H-H.T),initial=0)<1e-10
    return H

def projected_clock(B,s,Jmax,quad):
    assert B<12, 'Bright-block reduction here is only implemented below four-particle threshold.'
    blocks=[]; total_single_modes=0; full_pair_dim=0
    for j in range(2,Jmax+1):
        tj2=(j-1)*(j+2)//2
        ks=[k for k in range(1,math.floor(B)+1) if Fraction(k)+Fraction(tj2,k)<=B]
        if not ks:continue
        deg=2*(2*j+1); total_single_modes+=deg*len(ks)
        Mq,Mp,_=build_massive_clock_matrices(max(ks),math.sqrt(tj2)/(2*math.pi),*quad)
        ids=np.array(ks)-1; Mq=Mq[np.ix_(ids,ids)]; Mp=Mp[np.ix_(ids,ids)]
        A=(Mq+Mp)/2; K=(Mq-Mp)/2
        pairs=pair_basis([Fraction(k)+Fraction(tj2,k) for k in ks],B)
        if not pairs:continue
        block=pair_number_matrix(A,pairs)
        v=np.array([K[k,l]*(1/math.sqrt(2) if k==l else 1) for k,l in pairs])
        blocks.append((block,math.sqrt(deg)*v,j,pairs))
        full_pair_dim+=deg*len(pairs)
    dim=1+sum(len(v) for _,v,_,_ in blocks)
    H=np.zeros((dim,dim));off=1
    for block,v,_,_ in blocks:
        n=len(v);H[off:off+n,off:off+n]=block
        H[0,off:off+n]=v;H[off:off+n,0]=v;off+=n
    Xi=H/s
    ev,U=eigh(Xi);weights=U[0,:]**2
    def mgf(lam):return float(np.dot(weights,np.exp(lam*ev)))
    mean=mgf(1);variance=mgf(2)-mean**2
    m2=float(Xi[0]@Xi[:,0]);m4=float(np.linalg.norm(Xi@Xi[:,0])**2)
    eta=math.log(2);failure=float(np.sum(weights[ev>eta]))
    return {'B':B,'s':s,'J_max':Jmax,'single_particle_modes':total_single_modes,
            'vacuum_block_dimension_before_channel_symmetry':1+full_pair_dim,
            'exact_bright_block_dimension':dim,
            'vacuum_spectral_min':float(ev.min()),'vacuum_spectral_max':float(ev.max()),
            'mean_exp_Xi':mean,'variance_exp_Xi':max(0,variance),
            'Xi_second_moment':m2,'Xi_fourth_moment':m4,
            'exact_vacuum_failure_at_log2':failure,
            'second_moment_failure_bound':min(1,m2/eta**2),
            'fourth_moment_failure_bound':min(1,m4/eta**4),
            'Hermitian_residual':float(np.max(np.abs(Xi-Xi.T),initial=0))}

def independent_bounds(quad):
    q1=lambda x:math.sqrt(1-x**-2)+math.acosh(x)/x**2
    peak=minimize_scalar(lambda x:-q1(x),bounds=(1.00001,10),method='bounded')
    u,w,g2=quad
    g=np.sqrt(g2);gp=g*(1-2*u)/(4*u*u*(1-u)**2)
    int_gp=float(w@gp**2)
    x=sp.symbols('x',positive=True)
    f=sp.exp(2-1/(2*x*(1-x)))
    f4=sp.lambdify(x,sp.diff(f,x,4),'numpy')
    C4=float(w@np.abs(f4(u)))/(2*math.pi)**4
    ns=np.arange(1,81,dtype=float)
    coeff=np.cos(2*math.pi*ns[:,None]*u[None,:])@(w*g2)
    partial=float(np.sum((2*ns**2-ns)*np.abs(coeff)))
    # |c_n| <= C4/n^4 gives a genuine analytic tail expression. C4 and
    # the finite coefficients here are numerical quadratures, not intervals.
    tail=C4*(2*zeta(2,81)-zeta(3,81))
    uniform_floor=-.5*(partial+tail)
    table=[]
    for j in [2,5,10,20]:
        mu=math.sqrt((j-1)*(j+2)/2)/(2*math.pi)
        analytic_bound=-.5*C4*(mu*mu*zeta(3)*zeta(2)+zeta(2)-zeta(3)+2*mu*(zeta(3)-zeta(4)))
        for nr in [20,60,120,240]:
            Mq,Mp,_=build_massive_clock_matrices(nr,mu,*quad)
            A=(Mq+Mp)/2;K=(Mq-Mp)/2
            trace_lower=-.5*float(np.sum(np.abs(eigh(K,eigvals_only=True))))
            eq,Uq=eigh(Mq);root=(Uq*np.sqrt(np.maximum(eq,0)))@Uq.T
            freq2=eigh((root@Mp@root+(root@Mp@root).T)/2,eigvals_only=True)
            emin=.5*float(np.sum(np.sqrt(np.maximum(freq2,0))))-.5*float(np.trace(A))
            assert emin>=trace_lower-1e-6
            assert trace_lower>=analytic_bound-1e-6
            table.append({'j':j,'mu':mu,'radial_modes':nr,'numerical_Williamson_minimum':emin,
                          'finite_cutoff_trace_norm_lower_bound':trace_lower,
                          'analytic_radial_uniform_lower_bound_numerically_evaluated':analytic_bound,
                          'mass_and_radial_uniform_lower_bound_numerically_evaluated':uniform_floor,
                          'pair_Hilbert_Schmidt_squared':float(np.sum(K*K))})
            assert trace_lower>=uniform_floor-1e-6
    return {'FE_Q1_maximum':-peak.fun,'FE_Q1_maximum_x':peak.x,
            'int_g_prime_squared':int_gp,
            'FE_massless_T00_standard_normalization':-int_gp/(4*math.pi),
            'FE_mass_uniform_T00_envelope':1.2*(-int_gp/(4*math.pi)),
            'normalization_warning':'Standard timelike T00 theorem is not identified with the current spatial/canonical clock matrix.',
            'C4_numerical':C4,
            'uniform_mass_and_radial_bound':{
                'formula':'H_NO >= -K_g, K_g=1/2 sum_{n>=1}(2n^2-n)|c_n|',
                'partial_weighted_sum_through_80':partial,
                'analytic_tail_expression_evaluated_numerically':tail,
                'conservative_K_g_numerical':-uniform_floor,
                'infinite_series_K_g_estimate':partial/2,
                'rank_rule':'J <= floor(sqrt(s/log(3))); D(J)=2*((J+1)^2-4)',
                'Xi_lower_bound_for_all_areas_rank_rule_numerical':8*uniform_floor/(3*math.log(3)),
                'area_bound_multiplier_all_areas_upper_estimate':math.exp(-16*uniform_floor/(3*math.log(3))),
                'scope':'Analytic bound for the specified oscillator matrices, all mu>=0 and radial truncations. Not a scalar-field QEI attribution.'},
            'analytic_radial_bound':'-C4/2 * [mu^2 zeta(3) zeta(2) + zeta(2)-zeta(3) + 2 mu (zeta(3)-zeta(4))]',
            'radial_scans':table}

def full_fock_validation(B,quad):
    """Independent explicit occupation-state construction for two channels j=2.
    Checks the bright-sector reduction against the full budget sector.
    """
    ks=[k for k in range(1,math.floor(B)+1) if Fraction(k)+Fraction(2,k)<=B]
    cost=[Fraction(k)+Fraction(2,k) for k in ks]*2
    states=[]
    def rec(i,remaining,v):
        if i==len(cost):states.append(tuple(v));return
        for n in range(int(remaining//cost[i])+1):rec(i+1,remaining-n*cost[i],v+[n])
    rec(0,Fraction(B),[])
    index={v:i for i,v in enumerate(states)}
    Mq,Mp,_=build_massive_clock_matrices(max(ks),math.sqrt(2)/(2*math.pi),*quad)
    ids=np.array(ks)-1
    A=((Mq+Mp)/2)[np.ix_(ids,ids)];K=((Mq-Mp)/2)[np.ix_(ids,ids)]
    H=np.zeros((len(states),len(states)))
    def put(v,col,amp):
        row=index.get(tuple(v))
        if row is not None:H[row,col]+=amp
    for col,state in enumerate(states):
        for channel in range(2):
            off=channel*len(ks)
            for a in range(len(ks)):
                for b in range(len(ks)):
                    v=list(state)
                    if v[off+b]:
                        amp=math.sqrt(v[off+b]);v[off+b]-=1
                        amp*=math.sqrt(v[off+a]+1);v[off+a]+=1
                        put(v,col,A[a,b]*amp)
                    v=list(state);amp=math.sqrt(v[off+b]+1);v[off+b]+=1
                    amp*=math.sqrt(v[off+a]+1);v[off+a]+=1
                    put(v,col,.5*K[a,b]*amp)
                    v=list(state)
                    if v[off+b]:
                        amp=math.sqrt(v[off+b]);v[off+b]-=1
                        if v[off+a]:
                            amp*=math.sqrt(v[off+a]);v[off+a]-=1
                            put(v,col,.5*K[a,b]*amp)
    assert np.max(np.abs(H-H.T))<1e-12
    ev,U=eigh(H/(2*B));z=index[(0,)*len(cost)]
    m=float(np.dot(U[z,:]**2,np.exp(ev)))
    pairs=pair_basis(cost[:len(ks)],B);hn=pair_number_matrix(A,pairs)
    v=np.array([K[k,l]/math.sqrt(2) if k==l else K[k,l] for k,l in pairs])*math.sqrt(2)
    hb=np.block([[np.zeros((1,1)),v[None,:]],[v[:,None],hn]])/(2*B)
    ee,uu=eigh(hb);mb=float(np.dot(uu[0,:]**2,np.exp(ee)))
    assert abs(m-mb)<1e-12
    return {'B':B,'channels':2,'full_dimension':len(states),'bright_dimension':len(v)+1,'mgf_difference':abs(m-mb)}

if __name__=='__main__':
    OUT.mkdir(exist_ok=True,parents=True)
    quad=quadrature()
    result={'independent_lower_bounds':independent_bounds(quad),
            'occupation_projected_clock':[projected_clock(B,2*B,math.floor(math.sqrt(2*B/math.log(3))),quad) for B in [4,6,8,9,10,11]],
            'full_fock_cross_checks':[full_fock_validation(B,quad) for B in [4,6,8]],
            'scope':'Exact finite-budget vacuum statistics for B<12 with sufficient angular rank. Not large-area estimates.'}
    (OUT/'clock_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'clock':result['occupation_projected_clock'],'checks':result['full_fock_cross_checks'],
                      'Q1max':result['independent_lower_bounds']['FE_Q1_maximum']},indent=2))
