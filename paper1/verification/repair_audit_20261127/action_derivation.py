#!/usr/bin/env python3
"""Independent bulk-action checks. Does not edit the manuscript or legacy results.

Distinguishes the printed action, a minimal algebraic patch, and a replacement
derived from the Einstein-Hilbert spherical reduction. No full nonlinear fuzzy
Einstein equivalence is assumed or claimed.
"""
import json
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parents[2] / 'reports/repair-audit-2026-11-27'
t, r = sp.symbols('t r', real=True)
M, cosm, j = sp.symbols('M Lambda j', real=True)
X, N, F = (sp.Function(s)(t, r) for s in ('X', 'N', 'F'))

def simp(x):
    return sp.factor(sp.simplify(x.doit()))

def el_first(L, field):
    return simp(sp.diff(L, field) - sum(sp.diff(sp.diff(L, sp.diff(field, c)), c) for c in (t, r)))

def curvature_2d():
    coords = (t, r)
    metric = sp.Matrix([[-F, -N], [-N, 0]])
    inv = metric.inv()
    G = [[[simp(sum(inv[a, d]*(sp.diff(metric[d,c], coords[b]) + sp.diff(metric[d,b], coords[c]) - sp.diff(metric[b,c], coords[d])) for d in range(2))/2) for c in range(2)] for b in range(2)] for a in range(2)]
    Ric = sp.zeros(2)
    for b in range(2):
        for d in range(2):
            Ric[b,d] = simp(sum(sp.diff(G[a][b][d], coords[a])-sp.diff(G[a][a][d], coords[b]) + sum(G[a][a][e]*G[e][b][d]-G[a][b][e]*G[e][a][d] for e in range(2)) for a in range(2)))
    return simp(sum(inv[a,b]*Ric[a,b] for a in range(2) for b in range(2)))

def spherical():
    Rt = curvature_2d()
    expected_R = -sp.diff((sp.diff(F,r)-2*sp.diff(N,t))/N,r)/N
    assert simp(Rt-expected_R) == 0
    # Exact spherical EH reduction up to radial boundary derivatives:
    # sqrt(-h)[X R[h] + (grad X)^2/(2X) + 2 - 2 Lambda X].
    kinetic = -sp.diff(X,t)*sp.diff(X,r)/X + F*sp.diff(X,r)**2/(2*N*X)
    L_cov = N*X*Rt + kinetic + 2*N - 2*cosm*N*X
    L = sp.diff(X,r)*(sp.diff(F,r)-2*sp.diff(N,t))/N + kinetic + 2*N - 2*cosm*N*X
    boundary = -sp.diff(X*(sp.diff(F,r)-2*sp.diff(N,t))/N,r)
    assert simp(L_cov-L-boundary)==0
    eq = {name: el_first(L, field) for name,field in [('X',X),('N',N),('F',F)]}
    f = 1-2*M/r-cosm*r*r/3
    bg = {X:r*r, N:1, F:f}
    residuals = {name:simp(val.subs(bg)) for name,val in eq.items()}
    assert all(v==0 for v in residuals.values()), residuals
    # Gauge fixing h_rr=0 before variation omits one metric equation. Its
    # mass-conservation content must be supplied from the unreduced Einstein
    # equations; do not claim these three equations prove Birkhoff's theorem.
    mass_t=sp.Function('mass')(t)
    evolving = {X:r*r,N:1,F:1-2*mass_t/r-cosm*r*r/3}
    evolving_res = {name:simp(val.subs(evolving)) for name,val in eq.items()}
    # Restore the metric variation removed by h_rr=0 gauge fixing. These
    # covariant 2D equations follow by varying before choosing null gauge.
    fn=1-2*mass_t/r-cosm*r*r/3
    metric=sp.Matrix([[-fn,-1],[-1,0]]);inv=metric.inv();coords=(t,r)
    G=[[[simp(sum(inv[a,d]*(sp.diff(metric[d,c],coords[b])+sp.diff(metric[d,b],coords[c])-sp.diff(metric[b,c],coords[d])) for d in range(2))/2) for c in range(2)] for b in range(2)] for a in range(2)]
    dX=sp.Matrix([0,2*r])
    Hess=sp.Matrix(2,2,lambda a,b:sp.diff(r*r,coords[a],coords[b])-sum(G[c][a][b]*dX[c] for c in range(2)))
    boxX=sum(inv[a,b]*Hess[a,b] for a in range(2) for b in range(2))
    grad2=(dX.T*inv*dX)[0]
    Emetric=sp.Matrix(2,2,lambda a,b:simp(metric[a,b]*boxX-Hess[a,b]+dX[a]*dX[b]/(2*r*r)-metric[a,b]*grad2/(4*r*r)-metric[a,b]*(1-cosm*r*r)))
    assert simp(Emetric[0,0]+2*sp.diff(mass_t,t))==0
    assert Emetric[0,1]==Emetric[1,0]==Emetric[1,1]==0
    # Minimal reviewer patch: rho^2 N R_q -> 2 N, with all other spherical
    # terms unchanged. Integrate 2 X d_r(theta) by parts first.
    theta=sp.diff(X,t)/X-F*sp.diff(X,r)/(2*X)
    L_patch=-2*sp.diff(X,r)*theta + 2*N-2*cosm*N*X
    patch_eq={name:el_first(L_patch,field) for name,field in [('X',X),('N',N),('F',F)]}
    patch_res={name:simp(val.subs(bg)) for name,val in patch_eq.items()}
    assert patch_res['F']==4
    assert simp(patch_res['N']-(2-2*cosm*r*r))==0
    # Original B.11a-b scalar sector, affine gauge N=1, matches the new L.
    L_B11=-sp.diff(X,t)*sp.diff(X,r)/X + 2-2*cosm*X + sp.diff(X,r)*sp.diff(F,r) + F*sp.diff(X,r)**2/(2*X)
    assert simp(L.subs(N,1)-L_B11)==0
    return {'two_dimensional_curvature':str(Rt), 'replacement_bulk_density':str(L),
            'euler_lagrange':{k:str(v) for k,v in eq.items()},
            'schwarzschild_de_sitter_residuals':{k:str(v) for k,v in residuals.items()},
            'minimal_patch_residuals':{k:str(v) for k,v in patch_res.items()},
            'time_dependent_mass_gauge_fixed_residuals':{k:str(v) for k,v in evolving_res.items()},
            'restored_covariant_metric_equation_tt':str(Emetric[0,0]),
            'omitted_equation_warning':'Null gauge-fixed variation alone does not impose the missing mass-conservation equation.',
            'original_B11_scalar_match':True}

def axial():
    f=1-2*M/r
    ell=(j-1)*(j+2)
    h0,h1,h2,Lambda_g,Psi=(sp.Function(s)(t,r) for s in ['h0','h1','h2','gauge','Psi'])
    C=sp.diff(h1,t)-sp.diff(h0,r)+2*h0/r
    L=(C*C+ell*h0*h0/(r*r*f)-ell*f*h1*h1/(r*r))/2
    e0,e1=el_first(L,h0),el_first(L,h1)
    constraint=sp.diff(h0,t)-f*sp.diff(f*h1,r)
    identity=sp.diff(e0,t)+sp.diff(e1,r)+2*e1/r-ell*constraint/(r*r*f)
    assert simp(identity)==0
    # Eliminate the h0 constraint, then define Psi=f h1/r.
    h0t=f*sp.diff(f*h1,r)
    reduced=e1.subs(sp.diff(h0,t,r),sp.diff(h0t,r)).subs(sp.diff(h0,t),h0t)
    reduced=simp(reduced.subs(h1,r*Psi/f))
    V=f*(j*(j+1)/(r*r)-6*M/r**3)
    RW=-sp.diff(Psi,t,2)+f*sp.diff(f*sp.diff(Psi,r),r)-V*Psi
    assert simp(reduced-f**-1*r*RW)==0
    # Standard odd gauge completion, before Regge-Wheeler gauge h2=0.
    k0=h0-sp.diff(h2,t)/2
    k1=h1-sp.diff(h2,r)/2+h2/r
    trans={h0:h0-sp.diff(Lambda_g,t),h1:h1-sp.diff(Lambda_g,r)+2*Lambda_g/r,h2:h2-2*Lambda_g}
    assert simp(k0.subs(trans, simultaneous=True)-k0)==0
    assert simp(k1.subs(trans, simultaneous=True)-k1)==0
    return {'quadratic_density_up_to_positive_harmonic_prefactor':str(L),
            'constraint_identity_residual':'0','regge_wheeler_residual':'0',
            'derived_potential':str(V),'gauge_invariant_amplitudes':[str(k0),str(k1)],
            'scope':'Correct EH quadratic axial sector. Equivalence to the full proposed nonlinear matrix action is not proved.'}

if __name__=='__main__':
    OUT.mkdir(exist_ok=True,parents=True)
    result={'spherical':spherical(),'axial':axial()}
    (OUT/'action_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
