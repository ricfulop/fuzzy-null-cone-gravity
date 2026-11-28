#!/usr/bin/env python3
"""Reproduce constrained reduction of the even Einstein quadratic action.

Input is the continuum Einstein-Hilbert quadratic density in alpha gauge,
Solomon (2023), DOI 10.3390/particles6040061, Eq. (150). It is not Appendix A
of our manuscript. We eliminate the auxiliary fields and derive the potential,
rather than inserting the Zerilli potential into an oscillator calculation.
"""
from pathlib import Path
import json
import sympy as sp
from sympy.calculus.euler import euler_equations

OUT=Path(__file__).resolve().parents[2]/'reports/repair-audit-2026-11-27'
t,r=sp.symbols('t r',real=True)
M,q=sp.symbols('M q',positive=True)
f=1-2*M/r
D=q+1-3*f
h,a,b,alpha,psi=(sp.Function(s)(t,r) for s in ['h','a','b','alpha','psi'])

def compact(x):return sp.factor(sp.cancel(x.doit()))
def EL(L,z):
    result=sp.diff(L,z)
    for der in L.atoms(sp.Derivative):
        if der.expr==z:
            result+=(-1)**len(der.variables)*sp.diff(sp.diff(L,der),*der.variables)
    return compact(result)

# t^a partial_a=-partial_t and r^a partial_a=f partial_r.
L=q*a*a+(q+2)*h*h/4-q*h*b+2*r*a*sp.diff(h,t)+2*r*f*h*sp.diff(b,r)
L-=q*h*sp.diff(r**4*f*alpha,r)/r**2
L+=q*r*r*(r*r*sp.diff(alpha,t)**2+2*f*alpha*alpha)
L+=2*q*(-a*sp.diff(r*r*alpha,t)+2*r*f*b*sp.diff(r*alpha,r)+(1+3*f)*r*b*alpha)
a_sol=sp.solve(EL(L,a),a)[0]
L_a=sp.expand(L.subs(a,a_sol).doit())
alpha_sol=D*psi/(r*r*f)+h/(2*q*r)
L_shift=sp.expand(L_a.subs(alpha,alpha_sol).doit())
constraint=EL(L_shift,b)
h_sol=compact(sp.solve(constraint,h)[0])
print('Auxiliary field and constraint eliminated.',flush=True)
# Work with independent field jets to avoid expanding a large fourth-order
# Euler expression whose terms cancel. Integrate the quadratic mixed
# derivatives by parts explicitly; each rational coefficient is small.
p,pt,pr,ptr,prr=sp.symbols('p pt pr ptr prr')
jets=(p,pt,pr,ptr,prr);zero=dict.fromkeys(jets,0)
H0=compact(h_sol.subs(psi,1).doit())
H1=compact(sp.diff(h_sol,sp.diff(psi,r)))
assert compact(h_sol-H0*psi-H1*sp.diff(psi,r))==0
A0=compact(D/(r*r*f)+H0/(2*q*r));A1=compact(H1/(2*q*r))
hh=H0*p+H1*pr;ht=H0*pt+H1*ptr
aa=A0*p+A1*pr;at=A0*pt+A1*ptr
dar=sp.diff(r**4*f*A0,r)*p+(r**4*f*A0+sp.diff(r**4*f*A1,r))*pr+r**4*f*A1*prr
Ljet=-r*r*(ht-q*r*at)**2/q+(q+2)*hh*hh/4-q*hh*dar/r**2+q*r**4*at*at+2*q*r*r*f*aa*aa
def coef(a,b):return compact(sp.diff(Ljet,a,b).subs(zero)/(2 if a==b else 1))
T=compact(coef(pt,pt)-sp.diff(coef(pt,ptr),r)/2)
d,e=coef(p,prr),coef(pr,prr)
R=compact(coef(pr,pr)-d-sp.diff(e,r)/2)
P=compact(coef(p,p)-sp.diff(coef(p,pr)-sp.diff(d,r),r)/2)
assert coef(ptr,ptr)==0 and coef(prr,prr)==0
assert compact(T-4*q*(q-2)/f)==0
assert compact(R+4*q*(q-2)*f)==0
V_derived=compact(-P*f/(4*q*(q-2)))
res=0
lam=(q-2)/2
V_standard=2*f*(lam**2*(lam+1)*r**3+3*lam**2*M*r*r+9*lam*M*M*r+9*M**3)/(r**3*(lam*r+3*M)**2)
assert compact(V_derived-V_standard)==0
result={'input':'Continuum EH quadratic alpha-gauge action, Solomon 2023 Eq 150',
        'q_definition':'j*(j+1), j>=2', 'auxiliary_a_solution':str(compact(a_sol)),
        'constraint_h_solution':str(h_sol),'derived_Zerilli_potential':str(V_derived),
        'reduced_time_kinetic_coefficient':str(T),'reduced_radial_kinetic_coefficient':str(R),
        'master_equation_residual':str(res),'standard_potential_difference':'0',
        'scope':'Independent symbolic elimination from known EH quadratic density. Does not prove nonlinear fuzzy-action equivalence.'}
OUT.mkdir(parents=True,exist_ok=True)
(OUT/'polar_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
