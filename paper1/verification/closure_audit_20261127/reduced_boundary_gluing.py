#!/usr/bin/env python3
"""Exact finite gluing model after choosing a reduced boundary algebra.

The chosen algebra is an Abelian boost character (area generator) times
compact SU(2), with dual orientation on the second region. Boost angle is
a group parameter, not a self-adjoint canonical angle operator. This is a
constructive reduced-model example, not a derivation of the full gravitational
corner algebra or of the postulated area spectrum.
"""
from pathlib import Path
import json
import sympy as s

OUT=Path(__file__).resolve().parents[2]/'reports/nonlinear-consistency-2026-11-27'


def spin(j):
    dim=int(2*j+1);m=[j-a for a in range(dim)]
    plus=s.zeros(dim)
    for a in range(1,dim):plus[a-1,a]=s.sqrt((j-m[a])*(j+m[a]+1))
    return (plus+plus.T)/2,(plus-plus.T)/(2*s.I),s.diag(*m)


def main():
    # A/(8*pi*G*hbar) is dimensionless here; two selected, distinct characters.
    # The half/integer spins make this a generic algebra test, not a claim
    # about which corner multiplets the paper's physical model realizes.
    area=s.diag(1,1,2,2,2)
    Jhalf=spin(s.Rational(1,2));Jone=spin(s.S.One)
    J=[s.diag(Jhalf[i],Jone[i]) for i in range(3)]
    dim=5;eye=s.eye(dim)
    area_constraint=s.kronecker_product(area,eye)-s.kronecker_product(eye,area.T)
    rotation=[s.kronecker_product(j,eye)-s.kronecker_product(eye,j.T) for j in J]
    omega=[]
    for start,n in [(0,2),(2,3)]:
        v=s.zeros(dim*dim,1)
        for i in range(start,start+n):v[i*dim+i]=1/s.sqrt(n)
        omega.append(v)
    V=s.Matrix.hstack(*omega)
    P=V*V.conjugate().T
    assert V.conjugate().T*V==s.eye(2)
    assert P*P==P and P==P.conjugate().T
    assert area_constraint*P==s.zeros(dim*dim)
    for c in rotation:assert s.simplify(c*P)==s.zeros(dim*dim)
    # The intersection of the explicit constraints has exactly this dimension.
    stacked=s.Matrix.vstack(area_constraint,*rotation)
    assert len(stacked.nullspace())==2
    # Orientation matters: the same-sign area sum has no nonzero invariant state.
    wrong=s.kronecker_product(area,eye)+s.kronecker_product(eye,area.T)
    assert len(wrong.nullspace())==0
    # Associative sewing and evolution on the reduced two-dimensional space.
    H=s.Matrix([[1,s.Rational(1,3)],[s.Rational(1,3),2]])
    embedded=V*H*V.conjugate().T
    assert embedded*P==P*embedded==embedded
    # Exact rational Cayley evolution is unitary for this Hermitian H.
    dt=s.Rational(1,5)
    U=(s.eye(2)+s.I*dt*H/2).inv()*(s.eye(2)-s.I*dt*H/2)
    assert s.simplify(U.conjugate().T*U)==s.eye(2)
    K=V*U*V.conjugate().T
    assert s.simplify(K.conjugate().T*K-P)==s.zeros(dim*dim)
    assert s.simplify(K*K-V*(U*U)*V.conjugate().T)==s.zeros(dim*dim)
    result={'boundary_dimension_per_side':5,'tensor_product_dimension':25,
            'reduced_gluing_dimension':2,'projection_rank':P.rank(),
            'Hermiticity_residual':'0','idempotence_residual':'0',
            'area_constraint_residual':'0','rotation_constraint_residuals':['0','0','0'],
            'constraint_kernel_dimension':2,'same_orientation_kernel_dimension':0,
            'unitarity_on_reduced_space_residual':'0','discrete_sewing_residual':'0',
            'construction':'P=sum_n |Omega_n><Omega_n| with dual representation on the second boundary; boost parameter is not an operator.',
            'scope':'Exact for the explicitly selected reduced R x SU(2) boundary model. The full gravitational corner reduction, area spectrum, and coupling to the bulk are not derived.'}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'reduced_boundary_gluing.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
