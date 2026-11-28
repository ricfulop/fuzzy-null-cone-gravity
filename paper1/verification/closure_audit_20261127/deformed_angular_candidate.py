#!/usr/bin/env python3
"""An exact matrix Lie bracket, and its departure from the classical bracket.

This is a candidate only for the area-preserving angular subalgebra. It is
not the ADM hypersurface-deformation algebra or a symmetry proof for I_G,N.
"""
from itertools import permutations
from pathlib import Path
import json
import math
import sys
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'verification'))
from p2_frontier1_sun_fuzzy_action_and_saddles import su2_generators
from projected_sphere_algebra import poisson,x,y,z,r2


def sym_quantize(p,coords):
    N=coords[0].shape[0];out=np.zeros((N,N),complex)
    for powers,coeff in s.Poly(s.expand(p),x,y,z).terms():
        indices=[i for i,n in enumerate(powers) for _ in range(n)]
        orders=list(set(permutations(indices)))
        term=np.zeros((N,N),complex)
        for order in orders:
            product=np.eye(N,dtype=complex)
            for i in order:product=product@coords[i]
            term+=product
        out+=float(coeff)*term/len(orders)
    return out


def main():
    f,g,h=x*z,y*z,x*(5*z*z-r2)
    rows=[]
    for N in [4,8,16,32]:
        j=(N-1)/2;scale=math.sqrt(j*(j+1))
        coords=[J/scale for J in su2_generators(N-1)[:3]]
        def bracket(a,b):return -1j*scale*(a@b-b@a)
        F,G,H=[sym_quantize(v,coords) for v in [f,g,h]]
        jacobi=bracket(bracket(F,G),H)+bracket(bracket(G,H),F)+bracket(bracket(H,F),G)
        error=float(np.linalg.norm(jacobi,2))
        assert error<1e-10
        target=sym_quantize(poisson(f,g),coords)
        defect=float(np.linalg.norm(bracket(F,G)-target,2))
        assert defect>1e-7
        rows.append({'N':N,'matrix_Jacobi_residual':error,
                     'classical_bracket_defect_operator_norm':defect,
                     'Casimir_scaled_defect':defect*j*(j+1)})
    assert all(a['classical_bracket_defect_operator_norm']>b['classical_bracket_defect_operator_norm']
               for a,b in zip(rows,rows[1:]))
    result={'checks':rows,
            'algebraic_fact':'The scaled commutator closes and obeys Jacobi because matrix multiplication is associative.',
            'scope':'An explicit deformation of one angular subalgebra, with a low-polynomial comparison. It does not close the full gravitational constraints or establish invariance of the reconstructed bulk action.'}
    out=ROOT/'reports/nonlinear-consistency-2026-11-27';out.mkdir(parents=True,exist_ok=True)
    (out/'deformed_angular_candidate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
