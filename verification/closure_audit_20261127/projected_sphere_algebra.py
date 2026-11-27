#!/usr/bin/env python3
"""Exact audit of a candidate truncated angular gauge bracket.

This tests P_L{f,g} on scalar Hamiltonians for area-preserving sphere
diffeomorphisms. It is not the full canonical constraint calculation for the
Galerkin metric action. A failure rules out this particular naive bracket.
"""
from functools import lru_cache
from pathlib import Path
import json
import sympy as s

x, y, z = s.symbols('x y z', real=True)
xyz = (x, y, z)
r2 = x*x+y*y+z*z
OUT = Path(__file__).resolve().parents[2]/'reports/nonlinear-consistency-2026-11-27'


def lap(p):
    return s.expand(sum(s.diff(p, v, 2) for v in xyz))


def poisson(f, g):
    df = s.Matrix([s.diff(f, v) for v in xyz])
    dg = s.Matrix([s.diff(g, v) for v in xyz])
    return s.expand(s.Matrix(xyz).dot(df.cross(dg)))


@lru_cache(None)
def homogeneous_parts(p):
    parts = {}
    for powers, coeff in s.Poly(p, *xyz).terms():
        degree = sum(powers)
        parts[degree] = parts.get(degree, s.S.Zero) + coeff*s.prod(
            v**n for v, n in zip(xyz, powers))
    return tuple(sorted(parts.items()))


@lru_cache(None)
def fischer(p, degree):
    """P_d=sum_k r^(2k) H_(d-2k), all H harmonic, exactly."""
    if p == 0:
        return ()
    coeff, delta, harmonic = s.S.One, p, p
    for k in range(1, degree//2+1):
        coeff /= -2*k*(2*degree-2*k+1)
        delta = lap(delta)
        harmonic += coeff*r2**k*delta
    harmonic = s.expand(harmonic)
    assert lap(harmonic) == 0
    remainder = s.expand(p-harmonic)
    if remainder == 0:
        return ((degree, harmonic),)
    quotient, residual = s.div(remainder, r2, *xyz)
    assert residual == 0
    return ((degree, harmonic),)+fischer(s.expand(quotient), degree-2)


@lru_cache(None)
def project(p, cutoff):
    # Restriction to the unit sphere sets each power of r^2 to one.
    return s.expand(sum(h for degree, part in homogeneous_parts(s.expand(p))
                        for j, h in fischer(part, degree) if j <= cutoff))


def bracket(f, g, cutoff):
    return project(poisson(f, g), cutoff)


def main():
    assert poisson(x, y) == z
    for degree in range(1, 8):
        assert project(z**degree, degree) != 0
        assert project(project(z**degree, 3), 3) == project(z**degree, 3)
    f, g, h = x*z, y*z, x*(5*z*z-r2)
    continuum = s.expand(poisson(poisson(f,g),h)
                         +poisson(poisson(g,h),f)
                         +poisson(poisson(h,f),g))
    assert continuum == 0
    rows = []
    for L in [3,4,5]:
        j = s.factor(bracket(bracket(f,g,L),h,L)
                     +bracket(bracket(g,h,L),f,L)
                     +bracket(bracket(h,f,L),g,L))
        rows.append({'L':L,'Jacobiator':str(j),
                     'value_at_0_1_0':str(j.subs({x:0,y:1,z:0}))})
    assert rows[0]['value_at_0_1_0'] == '-15/7'
    assert rows[1]['Jacobiator'] == rows[2]['Jacobiator'] == '0'
    example={'f':str(f),'g':str(g),'h':str(h),'L':3,
             'Jacobiator':rows[0]['Jacobiator'],
             'continuum_Jacobiator':'0','L4_Jacobiator':'0',
             'lost_intermediate_harmonic':str(s.factor(project(poisson(g,h),4)-bracket(g,h,3)))}
    result={'counterexample':example,'cutoff_checks_for_this_triple':rows,
            'conclusion':'Naive projected angular Poisson bracket fails Jacobi at L=3.',
            'scope':'Necessary diagnostic for a proposed angular gauge bracket; not a derivation or universal no-go theorem for the full metric constraint algebra.'}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'projected_sphere_algebra.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
