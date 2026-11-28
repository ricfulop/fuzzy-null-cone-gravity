#!/usr/bin/env python3
"""Check the conformal-sign obstruction to a direct real action equivalence.

The HMC scalar statements are extracted unchanged from the archived driver
and evaluated with a small NumPy tensor adapter. No Modal/GPU job or sampler
is invoked. The probe has Z=0 and addresses direct off-shell identification;
it does not rule out a derived effective theory or complex-contour relation.
"""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json
import math
import sys
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/nonlinear-consistency-2026-11-27'


class Tensor(np.ndarray):
    def square(self):return np.square(self)
    def sum(self,dim=None,axis=None,**kwargs):
        return np.asarray(self).sum(axis=dim if axis is None else axis,**kwargs)


def tensor(x):return np.asarray(x,dtype=float).view(Tensor)


class TorchAdapter:
    Tensor=Tensor
    @staticmethod
    def roll(x,shifts,dims):return tensor(np.roll(np.asarray(x),shifts,axis=dims))
    @staticmethod
    def matmul(a,b):return tensor(np.asarray(a)@np.asarray(b))
    @staticmethod
    def zeros_like(a):return tensor(np.zeros_like(np.asarray(a)))


def scalar_action_from_archived_driver(N,mu2,lam4):
    source=(ROOT/'verification/modal_fsg_hmc.py').read_text()
    tree=ast.parse(source)
    engine=next(n for n in tree.body if isinstance(n,ast.FunctionDef)
                and n.name=='_create_gpu_fsg_engine')
    action=next(n for n in engine.body if isinstance(n,ast.FunctionDef)
                and n.name=='action_components')
    prefix=[]
    for node in action.body:
        prefix.append(node)
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='S_phi'
                                               for t in node.targets):break
    assert isinstance(prefix[-1],ast.Assign)
    action.body=prefix+[ast.Return(value=ast.Name(id='S_phi',ctx=ast.Load()))]
    module=ast.fix_missing_locations(ast.Module(body=[action],type_ignores=[]))
    j=(N-1)/2;m=j-np.arange(N)
    c=np.sqrt((j-m[1:])*(j+m[1:]+1))
    ns={'torch':TorchAdapter,'N':N,'mu2':mu2,'lam4':lam4,
        'diag_L2_phi':tensor(2*j*(j+1)-2*np.outer(m,m)),
        'off_L2_phi':tensor(np.outer(c,c)),
        'unpack_phi':lambda phi:(tensor(phi.real),tensor(phi.imag))}
    exec(compile(module,'archived_hmc_scalar_prefix','exec'),ns)
    # Reuse the archived diagonal harmonic basis without importing Modal.
    basis_def=next(n for n in tree.body if isinstance(n,ast.FunctionDef)
                   and n.name=='_build_scalar_diagonal_basis_np')
    basis_ns={'math':math}
    exec(compile(ast.Module(body=[basis_def],type_ignores=[]),'archived_scalar_basis','exec'),basis_ns)
    basis=basis_ns['_build_scalar_diagonal_basis_np'](N)
    T=np.diag(math.sqrt(N)*basis[(2,0)])
    return ns['action_components'],T,hashlib.sha256(source.encode()).hexdigest()


def main():
    # Independent curvature routine from the preceding action audit.
    path=ROOT/'verification/repair_audit_20261127/nonlinear_galerkin.py'
    spec=importlib.util.spec_from_file_location('galerkin_curvature',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    u,x,y,z=s.symbols('u x y z',real=True);A=s.Function('A')(u)
    R=mod.ricci_scalar(s.diag(A,A,A,A),(u,x,y,z))
    density=s.simplify(A*A*R)
    expected=-3*s.diff(A,u,2)+s.Rational(3,2)*s.diff(A,u)**2/A
    assert s.simplify(density-expected)==0
    eps=s.symbols('epsilon',real=True)
    conf=1+eps*s.cos(u)
    bulk=s.Rational(3,2)*s.diff(conf,u)**2/conf
    integrated_quad=s.integrate(s.expand(s.series(bulk,eps,0,3).removeO()),(u,0,2*s.pi))
    assert integrated_quad==3*s.pi*eps**2/2
    # Standard Euclidean convention S_E=-(16*pi*G)^(-1) int sqrt(g) R.
    # Per unit transverse coordinate volume and G=1:
    euclidean_quadratic=s.simplify(-integrated_quad/(16*s.pi*eps**2))
    assert euclidean_quadratic==-s.Rational(3,32)
    rows=[]
    for N in [3,4,6,8]:
        action,T,source_hash=scalar_action_from_archived_driver(N,0.25,0.1)
        assert abs(np.trace(T))<1e-12
        assert abs(np.trace(T@T)-N)<1e-11
        values=[]
        for e in [0.001,0.01,0.1]:
            phi=(e*T)[None,None,None,:,:].astype(complex)
            measured=float(action(phi,None,None)[0])
            exact=0.5*(6+0.25)*N*e*e+0.1/(4*N)*np.trace(T@T@T@T)*e**4
            assert abs(measured-exact)<1e-12
            values.append({'epsilon':e,'action':measured,'analytic_residual':float(abs(measured-exact))})
        rows.append({'N':N,'HMC_quadratic_coefficient':0.5*(6+0.25)*N,'evaluations':values})
    result={'direct_curvature_residual':'0','sqrt_g_R_density':str(density),
            'Euclidean_conformal_quadratic_per_unit_transverse_volume_G1':str(euclidean_quadratic),
            'HMC_scalar_checks':rows,'archived_HMC_sha256':source_hash,
            'execution':'Original scalar action statements and harmonic builder evaluated on CPU with a NumPy tensor adapter; no sampler run.',
            'conclusion':'The standard real Euclidean Einstein action has a negative conformal direction, while the archived HMC scalar quadratic form is positive. A direct real identification retaining this direction is impossible.',
            'scope':'This rules out the stated direct unreduced real action identification, not every possible gauge-reduced effective theory, analytic continuation, or continuum-universality relation. Those require separate derivations.'}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'euclidean_signature.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
