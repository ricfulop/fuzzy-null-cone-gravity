#!/usr/bin/env python3
"""Reproduce the first nonlinear-consistency checkpoint on CPU."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUT=ROOT/'reports/nonlinear-consistency-2026-11-27'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    runs=[]
    for name in ['projected_sphere_algebra','deformed_angular_candidate',
                 'boundary_compression','reduced_boundary_gluing','euclidean_signature']:
        start=time.monotonic()
        with (OUT/(name+'_run.log')).open('w') as log:
            subprocess.run([sys.executable,str(HERE/(name+'.py'))],cwd=ROOT,
                           env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),
                           stdout=log,stderr=subprocess.STDOUT,check=True)
        runs.append({'script':name+'.py','status':'passed','seconds':round(time.monotonic()-start,3)})
        print(name+': passed',flush=True)
    sources=list(HERE.glob('*.py'))+[ROOT/'verification/modal_fsg_hmc.py',
        ROOT/'verification/repair_audit_20261127/nonlinear_galerkin.py',
        ROOT/'sections/nonlinear_action.tex',ROOT/'paper.tex']
    result={'python':platform.python_version(),
            'dependencies':{name:importlib.metadata.version(name) for name in ['numpy','scipy','sympy']},
            'repair_base_commit':'ab303dd3f756c7f3f8d393ef090d4da35a7fdae2',
            'runs':runs,'GPU_jobs':0,
            'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            'goal_status':'First obstruction/candidate checkpoint complete; nonlinear gravitational completion remains open.'}
    (OUT/'reproducibility.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
