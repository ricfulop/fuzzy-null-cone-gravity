#!/usr/bin/env python3
"""Reproduce the current repairs in dependency order; never run GPU jobs."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / 'reports/repair-audit-2026-11-27'
JOBS = ['action_derivation', 'polar_derivation', 'matrix_realization',
        'nonlinear_galerkin', 'clock_derivation', 'verify_proof_bounds',
        'clock_large_budget']


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    runs = []
    for name in JOBS:
        start = time.monotonic()
        with (OUT / (name + '_run.log')).open('w') as log:
            subprocess.run([sys.executable, str(HERE / (name + '.py'))],
                           cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
                           check=True)
        runs.append({'script':name+'.py', 'status':'passed',
                     'seconds':round(time.monotonic()-start, 3)})
        print(name + ': passed', flush=True)
    sources = list(HERE.glob('*.py')) + list((ROOT/'sections').glob('*.tex'))
    sources += [ROOT/'paper.tex']
    manifest = {'python':platform.python_version(), 'platform':platform.platform(),
                'dependencies':{name:importlib.metadata.version(name)
                                for name in ['numpy','scipy','sympy']},
                'integration_base_commit':'49de5fa7956dad447414c95726f53d5a2597136a',
                'source_head_at_execution':subprocess.check_output(
                    ['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
                'runs':runs,
                'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(sources)}}
    (OUT/'reproducibility.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()
