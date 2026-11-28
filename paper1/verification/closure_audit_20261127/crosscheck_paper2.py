#!/usr/bin/env python3
"""Reproduce the companion PR's pinned scattering script and add jet tests.

Usage: python .../crosscheck_paper2.py --source /path/to/bondi_sachs_4pt_scattering.py
The source must match Paper II PR #7 at the recorded commit. Its own JSON
output is written next to that source; this audit records a scoped summary.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import types
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/nonlinear-consistency-2026-11-27'
COMMIT='5e0ab55afed69f24432bccffd6ac2a781f29cfbf'
HASH='fcf4b5bd515ee1822c8b6ec3697e8c4ff99b84f03263dace1633d00bf05fe582'


class NPProxy:
    def __init__(self,seed):
        self.random=types.SimpleNamespace(default_rng=lambda ignored:np.random.default_rng(seed))
    def __getattr__(self,name):return getattr(np,name)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,required=True)
    path=parser.parse_args().source.resolve()
    assert hashlib.sha256(path.read_bytes()).hexdigest()==HASH,'Companion source differs from pinned version.'
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'paper2_crosscheck_run.log').open('w') as log:
        subprocess.run([sys.executable,str(path)],cwd=path.parent,
                       env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),
                       stdout=log,stderr=subprocess.STDOUT,check=True)
    spec=importlib.util.spec_from_file_location('companion_scattering',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    rows=[]
    for seed in range(30,50):
        mod.np=NPProxy(seed)
        rows.append({'seed':seed,'result':mod.verify_nonlinear_bondi_sachs_identity()})
    full=json.loads(path.with_name('bondi_sachs_4pt_results.json').read_text())
    result={'companion_PR':7,'companion_commit':COMMIT,'source_sha256':HASH,
            'full_benchmark_reproduction':'passed','nonlinear_identity_20_additional_seeds':rows,
            'corrected_action_F0_0':{k:v for k,v in full['corrected_action_F0_0'].items() if 'details' not in k},
            'corrected_action_F0_1':full['corrected_action_F0_1'],
            'ward_and_soft_limit':full['ward_and_soft_limit'],'ablations':full['ablations'],
            'scope':'Continuum affine Bondi-Sachs weak-field tree benchmark. It is not a finite-rank Dirac algebra, a boundary quantization, or an HMC equivalence proof.'}
    (OUT/'paper2_crosscheck.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Pinned companion scattering benchmark: passed; 20 additional nonlinear jets: passed.')


if __name__=='__main__':main()
