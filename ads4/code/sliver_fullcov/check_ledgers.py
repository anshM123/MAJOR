#!/usr/bin/env python3
"""Fast reproducibility path: recompute analytic bounds and all archived ledger inequalities."""
from pathlib import Path
import subprocess,sys,os
root=Path(__file__).resolve().parent;env=dict(os.environ,OPENBLAS_NUM_THREADS='1')
def run(d,f):subprocess.run([sys.executable,f],cwd=root/d,env=env,check=True)
run('sliver_next_work','verify_new_signs.py')
run('physical_sliver_run','verify_sign_ledger.py')
run('sliver_merge','verify_previous.py')
run('sliver_merge','verify_merged.py')
run('sliver_merge','assemble_coverage.py')
