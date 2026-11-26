#!/usr/bin/env python3
"""Run the remaining certified a*(L) cells of the FULLPOS3 high-twist scan in parallel (e.g. on a many-core laptop).

    cd ads4/code/fullpos3
    python3 run_laptop_scan.py            # uses (cores - 2) workers
    python3 run_laptop_scan.py 20         # or an explicit worker count

Skips every cell already in results/ht_scan_done.jsonl or in cache/laptop_*.jsonl. Each cell runs as its own process
(`fp3_ht.py scan "ba,bb" laptop_<ba>_<bb>.jsonl`), single-threaded, and writes its own file, so the run can be stopped and
restarted at any time. Requirements: python3 with python-flint, numpy, mpmath, sympy.
At the end: python3 run_laptop_scan.py --summary
"""
import os, sys, json, glob, subprocess, multiprocessing as mp
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = [(0.39, 0.44), (0.44, 0.49), (0.49, 0.54), (0.54, 0.59), (0.59, 0.64), (0.64, 0.69), (0.69, 0.74), (0.74, 0.79), (0.79, 0.84),
         (0.84, 0.89), (0.89, 0.94), (0.94, 0.99), (0.99, 1.04), (1.04, 1.14), (1.14, 1.24), (1.24, 1.34), (1.34, 1.44), (1.44, 1.54),
         (1.54, 1.64), (1.64, 1.74), (1.74, 1.84), (1.84, 1.94), (1.94, 2.04), (2.04, 2.14), (2.14, 2.24), (2.24, 2.34), (2.34, 2.44),
         (2.44, 2.54), (2.54, 2.64), (2.64, 2.74), (2.74, 2.84), (2.84, 2.94), (2.94, 3.04), (3.04, 3.24), (3.24, 3.44), (3.44, 3.64),
         (3.64, 3.84), (3.84, 4.04), (4.04, 4.24), (4.24, 4.44), (4.44, 4.64), (4.64, 4.84), (4.84, 5.04), (5.04, 5.24), (5.24, 5.44),
         (5.44, 5.64), (5.64, 5.84), (5.84, 6.04), (6.04, 6.24), (6.24, 6.44), (6.44, 6.64), (6.64, 6.84), (6.84, 7.04), (7.04, 7.24),
         (7.24, 7.44), (7.44, 7.64), (7.64, 7.84), (7.84, 8.04), (8.04, 8.24), (8.24, 8.4)]


def key(ba, bb):
    return (str(Fr(str(ba))), str(Fr(str(bb))))


def done_cells():
    d = {}
    files = [os.path.join(HERE, 'results', 'ht_scan_done.jsonl')] + glob.glob(os.path.join(HERE, 'cache', 'laptop_*.jsonl')) \
        + glob.glob(os.path.join(HERE, 'cache', 'ht_scan*.jsonl'))
    for f in files:
        if os.path.exists(f):
            for l in open(f):
                try:
                    r = json.loads(l)
                except Exception:
                    continue
                d[key(r['ba'], r['bb'])] = r
    return d


def run_cell(c):
    ba, bb = c
    out = f'laptop_{ba}_{bb}.jsonl'
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with open(os.path.join(HERE, 'cache', f'laptop_{ba}_{bb}.log'), 'w') as log:
        p = subprocess.run([sys.executable, 'fp3_ht.py', 'scan', f'{ba},{bb}', out], cwd=HERE, env=env, stdout=log, stderr=subprocess.STDOUT)
    return c, p.returncode


def summary():
    d = done_cells()
    print(f'{len(d)}/{len(CELLS)} cells done')
    for ba, bb in CELLS:
        r = d.get(key(ba, bb))
        print(f'  [{ba},{bb}]  ' + (f"a*(b=0) = {r['a0_b0']}, a*(b=4) = {r['a0_b4']}" if r else 'pending'))


if __name__ == '__main__':
    if '--summary' in sys.argv:
        summary(); sys.exit(0)
    os.makedirs(os.path.join(HERE, 'cache'), exist_ok=True)
    n = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, (os.cpu_count() or 4) - 2)
    d = done_cells()
    todo = [c for c in CELLS if key(*c) not in d]
    print(f'{len(todo)} cells to run on {n} workers', flush=True)
    with mp.Pool(n) as pool:
        for c, rc in pool.imap_unordered(run_cell, todo):
            print(f'cell {c}: exit {rc}', flush=True)
    summary()
