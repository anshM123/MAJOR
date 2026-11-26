#!/usr/bin/env python3
"""Run the certified high-twist scan for LARGE impact parameter, beta in [8.4, 64] (L ~ 7 to 62.6), in parallel (FULLPOS5.md §4).

    cd ads4/code/fullpos3
    python3 run_large_scan.py              # (cores - 2) workers
    python3 run_large_scan.py 20           # explicit worker count
    python3 run_large_scan.py --summary    # table of what is done
    python3 run_large_scan.py --merge      # also write results/ht_large_done.jsonl (all finished cells, one line each)
    python3 run_large_scan.py --beyond 100 [n]   # optional extra cells on [64, 100] (after the main list)

Independent of run_laptop_scan.py: it runs fp3_ht_large.py (not fp3_ht.py), writes cache/large_<ba>_<bb>.jsonl / .log, and reads
only results/ht_large_done.jsonl and cache/large_*.jsonl.  It never touches results/ht_scan_done.jsonl, cache/laptop_* or
cache/ht_scan*.  Both launchers can run at the same time (they share fp3_ht.py read-only).
Each cell is one single-threaded process and its own output file, so the run can be stopped and restarted at any time.
Cells: width 0.1 on [8.4, 64] (556 cells), ascending, so the band below the contour tail (beta < 31.6, L < 30.16) finishes first.
Requirements: python3 with python-flint, numpy.
"""
import os, sys, json, glob, subprocess, multiprocessing as mp
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))


def _fmt(q):
    q = Fr(q)
    return str(q.numerator) if q.denominator == 1 else f'{float(q):.10g}'

CELLS = [(_fmt(Fr(84 + i, 10)), _fmt(Fr(85 + i, 10))) for i in range(556)]          # [8.4, 8.5], ..., [63.9, 64]
# optional: --beyond B adds width-0.1 cells on [64, B] (floor: LRINF's beta >= 64 tail formula, fp3_ht_large.fcert_tail).
# The certificate cannot reach beta = infinity (FULLPOS5 §3): G^r/F_cert ~ 8e-6 beta^2 alone reaches 1 near beta ~ 350.
for _i, _a in enumerate(sys.argv):
    if _a == '--beyond':
        _B = Fr(sys.argv[_i + 1]); _k = 640
        while Fr(_k + 1, 10) <= _B:
            CELLS.append((_fmt(Fr(_k, 10)), _fmt(Fr(_k + 1, 10)))); _k += 1
        del sys.argv[_i:_i + 2]
        break


def key(ba, bb):
    return (str(Fr(str(ba))), str(Fr(str(bb))))


def done_cells():
    d = {}
    files = [os.path.join(HERE, 'results', 'ht_large_done.jsonl')] + glob.glob(os.path.join(HERE, 'cache', 'large_*.jsonl'))
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
    out = f'large_{ba}_{bb}.jsonl'
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with open(os.path.join(HERE, 'cache', f'large_{ba}_{bb}.log'), 'w') as log:
        p = subprocess.run([sys.executable, 'fp3_ht_large.py', 'scan', f'{ba},{bb}', out], cwd=HERE, env=env, stdout=log,
                           stderr=subprocess.STDOUT)
    return c, p.returncode


def summary(merge=False):
    d = done_cells()
    n = sum(1 for c in CELLS if key(*c) in d)
    print(f'{n}/{len(CELLS)} cells done')
    worst = {}
    for ba, bb in CELLS:
        r = d.get(key(ba, bb))
        if r:
            e0 = r['rec0']['err']/r['rec0']['Fcert'] if r.get('rec0') else None
            e4 = r['rec4']['err']/r['rec4']['Fcert'] if r.get('rec4') else None
            print(f"  [{ba},{bb}]  a*(b=0) = {r['a0_b0']}, a*(b=4) = {r['a0_b4']}   err/F = {e0}, {e4}")
        else:
            print(f'  [{ba},{bb}]  pending')
    if merge:
        os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
        with open(os.path.join(HERE, 'results', 'ht_large_done.jsonl'), 'w') as f:
            for ba, bb in CELLS:
                r = d.get(key(ba, bb))
                if r:
                    f.write(json.dumps(r) + '\n')
        print('wrote results/ht_large_done.jsonl')


if __name__ == '__main__':
    if '--summary' in sys.argv or '--merge' in sys.argv:
        summary('--merge' in sys.argv); sys.exit(0)
    os.makedirs(os.path.join(HERE, 'cache'), exist_ok=True)
    n = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, (os.cpu_count() or 4) - 2)
    d = done_cells()
    todo = [c for c in CELLS if key(*c) not in d]
    print(f'{len(todo)} cells to run on {n} workers', flush=True)
    with mp.Pool(n) as pool:
        for c, rc in pool.imap_unordered(run_cell, todo):
            print(f'cell {c}: exit {rc}', flush=True)
    summary()
