#!/usr/bin/env python3
"""
ext_launch.py: resumable laptop launcher for the sliver_extension pipeline on L in [~1.5, 2.9] (ANALYTIC_SLIVER3 §11).
Every cell runs ext_cell.py (which calls the author's make_resid / check_norm / resid_budget / ellipse_extend unchanged).

  python3 ext_launch.py plan      [--lmin 1.5] [--lmax 2.9]          list the cells
  python3 ext_launch.py preflight [--lmin ..] [--lmax ..]            node counts per cell (cheap) -> cache/ext/preflight.jsonl
  python3 ext_launch.py estimate  [--rate SEC_PER_NODE_K]            core-hours per cell from the preflight and the measured rate
  python3 ext_launch.py run       [--procs 22] [--cells-parallel 2]  run all pending cells (skip cells whose ledger has a result)
  python3 ext_launch.py summary                                     per cell: status, E/F (b = 0, 4), wall, scalar evaluations

Cells: width 0.15 for L >= 2.3, 0.1 on [1.9, 2.3), 0.05 below (the d-range per cell, hence the norming degree, grows as L falls).
A FAIL in a cell is not a counterexample: split the cell (smaller d-range) or the frequency bands and rerun.
"""
import os, sys, json, math, argparse, subprocess, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
HERE = Path(__file__).resolve().parent
OUT = HERE / 'cache' / 'ext'; OUT.mkdir(parents=True, exist_ok=True)


def cells(lmin, lmax):
    out = []; L = lmax
    while L > lmin + 1e-9:
        w = 0.15 if L > 2.3 + 1e-9 else (0.1 if L > 1.9 + 1e-9 else 0.05)
        La = round(max(lmin, L - w), 4); out.append((La, round(L, 4))); L = La
    return out


def tag(c): return f'{float(c[0])}_{float(c[1])}'


def preflight(lmin, lmax):
    fo = open(OUT / 'preflight.jsonl', 'a')
    for c in cells(lmin, lmax):
        t0 = time.time()
        r = subprocess.run([sys.executable, str(HERE / 'ext_cell.py'), str(c[0]), str(c[1]), '--preflight'], capture_output=True, text=True)
        rows = []
        for l in r.stdout.splitlines():
            p = l.split()
            if len(p) >= 6 and p[0].endswith('.json'):
                rows.append(dict(file=p[0], name=p[1], p=int(p[2]), b=int(p[3]), nodes=int(p[4])))
        # K per box from the generated configs
        od = OUT / tag(c); K = {}
        for fn in od.glob('resid_low_*.json'): K.update({x['name']: x['K'] for x in json.load(open(fn))})
        for fn in od.glob('*coupled*.json'):
            if 'results' not in fn.name: K.update({x['name']: x['K'] for x in json.load(open(fn))})
        nk = sum(x['nodes']*K.get(x['name'], 150) for x in rows)
        rec = dict(cell=tag(c), nodes=sum(x['nodes'] for x in rows), nodesK=nk, boxes=len(rows), secs=time.time() - t0, err=r.stderr[-300:] if r.returncode else '')
        fo.write(json.dumps(rec) + '\n'); fo.flush(); print(json.dumps(rec), flush=True)


def estimate(rate, contour_s):
    rows = {}
    for l in open(OUT / 'preflight.jsonl'):
        r = json.loads(l); rows[r['cell']] = r
    tot = 0.0
    for t, r in sorted(rows.items(), key=lambda z: -float(z[0].split('_')[0])):
        ch = r['nodesK']*rate/3600 + 4*contour_s/3600
        tot += ch
        print(f"[{t}] {r['nodes']:>9d} scalar evaluations, sum nodes*K = {r['nodesK']:.3g}: ~{ch:.1f} core-h")
    print(f'total ~{tot:.0f} core-h')


def run_cell(c, jobs):
    t = tag(c); led = OUT / t / 'ledger.json'
    if led.exists() and 'result' in json.load(open(led)):
        return t, 'done'
    r = subprocess.run([sys.executable, str(HERE / 'ext_cell.py'), str(c[0]), str(c[1]), '--jobs', str(jobs)], capture_output=True, text=True)
    open(OUT / f'{t}.out', 'w').write(r.stdout + r.stderr)
    return t, ('ok' if r.returncode == 0 else 'error')


def run(lmin, lmax, procs, par):
    todo = cells(lmin, lmax)
    jobs = max(1, procs//par)
    with ThreadPoolExecutor(max_workers=par) as ex:
        for t, st in ex.map(lambda c: run_cell(c, jobs), todo):
            print(t, st, flush=True)


def summary(lmin, lmax):
    for c in cells(lmin, lmax):
        led = OUT / tag(c) / 'ledger.json'
        if not led.exists():
            print(f'[{tag(c)}] pending'); continue
        L = json.load(open(led))
        if 'result' not in L:
            print(f"[{tag(c)}] running/incomplete: steps {list(L['steps'])}"); continue
        rs = {r['b']: r['error_ratio_upper'] for r in L['result']}
        st = all(r['status'] == 'PASS_CONDITIONAL' for r in L['result'])
        print(f"[{tag(c)}] {'PASS' if st else 'FAIL'}  E/F b=0 {rs.get(0):.4f}  b=4 {rs.get(4):.4f}  wall {L['wall_total']/3600:.2f} h  "
              f"scalar evaluations {L['scalar_evaluations']}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd')
    ap.add_argument('--lmin', type=float, default=1.5); ap.add_argument('--lmax', type=float, default=2.9)
    ap.add_argument('--procs', type=int, default=22); ap.add_argument('--cells-parallel', type=int, default=2)
    ap.add_argument('--rate', type=float, default=None); ap.add_argument('--contour-s', type=float, default=120.0)
    A = ap.parse_args()
    if A.cmd == 'plan': print(cells(A.lmin, A.lmax))
    elif A.cmd == 'preflight': preflight(A.lmin, A.lmax)
    elif A.cmd == 'estimate': estimate(A.rate, A.contour_s)
    elif A.cmd == 'run': run(A.lmin, A.lmax, A.procs, A.cells_parallel)
    elif A.cmd == 'summary': summary(A.lmin, A.lmax)
