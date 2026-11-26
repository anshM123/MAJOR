#!/usr/bin/env python3
"""Resumable parallel launcher for the Lemma U cell certificates (lu_cell.py), one process per beta-cell.

    cd ads4/code/lemU_cert
    python3 run_lemU_scan.py --plan                 # write scan_plan.json (cells + parameters) and print the estimate
    python3 run_lemU_scan.py 20                     # run batch 1 (beta in [1.47, 5.1]) on 20 workers
    python3 run_lemU_scan.py 20 --batch 2           # batch 2 (beta in [0.486, 1.47])
    python3 run_lemU_scan.py --summary              # coverage, failures, margins, CPU time

Each cell runs `lu_cell.py u_lo u_hi prec I1 I0 D cache/scan/cell_<lo>_<hi>.jsonl rho_fac precC` single-threaded with an address-space
limit (default 1500 MB, --mem MB), writes its own JSONL (one record per channel b = 0, 4) and log, so the run can be stopped and
restarted at any time: cells with both channel records present are skipped.  A cell whose certificate fails (ok = false) is
reported by --summary; `--split-failed` replaces each failed cell by its two halves in scan_plan.json (rerun afterwards).
Requirements: python3 with python-flint (>= 0.6), numpy.
"""
import os, sys, json, glob, subprocess, resource, time, math
import multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__))
PLAN = os.path.join(HERE, 'scan_plan.json')
OUT = os.path.join(HERE, 'cache', 'scan')

BATCHES = {1: (1.47, 5.1), 2: (0.486, 1.47)}


def params(beta):
    """cell half-width (in u = beta_c) and Taylor-model parameters as a function of beta (see LEMU_CERT.md §5 for the sizing)."""
    L = -math.log(1 - math.tanh(beta/2)**2)
    if L >= 1.5:
        return dict(width=0.004, prec=1500, precC=900, I1=260, I0=60, D=60, rho=6.0)
    if L >= 0.5:
        return dict(width=0.004, prec=2500, precC=1200, I1=260, I0=60, D=60, rho=6.0)
    return dict(width=0.002, prec=4000, precC=2000, I1=260, I0=60, D=60, rho=6.0)


def make_plan():
    cells = []
    for bt, (lo, hi) in BATCHES.items():
        u = lo
        while u < hi - 1e-12:
            p = params(u)
            w = p['width']
            u2 = min(hi, round(u + w, 6))
            cells.append(dict(batch=bt, u_lo=round(u, 6), u_hi=u2, **{k: v for k, v in p.items() if k != 'width'}))
            u = u2
    json.dump(cells, open(PLAN, 'w'), indent=0)
    return cells


def cell_file(c):
    return os.path.join(OUT, f"cell_{c['u_lo']:.6f}_{c['u_hi']:.6f}.jsonl")


def records(c):
    f = cell_file(c)
    out = {}
    if os.path.exists(f):
        for l in open(f):
            try:
                r = json.loads(l)
                out[r['b']] = r
            except Exception:
                pass
    return out


def done(c):
    r = records(c)
    return 0 in r and 4 in r


def run_cell(args):
    c, mem = args
    def lim():
        resource.setrlimit(resource.RLIMIT_AS, (mem*2**20, mem*2**20))
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    f = cell_file(c)
    if os.path.exists(f):
        os.remove(f)                                  # partial (one channel) record: redo the cell
    cmd = [sys.executable, 'lu_cell.py', str(c['u_lo']), str(c['u_hi']), str(c['prec']), str(c['I1']), str(c['I0']), str(c['D']),
           f, str(c['rho']), str(c['precC'])]
    t = time.time()
    with open(f.replace('.jsonl', '.log'), 'w') as log:
        p = subprocess.run(cmd, cwd=HERE, env=env, stdout=log, stderr=subprocess.STDOUT, preexec_fn=lim)
    return c, p.returncode, time.time() - t


def summary(cells):
    tot = len(cells); nd = 0; fails = []; worst = None; cpu = 0.0
    for c in cells:
        r = records(c)
        if 0 in r and 4 in r:
            nd += 1
            for b in (0, 4):
                cpu += r[b].get('t_chan', 0)
                if not r[b]['ok']:
                    fails.append((c['u_lo'], c['u_hi'], b))
                elif worst is None or r[b]['margin_rel'] < worst[0]:
                    worst = (r[b]['margin_rel'], c['u_lo'], c['u_hi'], b)
            cpu += r[0].get('t_jets', 0) + r[0].get('t_yprime', 0)
    print(f'{nd}/{tot} cells done; {len(fails)} failed channel certificates; CPU in finished cells {cpu/3600:.1f} h')
    for bt, (lo, hi) in BATCHES.items():
        cs = [c for c in cells if c['batch'] == bt]
        ok = [c for c in cs if done(c) and all(records(c)[b]['ok'] for b in (0, 4))]
        cov = sum(c['u_hi'] - c['u_lo'] for c in ok)
        print(f'  batch {bt}: beta in [{lo}, {hi}]: {len(ok)}/{len(cs)} cells certified (both b), covered length {cov:.4f} of {hi - lo:.4f}')
        # contiguity check of the certified part
        gaps = []
        prev = lo
        for c in sorted(cs, key=lambda c: c['u_lo']):
            if c['u_lo'] > prev + 1e-12:
                gaps.append((prev, c['u_lo']))
            prev = max(prev, c['u_hi'])
        if gaps:
            print(f'    plan gaps: {gaps[:5]}')
    if worst:
        print(f'  smallest certified margin (lower/Phat00): {worst[0]:.4f} at [{worst[1]}, {worst[2]}], b = {worst[3]}')
    for f in fails[:20]:
        print(f'  FAILED: [{f[0]}, {f[1]}] b = {f[2]}')


if __name__ == '__main__':
    a = sys.argv[1:]
    os.makedirs(OUT, exist_ok=True)
    if '--plan' in a or not os.path.exists(PLAN):
        cells = make_plan()
        print(f'plan: {len(cells)} cells -> {PLAN}')
        if '--plan' in a:
            sys.exit(0)
    cells = json.load(open(PLAN))
    if '--summary' in a:
        summary(cells); sys.exit(0)
    if '--split-failed' in a:
        new = []
        for c in cells:
            r = records(c)
            if 0 in r and 4 in r and not (r[0]['ok'] and r[4]['ok']):
                m = round((c['u_lo'] + c['u_hi'])/2, 6)
                new += [dict(c, u_hi=m), dict(c, u_lo=m)]
            else:
                new.append(c)
        json.dump(new, open(PLAN, 'w'), indent=0); print(f'{len(new) - len(cells)} cells split'); sys.exit(0)
    batch = 1
    if '--batch' in a:
        batch = int(a[a.index('--batch') + 1])
    mem = int(a[a.index('--mem') + 1]) if '--mem' in a else 1500
    nums = [x for x in a if x.isdigit() and (a.index(x) == 0 or a[a.index(x) - 1] not in ('--batch', '--mem'))]
    n = int(nums[0]) if nums else max(1, (os.cpu_count() or 4) - 2)
    todo = [c for c in cells if c['batch'] == batch and not done(c)]
    print(f'{len(todo)} cells to run (batch {batch}) on {n} workers', flush=True)
    with mp.Pool(n) as pool:
        for c, rc, dt in pool.imap_unordered(run_cell, [(c, mem) for c in todo]):
            r = records(c)
            st = ' '.join(f"b{b}:{'ok' if r[b]['ok'] else 'FAIL'}" for b in sorted(r)) if r else 'no record'
            print(f"cell [{c['u_lo']}, {c['u_hi']}]: exit {rc}, {dt/60:.1f} min, {st}", flush=True)
    summary(cells)
