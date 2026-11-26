"""
tr_launch.py: resumable launcher for the signed-transfer boxes of Route T''' on L in [1.8, 3.65] (ANALYTIC_SLIVER3 §9.3).

  python3 tr_launch.py make   [--lmin 1.8] [--lmax 3.65] [--dl 0.15]    write cache/launch/jobs_La_Lb.jsonl for every L-cell
  python3 tr_launch.py run    [--procs 22] [--cells 3.2_3.35,...]        run every pending job (one checker process per job), resumable
  python3 tr_launch.py summary                                           per cell: done / failed / CPU hours / projected remaining
  python3 tr_launch.py estimate                                          core-hour estimate from the m^3 K cost model, before running

Each job = one (d-slab, frequency box, channel pair) box of code/signed_transfer/check_transfer.js; results are appended as JSONL
(cache/launch/results_La_Lb.jsonl, one line per job: id, status, P, Q, m, secs); a killed run restarts where it stopped (skip-PASS);
FAIL/ERROR jobs are retried with m + 6 and m + 12 (R >= 16).  The bound source is pluggable: any producer that writes the same
{id, status, P, Q} rows for these boxes (e.g. the author's Chebyshev-norming checker, ANGULAR_ROUTE §108) can replace the Taylor checker.
Piece selection per cell: residual boxes (nu in [500, 1060]) and the psi2 top line (t >= 400) for every cell; below L = 3.05 also the psi2
right edge and the psi3 top line / right edge (their termwise bounds fail there with the hump).  Assembly: tr_assemble.py La Lb (per cell).
"""
import os, sys, json, math, argparse, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tr_jobs as TJ
import transfer_kernel as TK
OUT = os.path.join(HERE, 'cache', 'launch')
os.makedirs(OUT, exist_ok=True)


def cells(lmin, lmax, dl):
    out = []; L = lmax
    while L > lmin + 1e-9:
        La = round(max(lmin, L - dl), 4); out.append((La, round(L, 4))); L = La
    return out


def nslab_for(La, Lb, width=0.0375):
    return max(1, int(math.ceil((Lb - La)/width - 1e-9)))


def make(lmin, lmax, dl):
    for (La, Lb) in cells(lmin, lmax, dl):
        pieces = ('res', 'top2') if La >= 3.05 else ('res', 'top2', 'right2', 'top3', 'right3')
        jobs = TJ.make(La, Lb, nslab_for(La, Lb), pieces)
        res = [j for j in jobs if j['kind'] == 'res']; oth = sorted([j for j in jobs if j['kind'] != 'res'], key=lambda j: -j['f'][0])
        with open(os.path.join(OUT, f'jobs_{La}_{Lb}.jsonl'), 'w') as f:
            for j in res + oth: f.write(json.dumps(j) + '\n')
        print(f'[{La}, {Lb}]: {len(jobs)} jobs ({nslab_for(La, Lb)} slabs)')


def all_jobs(sel=None):
    out = []
    for fn in sorted(os.listdir(OUT)):
        if fn.startswith('jobs_'):
            tag = fn[5:-6]
            if sel and tag not in sel: continue
            out += [(tag, json.loads(l)) for l in open(os.path.join(OUT, fn)) if l.strip()]
    return out


def run(procs, sel=None):
    jobs = all_jobs(sel)
    done = {}
    for tag in {t for t, _ in jobs}:
        for k, v in TJ.load_results(os.path.join(OUT, f'results_{tag}.jsonl')).items():
            if v['status'] == 'PASS': done[k] = v
    todo = [(t, j) for (t, j) in jobs if j['id'] not in done]
    print(f'{len(jobs)} jobs, {len(done)} PASS, {len(todo)} to run on {procs} processes', flush=True)
    def one(tj):
        t, j = tj
        r = TJ.run_one(j)
        with open(os.path.join(OUT, f'results_{t}.jsonl'), 'a') as fo:
            fo.write(json.dumps(r) + '\n')
        return t, r
    with ThreadPoolExecutor(max_workers=procs) as ex:
        for t, r in ex.map(one, todo):
            print(t, r['id'], r['status'], f"{r['secs']:.0f}s m={r.get('m')}", flush=True)


def summary(sel=None):
    jobs = all_jobs(sel)
    by = {}
    for t, j in jobs: by.setdefault(t, []).append(j)
    for t in sorted(by, key=lambda s: -float(s.split('_')[0])):
        res = TJ.load_results(os.path.join(OUT, f'results_{t}.jsonl'))
        ok = sum(1 for j in by[t] if res.get(j['id'], {}).get('status') == 'PASS')
        fail = sum(1 for j in by[t] if j['id'] in res and res[j['id']]['status'] != 'PASS')
        cpu = sum(r['secs'] for r in res.values())/3600
        avg = cpu*3600/max(1, len(res))
        print(f'[{t}] {len(by[t])} jobs: {ok} PASS, {fail} FAIL/ERROR, {len(by[t]) - ok - fail} pending; CPU {cpu:.2f} h; '
              f'projected {((len(by[t]) - ok - fail)*avg/3600):.2f} h')


COST_REF = (28, 200, 27.0)     # (m, K, seconds per channel) measured for the combo box (signed_transfer README)


def estimate(lmin, lmax, dl):
    tot = 0.0
    for (La, Lb) in cells(lmin, lmax, dl):
        L = (La + Lb)/2
        pieces = ('res', 'top2') if La >= 3.05 else ('res', 'top2', 'right2', 'top3', 'right3')
        jobs = TJ.make(La, Lb, nslab_for(La, Lb), pieces)
        s = 0.0
        for j in jobs:
            m, K = j['box']['m'], j['box']['K']
            s += 2*COST_REF[2]*(m/COST_REF[0])**3*(K/COST_REF[1])*1.3       # 2 channels, 30 % retry/overhead allowance
        tot += s
        print(f'[{La}, {Lb}]: {len(jobs)} jobs, ~{s/3600:.1f} core-h')
    print(f'total ~{tot/3600:.0f} core-h (model m^3 K from the README timing; uncontended; excludes assembly, ~0.5 core-h per cell)')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd'); ap.add_argument('--lmin', type=float, default=1.8); ap.add_argument('--lmax', type=float, default=3.65)
    ap.add_argument('--dl', type=float, default=0.15); ap.add_argument('--procs', type=int, default=22); ap.add_argument('--cells')
    A = ap.parse_args()
    sel = set(A.cells.split(',')) if A.cells else None
    if A.cmd == 'make': make(A.lmin, A.lmax, A.dl)
    elif A.cmd == 'run': run(A.procs, sel)
    elif A.cmd == 'summary': summary(sel)
    elif A.cmd == 'estimate': estimate(A.lmin, A.lmax, A.dl)
