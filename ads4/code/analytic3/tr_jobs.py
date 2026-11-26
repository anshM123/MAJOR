"""
tr_jobs.py: job lists for the signed-transfer route (ANALYTIC_SLIVER3 §9) and a resumable parallel runner.

  python3 tr_jobs.py make La Lb [--nslab 4] [--out jobs_La_Lb.jsonl]     write the job list for one L-cell
  python3 tr_jobs.py run JOBS.jsonl RESULTS.jsonl [--procs 3]            run pending jobs (skip-done), one checker process per job
  python3 tr_jobs.py summary JOBS.jsonl RESULTS.jsonl                    counts, failures, CPU time, projected remaining time

A job = one check_transfer.js box (one d-slab x one frequency box x its channels).  Kinds:
  res     real nu in [lo, lo + 20.1], channels (3,0), (3,4)       (real-axis residual, dressed F3)
  top2    nu = t + 6.25 i, t-box of width 10, channels (3,0), (3,4) (psi2 top line, Re alpha = 1/4 - 4)
  right2  nu = LAM_R + i y, y-box of width 1/4 (nuAxis imag), channels (3,0), (3,4)
  top3 / right3: the same for psi3 (p = 4, Y = 7.25), channels (4,0), (4,4).
A FAIL/ERROR job is not a counterexample: the assembler then falls back to the termwise bound for the cells it would have covered.
"""
import os, sys, json, math, time, subprocess, argparse
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import transfer_kernel as TK

LAM_R = 998.4998          # right edge at Re nu = LAM_R < Lam = 998.49987...; the real segment [LAM_R, Lam] is a separate tiny piece
NU_RES = (500.0, 1060.0)  # residual band handled by the transfer (termwise elsewhere)
TOP_T0 = 400.0            # top line handled by the transfer for t >= TOP_T0


def slabs(La, Lb, n):
    """d-slabs (d_lo, d_hi) covering [e^{-Lb}, e^{-La}], equal in log d."""
    out = []
    for i in range(n):
        lb = Lb - (Lb - La)*i/n; la = Lb - (Lb - La)*(i + 1)/n
        out.append((math.exp(-lb)*(1 - 1e-12), math.exp(-la)*(1 + 1e-12)))
    return out


def make(La, Lb, nslab=4, pieces=('res', 'top2', 'right2', 'top3', 'right3')):
    L = (La + Lb)/2
    jobs = []
    for si, (dlo, dhi) in enumerate(slabs(La, Lb, nslab)):
        tag = f'L{La}_{Lb}_s{si}'
        if 'res' in pieces:
            lo = NU_RES[0]
            while lo < NU_RES[1]:
                hi = lo + 20.1
                box = TK.make_box(f'{tag}_res_{lo:.0f}', dlo, dhi, lo, hi, TK.m_for(hi, L), channels=((3, 0), (3, 4)))
                jobs.append(dict(id=box['name'], kind='res', slab=si, d=[dlo, dhi], f=[lo, hi], box=box)); lo += 20.0
        for (kind, p, Y) in (('top2', 3, 6.25), ('top3', 4, 7.25)):
            if kind not in pieces: continue
            t = TOP_T0
            while t < 998.5:
                t1 = min(t + 10.0, 998.6)
                box = TK.make_box(f'{tag}_{kind}_{t:.0f}', dlo, dhi, t, t1, TK.m_for(t1, L), K=160, R=16,
                                  channels=((p, 0), (p, 4)), nu_imag=Y)
                jobs.append(dict(id=box['name'], kind=kind, slab=si, d=[dlo, dhi], f=[t, t1], Y=Y, box=box)); t = t1
        for (kind, p, Y) in (('right2', 3, 6.25), ('right3', 4, 7.25)):
            if kind not in pieces: continue
            n = int(math.ceil(Y/0.25 - 1e-9))
            for i in range(n):
                y0, y1 = i*0.25, min((i + 1)*0.25, Y)
                box = TK.make_box(f'{tag}_{kind}_{i}', dlo, dhi, y0, y1, TK.m_for(998.5, L) + 2, K=160, R=16,
                                  channels=((p, 0), (p, 4)), axis='imag')
                # imag axis: nu0 is the fixed real part, the box varies Im nu in [y0, y1]
                box['nuImag'] = TK._dec_up(TK.Fr(y0 + y1)/2, 4); box['nu0'] = f'{LAM_R:.4f}'; box['dnu'] = TK._dec_up(TK.Fr(y1 - y0)/2 + TK.Fr(1, 10**6), 6)
                jobs.append(dict(id=box['name'], kind=kind, slab=si, d=[dlo, dhi], f=[y0, y1], Y=Y, box=box))
    return jobs


def load_results(fn):
    out = {}
    if os.path.exists(fn):
        for line in open(fn):
            line = line.strip()
            if line:
                r = json.loads(line); out[r['id']] = r
    return out


def run_one(job, escal=(0, 6, 12)):
    """run a job, raising the Taylor order m (and R for the larger orders) on FAIL_REMAINDER/ERROR."""
    t0 = time.time()
    for dm in escal:
        rows = None; err = None
        box = dict(job['box']); box['m'] = job['box']['m'] + dm
        if dm >= 12: box['R'] = max(box['R'], 16)
        try:
            rows = TK.run_box(box)
        except Exception as e:
            err = str(e)[-500:]
        st = 'ERROR' if rows is None else ('PASS' if all(r['status'] == 'PASS' for r in rows.values()) else 'FAIL')
        if st == 'PASS':
            break
    res = dict(id=job['id'], status=st, secs=time.time() - t0, err=err, m=box['m'])
    if rows is not None:
        res['P'] = max(float(TK.Fr(r['P_bound'])) for r in rows.values() if r['status'] == 'PASS') if st == 'PASS' else None
        res['Q'] = max(float(TK.Fr(r['Q_bound'])) for r in rows.values() if r['status'] == 'PASS') if st == 'PASS' else None
        res['detail'] = {f'{k[0]},{k[1]}': dict(status=v['status'], P=v.get('P_bound', '')[:20], Q=v.get('Q_bound', '')[:20],
                                                rem=v.get('remainder', '')[:14], secs=v.get('seconds')) for k, v in rows.items()}
    return res


def run(jobs_fn, res_fn, procs=3):
    jobs = [json.loads(l) for l in open(jobs_fn) if l.strip()]
    done = {k: v for k, v in load_results(res_fn).items() if v['status'] == 'PASS'}
    todo = [j for j in jobs if j['id'] not in done]
    print(f'{len(jobs)} jobs, {len(done)} done, {len(todo)} to run with {procs} processes', flush=True)
    with ThreadPoolExecutor(max_workers=procs) as ex, open(res_fn, 'a') as fo:
        for res in ex.map(run_one, todo):
            fo.write(json.dumps(res) + '\n'); fo.flush()
            print(res['id'], res['status'], f"{res['secs']:.0f}s", res.get('P'), res.get('Q'), flush=True)


def summary(jobs_fn, res_fn):
    jobs = [json.loads(l) for l in open(jobs_fn) if l.strip()]
    done = load_results(res_fn)
    st = {}
    secs = 0.0
    for r in done.values():
        st[r['status']] = st.get(r['status'], 0) + 1; secs += r['secs']
    left = len([j for j in jobs if j['id'] not in done])
    avg = secs/max(1, len(done))
    print(json.dumps(dict(jobs=len(jobs), done=len(done), status=st, cpu_hours=secs/3600, avg_s=avg, left=left,
                          projected_cpu_hours_left=left*avg/3600,
                          failed=[k for k, r in done.items() if r['status'] != 'PASS'][:20]), indent=1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd'); ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('--nslab', type=int, default=4); ap.add_argument('--out'); ap.add_argument('--procs', type=int, default=3)
    ap.add_argument('--pieces', default='res,top2,right2,top3,right3')
    A = ap.parse_args()
    if A.cmd == 'make':
        La, Lb = float(A.a), float(A.b)
        jobs = make(La, Lb, A.nslab, tuple(A.pieces.split(',')))
        fn = A.out or os.path.join(HERE, 'cache', f'jobs_{La}_{Lb}.jsonl')
        with open(fn, 'w') as f:
            for j in jobs: f.write(json.dumps(j) + '\n')
        print(len(jobs), 'jobs ->', fn)
    elif A.cmd == 'run':
        run(A.a, A.b, A.procs)
    elif A.cmd == 'summary':
        summary(A.a, A.b)
