"""
y5 (EE_DG100.md): x = 0.05 node at Delta_g = 100 (nu ~ 4.9-6.5, below the fold) for the EE integrand, as the sum of the A+ and A-
Lefschetz thimbles (n4_run / n4_thimble, imported read-only; NEARGAP100_SMALLNU.md recipe): state + Delta_ref = 2e4 at the same nu on
common random numbers.  Kinds:
  even (Delta = 100 row, J = 4, 6, 8, 10):  'EE' (w0 integrand on the n1_pw point-wise backend), 'C2' (q^2, for comparison), 'C3H' (o_H q)
  odd  (Delta = 100 row J = 5, 7, 9; extra rows m = 1.05 ... 1.3 Lam, p8_k2 keys):  'odd:EE' (w0 make_build_odd_EE on p8's P8 module)
Only the graviton builders are patched in-process (n1_pw.P9.build_spin, p8_odd_pw.P8.build); nothing else is touched.
Files: cache/y5_raw.pkl (raw thimble records, key (mode, Delta, J, nu)), cache/y5_pts_<kind>_J<J>.pkl / cache/y5_odd_pts_<kind>_J<J>.pkl
(n1_run / p8 format, key (Delta, round(nu, 9))), claims in cache/y5_claims/ (several workers may run the list concurrently).
Usage: OMP_NUM_THREADS=1 python3 y5_n4.py list | run [even|odd|all]
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def patch():
    import n1_pw as PW
    import p8_odd_pw as P8M
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(PW, PW.S, PW.Cr, PW.A6)
    orig = PW.P9.build_spin
    if not getattr(PW.P9, '_y5', False):
        PW.P9.build_spin = lambda zl, nu, Delta, J, mons, sol, solt, kind: (bEE(zl, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                           orig(zl, nu, Delta, J, mons, sol, solt, kind))
        PW.P9._y5 = True
        P8M.P8.build = W0.make_build_odd_EE(P8M.P8)


def jobs():
    L = []
    for J in (4, 6, 8, 10):
        L.append(('even', J, 0.05*mphys(100.0, J), 100.0))
    for J in (5, 7, 9):
        L.append(('odd', J, 0.05*mphys(100.0, J), 100.0))
    for J in (5, 7, 9):
        R = pickle.load(open(os.path.join(C, 'p8_k2_C3H_J%d.pkl' % J), 'rb'))
        Ds = sorted(set(k[0] for k in R if k[2] == 1.0 and 100.0 < k[0] < 150.0))
        for D in Ds:
            nus = sorted(k[1] for k in R if k[0] == D and k[2] == 1.0)
            nu = [v for v in nus if abs(v/mphys(D, J) - 0.05) < 1e-6][0]
            L.append(('odd', J, float(nu), float(D)))
    # priority: binding / previously failing states first (odd J = 9 extra rows, then the Delta = 100 rows)
    pri = {('even', 4): 0, ('odd', 5): 1, ('odd', 9): 2, ('even', 6): 3, ('odd', 7): 4, ('even', 8): 5, ('even', 10): 6}
    L.sort(key=lambda j: (j[3] != 100.0 and not (j[0] == 'odd' and j[1] == 9), pri.get((j[0], j[1]), 9), j[3]))
    return L


def kinds_of(mode):
    return ('EE', 'C2', 'C3H') if mode == 'even' else ('odd:EE',)


def fn(kind, J, odd):
    return os.path.join(C, 'y5_%spts_%s_J%d.pkl' % ('odd_' if odd else '', kind, J))


def done(job):
    f = os.path.join(C, 'y5_raw.pkl')
    if not os.path.exists(f):
        return False
    db = pickle.load(open(f, 'rb'))
    mode, J, nu, D = job
    return all((mode, round(d, 9), J, round(nu, 9)) in db for d in (D, DREF))


def claim(job):
    d = os.path.join(C, 'y5_claims'); os.makedirs(d, exist_ok=True)
    f = os.path.join(d, '%s_%d_%.6f_%.4f' % job)
    try:
        fd = os.open(f, os.O_CREAT | os.O_EXCL | os.O_WRONLY); os.write(fd, str(os.getpid()).encode()); os.close(fd)
        return True
    except FileExistsError:
        return False


if __name__ == '__main__':
    if sys.argv[1] == 'list':
        for j in jobs():
            print(j, 'done' if done(j) else '')
        sys.exit()
    sel = sys.argv[2] if len(sys.argv) > 2 else 'all'
    patch()
    import n4_run as NRun
    L = jobs()
    for job in L:
        mode, J, nu, D = job
        if sel != 'all' and mode != sel:
            continue
        if done(job) or not claim(job):
            continue
        kinds = kinds_of(mode)
        t = time.time()
        res = None
        for nuc in (10.0, 9.0, 8.0):
            if nuc != NRun.NUC:
                NRun._TR.clear()
            NRun.NUC = nuc
            try:
                res = NRun.run_pair(J, nu, [D, DREF], kinds, more=[j[2] for j in L if j[1] == J and j[0] == mode])
                break
            except RuntimeError as ex:
                print(f"# job {job}: nu_c = {nuc}: {ex}; retrying lower", flush=True)
        if res is None:
            print(f"## job {job} FAILED", flush=True); continue
        for D_, raw in res.items():
            raw['nuc_used'] = NRun.NUC
            NRun._save(os.path.join(C, 'y5_raw.pkl'), (mode, round(D_, 9), J, round(nu, 9)), raw)
            for k in kinds:
                NRun._save(fn(k.replace('odd:', ''), J, mode == 'odd'), (D_, round(nu, 9)), NRun.combine_kind(raw, k))
        print(f'## job {job} done in {time.time() - t:.0f}s (nu_c {NRun.NUC})', flush=True)
