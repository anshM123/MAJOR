"""
n4_prod (NEARGAP100_SMALLNU.md): production list of the x = 0.05 graviton points (state + Delta_ref = 2e4 pair, A+ + A- thimbles,
common random numbers) and a 2-worker driver.
  even (C2, C3H; l9 Delta = 100 row):            J = 4, 6, 8, 10 at nu = 0.05 m(100, J)
  odd  (P8.build C2, C3H): Delta = 100 row      J = 5, 7, 9 at nu = 0.05 m(100, J)
                           extra m-rows         (p8_k2 keys: m = 1.05 (J = 5 only), 1.1, 1.2, 1.3 Lam) at their x = 0.05 node
Usage: OMP_NUM_THREADS=1 python3 n4_prod.py list | WORKER NWORKERS     (resumable: skips pairs already in cache/n4_raw.pkl)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def jobs():
    L = []
    for J in (4, 6, 8, 10):
        L.append(('even', J, 0.05*mphys(100.0, J), 100.0))
    for J in (5, 7, 9):
        L.append(('odd', J, 0.05*mphys(100.0, J), 100.0))
    for J in (5, 7, 9):
        R = pickle.load(open(os.path.join(HERE, 'cache', 'p8_k2_C2_J%d.pkl' % J), 'rb'))
        Ds = sorted(set(k[0] for k in R if k[2] == 1.0 and 100.0 < k[0] < 150.0))
        for D in Ds:
            nus = sorted(k[1] for k in R if k[0] == D and k[2] == 1.0)
            nu = [v for v in nus if abs(v/mphys(D, J) - 0.05) < 1e-6][0]
            L.append(('odd', J, float(nu), float(D)))
    return L


def done(job):
    f = os.path.join(HERE, 'cache', 'n4_raw.pkl')
    if not os.path.exists(f):
        return False
    db = pickle.load(open(f, 'rb'))
    mode, J, nu, D = job
    return all((mode, round(d, 9), J, round(nu, 9)) in db for d in (D, 20000.0))


if __name__ == '__main__':
    if sys.argv[1] == 'list':
        for j in jobs():
            print(j, 'done' if done(j) else '')
        sys.exit()
    w, nw = int(sys.argv[1]), int(sys.argv[2])
    import n4_run as NRun
    L = jobs()
    order = sorted(range(len(L)), key=lambda i: (L[i][0] == 'even', L[i][3] == 100.0, -L[i][1], abs(L[i][3] - 120.08) > 1)) \
        if os.environ.get('N4_ORDER', 'xrows') == 'xrows' else range(len(L))
    mine = [L[i] for k, i in enumerate(order) if k % nw == w]
    for job in mine:
        if done(job):
            continue
        mode, J, nu, D = job
        odd = mode == 'odd'
        kinds = ('odd:C2', 'odd:C3H') if odd else ('C2', 'C3H')
        t = time.time()
        res = NRun.run_pair(J, nu, [D, NRun.DREF], kinds, more=[j[2] for j in mine if j[1] == J and j[0] == mode])
        for D_, raw in res.items():
            NRun._save(os.path.join(HERE, 'cache', 'n4_raw.pkl'), (mode, round(D_, 9), J, round(nu, 9)), raw)
            for k in kinds:
                NRun._save(NRun.fn(k.replace('odd:', ''), J, odd), (D_, round(nu, 9)), NRun.combine_kind(raw, k))
        print(f'## job {job} done in {time.time() - t:.0f}s', flush=True)
