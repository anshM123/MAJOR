"""
q4 (S1 validation, S1_S4_FIXES.md §2.2): NON-PERTURBATIVE (direct RQMC, n1_direct contour) graviton crystal values for C^_2 / C^_3 at
large J on the Delta = 400 row (the row that carries the Delta_g = 300 near-gap states), compared with the K2 Laplace chain used by q3.
n1_run.point is used unchanged (read-only import); only its cache file is redirected to cache/q4_n1_<kind>_J<J>.pkl, and n1_analysis's
shadow table is extended by the q2 shadow caches.
Usage:  OMP_NUM_THREADS=1 python3 q4_direct.py run kind J Delta x1,x2,...      (state at Delta and reference at 2e4, same nu)
        python3 q4_direct.py analyse J [J ...]                                   -> cache/q4_direct_compare.txt
"""
import os, sys, glob, pickle, time
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
CACHE = os.path.join(HERE, 'cache')


def qfn(kind, J):
    return os.path.join(CACHE, 'q4_n1_%s_J%d.pkl' % (kind, J))


def qload(kind, J):
    f = qfn(kind, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def run(kind, J, Delta, xs):
    import n1_run as R
    m = R.mphys(Delta, J)
    for D in (Delta, 20000.0):
        for x in xs:
            nu = x*m
            key = (D, round(nu, 9))
            if key in qload(kind, J):
                continue
            t = time.time()
            r = R.point(kind, J, D, nu)
            db = qload(kind, J); db[key] = r
            pickle.dump(db, open(qfn(kind, J) + '.tmp', 'wb')); os.replace(qfn(kind, J) + '.tmp', qfn(kind, J))
            rel = np.max(np.abs(r['err']))/np.max(np.abs(r['E']))
            print(f"{kind} J={J} Delta={D:.0f} x={x} nu={nu:.3f}: max rel err {rel:.1e}, nmask {r['nmask']}, {time.time() - t:.0f}s", flush=True)


def analyse(Js, D0=400.0):
    import q3_analysis as Q3                     # applies the q3 loaders (q2 caches, phase strip) and de-shifts (C^_2', C^_3 none)
    import n1_analysis as N1
    for fn in glob.glob(os.path.join(CACHE, 'q2_shadow_k2*.pkl')):
        N1.SHK2.update(pickle.load(open(fn, 'rb')))
    N1.dload = qload
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# q4: direct (RQMC) vs K2 Laplace on the Delta = {D0:.0f} row, l9 units (Gamma = c2 units x LO shape); C^_3 shift: {Q3.C3SHIFT}")
    store = {}
    for J in Js:
        N1.DIRECT = True
        gd = N1.grid_direct(J, 2, ('C2', 'C3H'), Ds=[D0])
        N1.DIRECT = False
        g2 = N1.grid_direct(J, 2, ('C2', 'C3H'), Ds=[D0])
        N1.DIRECT = True
        if D0 not in gd:
            pr(f"  J={J}: no row"); continue
        r, r2 = gd[D0], g2[D0]
        store[J] = dict(direct=r, k2=r2)
        pr(f"  J={J} (m={r['m']:.1f}, J^2/m={J*J/r['m']:.2f}): direct points used {r['ndirect']}; a_h direct/K2: " +
           ", ".join(f"{h} {r['a'][h]/r2['a'][h]:.4f}" for h in r['a']))
        pr("    x       : " + " ".join(f"{x:7.2f}" for x in r['xs']))
        m = r['m']
        for kd in ('C2', 'C3H'):
            sg = {'++': 1.0, '+-': -1.0} if kd == 'C3H' else {'++': 1.0, '+-': 1.0}
            for e in ('++', '+-'):
                # direct re-expressed with the K2 per-state constant a_h (removes the congruence a_dir/a_K2):
                Dh = Q3.A.Dsh(e, J, r['xs'])
                ra = m*m*((r[(kd, e)]/m**2 + sg[e]*Dh)*r['a'][e]/r2['a'][e] - sg[e]*Dh)
                pr(f"    {kd:3s} {e:3s} direct (a_K2): " + " ".join(f"{v:7.0f}" for v in ra))
                pr(f"    {kd:3s} {e:3s} K2           : " + " ".join(f"{v:7.0f}" for v in r2[(kd, e)]))
                pr(f"    {kd:3s} {e:3s} RQMC err     : " + " ".join(f"{v:7.0f}" for v in r[(kd, e, 'err')]))
            pr(f"    {kd:3s} off direct       : " + " ".join(f"{v:7.2f}" for v in r[(kd, 'off')]*np.sqrt(r['a']['++']*r['a']['+-']/(r2['a']['++']*r2['a']['+-']))))
            pr(f"    {kd:3s} off K2           : " + " ".join(f"{v:7.2f}" for v in r2[(kd, 'off')]))
    if D0 == 400.0:
        pickle.dump(store, open(os.path.join(CACHE, 'q4_direct_grids.pkl'), 'wb'))
        open(os.path.join(CACHE, 'q4_direct_compare.txt'), 'w').write("\n".join(lines) + "\n")
    else:                                         # other rows: keyed (J, Delta) in a separate file
        fn = os.path.join(CACHE, 'q4_direct_grids_rows.pkl')
        old = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
        old.update({(J, D0): v for J, v in store.items()}); pickle.dump(old, open(fn, 'wb'))
        open(os.path.join(CACHE, 'q4_direct_compare_%d.txt' % D0), 'w').write("\n".join(lines) + "\n")


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), [float(v) for v in sys.argv[5].split(',')])
    else:
        analyse([int(a) for a in sys.argv[2:]], float(os.environ.get('Q4_ROW', '400')))
