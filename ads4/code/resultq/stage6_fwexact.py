"""Result Q, Stage 6 (RESULT_Q_STAGE6.md §1.3): exact CHMRS scalar-skeleton FORWARD moments ex4(0), Lam^2-free d ex4/d nu^2|_0, ex5(0) on
NEW heavy states outside Stage 3's grid, to validate the HYB forward model (stage6_model.py) where Stage 3 had no exact data:
   set 'far'  : J in {0, 2, 4, 5, 6, 8, 12, 16, 20, 24}, m/Lam in {6, 10, 20}        (m > 4 Lam, bhat < 20 for most)
   set 'bigJ' : J in {26, 28, 32, 40}, m/Lam in {2.6, 4, 8} (J = 40 only at 4)       (J > 24, bhat ~ 7 - 31: the overlap band)
Same numerics as stage3_fwexact.one (read-only import): N0 = 200, two independent cubic nu^2-fits.  Records keyed (round(Delta, 6), J) with
'mfl', 'Dg'.  Usage: OMP_NUM_THREADS=1 python3 stage6_fwexact.py Dg iworker nworkers -> cache/stage6_fw_<Dg>_Df3_w<i>.pkl (resumable)
"""
import os, sys, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import stage3_fwexact as S3F

CACHE = os.path.join(HERE, 'cache')


def states(Dg):
    L = S3F.lam(Dg); out = []
    for J in (0, 2, 4, 5, 6, 8, 12, 16, 20, 24):
        for f in (6.0, 10.0, 20.0):
            out.append((J, f))
    for J in (26, 28, 32):
        for f in (2.6, 4.0, 8.0):
            out.append((J, f))
    out.append((40, 4.0))
    res = []
    for J, f in out:
        m = f*L
        mg = np.sqrt((Dg - 2.0)*(Dg - 1.0 + 2*J))
        if m < mg:
            continue
        al = 0.5*(-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*m*m))
        res.append((float(al + J + 2), J, f))
    return res


def main():
    Dg = int(sys.argv[1]); iw = int(sys.argv[2]) if len(sys.argv) > 2 else 0; nw = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    fn = os.path.join(CACHE, 'stage6_fw_%d_Df%g_w%d.pkl' % (Dg, S3F.DF, iw))
    db = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    st = sorted(states(Dg), key=lambda s: (s[1], s[2]))[iw::nw]
    for (Delta, J, mfl) in st:
        key = (round(Delta, 6), J)
        if key in db:
            continue
        t = time.time()
        r = S3F.one(Delta, J); r['mfl'] = mfl; r['Dg'] = Dg
        db[key] = r
        pickle.dump(db, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
        print(f"Dg={Dg} J={J:2d} m/Lam={mfl:.4f} Delta={Delta:.3f}: ex4(0)={r[4]['v0']:.10f} (+-{r[4]['ev0']:.1e}) "
              f"d4={r[4]['d0']:.6e} (+-{r[4]['ed0']:.1e}) ex5(0)={r[5]['v0']:.10f} canc {r['canc']:.1e}  {time.time() - t:.1f}s", flush=True)


if __name__ == '__main__':
    main()
