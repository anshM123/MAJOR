"""
h13_hyb: hybrid (flat NF02 + leading AdS Regge curvature) D-safe action on the H13' region, all J (every integer J) at each mass.
F_hyb is the exact action minus the scalar-skeleton remainder S (h13_common docstring); S is computed separately (h13_exact.py).
Usage: python3 h13_hyb.py Dg [BMIN BMAX]   (default bhat in [20, 200]; masses mfl in [1.5, 10]; env H13_MFLS, H13_JSTRIDE)
Output (resumable, incremental): cache/h13_hyb_<Dg>_<BMIN>_<BMAX>.pkl  {(h, mfl): [records]}
"""
import os, sys, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from h13_common import flat_rows, hybrid, CACHE, Lam_of
from nf_common import state

MFLS = [round(v, 3) for v in np.concatenate([np.arange(1.5, 3.0, 0.1), np.arange(3.0, 5.0, 0.25), np.arange(5.0, 10.01, 0.5)])]


def main():
    Dg = int(sys.argv[1])
    BMIN = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
    BMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 200.0
    mfls = MFLS
    if os.environ.get('H13_MFLS'):
        mfls = [float(a) for a in os.environ['H13_MFLS'].split(',')]
    pk = os.path.join(CACHE, 'h13_hyb_%d_%g_%g.pkl' % (Dg, BMIN, BMAX))
    db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    Lam = Lam_of(Dg)
    for mfl in mfls:
        for h in ('++', '+-'):
            if (h, mfl) in db:
                continue
            t0 = time.time()
            Js = [J for J in range(0, int(BMAX*mfl) + 2) if BMIN <= (2*J + 1)/mfl <= BMAX
                  and not (h == '++' and J % 2) and not (h == '+-' and J < 4) and state(Dg, Lam, mfl, J) is not None]
            stride = int(os.environ.get('H13_JSTRIDE', '1'))
            if stride > 1:
                Js = Js[::stride] + ([Js[-1]] if (len(Js) - 1) % stride else [])
            if not Js:
                db[(h, mfl)] = []
                continue
            fl = flat_rows(h, 1/mfl**2, Js)
            recs = [hybrid(h, Dg, mfl, J, *fl[J]) for J in Js]
            db[(h, mfl)] = recs
            w = min(recs, key=lambda r: r['Fh']/r['Ns'])
            print(f"Dg={Dg} {h} m/Lam={mfl:5.2f}: {len(Js):4d} J in [{Js[0]},{Js[-1]}], min F_hyb/Ns = {w['Fh']/w['Ns']:+.4f} at J={w['J']} "
                  f"bhat={w['bh']:.1f};  min F_flat/N_flat = {min(r['Fflat']/r['Nflat'] for r in recs):+.4f}  ({time.time()-t0:.1f}s)", flush=True)
            pickle.dump(db, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    print('# done', flush=True)


if __name__ == '__main__':
    main()
