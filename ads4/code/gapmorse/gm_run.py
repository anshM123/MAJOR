"""
gm_run (GAP_MORSE.md §4): batch of certified Taylor-model data (gm_cert) on Hermitian balls -> cache/gm_cert_J<J>.json
(one record per (xi, rho, deg)).  Usage: OMP_NUM_THREADS=1 GC_DEG=8 python3 gm_run.py J1,J2 xi1,xi2 rho1,rho2
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gm_guard  # noqa: F401
import gm_cert as C

if __name__ == '__main__':
    Js = [int(a) for a in sys.argv[1].split(',')]
    xis = [float(a) for a in sys.argv[2].split(',')]
    rhos = [float(a) for a in sys.argv[3].split(',')]
    for J in Js:
        fn = os.path.join(HERE, 'cache', 'gm_cert_J%d.json' % J)
        for xi in xis:
            for rho in rhos:
                db = json.load(open(fn)) if os.path.exists(fn) else {}
                key = '%.4f_%.4f_%d' % (xi, rho, C.B.DEG)
                if key in db:
                    continue
                try:
                    o = C.cert(J, xi, rho)
                except Exception as e:          # e.g. a non-finite enclosure; recorded as a failed radius
                    o = dict(J=J, xi=xi, rho=rho, deg=C.B.DEG, bad=-1, err=repr(e)[:200])
                db = json.load(open(fn)) if os.path.exists(fn) else {}
                db[key] = o
                json.dump(db, open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
                print(f"J={J} xi={xi} rho={rho} deg={C.B.DEG}: q={o.get('q', -1):.3f} bad={o['bad']} remPhi={o.get('remPhi', -1):.3g} "
                      f"remh={o.get('remh', -1):.3g} t={o.get('t', 0):.1f}s", flush=True)
