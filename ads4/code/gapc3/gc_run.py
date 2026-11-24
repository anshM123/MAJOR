"""
gc_run (GAP_C3.md §4): batch of certified scalar nodes (gc_cert) -> cache/gc_cert_J<J>.json (one record per (xi, rho, deg)).
Usage: OMP_NUM_THREADS=1 GC_DEG=6 python3 gc_run.py J1,J2 xi1,xi2 rho1,rho2
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gc_guard  # noqa: F401
import gc_cert as C

DREF = 1002.0                                  # any Delta: at fixed xi the data are Delta-independent (Prop. S)

if __name__ == '__main__':
    Js = [int(a) for a in sys.argv[1].split(',')]
    xis = [float(a) for a in sys.argv[2].split(',')]
    rhos = [float(a) for a in sys.argv[3].split(',')]
    for J in Js:
        fn = os.path.join(HERE, 'cache', 'gc_cert_J%d.json' % J)
        for xi in xis:
            for rho in rhos:
                db = json.load(open(fn)) if os.path.exists(fn) else {}
                key = '%.4f_%.4f_%d' % (xi, rho, C.B.DEG)
                if key in db:
                    continue
                o, _ = C.cert(J, xi*DREF, DREF, rho)
                db = json.load(open(fn)) if os.path.exists(fn) else {}
                db[key] = o
                json.dump(db, open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
                print(f"J={J} xi={xi} rho={rho} deg={C.B.DEG}: q={o['q']:.3f} bad={o['bad']} S3={o['S3']:.3g} Mg={o['Mg']:.3g} "
                      f"remh={o['remh']:.3g} t={o['t']:.1f}s", flush=True)
