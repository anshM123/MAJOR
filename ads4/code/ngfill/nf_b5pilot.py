"""
nf_b5pilot (NEARGAP_FILL.md §4): the 1-node pilot of ext/jobs_batch5.txt -- ext/point.py run UNCHANGED (runpy) on one job, behind the
write guard, with its output directory in cache/b5pilot/<tag>/ (the layout run_ext.py uses).  Then the point is turned into Gamma (c2)
the way the analysis will (reference record present, replica errors, failure flag).
Usage: python3 nf_b5pilot.py TAG J DELTA NU
"""
import os, sys, time, runpy, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nf_guard  # noqa: F401
import numpy as np

if __name__ == '__main__':
    tag, J, D, nu = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    out = os.path.join(HERE, 'cache', 'b5pilot', tag, f'n1_J{J}_D{D:g}_nu{nu:.6f}_s4242')
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    sys.argv = ['point.py', 'n1', str(J), repr(D), repr(nu), '4242', out]
    runpy.run_path(os.path.join(HERE, '..', '..', 'ext', 'point.py'), run_name='__main__')
    print(f"# pilot wall time {time.time() - t0:.0f}s", flush=True)
    for f in sorted(os.listdir(out)):
        if f.endswith('.pkl'):
            P = pickle.load(open(os.path.join(out, f), 'rb'))
            for k, r in P.items():
                rel = float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300)))
                print(f"  {f} key {k}: max replica rel err {rel:.2e} ({'OK' if rel <= 0.2 else 'FAILED'})", flush=True)
