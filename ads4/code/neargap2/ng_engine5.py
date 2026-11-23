"""
ng_engine5 (NEARGAP2.md §1.3): the parity-odd analytic-radial engine (po_ee_odd, PARITY_ODD_CLOSE §2; read-only) for the EE k = 2
image (kind EE2, the Result-D kind) of the binding state J = 5 (+-, parity-odd) at its twist edge Delta = 1005 (Delta_g = 1000) and of
its reference Delta_ref = 2e4, on a small-nu grid that overlaps the own-mass K2 row (x = 0.0125 ... 0.05 <-> nu = 12.5 ... 50).
Truncations Kt = 1, 2, 3 are stored (the Kt-difference is the truncation estimate).
Usage: python3 ng_engine5.py Delta K depth   -> cache/ng_po_eeq2_D<Delta>_J5_K<K>.pkl
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF); sys.path.insert(0, HERE)
import ng_guard  # noqa: F401
os.chdir(SF)
import numpy as np
import po_ee_odd as PO

NGC = os.path.join(HERE, 'cache')
M5 = np.sqrt(998.0*1009.0)                 # m of (1005, 5)
NUS = np.array(sorted(set(np.round(np.arange(0.5, 20.001, 0.5), 6).tolist() +
                          [round(x*M5, 6) for x in (0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025, 0.0375, 0.05)])))

if __name__ == '__main__':
    D, K, depth = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    fn = os.path.join(NGC, 'ng_po_eeq2_D%g_J5_K%d.pkl' % (D, K))
    t0 = time.time()
    r = PO.images(D, 5, K, depth, nus=NUS, kinds=('EE2',))
    pickle.dump(r, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
    print(f'ng_engine5 D={D:g} K={K} depth={depth}: {time.time() - t0:.0f}s', flush=True)
