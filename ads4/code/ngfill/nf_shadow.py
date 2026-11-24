"""
nf_shadow (NEARGAP_FILL.md §1): the parity-odd <TTO> shadow coefficient of an odd-J own row with no K2 row (J = 9, 11 at Delta = 1000 + J),
needed to turn its direct points into chain records (hc_nodes.to_record) -- exactly the function nf_rows.py (= neargap2/ng_rows) calls
when it builds a K2 row (spinfunc/p2_shadow_odd.shadow_all at its defaults, cached in cache/nf_p2shadow.pkl), without the row itself.
Usage: python3 nf_shadow.py J Delta [J Delta ...]
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nf_rows as R
import p2_shadow_odd as P2

if __name__ == '__main__':
    R.setup_odd()
    a = sys.argv[1:]
    for J, D in zip(a[0::2], a[1::2]):
        t0 = time.time()
        b, logK, resid = P2.shadow_odd(float(D), int(J))
        print(f"shadow (Delta={float(D):g}, J={int(J)}): b_OO = {complex(__import__('numpy').asarray(b)[-1, -1]):.6e}, logK = {complex(logK):.6f}, "
              f"resid {resid} [{time.time() - t0:.0f}s]", flush=True)
