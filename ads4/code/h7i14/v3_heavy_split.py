"""v3 (§20.3): as v3_heavy_impact, but the spinning de-shift is applied only at nodes with nu >= NUC (scalar g below), to separate the
controlled window (nu >= 20, direct crystal data) from the nu-soft region.  Output: cache/v3_heavy_split.txt"""
import os, sys, io
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import v3_heavy_impact as H
buf = io.StringIO()
def pr(s=""):
    print(s, flush=True); buf.write(str(s) + "\n")
H.pr = pr
for NUC in (15.0, 20.0):
    def shift_c(J, h, k, xs, m, NUC=NUC):
        a = H.shift_spin(J, h, k, xs, m); b = H._orig_shift(J, h, k, xs, m)
        return np.where(np.asarray(xs)*m >= NUC, a, b) if k == 2 else b
    H.A.shift = shift_c
    H.run(f"spinning de-shift only at nu >= {NUC}; l9 Laplace rows", (100,), direct=False)
    H.run(f"spinning de-shift only at nu >= {NUC}; direct rows", (100,), direct=True)
open(os.path.join(HERE, 'cache', 'v3_heavy_split.txt'), 'w').write(buf.getvalue())
