"""
ee1000_w2rerun (EE_DG1000_LARGEJ.md §1): the EE_NEARGAP main chain (w2_chain / q3, all exact-skeleton + hybrid states) re-run with the
stable complex Jacobi de-shift (ee1000_common, EE1000_JACFIX=1) and the ee1000 row caches visible (ee1000_chain loaders); everything
else as w2 (env Q3_*, W2_*).  Output: cache/ee1000_w2rerun_<EE1000_TAG>.txt, cache/ee1000_w2rerun_<EE1000_TAG>.pkl (per-state records).
Usage: Q3_JEVEN=... Q3_JODD=... Q3_EXTRA_M=... W2_BASIS=ee EE1000_TAG=x python3 ee1000_w2rerun.py 1000 [300]
"""
import os, sys, io, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ee1000_chain as CH            # loaders + JACFIX (no states are run on import)
W, Q, C = CH.W, CH.Q, CH.C
buf = io.StringIO()
def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")
Dgs = [int(a) for a in sys.argv[1:]] or [1000]
t0 = time.time()
pr(f"# ee1000_w2rerun: JACFIX {CH.E.JACFIX}; JEVEN {Q.JEVEN}; JODD {Q.JODD}; extra m {os.environ.get('Q3_EXTRA_M', '')}; rows {os.environ.get('W2_ROWS', 'all')}")
grids, out = Q.run(Dgs, pr)
pr(f"# despiked nodes (J, Delta, x): {W.DESPIKE_LOG}")
pr(f"[{time.time() - t0:.0f}s]")
tag = os.environ.get('EE1000_TAG', 'x')
open(os.path.join(C, 'ee1000_w2rerun_%s.txt' % tag), 'w').write(buf.getvalue())
pickle.dump({Dg: W.perstate(out[Dg]['recs']) for Dg in Dgs}, open(os.path.join(C, 'ee1000_w2rerun_%s.pkl' % tag), 'wb'))
