"""
sm_chainfix (SMOOTH_BOUND.md §4): the project chains re-run with the curvature double count removed (sm_curvfix, SM_CURV=row) or
unchanged (SM_CURV=none, reproduction check).  New outputs only: cache/sm_chainfix_<mode>_<SM_TAG>.txt/.pkl.
Modes:
  w2  Dg [Dg]  : EE_NEARGAP main chain as re-run by ee1000_w2rerun.py (w2/q3 neargap over all exact-skeleton + hybrid states,
                 ee1000 loaders + Jacobi fix); env Q3_JEVEN, Q3_JODD, Q3_EXTRA_M, W2_* as there.
  ee1000       : ee1000_chain.run() (hybrid large-J states, shape model; env EE1000_*, W2_*); also variant (i) 'state':
                 lambda with base - (cv2 + cv3) and the UNPATCHED rows' dF (only meaningful with SM_CURV=none).
  z4  300      : z4_chain (EE_DG300_SMALLNU main chain; env Z4_*, Q3_*, W2_ROWS); SM_JACFIX=1 also applies the ee1000 Jacobi fix.
"""
import os, sys, io, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
os.environ.setdefault('W2_BASIS', 'ee')
os.environ.setdefault('W2_MODES', 'ee')
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

mode = sys.argv[1]
tag = os.environ.get('SM_TAG', 'x')
C = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


t0 = time.time()
import sm_curvfix as CF
pr(f"# sm_chainfix {mode} {sys.argv[2:]}: SM_CURV={CF.MODE}; env " +
   " ".join(f"{k}={v}" for k, v in sorted(os.environ.items()) if k.startswith(('Q3_', 'W2_', 'EE1000', 'Z4_', 'SM_'))))
out = {}
if mode == 'w2':
    import ee1000_chain as CH
    W, Q = CH.W, CH.Q
    Dgs = [int(a) for a in sys.argv[2:]] or [1000]
    grids, res = Q.run(Dgs, pr)
    out = {Dg: W.perstate(res[Dg]['recs']) for Dg in Dgs}
elif mode == 'ee1000':
    import ee1000_chain as CH
    recs, info = CH.run(pr)
    w = min(recs, key=lambda r: r['lam'])
    pr(f"# overall (these J): min lam/N_s = {w['lam']:+.4f} at J={w['J']}, m/Lam={w['mfl']:.4f}; states {len(recs)}")
    # variant (i): drop the state curvature from the base instead (base + dF - diag(cv)) -- valid only on unpatched rows
    for J in sorted(set(r['J'] for r in recs)):
        R = [r for r in recs if r['J'] == J]
        lam_i = []
        for r in R:
            cv = np.array(r['curv'])*r['Ns']
            if len(r['hs']) == 2:
                lam_i.append(CH.E.lam2(r['base'] - cv, r['Ns'], r['dF']))
            else:
                lam_i.append(float((r['base'][0] - cv[0] + r['dF'][0, 0])/r['Ns'][0]))
        i = int(np.argmin(lam_i))
        pr(f"   J={J}: variant (i) base-cv [meaningful for SM_CURV=none]: min {lam_i[i]:+.4f} at m/Lam={R[i]['mfl']:.4f}; "
           f"max cv/N_s {max(max(r['curv']) for r in R):.4f}")
    out = dict(recs=recs, info=info)
elif mode == 'z4':
    if os.environ.get('SM_JACFIX', '0') == '1':
        import ee1000_common            # noqa: F401  (Jacobi fix of the +- de-shift)
    import z4_chain as Z4
    W2, Q = Z4.W2, Z4.Q
    W2.set_mode(True)
    Dgs = [int(a) for a in sys.argv[2:]] or [300]
    grids, res = Q.run(Dgs, pr)
    Z4.node_table(pr)
    for s in Z4.LOG:
        pr("# " + s)
    out = {Dg: W2.perstate(res[Dg]['recs']) for Dg in Dgs}
else:
    raise SystemExit('mode?')
pr("# " + CF.summary())
pr(f"[{time.time() - t0:.0f}s]")
open(os.path.join(C, 'sm_chainfix_%s_%s.txt' % (mode, tag)), 'w').write(buf.getvalue())
pickle.dump(out, open(os.path.join(C, 'sm_chainfix_%s_%s.pkl' % (mode, tag)), 'wb'))
