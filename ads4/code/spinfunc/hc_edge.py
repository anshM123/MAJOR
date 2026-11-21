"""
hc_edge (HYPOTHESES_CLOSE.md, H_edge): variant table for the light-cone-edge window x in (0.8, Lam/m] with the direct hc edge nodes,
by running hc_chain.py (ee1000 mode) in subprocesses (one variant each; nothing shared on disk except the per-variant outputs
cache/hc_chain_ee1000_edge<tag>_*.txt/.pkl).
Variants (all with SM_CURV=row, odd despike, shape model K = 4 on x <= 0.8 as the protocol):
  base      : no injection, EE1000_XMAX = 0.8 (the protocol: c, s, r held beyond the last node)
  A         : direct nodes for x > 0.8 only (K2 kept at x = 0.8), edge model = constant-coefficient shape fit on x >= 0.8
  B         : direct nodes from x = 0.8 on (K2 at 0.8 replaced), same edge model           [main]
  B_lin / B_spl : B with linear / natural-spline interpolation of the edge nodes instead of the shape fit
  B_nd      : B, injected nodes exempt from the Hampel despike
  B_K1      : B with linear-in-x edge coefficients (needs >= 7 edge nodes)
  glob      : one shape fit (K = 4) over all nodes up to the last edge node
  jk        : B with replica i of the injected nodes left out, i = 0..7 -> jackknife error of every per-state lambda
Usage: python3 hc_edge.py J Delta [Dg] [variants]   (odd J: 1x1 +-)   -> cache/hc_edge_J<J>.txt
"""
import os, sys, subprocess, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); C = os.path.join(HERE, 'cache')
J = int(sys.argv[1]); D = float(sys.argv[2]); DG = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
VARS = sys.argv[4].split(',') if len(sys.argv) > 4 else ['base', 'A', 'B', 'B_lin', 'B_spl', 'B_nd', 'B_K1', 'glob', 'jk']
odd = J % 2 == 1
common = dict(W2_DMIN=str(1000 if DG == 1000 else 0), EE1000_DG=str(DG), EE1000_JEVEN='' if odd else str(J), EE1000_JODD=str(J) if odd else '',
              EE1000_MODEL='shape', EE1000_K='4', EE1000_XMIN='0.05', SM_CURV='row', OMP_NUM_THREADS='1', HC_NODES='0')
spec = {
    'base': dict(EE1000_XMAX='0.8', HC_INJ=''),
    'A': dict(EE1000_XMAX='0.99', HC_INJ='%d:%g' % (J, D), HC_IXMIN='0.81'),
    'B': dict(EE1000_XMAX='0.99', HC_INJ='%d:%g' % (J, D)),
    'B_lin': dict(EE1000_XMAX='0.99', HC_INJ='%d:%g' % (J, D), HC_EREP='linear'),
    'B_spl': dict(EE1000_XMAX='0.99', HC_INJ='%d:%g' % (J, D), HC_EREP='spline'),
    'B_nd': dict(EE1000_XMAX='0.99', HC_INJ='%d:%g' % (J, D), HC_DESPIKE_INJ='0'),
    'B_K1': dict(EE1000_XMAX='0.99', HC_INJ='%d:%g' % (J, D), HC_EK='1'),
    'glob': dict(EE1000_XMAX='0.99', HC_INJ='%d:%g' % (J, D), HC_EDGE_MODEL='global'),
}


def run(tag, extra):
    env = dict(os.environ); env.update(common); env.update(extra); env['HC_TAG'] = 'edge%d_%s' % (J, tag)
    p = subprocess.run([sys.executable, os.path.join(HERE, 'hc_chain.py'), 'ee1000'], env=env, capture_output=True, text=True, cwd=HERE)
    fn = os.path.join(C, 'hc_chain_ee1000_%s.pkl' % env['HC_TAG'])
    if p.returncode != 0 or not os.path.exists(fn):
        return None, p.stdout[-2000:] + p.stderr[-2000:]
    o = pickle.load(open(fn, 'rb'))
    txt = open(os.path.join(C, 'hc_chain_ee1000_%s.txt' % env['HC_TAG'])).read()
    return o['recs'], txt


lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


res = {}
for v in VARS:
    if v == 'jk':
        continue
    recs, txt = run(v, spec[v])
    if recs is None:
        pr(f"  {v}: FAILED\n{txt}"); continue
    res[v] = recs
    ef = [l for l in txt.splitlines() if l.startswith('# edge fit') or 'despiked' in l]
    w = min(recs, key=lambda r: r['lam'])
    g = recs[0]
    pr(f"  {v:6s}: gap m/Lam={g['mfl']:.4f}: lam/N_s {g['lam']:+.4f} (diag corr " + "/".join(f"{d:+.4f}" for d in g['diag']) + f"); "
       f"min over m {w['lam']:+.4f} at {w['mfl']:.3f}; near gap (<= 1.1): " +
       " ".join(f"{r['mfl']:.3f}:{r['lam']:+.4f}" for r in recs if r['mfl'] <= 1.1 + 1e-9))
    for l in ef:
        pr("          " + l[:260])
if 'jk' in VARS and 'B' in res:
    jk = []
    for i in range(8):
        recs, txt = run('jk%d' % i, dict(spec['B'], HC_JK=str(i)))
        if recs is not None:
            jk.append([r['lam'] for r in recs])
    if jk:
        jk = np.array(jk); n = len(jk)
        err = np.sqrt((n - 1)/n*np.sum((jk - jk.mean(0))**2, axis=0))
        pr("  B jackknife (8 replicas of the injected nodes): per m/Lam lambda +- err: " +
           " ".join(f"{r['mfl']:.3f}:{r['lam']:+.4f}+-{e:.4f}" for r, e in zip(res['B'], err)))
        dg = np.array([[0.0]])
open(os.path.join(C, 'hc_edge_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
