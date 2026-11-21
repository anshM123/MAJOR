"""oo_sum (HYPOTHESES_CLOSE2.md): summary of oo_t.sh outputs cache/oo_chain_<prefix>_<variant>.pkl: per J, min lambda_OO/N_s over the
near-gap states (m/Lam <= 1.3, incl. twist-edge masses) and over all states, per variant, and the jackknife error (jk0..7) at the
main variant's near-gap minimum and its maximum over near-gap states.  Usage: python3 oo_sum.py prefix [variant ...]"""
import os, sys, pickle
import numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
pre = sys.argv[1]
vs = sys.argv[2:] or ['best', 'dir', 'cv', 'nodesp', 'zero', 'lin', 'convA', 'c-1000', 'c-2000', 'c-4000', 'c-8000']


def load(v):
    f = os.path.join(C, 'oo_chain_%s_%s.pkl' % (pre, v))
    return pickle.load(open(f, 'rb'))['states'] if os.path.exists(f) else None


main = load('best')
Js = sorted(J for (Dg, J) in main)
jk = [load('jk%d' % i) for i in range(8)]
jk = [j for j in jk if j is not None]
for J in Js:
    rows = main[(300, J)]
    near = [t for t in rows if t['mfl'] <= 1.3 + 1e-9]
    w = min(near, key=lambda t: t['oo']); wa = min(rows, key=lambda t: t['oo'])
    line = f"J={J:2d} best: near-gap min {w['oo']:+.4f} at {w['mfl']:.4f} [{w['DF']}]"
    if jk:
        def val(s, m, DF):
            return next(t['oo'] for t in s[(300, J)] if abs(t['mfl'] - m) < 1e-9 and t['DF'] == DF)
        V = np.array([val(s, w['mfl'], w['DF']) for s in jk]); n = len(V)
        e = np.sqrt((n - 1)/n*((V - V.mean())**2).sum())
        A = np.array([[val(s, t['mfl'], t['DF']) for t in near] for s in jk])
        ea = np.sqrt((n - 1)/n*((A - A.mean(0))**2).sum(0))
        line += f" +- {e:.4f} (jk, {n}; max over near-gap states {ea.max():.4f})"
    line += f"; all m: {wa['oo']:+.4f} at {wa['mfl']:.2f}; EE ++ diag at that near-gap state {w['pp_ee']:+.4f}, EE 2x2 {w['lam_ee']:+.4f}"
    print(line)
    we = min(near, key=lambda t: t['lam_ee'])
    le = f"      EE block (run0) near-gap min {we['lam_ee']:+.4f} at {we['mfl']:.4f} [{we['DF']}]"
    if jk and not np.isnan(jk[0][(300, J)][0]['lam_ee']):
        E_ = np.array([next(t['lam_ee'] for t in s_[(300, J)] if abs(t['mfl'] - we['mfl']) < 1e-9 and t['DF'] == we['DF']) for s_ in jk])
        le += f" +- {np.sqrt((len(E_) - 1)/len(E_)*((E_ - E_.mean())**2).sum()):.4f} (jk)"
    le += f"; all m {min(t['lam_ee'] for t in rows):+.4f}"
    print(le)
    s = "      variants (near-gap min [m]): "
    for v in vs:
        o = load(v)
        if o is None or (300, J) not in o:
            continue
        nn = [t for t in o[(300, J)] if t['mfl'] <= 1.3 + 1e-9]
        u = min(nn, key=lambda t: t['oo'])
        ue = min(nn, key=lambda t: t['lam_ee'])
        s += f"{v} {u['oo']:+.3f} [{u['mfl']:.3f}]" + (f" (EE {ue['lam_ee']:+.3f})" if not np.isnan(ue['lam_ee']) else '') + "; "
    print(s)
