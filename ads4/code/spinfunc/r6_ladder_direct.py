"""
r6_ladder_direct (LARGEJ_OFFDIAG.md §2.3): the fixed-b-hat Regge ladder at Delta_g = 300 with NON-PERTURBATIVE (RQMC) values:
    (336, 18) [q2 tag _h K2 + q4 direct]   vs   (672, 36) [r1 K2 + r1 direct]       b-hat = 33.07 / 32.55, beta = 0.1107 / 0.1090
Real off-diagonal Gamma_off (C2) on the common nu-nodes 25 ... 67, through n1_analysis.grid_direct (same chain as q4: K2 helicity
basis, per-row a_h fit, Re, l9 units) restricted to the six direct nodes (so a_h comes from direct values only), and the K2 values on the
same nodes.  Predictions:  leading-Regge (r4):  Gamma_off(672,36)/Gamma_off(336,18) = (m beta^2)_36/(m beta^2)_18 at t = nu beta;
"O(1/m^2) at fixed (nu, beta)": ratio ~ 1/2.
Usage: python3 r6_ladder_direct.py -> cache/r6_ladder_direct.txt
"""
import os, sys, glob, pickle, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import r2_ladder                               # noqa: r1 K2 caches merged into the q3 loaders
import q3_analysis as Q
import n1_analysis as N1
import r1_run as R1

CACHE = os.path.join(HERE, 'cache')
for fn in glob.glob(os.path.join(CACHE, 'q2_shadow_k2*.pkl')) + [os.path.join(CACHE, 'r1_shadow_k2.pkl')]:
    if os.path.exists(fn):
        N1.SHK2.update(pickle.load(open(fn, 'rb')))
NUS = [25.049, 33.399, 41.749, 50.098, 58.448, 66.798]          # r1 keys (rounded); the (336, 18) row uses x * m exactly


def nus_for(J, D):
    if (J, D) == (18, 336.0):
        m = R1.mphys(336.0, 18)
        return [x*m for x in (0.075, 0.1, 0.125, 0.15, 0.175, 0.2)]
    return NUS


def dload(kind, J):
    db = {}
    for f in (os.path.join(CACHE, 'q4_n1_%s_J%d.pkl' % (kind, J)), os.path.join(CACHE, 'r1_n1_%s_J%d.pkl' % (kind, J))):
        if os.path.exists(f):
            db.update(pickle.load(open(f, 'rb')))
    return db


N1.dload = dload


def row(J, D, direct):
    m = R1.mphys(D, J)
    N1.DIRECT = direct
    g = N1.grid_direct(J, 2, ('C2', 'C3H'), XS=[nu/m for nu in nus_for(J, D)], Ds=[D])
    N1.DIRECT = True
    return g.get(D)


def main():
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    pr("# r6: fixed-b-hat ladder at Delta_g = 300, Gamma_off (l9 units) on nu = " + ", ".join(f"{v:.1f}" for v in NUS))
    res = {}
    for J, D in ((18, 336.0), (36, 672.0)):
        m = R1.mphys(D, J); be = np.arccosh(1 + (2*J + 1)**2/(2*m*m))
        for direct in (True, False):
            r = row(J, D, direct)
            if r is None:
                pr(f"  ({J},{D:.0f}) {'direct' if direct else 'K2'}: no row"); continue
            res[(J, direct)] = r
            pr(f"  ({J},{D:.0f}) m={m:.1f} beta={be:.4f} {'direct' if direct else 'K2    '} (direct nodes used {r['ndirect']}), a_h = "
               + ", ".join(f"{h} {v:.4f}" for h, v in r['a'].items()))
            for kd in ('C2', 'C3H'):
                pr(f"      {kd:3s} Gamma_off     : " + " ".join(f"{v:7.2f}" for v in r[(kd, 'off')]))
            pr(f"      C2  f/beta^2 = Gamma_off/(m beta^2): " + " ".join(f"{v:6.2f}" for v in r[('C2', 'off')]/(m*be*be)) +
               "   at t = " + " ".join(f"{nu*be:5.2f}" for nu in NUS))
            for h in ('++', '+-'):
                pr(f"      C2  {h} Gamma_hh   : " + " ".join(f"{v:7.0f}" for v in r[('C2', h)]) + "   (err " +
                   " ".join(f"{v:5.0f}" for v in r[('C2', h, 'err')]) + ")")
    for direct in (True, False):
        if (18, direct) in res and (36, direct) in res:
            a, b = res[(18, direct)][('C2', 'off')], res[(36, direct)][('C2', 'off')]
            if a @ b < 0:
                b = -b
            m18, m36 = R1.mphys(336, 18), R1.mphys(672, 36)
            b18, b36 = np.arccosh(1 + 37**2/(2*m18**2)), np.arccosh(1 + 73**2/(2*m36**2))
            pr(f"  {'direct' if direct else 'K2'}: Gamma_off(36)/Gamma_off(18) per node: " + " ".join(f"{v:5.2f}" for v in b/a) +
               f";  rms ratio {np.sqrt((b*b).mean()/(a*a).mean()):.3f}  (leading Regge m beta^2 ratio {m36*b36**2/(m18*b18**2):.3f}; "
               f"O(1/m^2) at fixed (nu,beta): {m18/m36:.3f})  => s = {np.log(np.sqrt((b*b).mean()/(a*a).mean()))/np.log(m36/m18):+.2f}")
    open(os.path.join(CACHE, 'r6_ladder_direct.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
