"""cert_qsummary: summary of the exact-p-integral certificates (cert_quad.py on top of cert_action.py).
Usage: python3 cert_qsummary.py Dg -> cache/certq_<Dg>_summary.txt"""
import os, sys, pickle, glob
HERE = os.path.dirname(os.path.abspath(__file__))
Dg = int(sys.argv[1])
out = []
def pr(s=""):
    print(s); out.append(s)
for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'certq_%d_Df*.pkl' % Dg))):
    d = pickle.load(open(fn, 'rb'))
    Df = os.path.basename(fn).split('_Df')[1][:-4]
    recs = [(k, h, v) for k, r in d.items() for h, v in r['res'].items()]
    ok = [x for x in recs if x[2]['cert']]; bad = [x for x in recs if not x[2]['cert']]
    pr(f"## Delta_g = {Dg}, Df = {Df}: {len(d)} states, {len(recs)} cases; exact p-integral certified F > 0: {len(ok)}; not: {len(bad)}")
    for lo, hi in ((1.0, 1.05), (1.05, 1.1), (1.1, 1.2), (1.2, 1.5), (1.5, 2.01), (2.01, 3.01)):
        for jl, jh in ((0, 10), (11, 40)):
            sel = [x for x in recs if lo <= x[0][0] < hi and jl <= x[0][1] <= jh]
            if not sel:
                continue
            w = min(sel, key=lambda x: x[2]['Fexact_lo']/x[2]['Ns'])
            e = max(sel, key=lambda x: x[2]['err']/x[2]['Ns'])
            pr(f"   m/Lam [{lo:.2f},{hi:.2f}) J {jl}-{jh}: {len(sel):3d} cases; min certified F_exact/Ns >= {w[2]['Fexact_lo']/w[2]['Ns']:+.4f} "
               f"({w[1]}, {w[0][0]:.2f}, J={w[0][1]}); max quadrature error bound/Ns {e[2]['err']/e[2]['Ns']:.1e} ({e[1]}, {e[0][0]:.2f}, J={e[0][1]})")
    for x in bad:
        pr(f"   NOT certified (exact integral): {x[1]} m/Lam={x[0][0]:.2f} J={x[0][1]}: F64_lo/Ns {x[2]['F_lo']/x[2]['Ns']:+.4f}, quad bound {x[2]['err']}")
open(os.path.join(HERE, 'cache', 'certq_%d_summary.txt' % Dg), 'w').write("\n".join(out) + "\n")
