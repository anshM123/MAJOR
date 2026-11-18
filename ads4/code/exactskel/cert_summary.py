"""cert_summary: summary of the certified D-safe action runs (cert_action.py).  Usage: python3 cert_summary.py Dg -> cache/cert_<Dg>_summary.txt"""
import os, sys, pickle, glob
HERE = os.path.dirname(os.path.abspath(__file__))
Dg = int(sys.argv[1])
out = []
def pr(s=""):
    print(s); out.append(s)
for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'cert_%d_Df*.pkl' % Dg))):
    d = pickle.load(open(fn, 'rb'))
    Df = os.path.basename(fn).split('_Df')[1][:-4]
    recs = [(k, h, v) for k, r in d.items() for h, v in r['res'].items()]
    ok = [x for x in recs if x[2]['cert']]
    bad = [x for x in recs if not x[2]['cert']]
    pr(f"## Delta_g = {Dg}, Df = {Df}: {len(d)} (m, J) states, {len(recs)} (state, channel) cases, certified F > 0: {len(ok)}, not certified: {len(bad)}")
    for lo, hi in ((1.0, 1.05), (1.05, 1.1), (1.1, 1.2), (1.2, 1.5), (1.5, 2.01), (2.01, 3.01)):
        sel = [x for x in recs if lo <= x[0][0] < hi]
        if not sel:
            continue
        w = min(sel, key=lambda x: x[2]['F_lo']/x[2]['Ns'])
        rr = max(x[2]['rad']/x[2]['Ns'] for x in sel)
        pr(f"   m/Lam in [{lo:.2f},{hi:.2f}): {len(sel):3d} cases; min certified lower bound F/Ns = {w[2]['F_lo']/w[2]['Ns']:+.4f} "
           f"({w[1]}, m/Lam={w[0][0]:.2f}, J={w[0][1]}); max radius/Ns = {rr:.2e}")
    for x in bad:
        pr(f"   NOT certified: {x[1]} m/Lam={x[0][0]:.2f} J={x[0][1]}: F/Ns in [{x[2]['F_lo']/x[2]['Ns']:+.4f}, {x[2]['F_hi']/x[2]['Ns']:+.4f}]")
    ms = sorted(set(k[0] for k in d)); Js = sorted(set(k[1] for k in d))
    pr(f"   m/Lam grid: {ms}")
    pr(f"   J values: {Js}")
open(os.path.join(HERE, 'cache', 'cert_%d_summary.txt' % Dg), 'w').write("\n".join(out) + "\n")
