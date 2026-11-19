"""
r4_collapse (LARGEJ_OFFDIAG.md §2.2): is the real helicity off-diagonal a LEADING-REGGE object of the form
    f(nu; m, J) := Gamma_off / m   (the F_+- density per unit psi-weight and per unit eo2, q3 units)  =  beta^2 g(nu beta)  ?
beta = arccosh(1 + (2J+1)^2/(2 m^2)) (the AdS impact parameter of the state; b-hat = Lam beta at small beta).
Collapse test: f/beta^2 against t = nu beta for every computed even-J row (K2; direct rows where q3 loads them with Q3_DIRECT=1).
Then the psi-integrated off-diagonal of a state (b-hat, Lam) in the small-beta form:
    X(b-hat, Lam) ~ int_0^1 dp psi2(p) eo2 beta^2 g(b-hat p)   (eo2 -> 2),    vs the flat base 1023.6 / b-hat^3 (NF02 scaling limit)
Usage: python3 r4_collapse.py -> cache/r4_collapse.txt, cache/r4_g.pkl
"""
import os, sys, io, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.interpolate import CubicSpline
import r2_ladder                                    # noqa: r1 caches merged
import q3_analysis as Q

CACHE = os.path.join(HERE, 'cache')
JS = [int(v) for v in os.environ.get('R4_JS', '6,8,10,12,14,16,18,20,24,28,32,36').split(',')]


def beta_of(m, J):
    return np.arccosh(1 + (2*J + 1)**2/(2*m*m))


def main():
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    rows = []
    for J in JS:
        try:
            g = Q.even_grid(J)
        except Exception as ex:
            pr(f"  J={J}: no grid ({ex})"); continue
        for D, r in sorted(g.items()):
            if D < 300:                       # Delta = 100, 200 rows: not converged (review S3); excluded from the collapse
                continue
            m = r['m']; be = beta_of(m, J)
            rows.append(dict(J=J, D=D, m=m, beta=be, t=r['xs']*m*be, f2=r[('C2', 'off')]/m/be**2, f3=r[('C3H', 'off')]/m/be**2,
                             nu=r['xs']*m, J2m=J*J/m))
    tg = np.arange(1.0, 12.01, 0.5)
    pr("# r4: collapse of f/beta^2 = Gamma_off/(m beta^2) vs t = nu beta (C2; sign of each row fixed by overlap with the reference row)")
    ref = max(rows, key=lambda r: (r['D'], r['J']))           # most asymptotic large-J row as sign reference
    def on(r, key):
        ok = (tg >= r['t'][0]) & (tg <= r['t'][-1])
        out = np.full(len(tg), np.nan)
        out[ok] = CubicSpline(r['t'], r[key])(tg[ok])
        return out
    gref = on(ref, 'f2')
    pr("  (J, Delta)   beta    J^2/m  nu-range     | t = " + " ".join(f"{t:5.1f}" for t in tg))
    G = []
    for r in sorted(rows, key=lambda r: (r['beta'], r['D'])):
        v = on(r, 'f2'); ok = np.isfinite(v) & np.isfinite(gref)
        if ok.sum() >= 3 and np.nansum(v[ok]*gref[ok]) < 0:
            v = -v; r['f2'] = -r['f2']; r['f3'] = -r['f3']
        G.append((r, v))
        pr(f"  ({r['J']:2d},{r['D']:5.0f})  {r['beta']:.4f}  {r['J2m']:5.2f}  [{r['nu'][0]:4.0f},{r['nu'][-1]:4.0f}] | " +
           " ".join("    ." if not np.isfinite(a) else f"{a:5.2f}" for a in v))
    # consensus g(t): median over rows with J^2/m <= 1 (K2 range) and D >= 400
    good = [v for r, v in G if r['J2m'] <= 1.0 and r['D'] >= 400]
    gmed = np.nanmedian(np.array(good), axis=0); gsd = np.nanstd(np.array(good), axis=0)
    pr("  median g(t) (J^2/m <= 1, Delta >= 400): " + " ".join(f"{a:5.2f}" for a in gmed))
    pr("  spread             : " + " ".join(f"{a:5.2f}" for a in gsd))
    pickle.dump(dict(rows=rows, tg=tg, gmed=gmed, gsd=gsd), open(os.path.join(CACHE, 'r4_g.pkl'), 'wb'))
    open(os.path.join(CACHE, 'r4_collapse.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
