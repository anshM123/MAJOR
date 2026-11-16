"""u4: positivity scan of the AdS heavy action of Omega^U (M02 primary; M01, F_g3 for comparison) from the exact scalar skeleton.

States: m/Lam in [1, 10] (grid of uniform/t3_scan.py), bhat <= 40, J <= 400, tau >= Dg, both helicities; NP Gauss nodes in p.
Per state (adsaction.py): hybrid Fh, signed skeleton remainder S (+- error Es), absolute remainder bound E, uncontrolled weight U.
Criteria [NE + M + I]:
  (abs)  fully controlled state: Fh - E - Sp > 0;  state with uncontrolled nodes: kappa = (Fh - E - Sp)/U  (positive if the relative AdS
         correction on uncontrolled nodes is below kappa);  Sp = Cs cosh(beta)/m^2 NU (spin-specific 1/m^2 envelope, Cs = 0, 40).
  (coh)  the same with the signed skeleton remainder: Fh + S - Es instead of Fh - E (cross-k cancellations of the skeleton kept).
  supported threshold n(crit): smallest grid m/Lam such that every scanned state with m >= n Lam is positive (controlled) or has kappa >= 1.
Usage: OMEGAU_DF=3 python3 u4_scan.py Dg [Dg ...]   -> cache/u4_<Dg>_Df<Df>.txt (+ skeleton cache pickle)
"""
import sys, os, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from common import load, FUNS, SIGMA, KS, Transforms, state, flat_unstripped_parts, flat_parts, gl01
from adsaction import action
import remk

NP, NPF = 64, 400
pn, pw = gl01(NP)
pf, wf = gl01(NPF)
CS = (0.0, 40.0)

def mfl_grid():
    return np.unique(np.round(np.concatenate([np.arange(1.0, 1.3, 0.01), np.arange(1.3, 3.0, 0.05), np.geomspace(3.0, 10.0, 14)]), 6))

def jlist(h, mfl, bmax=40.0):
    jmax = int(min(400, (bmax*mfl - 1)/2))
    js = list(range(0, min(jmax, 40) + 1)); j = 41
    while j <= jmax:
        js.append(j); j = int(j*1.06) + 1
    js = [x for x in js if x <= jmax]
    if h == '++':
        js = sorted(set(x + (x % 2) for x in js if x + (x % 2) <= max(jmax, 0)))
    else:
        js = [x for x in js if x >= 4]
    return np.array(js, int)

def run(Dg):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    Fs = {n: load(n) for n in FUNS}
    TRs = {n: Transforms(PH) for n, (F, PH) in Fs.items()}
    pk = os.path.join(HERE, 'cache', f'skel_u4_{Dg}_{remk.TAG}.pkl')
    cache = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    lines = []
    def pr(s):
        print(s, flush=True); lines.append(s)
    pr(f"# u4 Dg = {Dg}, Lam = {Lam:.4f}, skeleton {remk.TAG}; m/Lam in [1,10], bhat <= 40, J <= 400; NP = {NP}")
    recs = {(n, h): [] for n in FUNS for h in SIGMA}
    t0 = time.time()
    from common import Transforms as _T
    W0 = _T({k: np.array([1.0]) for k in KS})       # for W, P0 only
    for mfl in mfl_grid():
        u = 1/mfl**2
        for h in SIGMA:
            js = jlist(h, mfl)
            if js.size == 0:
                continue
            Jm = int(js.max())
            fl = {}
            for n, (F, PH) in Fs.items():
                Jsu, parts, FWu = flat_unstripped_parts(F, PH, h, u, Jm, pf, wf)
                Jss, B2, B3, FWs = flat_parts(F, h, u, Jm)
                fl[n] = (Jsu, parts, FWu, Jss, B2, B3, FWs)
            for J in js:
                st = state(Dg, Lam, mfl, int(J))
                if st is None:
                    continue
                key = (round(float(mfl), 6), int(J))
                if key not in cache:
                    _, W, P0 = W0.regge(Lam, st['beta'])
                    cache[key] = remk.skeleton_state(st, int(J), pn*np.sqrt(u), W, P0)
                sk = cache[key]
                for n, (F, PH) in Fs.items():
                    Jsu, parts, FWu, Jss, B2, B3, FWs = fl[n]
                    i = int(np.searchsorted(Jsu, J)); i2 = int(np.searchsorted(Jss, J))
                    r = action(F, PH, TRs[n], h, st, int(J), mfl, Lam, pn, pw, sk, parts[i], FWu[i], (B2[i2], B3[i2], FWs[i2]))
                    r.update(mfl=float(mfl), J=int(J), m=st['m'], flat=float(parts[i].sum() + FWu[i]), flatstrip=float(B2[i2] + B3[i2] + FWs[i2]))
                    recs[(n, h)].append(r)
        pickle.dump(cache, open(pk, 'wb')) if mfl in (1.1, 1.5, 3.0) else None
    pickle.dump(cache, open(pk, 'wb'))
    pickle.dump(recs, open(os.path.join(HERE, 'cache', f'recs_u4_{Dg}_{remk.TAG}.pkl'), 'wb'))
    summarize(recs, pr)
    pr(f"# time {time.time() - t0:.0f} s, {len(cache)} skeleton states")
    open(os.path.join(HERE, 'cache', f'u4_{Dg}_{remk.TAG}.txt'), 'w').write("\n".join(lines) + "\n")

BANDS = ((1.0, 1.1), (1.1, 1.2), (1.2, 1.3), (1.3, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 10.01))
NGRID = (1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0)

def crit(r, Cs, mode):
    """returns (controlled, value): controlled -> margin/Ns ; else kappa."""
    Sp = Cs*np.cosh(r['beta'])/r['m']**2*r['NU']
    if mode == 'abs':
        num = r['Fh'] - r['E'] - Sp
    else:
        num = r['Fh'] + r['S'] - r['Es'] - Sp
    if r['U'] == 0:
        return True, num/r['Ns']
    return False, num/r['U']

def summarize(recs, pr):
    pr("# max|flat(sum_k B' parts) - flat(stripped)|/Ns over states: %.2e" %
       max(abs(r['flat'] - r['flatstrip'])/r['Ns'] for v in recs.values() for r in v))
    for (n, h), v in recs.items():
        mF = min(v, key=lambda r: r['Fh']/r['Ns'])
        amp = max(v, key=lambda r: r['NU']/r['Ns'])
        pr(f"  {n:4s} {h}: {len(v)} states, {sum(r['U'] > 0 for r in v)} with uncontrolled nodes; min Fh/Ns = {mF['Fh']/mF['Ns']:+.4f} at "
           f"(m/Lam={mF['mfl']:.2f}, J={mF['J']});  max NU/Ns = {amp['NU']/amp['Ns']:.3g} at (m/Lam={amp['mfl']:.2f}, J={amp['J']})")
        for b in BANDS:
            vb = [r for r in v if b[0] <= r['mfl'] < b[1]]
            if not vb:
                continue
            s = f"     [{b[0]:.1f},{min(b[1], 10):.1f}) n={len(vb):4d} unc={sum(r['U'] > 0 for r in vb):4d}"
            for mode in ('abs', 'coh'):
                for Cs in CS:
                    c = [crit(r, Cs, mode) for r in vb]
                    lo = min([x for ok, x in c if ok], default=np.inf)
                    ka = min([x for ok, x in c if not ok], default=np.inf)
                    s += f" | {mode},Cs={Cs:.0f}: ctrl {lo:+.3f} kap {ka:.3g}"
            s += f" | max E/Ns {max(r['E']/r['Ns'] for r in vb):.3g} max U/Ns {max(r['U']/r['Ns'] for r in vb):.3g} max|S|/Ns {max(abs(r['S'])/r['Ns'] for r in vb):.3g}"
            pr(s)
        for mode in ('abs', 'coh'):
            for Cs in CS:
                nsup = None
                for n0 in NGRID:
                    ok = all((x > 0) if c else (x >= 1.0) for c, x in (crit(r, Cs, mode) for r in v if r['mfl'] >= n0 - 1e-9))
                    if ok:
                        nsup = n0; break
                pr(f"     supported threshold n [{mode}, Cs={Cs:.0f}, kappa>=1 on uncontrolled]: {nsup if nsup else '> 3'}")

if __name__ == '__main__':
    for Dg in [int(a) for a in sys.argv[1:]] or [100]:
        run(Dg)
