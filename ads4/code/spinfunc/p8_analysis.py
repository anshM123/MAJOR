"""
p8_analysis (PARITY_ODD_ROUTEA.md sect. 11): the p4 linear-form chain for the odd-J +- element (p4_odd_analysis.grid: heavy reference
Delta_ref = 2e4 at the same nu, per-row constant c0, a_{+-} from C2 at x <= 0.25 & nu >= 15 (else the 3 smallest x), CHMRS de-shift with
the scalar g, Re) with the Laplace state/reference ratio REPLACED by the non-perturbative (direct, p8_odd_pw) ratio wherever direct state
AND reference points exist with max relative RQMC error <= P8_RELMAX (default 0.05, as n1_analysis).  K2 records: the p3 caches (rows
100 ... 3200) plus the extra-row caches cache/p8_k2_* (rows m = 1.05, 1.1, 1.2, 1.3 Lam(100)).
  direct ratio: v = exp[(logN_s - logK_s) - (logN_r - logK_r) - c0] (E_s/b_s)/(E_r/b_r)          (b = beta_OO of p2, 1x1)
Errors: from the RQMC replicas of state and reference, in Gamma (c2) units.
Outputs: cache/p4_grid_direct.pkl (P8_DIRECT=1, default) or cache/p4_grid_k2x.pkl (P8_DIRECT=0: same rows, Laplace only), read by
p6_odd_neargap.py with P6_DIRECT=1 / P6_DIRECT=k2x; cache/p8_analysis.txt (direct vs K2 per point).
Usage: L1_DEG=6 python3 p8_analysis.py [J ...]
"""
import os, sys, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import p4_odd_analysis as A
import p2_shadow_odd as P2

DIRECT = os.environ.get('P8_DIRECT', '1') == '1'
RELMAX = float(os.environ.get('P8_RELMAX', '0.05'))
H = '+-'


def k2load(kind, J):
    R = A.load(kind, J)
    fn = os.path.join(HERE, 'cache', 'p8_k2_%s_J%d.pkl' % (kind, J))
    if os.path.exists(fn):
        R.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items() if k[2] == 1.0})
    return R


def dload(kind, J):
    fn = os.path.join(HERE, 'cache', 'p8_pts_%s_J%d.pkl' % (kind, J))
    return pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}


def dget(db, D, nu, relmax=None):
    relmax = RELMAX if relmax is None else relmax
    c = [k for k in db if abs(k[0] - D) < 1e-6 and abs(k[1] - nu) < 1e-7*nu]
    if not c:
        return None
    r = db[c[0]]
    return r if np.max(np.abs(r['err']))/np.max(np.abs(r['E'])) <= relmax else None


def drep(r, D, J):
    b, logK, _ = P2.shadow_odd(float(D), J)
    bi = 1.0/b[-1, -1]
    return r['logN'] - logK, r['E'][0, 0]*bi, np.array([rep[0, 0]*bi for rep in r['reps']])


def grid(J, K=2, kinds=('C2', 'C3H'), Ds=None, direct=None, dmax=3200.0):
    direct = DIRECT if direct is None else direct
    RC = {kd: k2load(kd, J) for kd in kinds}
    DB = {kd: dload(kd, J) for kd in kinds}
    Ds = Ds or sorted(set(k[0] for k in RC['C2'] if k[0] < A.DREF and k[0] <= dmax))
    out = {}
    for D in Ds:
        m = A.mphys(D, J); mr = A.mphys(A.DREF, J)
        xs = sorted(set(round(k[1]/m, 6) for k in RC['C2'] if k[0] == D))
        xs = [x for x in xs if all(A.getrec(RC[kd], D, x*m) is not None and A.getrec(RC[kd], A.DREF, x*m) is not None for kd in kinds)]
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        rows = {}; c0 = None; ndir = 0; used = {}
        for kd in kinds:
            k = 2 if kd == 'C2' else 3
            G, Ge, GK = [], [], []
            for x in xs:
                nu = x*m; st = A.getrec(RC[kd], D, nu); rf = A.getrec(RC[kd], A.DREF, nu); xr = nu/mr
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if c0 is None:
                    c0 = dl.real
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                sc = (A.SIG[H] if k == 3 else 1.0)*A.env(H, xr, k)*A.Dc(H, J, xr)
                vK = np.exp(dl)*A.pred([st[o][-1, -1] for o in ('M0', 'M2', 'M4')], [rf[o][-1, -1] for o in ('M0', 'M2', 'M4')], K)*fk
                ds, dr = (dget(DB[kd], D, nu), dget(DB[kd], A.DREF, nu)) if direct else (None, None)
                if ds is not None and dr is not None:
                    ndir += 1; used[(kd, round(x, 6))] = True
                    ls, vs, rs = drep(ds, D, J); lr, vr_, rr = drep(dr, A.DREF, J)
                    edl = np.exp(ls - lr - c0)
                    v = edl*vs/vr_*fk
                    n = min(len(rs), len(rr))
                    s1 = (edl*rs[:n]/vr_*fk).real; s2 = (edl*vs/rr[:n]*fk).real
                    sd = np.sqrt(np.var(s1, ddof=1)/n + np.var(s2, ddof=1)/n)
                else:
                    v = vK; sd = np.nan
                G.append(sc*v); Ge.append(abs(sc)*sd); GK.append(sc*vK)
            rows[kd] = (k, np.array(G), np.array(Ge), np.array(GK))
        sel = (xs <= A.XFIT + 1e-9) & (xs*m >= A.NUFIT)
        if sel.sum() < 3:
            sel = np.argsort(xs)[:3]
        def fit(g):
            Am = np.array([A.Dc(H, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(Am, (g/A.env(H, xs, 2))[sel], rcond=None)[0][0]
        a = fit(rows['C2'][1].real); a_raw = a
        if A.DESHIFT:
            sh = {k: A.shift(J, H, k, xs, m) for k in (2, 3)}
            for _ in range(3):
                a = fit(rows['C2'][1].real - a*sh[2])
            rows = {kd: (k, G - a*sh[k], Ge, GK - a*sh[k]) for kd, (k, G, Ge, GK) in rows.items()}
        rec = dict(m=m, xs=xs, a=a, a_raw=a_raw, c0=c0, ndirect=ndir, fitx=xs[sel])
        for kd, (k, G, Ge, GK) in rows.items():
            sg = A.SIG[H] if k == 3 else 1.0
            rec[kd] = m*m*(G.real/(a*A.env(H, xs, k)) - sg*A.Dc(H, J, xs))
            rec[(kd, 'err')] = m*m*Ge/abs(a*A.env(H, xs, k))
            rec[(kd, 'K2same_a')] = m*m*(GK.real/(a*A.env(H, xs, k)) - sg*A.Dc(H, J, xs))   # Laplace values, same a and de-shift
            rec[(kd, 'im')] = m*G.imag/(a*A.env(H, xs, k))
            rec[(kd, 'direct')] = np.array([(kd, round(x, 6)) in used for x in xs])
        out[D] = rec
    return out


if __name__ == '__main__':
    Js = [int(a) for a in sys.argv[1:]] or [5, 7, 9]
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    store = {}
    pr(f"# p8_analysis: odd-J +- Gamma (c2 units x LO shape), K2 linear chain, direct override {DIRECT} (relmax {RELMAX}), de-shift {A.DESHIFT}")
    for J in Js:
        g = grid(J)
        store[J] = g
        for D in sorted(g):
            r = g[D]
            if D > 200:
                continue
            pr(f"  J={J} Delta={D:.4f} m={r['m']:.2f} (m/Lam={r['m']/98.49873:.3f}) a={r['a']:.6g} (fit x {r['fitx']}) direct points {r['ndirect']}")
            pr(f"      x       : " + " ".join(f"{x:8.2f}" for x in r['xs']))
            for kd in ('C2', 'C3H'):
                pr(f"    {kd:3s} Gamma : " + " ".join(f"{v:+8.0f}" for v in r[kd]) )
                pr(f"    {kd:3s} +-err : " + " ".join(f"{v:8.0f}" if np.isfinite(v) else "       -" for v in r[(kd, 'err')]))
                pr(f"    {kd:3s} K2    : " + " ".join(f"{v:+8.0f}" for v in r[(kd, 'K2same_a')]))
    tag = 'direct' if DIRECT else 'k2x'
    pickle.dump(store, open(os.path.join(HERE, 'cache', 'p4_grid_%s%s.pkl' % (tag, '' if A.DESHIFT else '_noshift')), 'wb'))
    open(os.path.join(HERE, 'cache', 'p8_analysis_%s%s.txt' % (tag, '' if A.DESHIFT else '_noshift')), 'w').write("\n".join(lines) + "\n")
