"""
n1_analysis (doc 24): the l9 linear-form chain (l9_analysis.grid) with the Laplace-truncated state/reference ratio REPLACED by the
non-perturbative (direct, n1_run) ratio wherever direct points exist (Delta = 100 row, x >= 0.25), and the near-gap 2x2 test
(l9_neargap.main) re-run on these rows.  Everything else (helicity basis, per-row constant c0, a_h fit on x <= 0.25, de-shift with the
scalar g, Re prescription, base model, other rows) is exactly l9's.
  direct element:  v = exp[(logN_s - logK_s) - (logN_r - logK_r) - c0] (H_s)_{jj}/(H_r)_{jj},   H = Ti E inv(b)^T Ti^T
  off-diagonal:    v = exp[...] (H_s)_{+-,sym}/(H_r)_{++} x sq    (sq as l9: sqrt|ref0_++/ref0_+-| of the C2 reference at LO)
Errors: propagated from the RQMC replicas (reps) of state and reference, in Gamma (c2) units.
Usage: python3 n1_analysis.py            -> cache/n1_analysis.txt, cache/n1_neargap_*.txt
"""
import os, sys, pickle, types
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l9_analysis as A
import l8_analysis as L8

DIRECT = os.environ.get('N1_DIRECT', '1') == '1'
SHK2 = pickle.load(open(os.path.join(HERE, 'cache', 'l9_shadow_k2.pkl'), 'rb'))


def dload(kind, J):
    f = os.path.join(HERE, 'cache', 'n1_pts_%s_J%d.pkl' % (kind, J))
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


RELMAX = float(os.environ.get('N1_RELMAX', '0.05'))      # direct points with max|err|/max|E| above this are not used (nu <~ 20)


def dget(db, D, nu):
    c = [k for k in db if k[0] == float(D) and abs(k[1] - nu) < 1e-7*nu]
    if not c:
        return None
    r = db[c[0]]
    return r if np.max(np.abs(r['err']))/np.max(np.abs(r['E'])) <= RELMAX else None


def Hrep(r, B, D, J):
    """helicity-projected matrices per replica (l9 physical units without logC: E binv^T), and the log prefactor logN - logK."""
    b, logK = SHK2[(float(D), int(J))][0], SHK2[(float(D), int(J))][1]
    binv = np.linalg.inv(b); na = r['E'].shape[0]
    H = np.array([B['Ti'] @ (rep[:na] @ binv.T) @ B['Ti'].T for rep in r['reps']])
    lN = r['logN']
    if not r.get('ljc', False):                        # points stored before the log-Jacobian constant fix: add log prod c_i
        import l1_crystal as Cr
        c, _ = Cr.saddle(r['nu'] if 'nu' in r else NU[0], float(D))
        lN = lN + np.sum(np.log([c[i] for i in (1, 2, 3, 4)]))
    return H, lN - logK


NU = [None]


def grid_direct(J, K, kinds=('C2', 'C3H'), XS=None, Ds=None, pr=None, errs=None):
    """copy of l9_analysis.grid with the direct override (see module doc).  errs: optional dict filled with replica spreads."""
    RC = {k: A.Recs(k, J, K) for k in kinds}
    DB = {k: dload(k, J) for k in kinds}
    Rb = dict(L8.records(J)); Rb.update({kk: v for kk, v in A._load('C2', J, '_k2').items()}); Rb.update(A._load('C2', J, ''))
    Bref = L8.basis(J, A.DREF, Rb) if J >= 4 else dict(Ti=np.eye(1))
    hels = ('++', '+-') if J >= 4 else ('++',)
    Ds = Ds or sorted(set(k[0] for k in RC['C2'].R if k[0] < A.DREF))
    out = {}
    for D in Ds:
        m = A.mphys(D, J); mr = A.mphys(A.DREF, J)
        B = L8.basis(J, D, Rb) if J >= 4 else dict(Ti=np.eye(1))
        if B is None:
            continue
        xs = XS or sorted(set(round(k[1]/m, 6) for k in RC['C2'].R if k[0] == D))
        xs = [x for x in xs if all(RC[kd].get(D, x*m) is not None and RC[kd].get(A.DREF, x*m) is not None for kd in kinds)]
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        rows = {}
        c0 = None
        ndir = 0
        for kd in kinds:
            k = 2 if kd == 'C2' else 3
            G = {h: [] for h in hels}; Go = []; Gerr = {h: [] for h in hels}
            for x in xs:
                nu = x*m; st = RC[kd].get(D, nu); rf = RC[kd].get(A.DREF, nu); xr = nu/mr
                Ms, Mr = A.helmat(st, B, K), A.helmat(rf, Bref, K)
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if c0 is None:
                    c0 = dl.real
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                ds, dr = (dget(DB[kd], D, nu), dget(DB[kd], A.DREF, nu)) if DIRECT else (None, None)
                NU[0] = nu
                use = ds is not None and dr is not None and J >= 4
                if use:
                    ndir += 1
                    Hs, ls = Hrep(ds, B, D, J); Hr, lr = Hrep(dr, Bref, A.DREF, J)
                    edl = np.exp(ls - lr - c0)
                    Hsm, Hrm = Hs.mean(0), Hr.mean(0)
                for j, h in enumerate(hels):
                    if use:
                        v = edl*Hsm[j, j]/Hrm[j, j]*fk
                        n = min(len(Hs), len(Hr))
                        vr = edl*Hs[:n, j, j]/Hrm[j, j]*fk          # state replica spread (dominant)
                        vr2 = edl*Hsm[j, j]/Hr[:n, j, j]*fk
                        sd = np.sqrt(np.var(vr.real, ddof=1)/n + np.var(vr2.real, ddof=1)/n)
                    else:
                        v = np.exp(dl)*A.pred([Mo[j, j] for Mo in Ms], [Mo[j, j] for Mo in Mr], K)*fk
                        sd = np.nan
                    sc = (A.SIG[h] if k == 3 else 1.0)*A.env(h, xr, k)*A.Dsh(h, J, xr)
                    G[h].append(sc*v); Gerr[h].append(abs(sc)*sd)
                if len(hels) == 2:
                    r2 = Mr if kd == 'C2' else A.helmat(RC['C2'].get(A.DREF, nu), Bref, K)
                    sq = np.sqrt(abs((r2[0][0, 0]/r2[0][1, 1]).real))
                    if use:
                        if kd == 'C2':
                            Hr2m = Hrm
                        else:
                            Hr2, _ = Hrep(dget(DB['C2'], A.DREF, nu), Bref, A.DREF, J); Hr2m = Hr2.mean(0)
                            NU[0] = nu
                        v = edl*0.5*(Hsm[0, 1] + Hsm[1, 0])/Hrm[0, 0]*fk*sq
                    else:
                        off = 0.5*np.array([Mo[0, 1] + Mo[1, 0] for Mo in Ms])
                        v = np.exp(dl)*A.pred(off, [Mo[0, 0] for Mo in Mr], K)*fk*sq
                    Go.append(v*np.sqrt(A.env('++', xr, k)*A.env('+-', xr, k)*A.Dsh('++', J, xr)*A.Dsh('+-', J, xr)))
            rows[kd] = (k, {h: np.array(G[h]) for h in hels}, np.array(Go), {h: np.array(Gerr[h]) for h in hels})
        a = {}
        sel = (xs <= A.XFIT + 1e-9) & (xs*m >= A.NUFIT)
        if sel.sum() < 3:
            sel = np.argsort(xs)[:3]
        def fit(h, g):
            Am = np.array([A.Dsh(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(Am, (g/A.env(h, xs, 2))[sel], rcond=None)[0][0]
        for h in hels:
            a[h] = fit(h, rows['C2'][1][h].real)
        if A.DESHIFT:
            sh = {(k, h): A.shift(J, h, k, xs, m) for k in (2, 3) for h in hels}
            for _ in range(3):
                for h in hels:
                    a[h] = fit(h, rows['C2'][1][h].real - a[h]*sh[(2, h)])
            for kd, (k, G, Go, Ge) in rows.items():
                for h in hels:
                    G[h] = G[h] - a[h]*sh[(k, h)]
        rec = dict(m=m, xs=xs, a=a, ndirect=ndir)
        for kd, (k, G, Go, Ge) in rows.items():
            for h in hels:
                sg = A.SIG[h] if k == 3 else 1.0
                rec[(kd, h)] = m*m*(G[h].real/(a[h]*A.env(h, xs, k)) - sg*A.Dsh(h, J, xs))
                rec[(kd, h, 'err')] = m*m*Ge[h]/abs(a[h]*A.env(h, xs, k))
                rec[(kd, h, 'im')] = m*G[h].imag/(a[h]*A.env(h, xs, k))
            if len(hels) == 2:
                eo = np.sqrt(A.env('++', xs, k)*A.env('+-', xs, k))
                rec[(kd, 'off')] = m*Go.real/np.sqrt(a['++']*a['+-'])/eo
                rec[(kd, 'off', 'im')] = m*Go.imag/np.sqrt(a['++']*a['+-'])/eo
        out[D] = rec
    return out


def compare(pr, Js=(4, 6, 8, 10), kinds=('C2', 'C3H')):
    """Gamma at the Delta = 100 row: direct vs K2 (vs K3 for C2 where available), per x."""
    for J in Js:
        gd = grid_direct(J, 2, kinds, Ds=[100.0])
        global DIRECT
        sv = DIRECT; DIRECT = False
        g2 = grid_direct(J, 2, kinds, Ds=[100.0])
        try:
            g3 = A.grid(J, 3, ('C2',), Ds=[100.0])
        except Exception:
            g3 = {}
        DIRECT = sv
        if 100.0 not in gd:
            continue
        r, r2 = gd[100.0], g2[100.0]
        pr(f"## J = {J}, Delta = 100 (m = {r['m']:.2f}); {r['ndirect']} direct (x, kind) points; a_h direct {r['a']} / K2 {r2['a']}")
        hels = ('++', '+-') if J >= 4 else ('++',)
        for kd in kinds:
            for h in hels + ('off',):
                key = (kd, h)
                line = f"   {kd:3s} {h:3s} x: " + " ".join(f"{x:.2f}" for x in r['xs'])
                pr(line)
                pr(f"      direct : " + " ".join(f"{v:7.0f}" for v in r[key]))
                if h != 'off':
                    pr(f"      +- err : " + " ".join(f"{v:7.0f}" if np.isfinite(v) else "      -" for v in r[(kd, h, 'err')]))
                pr(f"      K2     : " + " ".join(f"{v:7.0f}" for v in r2[key]))
                if kd == 'C2' and 100.0 in g3 and key in g3[100.0] and len(g3[100.0]['xs']) == len(r['xs']):
                    pr(f"      K3     : " + " ".join(f"{v:7.0f}" for v in g3[100.0][key]))


def neargap(pr, variant='C23H_K2', Dgs=(100,)):
    import l9_neargap as NG
    tmp = os.path.join(HERE, 'cache', 'n1_tmp'); os.makedirs(os.path.join(tmp, 'cache'), exist_ok=True)
    NG.HERE = tmp
    orig = NG.A.grid
    try:
        for lab, use in (('l9 (Laplace) rows', False), ('direct rows (x >= 0.25 at Delta = 100)', True)):
            global DIRECT
            DIRECT = use
            NG.A.grid = lambda J, K, kinds=('C2', 'C3H'), **kw: grid_direct(J, K, kinds, **kw)
            pr(f"# near-gap 2x2 [{variant}] with {lab}")
            NG.main(variant, list(Dgs))
            pr(open(os.path.join(tmp, 'cache', 'l9_neargap_%s_%s.txt' % (variant, '_'.join(map(str, Dgs))))).read())
    finally:
        NG.A.grid = orig
        DIRECT = True


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    compare(pr)
    for v in sys.argv[1:] or ['C23H_K2']:
        neargap(pr, v)
    open(os.path.join(HERE, 'cache', 'n1_analysis.txt'), 'w').write("\n".join(lines) + "\n")
