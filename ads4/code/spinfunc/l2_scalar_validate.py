"""
l2 (M6 step 1 validation): Laplace expansion of the CLS-D crystal integral (l1_crystal) vs the EXACT scalar skeleton (exactskel e1).

Exact side.  CLS-D Psi_{2,nu} with f = (v'-u')/(u'v') is (up to a nu-only factor N(nu)) the J = -1 Lorentzian inversion
Phi_{-1,nu}[(v'-u')G/(u'v')], which by CHMRS's magic block identity equals
      X(nu) = g(nu) C2(nu + i) + g(-nu) C2(nu - i),   g(nu) = ((2-d)/2 + i nu)/(i nu) * (gamma_7(nu + i)/gamma_7(nu))^2,
with C2 = CHMRS C_{2,nu} = ex_2/m^4 (e1 exact, Delta_phi = 3).  C2(nu +- i) at fixed m: x -> x +- i/m; exact C2(nu) at the node plus
the shift terms from the leading form C0 = E_2(x) P_J(1-2x^2)/m^4 (error O(m^-4) relative).
Laplace side.  Y = exp(logY0) (1 + K1 + K2)  (l1_crystal.laplace; includes 1/beta_{Delta,J} exactly; sin^2 stripped on both sides).
Normalisation: e1 is the heavy DENSITY C2[G]/(q_{Delta,J} 2 sin^2) (CHMRS eq. definitionofheavydensity), so Y -> Y/q_{Delta,J}.
Test: |Y/(q X)| = N(nu) n(J) for all states (m, J) at the same node nu  (n(J): J-only structure normalisation).  Fit log|Y/X| = a_nu + b_J
and convert residuals eps to scalar-remainder units:  Delta c2 = m^2 eps ex/E_2  (c2 = m^2 rho of m1/e1).
Usage: python3 l2_scalar_validate.py Dg [Dg ...]   -> cache/l2_<Dg>.pkl (incremental), cache/l2_scalar_validate.txt
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, mpmath as mp
from scipy.special import eval_legendre
import l1_crystal as C
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
from exactmom import log_q                      # read-only import: CHMRS q_{Delta,J} (heavy-density normalisation)

EC = os.path.join(HERE, '..', 'exactskel', 'cache')
OUT = os.path.join(HERE, 'cache')
XW = (0.3, 0.7)
JS = (0, 2, 4, 6)
MFL = (1.0, 1.05, 1.1, 1.2, 1.3, 1.5)


def E2(x):
    return (2 - x*x)/(1 - x*x)**2


def C0x(x, J):
    """E_2(x) P_J(1-2x^2) and its first three x-derivatives (numerically exact via mpmath)."""
    f = lambda t: (2 - t*t)/(1 - t*t)**2*mp.legendre(J, 1 - 2*t*t)
    return [complex(mp.diff(f, mp.mpf(x), n)) for n in range(4)]


def lgam7(nu):
    a = 7; d = 3
    return mp.loggamma((1 + a - d/2 - 1j*nu)/2) + mp.loggamma((1 + a - d/2 + 1j*nu)/2)


def gfac(nu):
    nu = mp.mpf(nu)
    gp = ((2 - 3)/2 + 1j*nu)/(1j*nu)*mp.exp(2*(lgam7(nu + 1j) - lgam7(nu)))
    gm = ((2 - 3)/2 + 1j*(-nu))/(1j*(-nu))*mp.exp(2*(lgam7(-nu + 1j) - lgam7(-nu)))     # g(-nu)
    return complex(gp), complex(gm)


def Xexact(nu, m, J, ex2):
    x = nu/m
    c = [v/m**4 for v in C0x(x, J)]
    d = 1j/m
    Cn = ex2/m**4
    gp, gm = gfac(nu)
    Xp = Cn + d*c[1] + d*d*c[2]/2 + d**3*c[3]/6
    Xm = Cn - d*c[1] + d*d*c[2]/2 - d**3*c[3]/6
    return gp*Xp + gm*Xm, (gp + gm)*c[0]


def run(Dg):
    db = pickle.load(open(os.path.join(EC, 'e1_%d_Df3.pkl' % Dg), 'rb'))
    pk = os.path.join(OUT, 'l2_%d.pkl' % Dg)
    res = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    t0 = time.time(); n = 0
    for key in sorted(db):
        mfl, J = key
        if J not in JS or round(mfl, 3) not in MFL:
            continue
        r = db[key]
        for i, (x, nu) in enumerate(zip(r['x'], r['nus'])):
            if not (XW[0] <= x <= XW[1]) or (key, i) in res:
                continue
            L = C.laplace(nu, r['Delta'], J)
            X, Xlead = Xexact(nu, r['m'], J, r['ex'][2][i])
            res[(key, i)] = dict(mfl=mfl, J=J, m=r['m'], Delta=r['Delta'], x=x, nu=nu, ex2=r['ex'][2][i], rho2=r['rho'][2][i],
                                 logY0=L['logY0'], K1=L['K1'], K2=L['K2'], X=X, Xlead=Xlead)
            n += 1
            if n % 20 == 0:
                pickle.dump(res, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
                print(f"  Dg={Dg}: {n} new nodes, {time.time() - t0:.0f}s", flush=True)
    pickle.dump(res, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    return res


def analyse(Dg, res, pr):
    keys = sorted(res)
    for k in keys:                       # (re)compute the exact side from the stored exact data (cheap)
        q = res[k]; q['X'], q['Xlead'] = Xexact(q['nu'], q['m'], q['J'], q['ex2'])
        if 'logq' not in res[k]:
            res[k]['logq'] = log_q(res[k]['Delta'], res[k]['J'], 3.0)
    nus = sorted(set(round(res[k]['nu'], 9) for k in keys))
    pr(f"## Dg = {Dg}: {len(keys)} (state, node) points, {len(nus)} nodes nu, states with J in {JS}, m/Lam in {MFL}")
    for lab, fy in (('Y0 only', lambda q: 0.0), ('Y0(1+K1)', lambda q: q['K1']), ('Y0(1+K1+K2)', lambda q: q['K1'] + q['K2'])):
        rows = []
        for k in keys:
            q = res[k]
            ly = q['logY0'].real - q['logq'] + np.log(abs(1 + fy(q)))
            lr = ly - np.log(abs(q['X']))
            rows.append((round(q['nu'], 9), q['J'], lr, q))
        # least squares: lr = a_nu + b_J (b_0 = 0)
        iv = {v: j for j, v in enumerate(nus)}; Jl = [J for J in JS if J != 0]
        Acols = len(nus) + len(Jl)
        A = np.zeros((len(rows), Acols)); b = np.zeros(len(rows))
        for r_, (v, J, lr, q) in enumerate(rows):
            A[r_, iv[v]] = 1
            if J != 0:
                A[r_, len(nus) + Jl.index(J)] = 1
            b[r_] = lr
        wts = np.array([q['m']**2*abs(q['ex2'])/E2(q['x']) for (v, J, lr, q) in rows])   # residuals weighted in c2 units
        sol, *_ = np.linalg.lstsq(A*wts[:, None], b*wts, rcond=None)
        eps = b - A @ sol
        dc2 = np.array([q['m']**2*e*q['ex2']/E2(q['x']) for e, (v, J, lr, q) in zip(eps, rows)])
        c2 = np.array([q['m']**2*q['rho2'] for (v, J, lr, q) in rows])
        pr(f"  {lab:12s}: rms eps = {np.sqrt(np.mean(eps**2)):.2e}, max |eps| = {np.abs(eps).max():.2e};  in c2 units: rms |Delta c2| = "
           f"{np.sqrt(np.mean(dc2**2)):.3g}, max = {np.abs(dc2).max():.3g}  (exact |c2| median {np.median(np.abs(c2)):.3g}, max {np.abs(c2).max():.3g}); "
           f"n(J) = {', '.join('%d:%.4f' % (J, np.exp(sol[len(nus) + j])) for j, J in enumerate(Jl))}")
        if lab == 'Y0(1+K1+K2)':
            for J in JS:
                sel = [j for j, (v, JJ, lr, q) in enumerate(rows) if JJ == J]
                rel = np.abs(dc2[sel])/np.maximum(np.abs(c2[sel]), 1.0)
                pr(f"     J={J}: {len(sel)} pts, max |Delta c2| = {np.abs(dc2[sel]).max():.3g}, median |Delta c2|/|c2| = {np.median(rel):.2%}, "
                   f"90th pct = {np.percentile(rel, 90):.2%}")
            ph = np.array([np.angle(np.exp(1j*(q['logY0'].imag)) * (1 + q['K1'] + q['K2'])/q['X']) for (v, J, lr, q) in rows])
            pr(f"     phase of Y/X: spread over states at fixed nu (max over nu) = "
               f"{max(np.ptp(np.unwrap(ph[[j for j, rw in enumerate(rows) if rw[0] == v]])) for v in nus if sum(rw[0] == v for rw in rows) > 1):.2e} rad")
    return


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    for Dg in [int(a) for a in sys.argv[1:]] or [1000]:
        res = run(Dg)
        analyse(Dg, res, pr)
    open(os.path.join(OUT, 'l2_scalar_validate_%s.txt' % '_'.join(sys.argv[1:] or ['1000'])), 'w').write("\n".join(out) + "\n")


def residual_table(Dg, res, pr, JJ=(0,)):
    """eps per point (same fit as analyse, full Y) listed by x for given J, in c2 units."""
    keys = sorted(res)
    nus = sorted(set(round(res[k]['nu'], 9) for k in keys))
    iv = {v: j for j, v in enumerate(nus)}; Jl = [J for J in JS if J != 0]
    rows = []
    for k in keys:
        q = res[k]
        lr = q['logY0'].real - q['logq'] + np.log(abs(1 + q['K1'] + q['K2'])) - np.log(abs(q['X']))
        rows.append((round(q['nu'], 9), q['J'], lr, q))
    A = np.zeros((len(rows), len(nus) + len(Jl))); b = np.array([r[2] for r in rows])
    for r_, (v, J, lr, q) in enumerate(rows):
        A[r_, iv[v]] = 1
        if J != 0:
            A[r_, len(nus) + Jl.index(J)] = 1
    wts = np.array([r[3]['m']**2*abs(r[3]['ex2'])/E2(r[3]['x']) for r in rows])
    sol, *_ = np.linalg.lstsq(A*wts[:, None], b*wts, rcond=None); eps = b - A @ sol
    return rows, eps


def richardson(pr, gaps=(100, 300, 1000)):
    """Per matched point (m/Lam, J, node) present at all gaps: Delta c2(Dg) = Laplace - exact (weighted fit, c2 units);
    fit Delta c2 = d0 + d1/m + d2/m^2 over the gaps; d0 = mismatch of the 1/m^2 coefficient itself (truncation K3/m^3 removed)."""
    D = {}
    for Dg in gaps:
        res = pickle.load(open(os.path.join(OUT, 'l2_%d.pkl' % Dg), 'rb'))
        for k in res:
            q = res[k]; q['X'], q['Xlead'] = Xexact(q['nu'], q['m'], q['J'], q['ex2'])
            q['logq'] = log_q(q['Delta'], q['J'], 3.0)
        rows, eps = residual_table(Dg, res, pr)
        for (v, J, lr, q), e in zip(rows, eps):
            key = (round(q['mfl'], 6), J, round(q['x'], 9))
            D.setdefault(key, {})[Dg] = (q['m'], e*q['m']**2*q['ex2']/E2(q['x']), q['m']**2*q['rho2'])
    out = []
    for key, dd in D.items():
        if len(dd) < len(gaps):
            continue
        ms = np.array([dd[g][0] for g in gaps]); dc = np.array([dd[g][1] for g in gaps])
        A = np.c_[np.ones(len(gaps)), 1/ms, 1/ms**2][:, :len(gaps)]
        d = np.linalg.solve(A, dc) if len(gaps) == 3 else np.linalg.lstsq(A, dc, rcond=None)[0]
        out.append((key, d, dd[gaps[-1]][2], dc))
    d0 = np.array([o[1][0] for o in out]); c2 = np.array([o[2] for o in out]); d1 = np.array([o[1][1] for o in out])
    rel = np.abs(d0)/np.abs(c2)
    big = np.abs(c2) > 10
    pr(f"## Richardson over Dg = {gaps}: {len(out)} matched points (x in {XW}, J in {JS})")
    pr(f"   raw |Delta c2| at Dg = {gaps}: rms " + ", ".join('%.2f' % np.sqrt(np.mean(np.array([o[3][i] for o in out])**2)) for i in range(len(gaps))))
    pr(f"   extrapolated |d0| (the 1/m^2 coefficient mismatch): rms {np.sqrt(np.mean(d0**2)):.3f}, max {np.abs(d0).max():.3f}; "
       f"|d0|/|c2|: median {np.median(rel):.2%}, 90th pct {np.percentile(rel, 90):.2%}; on |c2| > 10 ({big.sum()} pts): median "
       f"{np.median(rel[big]):.2%}, 90th pct {np.percentile(rel[big], 90):.2%}, max {rel[big].max():.2%}")
    pr(f"   1/m^3 truncation coefficient d1 (= K3-type, c2 units x m): median |d1| {np.median(np.abs(d1)):.0f}, max {np.abs(d1).max():.0f}")
    for J in JS:
        s = np.array([o[0][1] == J for o in out])
        if s.any():
            pr(f"     J={J}: {s.sum()} pts, rms d0 {np.sqrt(np.mean(d0[s]**2)):.3f}, max |d0| {np.abs(d0[s]).max():.3f}, "
               f"median |d0|/|c2| {np.median(rel[s]):.2%}  (|c2| range {np.abs(c2[s]).min():.1f}-{np.abs(c2[s]).max():.1f})")
    return out


def summary(pr):
    from scipy.special import eval_legendre
    out = richardson(pr)
    pj = np.array([abs(eval_legendre(o[0][1], 1 - 2*o[0][2]**2)) for o in out])
    d0 = np.array([o[1][0] for o in out]); c2 = np.array([o[2] for o in out])
    for cut in (0.05,):
        s = pj > cut
        r = np.abs(d0[s])/np.maximum(np.abs(c2[s]), 1e-12); b = s & (np.abs(c2) > 5)
        pr(f"   excluding the near-zeros of P_J(1-2x^2) (|P_J| <= {cut}, where log Q_J is singular; {np.sum(~s)} pts): {s.sum()} pts, "
           f"max |d0| = {np.abs(d0[s]).max():.3f} (c2 units), rms {np.sqrt(np.mean(d0[s]**2)):.3f}; |d0|/|c2| on |c2| > 5 ({b.sum()} pts): "
           f"median {np.median(np.abs(d0[b])/np.abs(c2[b])):.3%}, max {np.max(np.abs(d0[b])/np.abs(c2[b])):.2%}")
    return out


if __name__ == '__main__' and len(sys.argv) == 1:
    pass
