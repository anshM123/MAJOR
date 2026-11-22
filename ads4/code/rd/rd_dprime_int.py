"""
rd_dprime_int (RD_CLOSE.md section 2.4): INTEGRATED far-region weight on the truncated real section D',
    W_kind(r, Delta) = int_{D', S_x >= r or S_y >= r} |Phi_kind| J  /  int_{D'} |Phi_kind| J ,   S = a(X1 + X2) (resp. Y),
by importance sampling (near: Gamma radial proposal; far: S - r ~ exponential with rate kappa Delta, kappa = d/dr[2 artanh(r/2)]).
|Phi| is nu-independent on D' (A > 0), so W bounds the far contribution to |I(nu)| up to the cancellation factor of the near region.
Fit log W = -c Delta + p log Delta + b over Delta = 30, 60, 100, 200 and compare c with 2 artanh(r/2).
Usage: python3 rd_dprime_int.py -> cache/rd_dprime_int.txt
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.stats import gamma as GD, beta as BD
import rd_core as C
from rd_dlib import logw, inside
np.seterr(all='ignore')

lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)


def near_total(Delta, kind, n, rng, sh=(4.0, 3.0, 4.0, 3.5), sc=1.6):
    if kind == 'EE':
        sh = (9.0, 5.0, 4.0, 3.5)
    X1, X2, Y1, Y2 = (GD.rvs(k, size=n, random_state=rng)*sc for k in sh)
    bs = BD.rvs(2.0, 2.0, size=n, random_state=rng); s = 0.97*(2*bs - 1)
    bp = BD.rvs(2.0, 2.0, size=n, random_state=rng); phi = 1.55*(2*bp - 1)
    ba = BD.rvs(0.7, 0.7, size=n, random_state=rng); al = np.pi*ba
    lq = (sum(GD.logpdf(v/sc, k) - np.log(sc) for v, k in zip((X1, X2, Y1, Y2), sh)) + BD.logpdf(bs, 2, 2) - np.log(2*0.97)
          + BD.logpdf(bp, 2, 2) - np.log(2*1.55) + BD.logpdf(ba, 0.7, 0.7) - np.log(np.pi))
    ok = inside(X1, X2, s, Y1, Y2, phi, Delta)
    L, _ = logw(X1[ok], X2[ok], s[ok], Y1[ok], Y2[ok], phi[ok], al[ok], Delta, kind)
    w = np.zeros(n); w[ok] = np.exp(L - lq[ok])
    bad = ~np.isfinite(w); NBAD[0] += int(bad.sum()); w[bad] = 0.0
    return w.mean(), w.std()/np.sqrt(n), L[np.isfinite(L)].max()


def far_part(Delta, kind, r, n, rng, side='x', sh=(4.0, 3.5), sc=1.6):
    kap = 1.0/(1 - r*r/4)
    t = rng.exponential(1.0/(kap*Delta), size=n)
    S = r + t
    f = BD.rvs(1.5, 1.5, size=n, random_state=rng)
    Ya, Yb = (GD.rvs(k, size=n, random_state=rng)*sc for k in sh)
    bs = BD.rvs(2.0, 2.0, size=n, random_state=rng); s = 0.97*(2*bs - 1)
    bp = BD.rvs(2.0, 2.0, size=n, random_state=rng); phi = 1.55*(2*bp - 1)
    ba = BD.rvs(0.7, 0.7, size=n, random_state=rng); al = np.pi*ba
    Xa = S*f*Delta; Xb = S*(1 - f)*Delta
    jac = np.log(S*Delta*Delta)                                   # dX1 dX2 = Delta^2 S dS df
    lq = (np.log(kap*Delta) - kap*Delta*t + BD.logpdf(f, 1.5, 1.5) + sum(GD.logpdf(v/sc, k) - np.log(sc) for v, k in zip((Ya, Yb), sh))
          + BD.logpdf(bs, 2, 2) - np.log(2*0.97) + BD.logpdf(bp, 2, 2) - np.log(2*1.55) + BD.logpdf(ba, 0.7, 0.7) - np.log(np.pi) - jac)
    if side == 'x':
        X1, X2, Y1, Y2 = Xa, Xb, Ya, Yb
    else:
        X1, X2, Y1, Y2 = Ya, Yb, Xa, Xb
    ok = inside(X1, X2, s, Y1, Y2, phi, Delta) & (S < 1)
    L, _ = logw(X1[ok], X2[ok], s[ok], Y1[ok], Y2[ok], phi[ok], al[ok], Delta, kind)
    w = np.zeros(n); w[ok] = np.exp(L - lq[ok])
    bad = ~np.isfinite(w); NBAD[0] += int(bad.sum()); w[bad] = 0.0
    return w.mean(), w.std()/np.sqrt(n)


NBAD = [0]
rng = np.random.default_rng(21)
RS = (0.2, 0.3, 0.5)
DS = (30.0, 60.0, 100.0, 200.0)
res = {}
pr('# rd_dprime_int: integrated far-region weight on D\' (relative to int_{D\'} |Phi| J); n = 2^18 per estimate')
t0 = time.time()
for kind in ('S', 'EE'):
    for Dl in DS:
        T, eT, _ = near_total(Dl, kind, 1 << 18, rng)
        row = f'  {kind:2s} Delta = {Dl:5g}: total {T:.4e} (+-{eT/T:.1%}); W(r) ='
        for r in RS:
            Fx, ex = far_part(Dl, kind, r, 1 << 18, rng, 'x')
            Fy, ey = far_part(Dl, kind, r, 1 << 18, rng, 'y')
            W = (Fx + Fy)/T; eW = np.hypot(ex, ey)/T
            res[(kind, Dl, r)] = (W, eW)
            row += f'  r={r}: {W:.3e} (+-{100*eW/max(W, 1e-300):.1f}%)'
        pr(row + f'   ({time.time() - t0:.0f}s; non-finite samples dropped so far: {NBAD[0]})')
pr('')
pr('# fit log W = -c Delta + p log Delta + b (points with relative error < 50 %); analytic c = 2 artanh(r/2)')
for kind in ('S', 'EE'):
    for r in RS:
        Ds = []; y = []
        for Dl in DS:
            W, eW = res[(kind, Dl, r)]
            if W > 0 and eW < 0.5*W:
                Ds.append(Dl); y.append(np.log(W))
        if len(Ds) >= 3:
            Ds = np.array(Ds); y = np.array(y)
            A = np.stack([-Ds, np.log(Ds), np.ones_like(Ds)], axis=1)
            coef = np.linalg.lstsq(A, y, rcond=None)[0]
            pr(f'   {kind:2s} r = {r}: c = {coef[0]:.4f} (analytic {2*np.arctanh(r/2):.4f}), p = {coef[1]:+.2f}; log W = {np.round(y, 1)} at Delta = {Ds}')
        else:
            pr(f'   {kind:2s} r = {r}: too few resolved points ({len(Ds)})')
open(os.path.join(HERE, 'cache', 'rd_dprime_int.txt'), 'w').write('\n'.join(lines) + '\n')
