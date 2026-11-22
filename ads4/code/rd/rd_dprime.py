"""
rd_dprime (RD_CLOSE.md section 2): the truncated real section
    D' = { D-point : 1 - s_x + q_x > 0, 1 - s_y + q_y > 0 (component of the origin) },   s = u + v, q = u v - p^2 (raw),
i.e. the real future-cone section cut at the zeros of B (m45 = 0, m35 = 0: the light cones of x5 = e from x4 and x3).

 (A) admissibility: on uniform samples of D' every factor of the scalar integrand with a non-integer exponent keeps one sign, A > 0,
     FP1 >= 0.5 (Bx_num + By_num) (so m12 = 0 is reached only at the corner where B^{Delta/2} vanishes), FP2 >= 1, FP3 > 0;
     the full real section D is NOT admissible: beyond the B-zero it crosses FP1 = 0 (m12 = 0), a pole-type branch divisor
     (exponent -11/2 scalar, -19/2 EE).
 (B) the B bound |B|^{Delta/2} <= exp(-2 Delta [artanh(s_x/2) + artanh(s_y/2)]) on D' (q <= s^2/4): checked on the samples.
 (C) far-region weight on D': for kind in (scalar, EE), M(r) = max over D' with max(s_x, s_y) >= r of log(|Phi_kind| J_D), minus the
     same max over all of D' (the near peak); multi-start L-BFGS-B in 7 real variables (|s| <= 0.97, where the double-precision EE
     contraction is reliable).  The exponent c_D'(r) = -[M(r) - M(0)]/Delta is fitted against Delta = 30, 100, 300, 1000 and compared
     with the analytic value 2 artanh(r/2) (+ O(log Delta/Delta)).
Usage: python3 rd_dprime.py -> cache/rd_dprime.txt (+ .pkl)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.optimize import minimize
import rd_core as C
np.seterr(all='ignore')

lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)


def raw_parts(Sx, fx, s, Sy, fy, phi, al, Delta):
    """D' coordinates: S = a(X1+X2) in (0,1), f = X1/(X1+X2) -> D variables."""
    X1 = Sx*fx*Delta; X2 = Sx*(1 - fx)*Delta; Y1 = Sy*fy*Delta; Y2 = Sy*(1 - fy)*Delta
    return X1, X2, Y1, Y2


def logw(X1, X2, s, Y1, Y2, phi, al, Delta, kind):
    """log(|Phi_kind| x D-Jacobian) at real D points (arrays); kind 'S' or 'EE'. Also returns geometry values."""
    z, y0 = C.raw_from_D(X1, X2, s, Y1, Y2, phi, al, Delta)
    a = 1.0/Delta
    g = C.geometry([z[i] for i in range(7)], y0val=y0.astype(complex), xsqval=(-a*a*X1*X2*(1 - s)*(1 + s)).astype(complex))
    F = C.Flow(0.0, Delta, np.ones(7))
    vals = {k: v for k, v in g.items() if not isinstance(v, list)}
    L = F.logabs_parts(vals)
    psi = al - phi
    L = L + np.log(np.sqrt(X1*X2)*np.sqrt(Y1*Y2)*np.cos(phi)**2*np.abs(np.sin(psi)))
    if kind == 'EE':
        L = L + np.log(np.abs(C.Qratio(g, Delta))) + np.log(a)
    return L, g


def inside(X1, X2, s, Y1, Y2, phi, Delta):
    a = 1.0/Delta
    ux, vx, uy, vy = a*X1, a*X2, a*Y1, a*Y2
    qx = ux*vx*(1 - s*s); qy = uy*vy*np.cos(phi)**2
    return (1 - (ux + vx) + qx > 0) & (1 - (uy + vy) + qy > 0) & (ux < 1) & (vx < 1) & (uy < 1) & (vy < 1)


# ---------------------------------------------------------------------------------------------------- (A), (B)
pr('# rd_dprime: the truncated real section D\' (cut at the B zeros m45 = 0, m35 = 0)')
pr('## (A) admissibility and (B) the B bound, uniform samples of D\' (raw light-cone components u, v in (0,1))')
rng = np.random.default_rng(11)
for Dl in (30.0, 300.0):
    n = 400000
    ux, vx, uy, vy = (rng.uniform(0, 1, n) for _ in range(4))
    s = rng.uniform(-0.999, 0.999, n); phi = rng.uniform(-1.55, 1.55, n); al = rng.uniform(0.01, np.pi - 0.01, n)
    X1, X2, Y1, Y2 = ux*Dl, vx*Dl, uy*Dl, vy*Dl
    m = inside(X1, X2, s, Y1, Y2, phi, Dl)
    for t in np.linspace(0.05, 0.95, 10):                        # component of the origin
        m &= inside(t*X1, t*X2, s, t*Y1, t*Y2, phi, Dl)
    X1, X2, Y1, Y2, s, phi, al = X1[m], X2[m], Y1[m], Y2[m], s[m], phi[m], al[m]
    z, y0 = C.raw_from_D(X1, X2, s, Y1, Y2, phi, al, Dl)
    a = 1.0/Dl
    g = C.geometry([z[i] for i in range(7)], y0val=y0.astype(complex), xsqval=(-a*a*X1*X2*(1 - s)*(1 + s)).astype(complex))
    F = C.Flow(1.0, Dl, np.ones(7))
    bad = []
    for k, e in F.fac:
        v = g[k]
        if np.max(np.abs(np.imag(v))) > 1e-9*np.max(np.abs(v)):
            bad.append(f'{k}: complex')
            continue
        if not float(np.real(e)).is_integer() or np.imag(e) != 0:
            sg = np.sign(np.real(v))
            if not (np.all(sg == sg[0])):
                bad.append(f'{k}: sign changes ({np.mean(sg > 0):.3f} positive)')
    A = np.real(g['Nn']**2*g['c13']/(g['m24']*g['m40']*g['m20']*g['c10']*g['c30']))
    Bxn = 1 - (z[1] + z[2]).real + (z[1]*z[2] - z[5]**2).real; Byn = 1 - (z[3] + z[4]).real + (z[3]*z[4] - z[6]**2).real
    sx = (z[1] + z[2]).real; sy = (z[3] + z[4]).real
    lB = 0.5*Dl*np.log(np.abs(np.real(np.exp(F.logB({k: v for k, v in g.items() if not isinstance(v, list)})))))
    bound = -2*Dl*(np.arctanh(sx/2) + np.arctanh(sy/2))
    Q = C.Qratio(g, Dl)
    pr(f'  Delta = {Dl:g}: {m.sum()} samples in D\'; non-integer-power factors with a sign change or complex value: {bad if bad else "none"}; '
       f'min A = {A.min():.3e} (A > 0: {np.all(A > 0)})')
    pr(f'      min FP1/(Bx_num + By_num) = {(np.real(g["FP1"])/(Bxn + Byn)).min():.3f};  min FP2 = {np.real(g["FP2"]).min():.3f};  '
       f'min FP3 = {np.real(g["FP3"]).min():.3f};  min |m15| = {np.abs(g["m15"]).min():.3f}, min |m25| = {np.abs(g["m25"]).min():.3f}')
    pr(f'      B bound: max[(Delta/2) log|B| - (-2 Delta (artanh(s_x/2) + artanh(s_y/2)))] = {np.max(lB - bound):.2e} (<= 0 required); '
       f'Q real: max|Im Q/Q| = {np.max(np.abs(Q.imag/Q.real)):.1e}; Q sign changes on D\': {not np.all(np.sign(Q.real) == np.sign(Q.real[0]))}')

# ---------------------------------------------------------------------------------------------------- (C)
pr('')
pr('## (C) far-region weight on D\': c_D\'(r) = -[max_{D\', max(s_x,s_y) >= r} log(|Phi| J) - max_{D\'} log(|Phi| J)]/Delta')
RS = (0.1, 0.2, 0.3, 0.5, 0.7)
SMAXABS = 0.97


def unpack(p, Delta, mode, r):
    """p (7,) unconstrained -> D' point with the escaping sum S >= r (mode 'x' or 'y') or unconstrained ('near')."""
    sig = lambda t: 1/(1 + np.exp(-t))
    if mode == 'near':
        Sx = 0.98*sig(p[0]); Sy = 0.98*sig(p[3])
    elif mode == 'x':
        Sx = r + (0.999 - r)*sig(p[0]); Sy = 0.98*sig(p[3])
    else:
        Sx = 0.98*sig(p[0]); Sy = r + (0.999 - r)*sig(p[3])
    fx = sig(p[1]); fy = sig(p[4])
    s = SMAXABS*np.tanh(p[2]); phi = 1.55*np.tanh(p[5]); al = np.pi*sig(p[6])
    X1, X2, Y1, Y2 = Sx*fx*Delta, Sx*(1 - fx)*Delta, Sy*fy*Delta, Sy*(1 - fy)*Delta
    return X1, X2, s, Y1, Y2, phi, al


def objective(p, Delta, kind, mode, r):
    X1, X2, s, Y1, Y2, phi, al = unpack(p, Delta, mode, r)
    args = [np.atleast_1d(v) for v in (X1, X2, s, Y1, Y2, phi, al)]
    if not inside(*args[:6], Delta)[0]:
        return 1e6
    L, _ = logw(*args, Delta, kind)
    v = -float(L[0])
    return v if np.isfinite(v) else 1e6


def maximise(Delta, kind, mode, r, nbatch=60000, nloc=6, seed=0):
    """vectorised random search over the unconstrained 7-vector, then Nelder-Mead from the best nloc points."""
    rg = np.random.default_rng(seed)
    P = rg.normal(size=(nbatch, 7))*1.5
    if mode == 'near':
        P[:, 0] = rg.uniform(-9, -1, nbatch); P[:, 3] = rg.uniform(-9, -1, nbatch)
    else:
        P[:, 0 if mode == 'x' else 3] = rg.uniform(-8, 3, nbatch)
        P[:, 3 if mode == 'x' else 0] = rg.uniform(-9, -1, nbatch)
    X1, X2, s, Y1, Y2, phi, al = unpack(P.T, Delta, mode, r)
    ok = inside(X1, X2, s, Y1, Y2, phi, Delta)
    L = np.full(nbatch, -np.inf)
    Lv, _ = logw(X1[ok], X2[ok], s[ok], Y1[ok], Y2[ok], phi[ok], al[ok], Delta, kind)
    L[ok] = np.where(np.isfinite(Lv), Lv, -np.inf)
    best = (np.inf, None)
    for j in np.argsort(-L)[:nloc]:
        res = minimize(objective, P[j], args=(Delta, kind, mode, r), method='Nelder-Mead',
                       options=dict(maxiter=3000, xatol=1e-7, fatol=1e-10))
        if res.fun < best[0]:
            best = (res.fun, res.x)
    return -best[0], best[1]


out = {}
t0 = time.time()
for kind in ('S', 'EE'):
    for Dl in (30.0, 100.0, 300.0, 1000.0):
        M0, p0 = maximise(Dl, kind, 'near', 0.0)
        row = [f'  {kind:2s} Delta = {Dl:6g}: near peak log(|Phi|J) = {M0:9.2f};  c(r) for r =']
        for r in RS:
            Mx, _ = maximise(Dl, kind, 'x', r, seed=1)
            My, _ = maximise(Dl, kind, 'y', r, seed=2)
            Mf = max(Mx, My)
            out[(kind, Dl, r)] = (M0, Mf, Mx, My)
            row.append(f' {r}: {-(Mf - M0)/Dl:.4f}')
        pr(''.join(row) + f'   ({time.time() - t0:.0f}s)')
pr('')
pr('  analytic (Delta -> oo):        ' + ''.join(f' {r}: {2*np.arctanh(r/2):.4f}' for r in RS) + '   [2 artanh(r/2): u = v = r/2, w = 1 maximises B on the shell]')
pr('  fit  -(M(r) - M(0)) = c_inf(r) Delta - p log Delta + b  over Delta = 30 ... 1000:')
for kind in ('S', 'EE'):
    for r in RS:
        Ds = np.array([30.0, 100.0, 300.0, 1000.0])
        y = np.array([-(out[(kind, D, r)][1] - out[(kind, D, r)][0]) for D in Ds])
        Amat = np.stack([Ds, -np.log(Ds), np.ones_like(Ds)], axis=1)
        coef, *_ = np.linalg.lstsq(Amat, y, rcond=None)
        resid = y - Amat @ coef
        pr(f'    {kind:2s} r = {r}: c_inf = {coef[0]:.4f} (analytic {2*np.arctanh(r/2):.4f}), p = {coef[1]:+.2f}, b = {coef[2]:+.2f}; '
           f'max residual {np.abs(resid).max():.2f}; drop at Delta = 300 / 1000: {y[2]:.1f} / {y[3]:.1f}')
pickle.dump(out, open(os.path.join(HERE, 'cache', 'rd_dprime.pkl'), 'wb'))
open(os.path.join(HERE, 'cache', 'rd_dprime.txt'), 'w').write('\n'.join(lines) + '\n')
