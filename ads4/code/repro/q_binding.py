"""repro/q_binding.py -- Result Q's binding state (++ J = 0 at the gap, Dg = 300, 1000): FF/N of the delta = 0.2 functional.

What is recomputed here with own code:
  (1) flat Stage-0 action of Q3-full-delta0.2 on a grid of flat states: min FF/N (should be delta = 0.2) and its location;
  (2) the HYB (flat exact + leading-Regge AdS curvature) action FF, N at ++ J = 0, u = 1 (Legendre functions P_{-1/2+i nu}(cosh b)
      by own hypergeometric series, D = -d_{nu^2} P|_0 by own series, J0 by scipy);
  (3) s_d = 0.04 sum_nu chi_d(nu) Re Gamma_K(nu) and s_v from the ENGINE RATIO rho(nu) = I_st/I_ref (read from the production cache
      spinfunc/cache/bg_qan_D{300,1000}_J0.pkl, K = 4): own T(nu) = sum_j R_j P_{-1/2+i(nu-ij)}(cosh beta) for state and reference,
      own Gamma_exact = m^2 (c T_ref rho - T_st), Gamma_K = Gamma_exact - 0.75 sum_j R_j (nu-ij)^2 - const N^_4.
      R_j: the exact rational functions stored (as sympy srepr strings) in resultq/cache/stage4_hj.pkl ['EE3'] -- read as data.
      chi_d, chi_v (H0 min-norm gauge, dnu = 0.04): resultq/cache/stage6_gauge.pkl -- read as data.
  (4) FF/N = (FF - 2 a u^3 s_d - c (2u^3/Lam^2) s_v)/N.
Definitions (Gamma_K, the HYB model, the s_d/FF coupling) were READ from spinfunc/bg_qan.py, bf_qan.py and resultq/stage6_model.py
(nothing imported).  NOT reproduced: the engine ratio rho(nu) itself (BELOW_FOLD2's 7-variable analytic-radial 1/Delta engine), and the
de-shift weights chi (Stage 4/6/7 Abel solutions).
Usage: python3 q_binding.py  -> stdout, cache/q_binding.txt
"""
import json, os, pickle
import numpy as np
import sympy as sp
from scipy.special import j0, roots_legendre

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
out = []


def pr(s=''):
    print(s, flush=True); out.append(s)


Q = json.load(open(os.path.join(ROOT, 'resultq', 'functionals_q3.json')))['Q3-full-delta0.2']
c2 = np.array(Q['psi2']); c3 = np.array(Q['psi3']); a, c = Q['a'], Q['c']; delta = Q['delta']
psi2 = lambda p: np.polynomial.polynomial.polyval(p, c2)
psi3 = lambda p: np.polynomial.polynomial.polyval(p, c3)


def gl01(n):
    x, w = roots_legendre(n); return (x + 1)/2, w/2


def jac_rows(nmax, al, be, x):
    """P^{(al,be)}_n(x), n = 0..nmax (standard recurrence; own code)."""
    P = np.zeros((nmax + 1,) + x.shape); P[0] = 1
    if nmax >= 1:
        P[1] = (al - be)/2 + (al + be + 2)*x/2
    for n in range(2, nmax + 1):
        A = 2*n*(n + al + be)*(2*n + al + be - 2)
        B = (2*n + al + be - 1)*((2*n + al + be)*(2*n + al + be - 2)*x + al*al - be*be)
        Cc = 2*(n + al - 1)*(n + be - 1)*(2*n + al + be)
        P[n] = (B*P[n - 1] - Cc*P[n - 2])/A
    return P


def flat_scan():
    """flat Stage-0 action: FF = int psi2 (2-y) D + s u int psi3 D + a (2JJ - c4 + s) u^3 - 2 c u^2, y = u p^2,
    D_{++} = P_J(1-2y), D_{+-} = P^{(0,8)}_{J-4}(1-2y);  N = |B2| + |B3| + |FW| + 2 c u^2."""
    best = []
    us0 = np.unique(np.concatenate([np.logspace(-4, 0, 400), np.linspace(0.9, 1, 101)]))
    for h, s, c4, al, be, J0 in (('++', 1, 1, 0, 0, 0), ('+-', -1, 41, 0, 8, 4)):
      mn = (9, None)
      for us in (us0, 'refine'):
        if isinstance(us, str):
            uc = mn[1][0]; us = np.linspace(uc*0.9, min(1.0, uc*1.1), 201)
        for u in us:
            Jmax = int(min(600, max(12, 60/np.sqrt(u))))
            p, w = gl01(Jmax + 40); y = u*p*p
            D = jac_rows(Jmax - J0, al, be, 1 - 2*y)
            B2 = D @ (w*psi2(p)*(2 - y)); B3 = D @ (w*s*u*psi3(p))
            J = np.arange(J0, Jmax + 1); JJ = J*(J + 1.0)
            FW = a*(2*JJ - c4 + s)*u**3; CB = 2*c*u*u
            FF = B2 + B3 + FW - CB; N = np.abs(B2) + np.abs(B3) + np.abs(FW) + CB
            r = FF/N
            if h == '++':
                r = np.where(J % 2 == 0, r, np.inf)
            i = int(np.argmin(r))
            if r[i] < mn[0]:
                mn = (r[i], (u, int(J[i])))
      best.append((h, mn))
    return best


def legP_series(mu, eta, nmax=200):
    """P_{-1/2+i mu}(eta) = 2F1(1/2+i mu, 1/2-i mu; 1; (1-eta)/2), mu complex array; own series (|1-eta| << 1 here)."""
    mu = np.asarray(mu, complex); w = (1 - eta)/2
    term = np.ones_like(mu); tot = term.copy()
    for n in range(1, nmax):
        term = term*((n - 0.5)**2 + mu*mu)*w/(n*n)
        tot = tot + term
        if np.max(np.abs(term)) < 1e-18*np.max(np.abs(tot)):
            break
    return tot


def Dder(eta, nmax=200):
    """D = -d/d(nu^2) P_{-1/2+i nu}(eta) at nu = 0 (own series)."""
    w = (1 - eta)/2; tot = 0.0; prod = 1.0; harm = 0.0
    for n in range(1, nmax):
        prod *= (n - 0.5)**2/(n*n); harm += 1/(n - 0.5)**2
        tot += prod*harm*w**n
    return -tot


def hyb_J0(Dg):
    """HYB action of ++ J = 0 at the gap (u = 1, bhat = 1)."""
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); u = 1.0; bh = 1.0
    eta = 1 + bh*bh/(2*Lam*Lam)
    p, w = gl01(200)
    P = legP_series(Lam*p, eta).real; Jf = j0(bh*p)
    dI20 = np.sum(w*psi2(p)*(P - Jf)); dI22 = np.sum(w*psi2(p)*p*p*(P - Jf)); dI3 = np.sum(w*psi3(p)*(P - Jf))
    P0 = legP_series(np.array([0.0]), eta).real[0]; LD = Lam*Lam*Dder(eta)
    dP0 = P0 - 1; dLD = LD - bh*bh/4
    B2f = np.sum(w*psi2(p)*(2 - p*p)); B3f = np.sum(w*psi3(p))
    B2 = B2f + 2*dI20 - u*dI22; B3 = B3f + u*dI3
    B4 = 2*u*u*P0; dB4 = (1 - 0)*u**3 - 2*u*u*dLD + 1*u**3*dP0; B5 = u**3*P0
    FW = -a*(dB4 - B5); FF = B2 + B3 + FW - c*B4; N = abs(B2) + abs(B3) + abs(FW) + c*abs(B4)
    return dict(Lam=Lam, FF=FF, N=N, B2=B2, B3=B3, FW=FW, B4=B4, flat=(B2f + B3f - 1)/(abs(B2f) + abs(B3f) + 1),
                dcurv=(dI20, dI22, dI3, dP0, dLD))


def Rj_funcs():
    D = pickle.load(open(os.path.join(ROOT, 'resultq', 'cache', 'stage4_hj.pkl'), 'rb'))['EE3']
    nu = sp.Symbol('nu')
    return {j: sp.lambdify(nu, sp.sympify(v), 'numpy') for j, v in D.items()}


def sd_from_ratio(Dg, Rj, W):
    db = pickle.load(open(os.path.join(ROOT, 'spinfunc', 'cache', 'bg_qan_D%d_J0.pkl' % Dg), 'rb'))
    rec = db['rs'][4]
    nus = rec['nus']; rho = rec['ratio']
    def mass(D, J):
        al = D - J - 2.0; m2 = al*(al + 2*J + 1); return np.sqrt(m2), 1 + (2*J + 1)**2/(2*m2)
    m, eta = mass(float(Dg), 0); mr, etar = mass(2e4, 0)
    R = {j: np.array(f(nus.astype(complex)), complex)*np.ones(len(nus)) for j, f in Rj.items()}
    T = lambda e: sum(R[j]*legP_series(nus - 1j*j, e) for j in R)
    Tst, Tref = T(eta), T(etar)
    NH = sum(R.values())
    tot = Tref*rho; cc = Tst[0]/tot[0]
    Gex = m*m*(cc*tot - Tst)
    GK = Gex - 0.75*sum(R[j]*(nus - 1j*j)**2 for j in R)
    GK = GK - (GK[0]/NH[0])*NH
    nuw, chiv, chid = W
    k = min(len(nus), len(chid)); assert np.allclose(nus[:k], nuw[:k])
    sel = nus[:k] < 41
    sd = 0.04*np.sum(chid[:k][sel]*GK[:k][sel].real); sv = 0.04*np.sum(chiv[:k][sel]*GK[:k][sel].real)
    # window decomposition of s_d and weight shares
    sabs = sum(np.abs(R[j]) for j in R)
    wts = np.abs(chid[:k])*sabs[:k]; win = [0, 1, 2, 3, 4, 8, 41]
    parts = [(lo, hi, 0.04*np.sum((chid[:k]*GK[:k].real)[(nus[:k] >= lo) & (nus[:k] < hi)])) for lo, hi in zip(win[:-1], win[1:])]
    share = [(lo, hi, wts[(nus[:k] >= lo) & (nus[:k] < hi)].sum()/wts[sel].sum()) for lo, hi in zip([0, 1, 2, 3], [1, 2, 3, 4])]
    cmp = dict(Tst=np.max(np.abs(Tst - rec['Tst'])/np.abs(rec['Tst'])), Tref=np.max(np.abs(Tref - rec['Tref'])/np.abs(rec['Tref'])),
               GK=np.max(np.abs(GK - rec['GK'])/sabs), GKrel=np.max(np.abs(GK - rec['GK'])/np.maximum(np.abs(rec['GK']), 1e-300)),
               NH=np.max(np.abs(NH - rec['NH'])/np.abs(rec['NH'])))
    # uniform per-harmonic error dG (worst sign) on a window that uses up s_fail - s_d:  0.04 sum |chi_d| sum_j|R_j| dG
    tolw = [(lo, hi, 0.04*wts[(nus[:k] >= lo) & (nus[:k] < hi)].sum()) for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 4))]
    return dict(sd=sd, sv=sv, parts=parts, share=share, cmp=cmp, tolw=tolw, stored=db['P'][4]['H0'], m=m, GK=GK, nus=nus)


if __name__ == '__main__':
    pr('# repro/q_binding.py  (Q3-full-delta0.2: a = %.6f, c = %.2f, delta = %.2f)' % (a, c, delta))
    fs = flat_scan()
    for h, (r, at) in fs:
        pr(f'(1) flat scan, {h}: min FF/N = {r:.6f} at (u, J) = ({at[0]:.4f}, {at[1]})   [LP margin delta = 0.2]')
    Rj = Rj_funcs()
    G = pickle.load(open(os.path.join(ROOT, 'resultq', 'cache', 'stage6_gauge.pkl'), 'rb'))
    W = (G[('C4', 0.04, 'nu')], G[('C4v', 0.04)]['H0 min-norm']['chi'], G[('C4d', 0.04)]['H0 min-norm']['chi'])
    claimed = {300: 0.1427, 1000: 0.1428}
    sdK5 = {300: 0.5096, 1000: 0.5077}
    for Dg in (300, 1000):
        H = hyb_J0(Dg)
        pr(f'\nDg = {Dg} (Lam = {H["Lam"]:.3f}):')
        pr(f'(2) HYB ++ J=0 gap: FF = {H["FF"]:.6f}, N = {H["N"]:.6f}, FF/N (s_d = 0) = {H["FF"]/H["N"]:.6f}; flat value {H["flat"]:.6f}; '
           f'curvature pieces dI20, dI22, dI3, dP0, dLD = ' + ', '.join(f'{x:+.2e}' for x in H['dcurv']))
        slope = 2*a/H['N']
        pr(f'    slope dFF/N / ds_d = -2a/N = {-slope:.4f} [-0.1126];  s_fail = FF/(2a) = {H["FF"]/(2*a):.4f} [1.78]')
        S = sd_from_ratio(Dg, Rj, W)
        pr(f'(3) own T, Gamma_K, s_d from the cached engine ratio (K = 4): s_d = {S["sd"]:+.5f}, s_v = {S["sv"]:+.4f};  '
           f'stored production (bg_qan, K = 4, H0): s_d = {S["stored"][0]:+.5f}, s_v = {S["stored"][1]:+.4f}')
        pr(f'    max rel diff vs stored: T_st {S["cmp"]["Tst"]:.1e}, T_ref {S["cmp"]["Tref"]:.1e}, N^_4 {S["cmp"]["NH"]:.1e}; '
           f'max |Gamma_K - stored| = {S["cmp"]["GK"]:.1e} c2 per harmonic')
        pr('    per-harmonic tolerance (uniform worst-sign |dG| that uses up s_fail - s_d, delta = 0.2): ' +
           ', '.join(f'[{lo},{hi}): {(H["FF"]/(2*a) - S["sd"])/v:.1f}' for lo, hi, v in S['tolw']) + '   [30.2, 16.9, 80.2, 65.2]')
        pr('    s_d by window: ' + ', '.join(f'[{lo},{hi}): {v:+.3f}' for lo, hi, v in S['parts']) +
           '   [HRS: +0.10, +0.57, -0.05, -0.10, -0.01 (4-8), ~0]')
        pr('    weight share |chi_d| sum|R_j|: ' + ', '.join(f'[{lo},{hi}): {100*v:.1f} %' for lo, hi, v in S['share']) +
           '   [24.5, 43.9, 9.2, 11.4 %]')
        for lab, sd, sv in (('K=4 own', S['sd'], S['sv']), ('K=5 quoted s_d (HRS_CLOSE)', sdK5[Dg], S['sv'])):
            ffn = (H['FF'] - 2*a*sd - c*(2/H['Lam']**2)*sv)/H['N']
            pr(f'(4) FF/N with {lab}: s_d = {sd:+.5f}  ->  FF/N = {ffn:+.5f}   [claimed {claimed[Dg]:+.4f}]')
    open(os.path.join(HERE, 'cache', 'q_binding.txt'), 'w').write('\n'.join(out) + '\n')
