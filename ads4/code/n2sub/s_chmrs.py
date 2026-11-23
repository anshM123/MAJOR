"""
s_chmrs.py (N2_SUBREGGE.md §3.2): CHMRS (2106.10274) eq. (regge subleading) -- the first 1/m^2 correction of the scalar C_{k,nu} at fixed
beta -- checked against the project's EXACT scalar skeleton (spinfunc/cache/sm_exact_<Dg>.pkl, read-only: ex[2] = m^4 C_{2,nu}, Df = 3).

   (m^{2k}/2) C_{k,nu} = P(eta) + delta(nu; eta)/m^2 + O(1/m^4),
   delta = (1/6){ P'(eta) [(eta^2-1)(nu^2 - 3a^2 + 6a(2Dphi-d) - 2) - (d-4)(d-2)]
                 + eta P(eta) [(nu^2 - a(a+d-2))(7d/2 + 3a - 1 - 6Dphi) - 2(Dphi-d-1)_3 - 2(a+(d-2)/2)(1 - a^2 - a(d-2)/2)] },
   a = 2 Dphi + k - d/2,  P = P_{(2-d)/2+i nu},  m^2 = (Delta-J-d+1)(Delta+J-1)  (= the project's m^2 in d = 3).
Prints, per state, on the nodes nu < 40: exact m^2 (ex2/2 - P)  vs  delta,  and the residual (expected O(nu^4/m^2)).
Usage: python3 s_chmrs.py -> cache/s_chmrs.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
import numpy as np
import n2_common as N

d, k = 3.0, 2.0


def poch(x, n):
    out = 1.0
    for i in range(n):
        out *= (x + i)
    return out


def delta_chmrs(nu, beta, Dphi=3.0, kk=2.0):
    a = 2*Dphi + kk - d/2
    eta = np.cosh(beta)
    P = N.Pnu(beta, nu).real
    h = 1e-5
    Pp = ((N.Pnu(beta + h, nu) - N.Pnu(beta - h, nu)).real/(2*h))/np.sinh(beta)       # dP/deta
    c1 = (eta**2 - 1)*(nu**2 - 3*a*a + 6*a*(2*Dphi - d) - 2) - (d - 4)*(d - 2)
    c2 = (nu**2 - a*(a + d - 2))*(3.5*d + 3*a - 1 - 6*Dphi) - 2*poch(Dphi - d - 1, 3) - 2*(a + (d - 2)/2)*(1 - a*a - a*(d - 2)/2)
    return (Pp*c1 + eta*P*c2)/6.0, P, Pp, c1, c2


if __name__ == '__main__':
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    for Dg in (1000, 300):
        fn = os.path.join(N.ADS4, 'code', 'spinfunc', 'cache', 'sm_exact_%d.pkl' % Dg)
        if not os.path.exists(fn):
            continue
        db = pickle.load(open(fn, 'rb'))
        for key in sorted(db, key=lambda q: (q[0], len(q))):
            r = db[key]
            if 2 not in r['ex']:
                continue
            st = r['st']; J = r['J']; m = st['m']; beta = st['beta']; Df = r.get('Df', 3.0)
            nu = np.asarray(r['nus']); ex2 = np.asarray(r['ex'][2])
            sel = nu < 40
            if not np.all(np.isfinite(ex2[sel])):
                pr(f"Dg={Dg} {key}: non-finite ex2, skipped"); continue
            dl, P, Pp, c1, c2 = delta_chmrs(nu[sel], beta, Df)
            dex = m*m*(ex2[sel]/2 - P)
            res = dex - dl
            pr(f"Dg={Dg} key={key} J={J} Delta={st['Delta']:g} m={m:.1f} beta={beta:.4f} Dphi={Df}: c2(nu=0) = {c2[0]:.2f}, "
               f"(eta^2-1) = {np.cosh(beta)**2 - 1:.4f}")
            idx = np.unique(np.searchsorted(nu[sel], [0.0, 2, 5, 10, 16, 25, 35]).clip(0, sel.sum() - 1))
            for i in idx:
                pr(f"    nu={nu[sel][i]:6.2f}: exact m^2(ex/2-P) = {dex[i]:+10.3f}   CHMRS delta = {dl[i]:+10.3f}   resid = {res[i]:+9.3f}   "
                   f"resid*m^2/(1+nu^2)^2 = {res[i]*m*m/(1 + nu[sel][i]**2)**2:+.3g}")
            w = nu[sel] < 16
            cP = float((res[w]*P[w]).sum()/(P[w]*P[w]).sum())
            r2 = res - cP*P
            pr(f"    fit resid = c P on nu<16: c = {cP:+.3f};  after it: max|resid| nu<16 {np.max(np.abs(r2[w])):.3f}, nu<40 {np.max(np.abs(r2)):.3f}")
            pr(f"    max|exact| = {np.max(np.abs(dex)):.2f};  max|resid| (nu<16) = {np.max(np.abs(res[nu[sel] < 16])):.3f};  (nu<40) = {np.max(np.abs(res)):.3f}")
    open(os.path.join(N.CACHE, 's_chmrs.txt'), 'w').write('\n'.join(lines) + '\n')
