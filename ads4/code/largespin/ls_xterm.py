"""
largespin/ls_xterm.py (LARGESPIN_ANALYTIC.md §3.4): the de-shift (min-norm Abel weight X) correction to the structure tolerances.

The k = 2 part of Omega^EE acts on RAW images with the weight X (prl/data/d_abel_X_1000.txt, read-only), not psi2; X o raw = psi2 o stripped
holds exactly for the leading-Regge raw image (de-shift identity, LRC*), hence by linearity for any correction proportional to it.
For a correction of unknown raw form we follow HLJ_ANALYTIC §2 and add the sign-free real-axis term
    cost_X(S) = 2 * sum_{nu_k < 30} |X(nu_k) - psi2(nu_k)| |S(nu_k)| dnu / Lam
(X - psi2 is <= 3e-3 relative beyond nu = 25).  Printed: a_tol for the structures F, G (smooth; X-term always applies) and, as a stress
test, for A, D, E with the X-term (applies only if the correction is NOT raw-proportional to the leading-Regge image).
Reads cache/struct_<Dg>.json (ls_struct.py).  Output: cache/xterm_<Dg>.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
import numpy as np
import ls_core as C

Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
L = C.lam(Dg)
d = np.loadtxt(os.path.join(C.HERE, '..', '..', 'prl', 'data', 'd_abel_X_%d.txt' % Dg), comments='#')
nu, X, ps = d[:, 0], d[:, 1] + 1j*d[:, 2], d[:, 5]
sel = nu < 30
nu, dX = nu[sel], np.abs(X[sel] - ps[sel])
dnu = 0.04
F = C.NFw(Dg)
chk = np.max(np.abs(ps[sel] - F.psi2(nu/L)))
rows = json.load(open(os.path.join(C.CACHE, 'struct_%d.json' % Dg)))
out = ['Delta_g = %d; psi2 column check vs nf_common: %.2e; sum |X - psi2| dnu (nu < 30) = %.4f' % (Dg, chk, (dX*dnu).sum())]


def Pnu(beta, nus, shift=False):
    ph, w = C.gl01(400); ph = ph*np.pi/2; w = w*np.pi/2
    s = 2*np.arcsinh(np.sinh(beta/2)*np.sin(ph))
    wt = w/np.cosh(s/2)*(np.cosh(s) if shift else 1.0)
    return (2/np.pi)*(np.cos(np.outer(nus, s))*wt).sum(1)


res = []
for r in rows:
    beta, m = r['beta'], r['m']
    fl = r['floor_rel']*r['Ns']
    P0 = Pnu(beta, nu); P1 = Pnu(beta, nu, True)
    S = dict(A=P0/m, D=P1/m, E=P0*0.25/(nu**2 + 0.25)/m, F=np.ones_like(nu)/m**2, G=np.exp(-(nu*beta/3)**2)/m**2)
    cx = {k: 2*(dX*np.abs(v)).sum()*dnu/L for k, v in S.items()}
    at = {}
    for k in S:
        base = fl/r['atol'][k]                         # psi2 cost per unit a
        at[k] = fl/(base + cx[k])
    res.append((r, at))
for reg in ('LJ', 'R'):
    rr = [(r, a) for r, a in res if r['reg'] == reg]
    out.append('\nregion %s:' % reg)
    for k in 'FGADE':
        v = np.array([a[k] for r, a in rr]); v0 = np.array([r['atol'][k] for r, a in rr])
        i = int(np.argmin(v))
        out.append('  a_tol[%s] with X-term: min %.4g (J=%d, m/Lam=%.3f, bhat=%.1f; without X-term there %.4g); max ratio (without/with) %.2f'
                   % (k, v[i], rr[i][0]['J'], rr[i][0]['mfl'], rr[i][0]['bhat'], v0[i], (v0/v).max()))
    out.append('  per state (J, m/Lam, bhat): a_tol[G] psi2-only -> with X;  a_tol[F] psi2-only -> with X')
    for r, a in rr:
        if r['bhat'] >= 40:
            out.append('    %5d %6.3f %8.1f : G %9.4g -> %9.4g ;  F %9.4g -> %9.4g ;  A %8.4g -> %8.4g' %
                       (r['J'], r['mfl'], r['bhat'], r['atol']['G'], a['G'], r['atol']['F'], a['F'], r['atol']['A'], a['A']))
txt = '\n'.join(out)
print(txt)
open(os.path.join(C.CACHE, 'xterm_%d.txt' % Dg), 'w').write(txt + '\n')
