"""
h2 (M3b): exact near-gap shift of the route-A' functional from the exact eta-inversion of psi_2 (h1 pickles), replacing the large-nu symbol
of f6.  Per state (all e1 states, both Df):  F_inv = F_base + dF,
   F_base = stripped B2 + B3 + psi_2 2(P - J0) + psi_3 u sigma (P - J0) + remainder (stripped envelope)          [= f6, 'strip', no inversion]
   dF     = int_0^1 dp (psit_2 - psi_2)(p) [ (2-y) D~_h(y) + 2 (P_{-1/2+i Lam p}(cosh beta) - J0((2J+1)x)) ]   (fine nu grid of h1)
            + sum_nodes (psit_2 - psi_2)(p_i) w_i env_2 fac rho_2(p_i)                                          (64 e1 nodes)
P_{-1/2+i nu}(cosh beta) at arbitrary nu by the Mehler-Dirichlet integral (vectorised; checked against remk.Pders_vec).
The tails of psit_2 beyond p = 1 (mass ~1e-9 of psi_2) are not included (their heavy action at nu > m is not given by the bulk-point formula) [I].
Usage: python3 h2_neargap_exact.py [Dg ...] -> cache/h2_neargap_exact_<Dgs>.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
from spf_common import NearGap, load_states
from common import gl01, Dh_one, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk

out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n); return 0.5*(b - a)*x + 0.5*(b + a), 0.5*(b - a)*w

def Pcone(nus, eta):
    """P_{-1/2+i nu}(eta), eta = cosh(beta) >= 1, Mehler-Dirichlet with t = beta(1 - s^2)."""
    beta = np.arccosh(max(eta, 1.0))
    if beta < 1e-12:
        return np.ones_like(nus)
    n = int(max(nus)*beta/2) + 200
    s, ws = gl(0.0, 1.0, n); tt = beta*(1 - s*s)
    den = 2*np.sinh(0.5*(beta + tt))*np.sinh(0.5*(beta - tt))
    kern = ws*2*beta*s/np.sqrt(den)
    res = np.empty_like(nus)
    ch = max(1, int(3e6/n))
    for i in range(0, len(nus), ch):
        res[i:i + ch] = (np.sqrt(2)/np.pi)*(np.cos(np.outer(nus[i:i + ch], tt)) @ kern)
    return res

def phit_at(d, nu):
    t, Phit = d['t'], d['Phit']; dt = t[1] - t[0]; wt = np.full_like(t, dt); wt[0] *= 0.5
    return (2/np.pi)*(np.cos(np.outer(nu, t)) @ (wt*Phit))

def main(Dgs):
    pr("# h2: exact eta-inversion shift of the route-A' near-gap action (psi_2 only; psi_3 has j = 0)")
    for Dg in Dgs:
        d = pickle.load(open(os.path.join(HERE, 'cache', 'h1_%d.pkl' % Dg), 'rb'))
        NG = NearGap(Dg); F = NG.F; pn, pw, L = NG.pn, NG.pw, NG.Lam
        assert abs(d['Lam'] - L) < 1e-9
        m = d['nus'] <= L
        nuf, wf = d['nus'][m], d['wns'][m]/L; pfine = nuf/L; dpf = (d['phit'] - d['phi'])[m]
        dpn = phit_at(d, L*pn) - F.psi2(pn)                      # at the 64 nodes (remainder part)
        # check Pcone
        chk = Pcone(np.array([0.5, 7.0, 30.0]), 1.05); ref = remk.Pders_vec(np.array([0.5, 7.0, 30.0]), 1.05, nmax=0)[0]
        pr(f"  Delta_g = {Dg}: Pcone check vs remk at eta = 1.05: max rel diff {np.max(np.abs(chk - ref)/np.abs(ref)):.1e}")
        recs = []
        for h, DF, r0, ref7 in load_states(Dg):
            mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]
            x = pn*np.sqrt(u); y = x*x
            P = remk.Pders_vec(L*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
            pf2, wf2 = gl01(max(400, 2*J + 100))
            Js, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
            p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
            cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
            fac = (1 - y)**2 if h == '++' else (1 - y)**-2
            rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
            Ns = abs(B2 + cv2) + abs(B3 + cv3)
            base = B2 + B3 + cv2 + cv3 + rem
            # exact shift on the fine grid
            xf = pfine*np.sqrt(u); yf = xf*xf
            Dh = Dh_one(h, J, yf); Dt = Dh if h == '++' else Dh/(1 - yf)**4
            Pf = Pcone(nuf, r0['eta'])
            T = (2 - yf)*Dt + 2*(Pf - j0((2*J + 1)*xf))
            dF = np.sum(wf*dpf*T) + np.sum(dpn*pw*envk(2, x)*fac*r0['rho'][2])
            recs.append(dict(h=h, Df=DF, mfl=mfl, J=J, Ns=Ns, base=base, dF=dF))
        grid = sorted(set(r['mfl'] for r in recs))
        for lab, Fv in (('no inversion', lambda r: r['base']), ('exact inversion', lambda r: r['base'] + r['dF'])):
            for delta in (0.0, 0.03):
                n = next((n0 for n0 in grid if all(Fv(r) - delta*r['Ns'] > 0 for r in recs if r['mfl'] >= n0 - 1e-9)), None)
                w = min(recs, key=lambda r: Fv(r)/r['Ns'])
                pr(f"    {lab:16s} delta={delta}: n = {n}; min F/N_s = {Fv(w)/w['Ns']:+.4f} [{w['h']} m/Lam={w['mfl']} J={w['J']} Df={w['Df']}]")
        wmax = max(recs, key=lambda r: abs(r['dF'])/r['Ns'])
        pr(f"    max |dF|/N_s (exact) = {abs(wmax['dF'])/wmax['Ns']:.3e} [{wmax['h']} m/Lam={wmax['mfl']} J={wmax['J']}]; "
           f"median {np.median([abs(r['dF'])/r['Ns'] for r in recs]):.2e}")
    open(os.path.join(HERE, 'cache', 'h2_neargap_exact_' + '_'.join(map(str, Dgs)) + '.txt'), 'w').write("\n".join(out) + "\n")

if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or [100, 300, 1000])
