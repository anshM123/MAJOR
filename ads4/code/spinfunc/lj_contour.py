"""
lj_contour (HLJ_ANALYTIC.md §2): explicit constants of the contour/endpoint lemma for the psi_2-weighted beyond-model correction
on the hybrid near-gap states (k = 2 part and conventions exactly as sm_bound.py / SUBLEADING_REGGE §0.1: Gamma in c2 units,
correction dF = sum p2 env_h Gamma / m^2, tolerance = 0.9 lam0 / weight, worse helicity).

Lemma (HLJ_ANALYTIC §2): split [0, 1] at p_B = nu_B/Lam, nu_B = c m/(2J+1) (the Bessel/Regge window: l theta <~ c).
If on p >= p_B the correction is Gamma = Gamma_+ + Gamma_- with Gamma_+- analytic in the half-strips +-Im p in [0, T] and
|Gamma_+-(p +- i y/Lam)| <= M e^{-beta_* y} (class AG), then
   |dF|/N_s <= W_B * max_{nu < nu_B} |Gamma| + E_bd * M + E_top * M + E_edge * M,
   W_B   = sum_{nu < nu_B} p2 env / (m^2 N_s)                                      (Bessel window, sign-free)
   E_bd  = 2 |psi2(p_B)| env(x_B) / (m^2 Lam beta_* N_s) * (1 + 1/(nu_B beta_*))     (vertical sides at p_B, both waves;
                                                                                    the 2nd term bounds the growth of psi2 ~ p up the side)
   E_top = 2 e^{-beta_* y_T} max_top |psi2 env| / (m^2 N_s)                           (top side at height y_T = 30/beta_*)
   E_edge= 2 * 5! |X(1)| env(x_max) / (m^2 phi'^6 N_s), phi' = Lam beta_*/sqrt(1-x_max^2) (p = 1 side; psi2 ~ X(1)(1-p)^5)
beta_* = 0.97 beta (lj_scan: min Re kappa~ >= 1.002 beta on the real axis, decay >= 0.97 beta on vertical lines up to y = 6 nu_r).
Tolerance on max|Gamma| (taking M = max|Gamma| in the window as the envelope at p_B): T_B = 0.9 lam0/(W_B + E_bd + E_top + E_edge).
Compared with SMOOTH_BOUND's uniform (all nu) and fold (nu < 30) tolerances (same code path as sm_bound).
Usage: python3 lj_contour.py Dg J1,J2,... [npn]  -> stdout, appended to cache/lj_contour.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import ee1000_common as E

CS = (3.0, 5.0, 8.0)
NUX = 25.0
_XC = {}


def _xcache(Dg):
    if Dg not in _XC:
        import pickle
        d = pickle.load(open(os.path.join(HERE, 'cache', 'u3_chi_EE_%d.pkl' % Dg), 'rb'))
        nu = np.asarray(d['nu']); X = np.asarray(d['X']); ph = np.asarray(d['phi'])
        _XC[Dg] = (nu, np.abs(X - ph))
    return _XC[Dg]
BSTAR = 0.97


def run(Dg, J, npn=1024):
    F = E.Q.NearGap(Dg).F
    X1 = float(np.polynomial.polynomial.polyval(1.0, F.X))
    Lam = E.lam_of(Dg); mfl = E.mgap(Dg, J); m = mfl*Lam
    beta = 2*np.arcsinh((J + .5)/m); bs = BSTAR*beta
    rows = []
    for h in ('++', '+-'):
        st = E.hyb(Dg, J, mfl, h, npn)
        x, Ns = st['x'], st['Ns']; nu = x*m
        env = E.Q.A.env(h, x, 2)
        w = st['p2']*env/m**2/Ns
        lam0 = st['base']/Ns
        r = dict(h=h, lam0=lam0, Ns=Ns, tol_all=0.9*lam0/w.sum(), tol30=0.9*lam0/w[nu < 30].sum())
        xm = 1/mfl
        envm = float(E.Q.A.env(h, np.array([xm*0.999999]), 2)[0])
        php = Lam*bs/np.sqrt(1 - xm**2)
        Eedge = 2*120*abs(X1)*envm/(m**2*php**6*Ns)
        yT = 30/bs
        pt = np.linspace(0, 1, 401) + 1j*yT/Lam
        xt = pt/mfl
        Etop = 2*np.exp(-bs*yT)*np.max(np.abs(F.psi2(pt)*E.Q.A.env(h, xt, 2)))/(m**2*Ns)
        w3 = np.abs(st['p3']*st['u']*E.Q.A.env(h, x, 3))/m**2/Ns          # C3H (psi_3) weight, as ee1000_common.dF_state
        # non-analytic part of the EE weight: the min-norm Abel solution X (u3_chi_EE cache) differs from psi_2 by O(1) at nu <~ 10,
        # 2-4 % at 12-18 and by the analytic local symbol beyond ~20; |X - psi2| on nu in [nu_B, NUX] is a real-axis (window) term.
        XC = _xcache(Dg)
        dX = np.interp(nu, XC[0], XC[1])                                     # |X - psi2| (p-measure, same normalisation as psi2)
        wX = dX*(np.asarray(st['p2'])/np.maximum(F.psi2(st['pn']), 1e-300))*env/m**2/Ns   # GL weight * |X - psi2| * env/(m^2 N_s)
        for c in CS:
            nuB = c*m/(2*J + 1); pB = nuB/Lam; xB = pB/mfl
            WB = w[nu < nuB].sum(); WB3 = w3[nu < nuB].sum()
            envB = float(E.Q.A.env(h, np.array([xB]), 2)[0]); envB3 = float(E.Q.A.env(h, np.array([xB]), 3)[0])
            Ebd = 2*abs(F.psi2(np.array([pB]))[0])*envB/(m**2*Lam*bs*Ns)*(1 + 1/(nuB*bs))
            Ebd3 = 2*abs(F.psi3(np.array([pB]))[0])*st['u']*envB3/(m**2*Lam*bs*Ns)*(1 + 1/(nuB*bs))
            WX = wX[(nu >= nuB) & (nu < NUX)].sum()
            r[c] = dict(nuB=nuB, WB=WB, Ebd=Ebd, WX=WX, TB=0.9*lam0/(WB + Ebd + Etop + Eedge), TBX=0.9*lam0/(WB + Ebd + WX + Etop + Eedge),
                        TB3=0.9*lam0/(WB3 + Ebd3 + Etop + Eedge), TBc=0.9*lam0/(WB + Ebd + WB3 + Ebd3 + 2*(Etop + Eedge)))
        r.update(Etop=Etop, Eedge=Eedge)
        rows.append(r)
    return dict(Dg=Dg, J=J, mfl=mfl, m=m, beta=beta, bhat=(2*J + 1)/mfl, rows=rows)


if __name__ == '__main__':
    Dg = int(sys.argv[1]); Js = [int(v) for v in sys.argv[2].split(',')]
    npn = int(sys.argv[3]) if len(sys.argv) > 3 else 1024
    f = open(os.path.join(HERE, 'cache', 'lj_contour.txt'), 'a')
    def pr(s):
        print(s); f.write(s + '\n'); f.flush()
    for J in Js:
        R = run(Dg, J, npn)
        pr(f"Dg={Dg} J={J} m/L={R['mfl']:.4f} bhat={R['bhat']:.1f} beta={R['beta']:.4f}")
        for r in R['rows']:
            s = (f"   {r['h']}: lam0 {r['lam0']:+.3f} | SMOOTH_BOUND tol: uniform {r['tol_all']:.2f}, nu<30 {r['tol30']:.0f} c2 | "
                 f"E_top {r['Etop']:.1e} E_edge {r['Eedge']:.1e} /c2 | ")
            s += ' ; '.join(f"c={c:g}: nu_B={r[c]['nuB']:.1f} W_B={r[c]['WB']:.2e} E_bd={r[c]['Ebd']:.2e} T_B(C2)={r[c]['TB']:.0f} "
                            f"T_B(C3H)={r[c]['TB3']:.0f} T_B(both)={r[c]['TBc']:.0f} W_X={r[c]['WX']:.2e} T_B(C2,+X)={r[c]['TBX']:.0f}" for c in CS)
            pr(s)
    f.close()
