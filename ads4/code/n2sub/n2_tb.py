"""
n2_tb.py (N2_SUBREGGE.md §2): HLJ_ANALYTIC §2's window tolerance T_B recomputed with the weight that actually multiplies the rows' Gamma.

T_B = 0.9 lam0 / (W_win + E_bd + W_X[nu_B, 25] + E_top + E_edge)  (lj_contour.py conventions; floor 0.9 lam0, worse of the two numbers
per helicity is NOT taken here -- both helicities printed).  Window weight variants (all on the 0.04 grid of d_abel_X_1000.txt):
   psi  : W = sum_{nu<nu_B} psi2 env dnu/(Lam m^2 N_s)           (HLJ, as published; tolerance on sup|Gamma - kappa|)
   |X|  : W = sum_{nu<nu_B} |X| env ...                          (exact weight, sign-free; tolerance on sup|Gamma_c| -- Gamma_c is the
                                                                   forward-normalised raw correction, singular at nu = 1.171, 6.643)
   |Y|  : W = sum_{nu<nu_B} |Y| env ...                          (tolerance on sup|a Gamma_c|, a = (R_-+R_+)/n, |a| <= 1: the regular,
                                                                   n-normalised raw correction; X Gamma_c = Y (a Gamma_c))
E_bd, E_top, E_edge, W_X exactly as lj_contour: the contour on nu >= nu_B deforms the analytic local weight psi2, and X - psi2 there is the
real-axis term W_X (HLJ_ANALYTIC sec. 2); only the window weight changes.
Also: the window-bump structure G = a e^{-(nu beta/3)^2}/m^2 (LARGESPIN §3.3) with |X| on nu < 30 instead of psi2 + |X - psi2|,
and the ring-amplitude class c0 D_h (HWIN §4.2) under the two readings (R1: raw = A[env c0 D_h]: exact psi2 action; R2: raw = env c0 D_h).
Usage: python3 n2_tb.py  -> cache/n2_tb.txt, cache/n2_tb.json
"""
import os, sys, json
sys.dont_write_bytecode = True
import numpy as np
import n2_common as N
sys.path.insert(0, os.path.join(N.ADS4, 'code', 'spinfunc'))
import ee1000_common as E          # read-only: hybrid state pieces exactly as lj_contour

Dg = 1000
L = N.lam(Dg)
XD = N.xdata(Dg)
nuG, XG, YG, psG, dn = XD['nu'], XD['X'], XD['Y'], XD['psi2'], XD['dnu']
NUX = 25.0
STATES = [(J, 1.0) for J in (48, 50, 64, 100, 150, 200, 250, 300, 344)] + [(100, 1.1), (200, 1.1), (300, 1.025)]
CS = (3.0, 5.0)


def run(J, mf):
    F = E.Q.NearGap(Dg).F
    X1 = float(np.polynomial.polynomial.polyval(1.0, F.X))
    mfl = E.mgap(Dg, J)*mf; m = mfl*L
    if mfl > 1.3 + 1e-9:
        return None
    beta = 2*np.arcsinh((J + .5)/m); bs = 0.97*beta
    res = dict(J=J, mfl=mfl, m=m, beta=beta, bhat=(2*J + 1)/mfl, h={})
    for h in ('++', '+-'):
        st = E.hyb(Dg, J, mfl, h, 1024)
        Ns = st['Ns']; lam0 = st['base']/Ns
        x = st['x']; nu = x*m
        envGL = E.Q.A.env(h, x, 2)
        wGL = st['p2']*envGL/m**2/Ns
        xm = 1/mfl
        envm = float(E.Q.A.env(h, np.array([xm*0.999999]), 2)[0])
        php = L*bs/np.sqrt(1 - xm**2)
        Eedge = 2*120*abs(X1)*envm/(m**2*php**6*Ns)
        yT = 30/bs
        pt = np.linspace(0, 1, 401) + 1j*yT/L
        Etop = 2*np.exp(-bs*yT)*np.max(np.abs(F.psi2(pt)*E.Q.A.env(h, pt/mfl, 2)))/(m**2*Ns)
        envg = N.env(h, nuG/m, 2).real
        cg = envg*dn/(L*m**2*Ns)                          # grid measure
        r = dict(lam0=lam0, Ns=Ns)
        for c in CS:
            nuB = c*m/(2*J + 1); pB = nuB/L; xB = pB/mfl
            envB = float(E.Q.A.env(h, np.array([xB]), 2)[0])
            Ebd = 2*abs(F.psi2(np.array([pB]))[0])*envB/(m**2*L*bs*Ns)*(1 + 1/(nuB*bs))
            rXB = np.interp(nuB, nuG, np.abs(XG))/np.interp(nuB, nuG, psG)
            win = nuG < nuB
            xr = (nuG >= nuB) & (nuG < NUX)
            WX = (np.abs(XG - psG)*cg)[xr].sum()
            Wps = (psG*cg)[win].sum(); WaX = (np.abs(XG)*cg)[win].sum(); WaY = (np.abs(YG)*cg)[win].sum()
            WpsGL = wGL[nu < nuB].sum()
            den0 = Ebd + WX + Etop + Eedge
            r[c] = dict(nuB=nuB, WpsGL=WpsGL, Wps=Wps, WaX=WaX, WaY=WaY, Ebd=Ebd, WX=WX,
                        TB_psi=0.9*lam0/(WpsGL + Ebd + Etop + Eedge),          # = lj_contour T_B(C2) (no W_X)
                        TB_psiX=0.9*lam0/(Wps + den0),                        # = lj_contour T_B(C2,+X)
                        TB_X=0.9*lam0/(WaX + Ebd + WX + Etop + Eedge),
                        TB_Y=0.9*lam0/(WaY + Ebd + WX + Etop + Eedge))
        # window bump G (LARGESPIN structure G): a e^{-(nu beta/3)^2}/m^2, nu < 30 window part with |X|, psi2 elsewhere
        S = np.exp(-(nuG*beta/3)**2)
        lo = nuG < 30
        cG_psi = (psG*S*cg).sum()
        cG_X = (np.abs(XG)*S*cg)[lo].sum() + (psG*S*cg)[~lo].sum()
        cG_psiX = cG_psi + (np.abs(XG - psG)*S*cg)[lo].sum()
        r['G'] = dict(psi=0.9*lam0/cG_psi, X=0.9*lam0/cG_X, psiX=0.9*lam0/cG_psiX)
        # D_h class: Gamma = c0 D_h(x) on all nu (fit class of fw_gamma); readings R1 (raw = A[env D]) and R2 (raw = env D)
        sel = nuG < m*min(1.0/mfl, 1.0)*0.999
        Dg_ = N.Dsh(h, J, nuG[sel]/m).real
        ring = (psG[sel]*Dg_*cg[sel]).sum()                 # R1: exact psi2 action per unit c0 (= de-shift identity)
        R2 = ((XG[sel]*Dg_).real*cg[sel]).sum()             # R2: Re sum X env D_h
        R2w = (((XG - psG)[sel]*Dg_).real*cg[sel])[nuG[sel] < 30].sum()
        r['Dh'] = dict(ring=ring, R2=R2, R2win=R2w, c0tol_R1=0.9*lam0/abs(ring), c0tol_R2=0.9*lam0/abs(R2))
        res['h'][h] = r
    return res


if __name__ == '__main__':
    out, lines = [], []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    pr('n2_tb: T_B (0.9 lam0 / weight, c2 units) with the window weight psi2 (HLJ) vs |X|, |Y|; G-bump; D_h class.  Dg = 1000')
    for J, mf in STATES:
        R = run(J, mf)
        if R is None:
            continue
        out.append(R)
        pr(f"J={J} m/L={R['mfl']:.4f} m={R['m']:.1f} beta={R['beta']:.4f} bhat={R['bhat']:.1f}")
        for h, r in R['h'].items():
            pr(f"  {h}: lam0={r['lam0']:+.4f} Ns={r['Ns']:.3e}")
            for c in CS:
                q = r[c]
                pr(f"     c={c:g} nu_B={q['nuB']:.2f}: W_psi(GL)={q['WpsGL']:.3e} W_psi(grid)={q['Wps']:.3e} W_|X|={q['WaX']:.3e} "
                   f"W_|Y|={q['WaY']:.3e} E_bd={q['Ebd']:.2e} W_X[nuB,25]={q['WX']:.2e} | T_B psi={q['TB_psi']:.0f} psi+X={q['TB_psiX']:.0f} "
                   f"|X|={q['TB_X']:.0f} |Y|={q['TB_Y']:.0f}")
            g = r['G']; d = r['Dh']
            pr(f"     G-bump a_tol (c2): psi2 {g['psi']:.0f}; psi2+|X-psi2| {g['psiX']:.0f}; |X| {g['X']:.0f}")
            pr(f"     D_h class per unit c0 (in N_s): R1 (A-form, exact = ring) {d['ring']:+.3e}; R2 (raw = D_h) {d['R2']:+.3e} "
               f"(window part of X-psi2: {d['R2win']:+.3e});  c0_tol R1 {d['c0tol_R1']:.3g}  R2 {d['c0tol_R2']:.3g}")
    open(os.path.join(N.CACHE, 'n2_tb.txt'), 'w').write('\n'.join(lines) + '\n')
    json.dump(out, open(os.path.join(N.CACHE, 'n2_tb.json'), 'w'), default=float, indent=0)
