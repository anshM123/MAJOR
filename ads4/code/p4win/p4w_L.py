"""p4win/p4w_L.py (P4_WINDOW.md §5): Lemma P4W along Theorem L''s range (J >= J*(m) = ceil((0.37211 m - 1)/2), gap and other masses,
J up to 10^6).  No computed spinning O(1/m^2) exists there, so the window O(1/m^2) part is Theorem L''s premise P_L'(i) in the window:
Lam_amp x the scalar CHMRS beyond-hybrid O(1/m^2) term B2 = (T_1 - nu^2 g_1)/m^2 (sign-free), Lam_amp = 3 (and 10).
Budget: Theorem L''s per-entry half floor, F/2, with F = F_cert (LRINF_TIGHT §5 [certified]) and F_num (numerical min F_LR), quoted.
Window nu < 3/beta; k = 2 (weight 2 psi2) and k = 3 (weight u |psi3|); O(1/m^4) scalar seed S4 and the c4w unit; boundary term with
M_w = 2 sup_{[0.7,1.3] nu_w} |Lam_amp B2 + S4/m^4|; edge with the scalar envelope e_sc (+- helicity factor (1 - x^2)^-2), local rate.
-> cache/L.txt, cache/L.json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
import p4w_common as C, p4w_contour as K

OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))
Lg = C.L
# LRINF_TIGHT §5 table: (label, J, m/Lam, F_num, F_cert)
ROWS = [('gap', 223, 1.203, 3.41e-5, 1.92e-5), ('gap', 250, 1.225, 2.91e-5, 1.89e-5), ('gap', 300, 1.265, 2.37e-5, 1.78e-5),
        ('gap', 344, 1.299, 2.10e-5, 1.69e-5), ('gap', 400, 1.342, 1.89e-5, 1.60e-5), ('gap', 500, 1.415, 1.69e-5, 1.49e-5),
        ('gap', 700, 1.550, 1.51e-5, 1.36e-5), ('gap', 857, 1.648, 1.44e-5, 1.30e-5), ('gap', 10**4, 4.585, 6.86e-6, 5.53e-6),
        ('gap', 10**5, 14.19, 2.08e-6, 1.68e-6), ('gap', 10**6, 44.76, 5.51e-7, 4.13e-7), ('1.1 m_gap', 250, 1.347, 3.43e-5, 1.92e-5),
        ('1.3 m_gap', 307, 1.652, 3.46e-5, 1.92e-5), ('1.5 Lam', 279, 1.5, 3.44e-5, 1.92e-5), ('2 Lam', 372, 2.0, 3.48e-5, 1.92e-5),
        ('5 Lam', 929, 5.0, 3.54e-5, 1.92e-5)]


def mfl_exact(label, J, mfl_tab):
    if 'gap' in label:
        mg = np.sqrt((1000 - 2.0)*(1000 - 1.0 + 2*J))/Lg
        return mg*(1.1 if label.startswith('1.1') else 1.3 if label.startswith('1.3') else 1.0)
    return mfl_tab


def main():
    pr(__doc__.split('->')[0])
    pr(' state       J      m/L    beta   nu_w  | budget F_cert/2 | window O(1/m^2) (Lam_amp=3) | O(1/m^4) | bd    | edge+top (+-) | '
       'd/budget (cert) [Lam_amp = 3 / 10] | d/budget (num) | K4w tol | edge amplitude tol vs scalar (c2)')
    res = []
    for lab, J, mtab, Fn, Fc in ROWS:
        mfl = mfl_exact(lab, J, mtab); m = mfl*Lg; beta = C.beta_of(J, mfl); u = 1/mfl**2
        nuw = 3.0/beta; pw = nuw/Lg
        nu = np.linspace(0, nuw, 4001)[1:]; dn = nu[1] - nu[0]
        P, dPb = C.Pd(beta, nu); P = P.real
        bb = (2*J + 1)/mfl*nu/Lg
        B2 = (C.T(1, beta, nu, P, dPb)[0] - nu**2*C.gflat(1, bb))/m**2
        S4 = (C.T(2, beta, nu, P, dPb)[0] - nu**4*(C.gflat(2, bb) + (3/16)*j0(bb)))/m**4
        ps2 = C.psi2(nu/Lg); ps3 = C.psi3(nu/Lg)
        act2 = lambda f: 2*float((ps2*f).sum())*dn/Lg
        w2 = abs(act2(B2)) + u*float((np.abs(ps3)*np.abs(B2)).sum())*dn/Lg
        w4 = abs(act2(S4)) + u*float((np.abs(ps3)*np.abs(S4)).sum())*dn/Lg
        c4w = (abs(act2(np.cosh(beta)**2*nu**4*P)) + abs(act2(np.cosh(beta)*np.sinh(beta)*nu**4*dPb.real)))/m**4
        nn = np.linspace(0.7*nuw, 1.3*nuw, 41); Pn, dPn = C.Pd(beta, nn); bn = (2*J + 1)/mfl*nn/Lg
        B2n = (C.T(1, beta, nn, Pn.real, dPn)[0] - nn**2*C.gflat(1, bn))/m**2
        S4n = (C.T(2, beta, nn, Pn.real, dPn)[0] - nn**4*(C.gflat(2, bn) + (3/16)*j0(bn)))/m**4
        row = dict(label=lab, J=J, mfl=mfl, beta=beta, nuw=nuw, Fc=Fc, Fn=Fn, c4w=c4w)
        for La in (3.0, 10.0):
            Mw = 2*float(np.max(np.abs(La*B2n + S4n)))
            env = K.e_sc(beta, mfl, hel=-2.0)
            bd2, ed2, tp2 = K.contour_bound(beta, mfl, pw, Mw, env, 'local', k=2, Aunif=K.A_curv(beta, J, mfl))
            bd3, ed3, tp3 = K.contour_bound(beta, mfl, pw, Mw, env, 'local', k=3, u=u, Aunif=K.A_curv(beta, J, mfl))
            d = La*w2 + w4 + bd2 + bd3 + ed2 + tp2 + ed3 + tp3
            row[f'L{La:g}'] = dict(Mw=Mw, w2=La*w2, w4=w4, bd=bd2 + bd3, ag=ed2 + tp2 + ed3 + tp3, d=d,
                                   K4w_tol=(Fc/2 - (d - w4))/c4w, Medge_tol_c2=(Fc/2 - (d - ed2 - ed3))/((ed2 + ed3)/float(env(np.array([1.0]))[0]))*m**2,
                                   Msc_c2=float(env(np.array([1.0]))[0])*m**2)
        r3 = row['L3']; r10 = row['L10']
        pr(f" {lab:10s} {J:7d} {mfl:7.3f} {beta:6.3f} {nuw:6.2f} | {Fc/2:.2e} | {r3['w2']/(Fc/2):.1e} | {w4/(Fc/2):.1e} | {r3['bd']/(Fc/2):.1e} | "
           f"{r3['ag']/(Fc/2):.1e} | {r3['d']/(Fc/2):.2e} / {r10['d']/(Fc/2):.2e} | {r3['d']/(Fn/2):.2e} | {r3['K4w_tol']:.2e} | "
           f"{r3['Medge_tol_c2']:.1e} vs {r3['Msc_c2']:.1e}")
        res.append(row)
    pr('max d/budget (certified floor, Lam_amp = 3): %.3g; (Lam_amp = 10): %.3g' % (max(r['L3']['d']/(r['Fc']/2) for r in res),
                                                                                   max(r['L10']['d']/(r['Fc']/2) for r in res)))
    json.dump(res, open(os.path.join(C.CACHE, 'L.json'), 'w'), default=float)
    open(os.path.join(C.CACHE, 'L.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
