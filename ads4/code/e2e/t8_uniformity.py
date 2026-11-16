"""
Task 5 (uniformity) from the stored mode data (data_modes/):
 (a) nu growing with m: nu = x m at fixed beta; ratios exact/leading for the helicity-preserving and flip components;
 (b) small impact parameter: J = 2, 3, 5 (beta ~ (2J+1)/m) at m ~ 20, 40, 80, including beta ~ 1/nu and below;
     leading formula, and the same with the exact normalisation polynomials of the conserved structures restored:
         l = 2:  (J+1/2)^4 -> (J-1)J(J+1)(J+2),  (Delta-3/2)^4 -> Delta(Delta-1)(Delta-2)(Delta-3)
         l = 1:  (J+1/2)^2 -> J(J+1),            (Delta-3/2)^2 -> (Delta-1)(Delta-2)
     (these are the exact prefactors L_0(Delta,J) L_0(3-Delta,-1-J) of the conserved structures, see SPINNING_E2E.md sect. 4.4).
Numerical accuracy: for large nu the moment is suppressed like gamma_{a_x} gamma_{a_y} ~ e^{-pi nu} against an integrand of order one;
the column 'supp' gives |moment|/sum|integrand| (results with supp < 1e-6 are not reliable in double precision and are flagged).
Usage: python3 t8_uniformity.py l
"""
import sys, os, numpy as np, mpmath as mp
from t4_spinning_e2e import main_cases, small_cases, get_moments, fname, DPHI
from t7_corrections import measured_components
from predict import scalar_prefactor, components
from formulas import Fnu_np, Fnu_mp
from cases import info
mp.mp.dps = 25

def moments_with_supp(mom, k, nu, modes):
    """returns dict mode -> (-1)^J value (mp) and the suppression ratio |sum|/sum|.| of the least suppressed mode"""
    g = mom.g; U, S = g['U'], g['S']; xi = -np.log(g['us'])
    res = {}; supp = 1.0
    rho, rhob = S**2, (S*U)**2
    with np.errstate(divide='ignore', invalid='ignore'):
        logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k-4)*np.log(rho*rhob) + mom.logscale
    logw[~g['ok']] = -np.inf
    sh = np.max(logw)
    if not hasattr(mom, '_Fcache'): mom._Fcache = {}
    if nu not in mom._Fcache:
        mom._Fcache[nu] = Fnu_np(nu, list(range(0, 2*mom.l + 1)), xi)
    Fs = mom._Fcache[nu]
    for (q1, q3) in modes:
        integ = np.exp(logw - sh)*Fs[abs(q1+q3)][:, None]*mom.modes[(q1, q3)]*g['W']
        val = integ.sum()
        supp = min(supp, abs(val)/np.abs(integ).sum())
        res[(k, nu, q1, q3)] = (-1)**mom.J*mp.mpc(complex(val))*mp.exp(mp.mpf(sh))
    return res, supp

if __name__ == "__main__":
    l = int(sys.argv[1])
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    allmodes = [(q1, q3) for q1 in range(-l, l+1) for q3 in range(-l, l+1)]
    from moments import ScalarMoments
    from t3_scalar_mixed import prediction as scalar_prediction
    pr(f"# l = {l}.  (a) nu = x m at fixed beta.  R_pres, R_flip: measured/leading prediction for the two helicity components (two-mode extraction")
    pr("#     of t7);  R_scal: the same ratio for the scalar skeleton with dimensions (l+1, Dphi) [exact scalar block / scalar leading formula];")
    pr("#     zeroP = |P_{-1/2+i nu}(cosh beta)| / envelope: the ratios R are ill-conditioned near the zeros of the harmonic function (zeroP << 1);")
    pr("#     supp = |moment| / sum |integrand|: double precision fails for supp < 1e-6.")
    pr("#     corr = [exact normalisation polynomials / leading] * P(eta_AdS - 2 l^2/m^2)/P(eta_AdS): the spin-dependent improvement of sect. 6.3")
    pr("#            (for l = 2 the additional measured piece of sect. 6.3 is NOT included)")
    for beta_t in (0.4, 1.0, 2.0):
        for k in (4, 6):
            pr(f"# beta ~ {beta_t}, k = {k}")
            pr("#     m       nu     x=nu/m    R_pres       R_flip      R_scal     R_pres/R_scal  R_pres/(R_scal*corr)  (R_scal-1)/x^2   zeroP    supp")
            for c in main_cases():
                if not os.path.exists(fname(l, *c)): continue
                m, b_, eta = info(*c)
                if abs(float(b_) - beta_t) > 0.06: continue
                mom = get_moments(l, *c)
                sm = ScalarMoments(c[0], c[1], float(l + 1), DPHI, int(min(14000, max(800, 60*c[0]))))
                for x in (0.05, 0.1, 0.2, 0.3, 0.5, 0.8):
                    nu = round(x*float(m), 3)
                    if nu > 12.5: continue
                    res, supp = moments_with_supp(mom, k, nu, allmodes)
                    Rp, Rf, worst, bb = measured_components(l, c[0], c[1], DPHI, k, nu, res)
                    sg, lg = sm.moment(k, nu)
                    Rs = sg*mp.exp(mp.mpf(lg))/scalar_prediction(c[0], c[1], float(l + 1), DPHI, k, nu).real
                    P0 = Fnu_mp(nu, 0, bb); dP = mp.diff(lambda t: Fnu_mp(nu, 0, t), bb)
                    zeroP = abs(P0)/mp.sqrt(abs(P0)**2 + abs(dP)**2/(nu*nu + 0.25))
                    flag = "" if supp > 1e-6 else "  (unreliable)"
                    # spin-dependent improvement of sect. 6.3: exact normalisation polynomials and eta_eff = eta - 2 l^2/m^2
                    D_, J_ = mp.mpf(c[0]), c[1]
                    Am, Ap = J_ + mp.mpf(1)/2, D_ - mp.mpf(3)/2
                    fac = mp.mpf(1)
                    for i in range(1, l + 1):
                        fac *= (Am**2 - (i - mp.mpf(1)/2)**2)/Am**2*(Ap**2 - (i - mp.mpf(1)/2)**2)/Ap**2
                    Pe = Fnu_mp(nu, 0, mp.acosh(eta - 2*l*l/(m*m)))
                    corr = fac*Pe/P0
                    pr(f"  {float(m):7.2f}  {nu:7.3f}  {x:5.2f}   {float(mp.re(Rp)):10.6f}  {float(mp.re(Rf)):10.6f}  {float(Rs):10.6f}   {float(mp.re(Rp)/Rs):10.6f}   {float(mp.re(Rp/(Rs*corr))):10.6f}     {float((Rs-1)/x**2):9.3f}      {float(zeroP):5.2f}   {supp:.0e}{flag}")
    pr("# (b) small impact parameter (fixed J, beta ~ (2J+1)/m):  exact/prediction for the modes (l,-l) and (l,l)")
    pr("#     lead : the leading formula (predict.py)")
    pr("#     impr : the same with the exact normalisation polynomials of the conserved structures restored,")
    pr("#              (J+1/2)^{2l} -> prod_{i=1..l} [(J+1/2)^2 - (i-1/2)^2],   (Delta-3/2)^{2l} -> prod_{i=1..l} [(Delta-3/2)^2 - (i-1/2)^2]")
    pr("#     impr2: in addition the helicity-preserving harmonic function evaluated at  eta_eff = eta_AdS - 2 l^2/m^2")
    pr("#              (i.e. (J+1/2)^2 -> (J+1/2)^2 - l^2 and (Delta-3/2)^2 -> (Delta-3/2)^2 - l^2 inside P; see SPINNING_E2E.md sect. 6.3)")
    pr("#     impr3: in addition the flip harmonic function F_{2l}(beta) ~ ((J+1/2) nu/m)^{2l}/(2l)! multiplied by prod[(J+1/2)^2-(i-1/2)^2]/(J+1/2)^{2l}")
    pr("#              (the flat-space ratio d^J_{l,-l}/d^J_{l,l} ~ [(J+l)!/(J-l)!] (theta/2)^{2l}/(2l)!)")
    pr("#     nu*beta is the argument of the flat-space Bessel functions;  supp as in (a)")
    for k in (4, 6):
        pr(f"# k = {k}")
        pr("#     m     J    beta     nu   nu*beta |  (l,-l):  lead      impr     impr2     impr3  |  (l,l):  lead      impr     impr2     impr3  |  supp")
        for c in small_cases(l) + [cc for cc in main_cases() if abs(float(info(*cc)[1]) - 0.4) < 0.06]:
            if not os.path.exists(fname(l, *c)): continue
            m, b_, eta = info(*c)
            mom = get_moments(l, *c)
            D, J = mp.mpf(c[0]), c[1]
            Am, Ap = J + mp.mpf(1)/2, D - mp.mpf(3)/2
            fac = mp.mpf(1)
            for i in range(1, l + 1):
                fac *= (Am**2 - (i - mp.mpf(1)/2)**2)/Am**2*(Ap**2 - (i - mp.mpf(1)/2)**2)/Ap**2
            facJ = mp.mpf(1)
            for i in range(1, l + 1):
                facJ *= (Am**2 - (i - mp.mpf(1)/2)**2)/Am**2
            eta_eff = eta - 2*l*l/(m*m)
            for nu in (0.0, 1.0, 3.0, 6.0):
                res, supp = moments_with_supp(mom, k, nu, [(l, -l), (l, l)])
                row = []
                for (q1, q3) in ((l, -l), (l, l)):
                    pref, beta = scalar_prefactor(l, c[0], c[1], DPHI, k, nu)
                    A, B = components(l, k, nu, q1, q3)
                    P0 = Fnu_mp(nu, 0, beta); F0 = Fnu_mp(nu, 2*l, beta)
                    pd = pref*(A*P0 + B*F0)
                    r = mp.re(res[(k, nu, q1, q3)]/pd)
                    if eta_eff > 1:
                        Pe = Fnu_mp(nu, 0, mp.acosh(eta_eff))
                    else:
                        Pe = mp.hyp2f1(mp.mpf(1)/2 - 1j*nu, mp.mpf(1)/2 + 1j*nu, 1, (1 - eta_eff)/2)
                    pd2 = pref*fac*(A*Pe + B*F0)
                    r2 = mp.re(res[(k, nu, q1, q3)]/pd2)
                    pd3 = pref*fac*(A*Pe + B*F0*facJ)
                    r3 = mp.re(res[(k, nu, q1, q3)]/pd3)
                    row.append(f"{float(r):9.5f} {float(r/fac):9.5f} {float(r2):9.5f} {float(r3):9.5f}")
                pr(f"  {float(m):6.2f}  {J:2d}  {float(b_):6.4f}  {nu:4.1f}  {float(nu*b_):6.3f}  |  " + "   |  ".join(row) + f"   |  {supp:.0e}")
    open(f'out_t8_uniformity_l{l}.txt', 'w').write("\n".join(out) + "\n")
