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
    pr(f"# l = {l}.  (a) nu = x m at fixed beta:  R_pres, R_flip = measured/leading prediction  (two-mode extraction of t7)")
    for beta_t in (0.4, 1.0, 2.0):
        for k in (4, 6):
            pr(f"# beta ~ {beta_t}, k = {k}")
            pr("#     m       nu     x=nu/m    R_pres       R_flip      supp      (R_pres-1)/x^2   (R_flip-1)/x^2")
            for c in main_cases():
                if not os.path.exists(fname(l, *c)): continue
                m, b_, eta = info(*c)
                if abs(float(b_) - beta_t) > 0.06: continue
                mom = get_moments(l, *c)
                for x in (0.05, 0.1, 0.2, 0.3, 0.5, 0.8):
                    nu = round(x*float(m), 3)
                    if nu > 12.5: continue
                    res, supp = moments_with_supp(mom, k, nu, allmodes)
                    Rp, Rf, worst, bb = measured_components(l, c[0], c[1], DPHI, k, nu, res)
                    flag = "" if supp > 1e-6 else "  (unreliable)"
                    pr(f"  {float(m):7.2f}  {nu:7.3f}  {x:5.2f}   {float(mp.re(Rp)):10.6f}  {float(mp.re(Rf)):10.6f}   {supp:.1e}   {float((mp.re(Rp)-1)/x**2):10.4f}  {float((mp.re(Rf)-1)/x**2):10.4f}{flag}")
    pr("# (b) small impact parameter (J = 1, 2, 3, 5):  exact/leading for the modes (l,-l) and (l,l), and the same with the exact")
    pr("#     normalisation polynomials restored ('impr');  nu*beta is the argument of the flat-space Bessel functions")
    for k in (4, 6):
        pr(f"# k = {k}")
        pr("#     m     J    beta     nu   nu*beta |  (l,-l): leading    impr   |  (l,l): leading     impr   |   supp")
        for c in small_cases(l):
            if not os.path.exists(fname(l, *c)): continue
            m, b_, eta = info(*c)
            mom = get_moments(l, *c)
            D, J = mp.mpf(c[0]), c[1]
            if l == 2:
                fac = ((J-1)*J*(J+1)*(J+2))/(J + mp.mpf(1)/2)**4*(D*(D-1)*(D-2)*(D-3))/(D - mp.mpf(3)/2)**4
            else:
                fac = (J*(J+1))/(J + mp.mpf(1)/2)**2*((D-1)*(D-2))/(D - mp.mpf(3)/2)**2
            for nu in (0.0, 1.0, 3.0, 6.0):
                res, supp = moments_with_supp(mom, k, nu, [(l, -l), (l, l)])
                row = []
                for (q1, q3) in ((l, -l), (l, l)):
                    pref, beta = scalar_prefactor(l, c[0], c[1], DPHI, k, nu)
                    A, B = components(l, k, nu, q1, q3)
                    pd = pref*(A*Fnu_mp(nu, 0, beta) + B*Fnu_mp(nu, 2*l, beta))
                    r = mp.re(res[(k, nu, q1, q3)]/pd)
                    row.append(f"{float(r):9.5f}  {float(r/fac):9.5f}")
                pr(f"  {float(m):6.2f}  {J:2d}  {float(b_):6.4f}  {nu:4.1f}  {float(nu*b_):6.3f}  |  " + "   |  ".join(row) + f"   |  {supp:.0e}")
    open(f'out_t8_uniformity_l{l}.txt', 'w').write("\n".join(out) + "\n")
