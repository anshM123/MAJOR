"""
lsl/lsl_window.py (LARGESPIN_L.md §4): the Regge-window tolerance of Omega^EE along the near-gap large-spin line, J = 26 ... 10^6,
in two conventions.

For a state (J, m) at Delta_g = 1000 (beta, u = Lam^2/m^2, worse helicity sigma):
  floor  F_LR = 2 Omega_2(beta) + u sigma Omega_3(beta)          (leading-Regge per-state action, lsl_core.Omega)
  E(nu)  = local amplitude of the de-shifted LR image P_{-1/2+i nu}(cosh beta)      (running max of |P| over |nu'-nu| <= pi/beta)
  Er(nu) = local amplitude of the RAW LR image A[P](nu) = g+ P(nu+i) + g- P(nu-i)    (n2_common.Aop; running max likewise)
  * de-shift convention (R1 holds: the correction is A[f] with f strip-analytic; its k = 2 cost is the psi_2 cost of f):
        W_ds(nu_W) = int_0^{nu_W} (2 psi_2 + u |psi_3|)(nu/Lam) E(nu) dnu/Lam,        eps_ds = kappa F_LR / W_ds
  * raw convention (sign-free, any nu-shape on the raw row, X weight, N2):
        W_raw(nu_W) = int_0^{nu_W} [2 |chi(nu)| Er(nu) + u |psi_3(nu/Lam)| E(nu)] dnu/Lam,  eps_raw = kappa F_LR / W_raw
  eps = tolerated |relative correction| (relative to the local LR amplitude) on nu < nu_W; kappa = 1/2.
  * absolute (HLJ-style, sup|Gamma| in c2 on the window, de-shift weights): Gam_tol = kappa F_LR m^2 / int_0^{nu_W}(2psi_2 + u|psi_3|) dnu/Lam
Windows: nu_B = 3m/(2J+1) (HLJ, c = 3), nu_W = 10 (fold region, J-independent), nu_W = 20.
Output: cache/window_map.txt, cache/window_map.json.
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'n2sub'))
import numpy as np
import lsl_core as C
import n2_common as N

DG = 1000
L = C.lam(DG)
KAP = 0.5


def runmax(y, nu, half):
    out = np.empty_like(y)
    j0 = 0; j1 = 0; n = len(nu)
    for i in range(n):
        while nu[j0] < nu[i] - half:
            j0 += 1
        while j1 < n - 1 and nu[j1 + 1] <= nu[i] + half:
            j1 += 1
        out[i] = y[j0:j1 + 1].max()
    return out


def state(J, mfac=1.0):
    m = C.mgap(DG, J)*mfac
    beta = C.beta_of(J, m); u = (L/m)**2
    nuB = 3*m/(2*J + 1)
    numax = max(25.0, nuB) + 2
    dnu = min(0.01, 0.05/beta)
    nu = np.arange(dnu/2, numax + dnu, dnu)
    P, dP = C.P_real(beta, nu, nphi=max(1500, int(400*beta*numax/50) + 1500))
    half = np.pi/beta
    E = runmax(np.abs(P), nu, half)
    xd = N.xdata(DG)
    Rm, Rp = N.Rk()
    # chi = X/(R_- + R_+): the weight of the UN-normalised raw image R_- P(nu+i) + R_+ P(nu-i) (X A[P] = chi (R_- P_+ + R_+ P_-));
    # regular at the real zeros of R_- + R_+ (nu = 1.1714, 6.6428), where it is interpolated across +-0.03
    gx = xd['nu']; chi = np.abs(xd['X']/(Rm(gx + 0j) + Rp(gx + 0j)))
    bad = (np.abs(gx - 1.1714) < 0.03) | (np.abs(gx - 6.6428) < 0.03)
    chi = np.interp(gx, gx[~bad], chi[~bad])
    Xabs = np.interp(nu, gx, chi)
    Pf = lambda z: N.Pnu(beta, z, nphi=max(1500, int(400*beta*numax/50) + 1500))
    Ar = Rm(nu + 0j)*Pf(nu + 1j) + Rp(nu + 0j)*Pf(nu - 1j)
    Er = runmax(np.abs(Ar), nu, half)
    p = nu/L
    w2 = C.psi(2, p); w3 = np.abs(C.psi(3, p))
    fl = {s: C.LR_action(beta, u, C.SIGMA[s], DG) for s in ('++', '+-')}
    sw = min(fl, key=fl.get); F = fl[sw]
    res = dict(J=J, mfac=mfac, m=m, beta=beta, bhat=C.bhat_of(J, m, DG), u=u, nuB=nuB, F_LR=F, helicity=sw,
               F_LR_pp=fl['++'], F_LR_pm=fl['+-'])
    for tag, nW in (('B', nuB), ('10', 10.0), ('20', 20.0)):
        sel = nu <= nW
        Wds = np.trapezoid(((2*w2 + u*w3)*E)[sel], nu[sel])/L
        Wraw = np.trapezoid((2*Xabs*Er + u*w3*E)[sel], nu[sel])/L
        Wabs = np.trapezoid((2*w2 + u*w3)[sel], nu[sel])/L
        res['eps_ds_' + tag] = KAP*F/Wds
        res['eps_raw_' + tag] = KAP*F/Wraw
        res['Gtol_' + tag] = KAP*F*m*m/Wabs
        res['E0'] = float(E[0]); res['Er0'] = float(Er[0])
    # where the floor is decided: nu below which half of |int psi2 P| accumulates (signed integral of the k=2 LR action)
    dens = 2*w2*P/L
    cum = np.cumsum(dens)*dnu
    res['F2_partial_25'] = float(cum[-1])
    res['Er_over_E_at_1'] = float(np.interp(1.0, nu, Er/np.maximum(E, 1e-300)))
    res['chiEr_over_psiE_at_1'] = float(np.interp(1.0, nu, 2*Xabs*Er/np.maximum(2*w2*E, 1e-300)))
    return res


def main():
    Js = [26, 50, 100, 150, 200, 344, 500, 1000, 2000, 3000, 5000, 10000, 30000, 100000, 300000, 1000000]
    out = []
    t0 = time.time()
    for mf in (1.0, 1.3):
        for J in Js:
            r = state(J, mf)
            out.append(r)
            print(f"J={J:8d} m/mgap={mf} beta={r['beta']:.4f} bhat={r['bhat']:.1f} F_LR={r['F_LR']:.4e}({r['helicity']}) "
                  f"eps_ds[B,10,20]={r['eps_ds_B']:.3e},{r['eps_ds_10']:.3e},{r['eps_ds_20']:.3e} "
                  f"eps_raw[B,10,20]={r['eps_raw_B']:.3e},{r['eps_raw_10']:.3e},{r['eps_raw_20']:.3e} "
                  f"Gtol[B,10,20]={r['Gtol_B']:.3g},{r['Gtol_10']:.3g},{r['Gtol_20']:.3g} nuB={r['nuB']:.2f} "
                  f"rawcost/dscost(nu=1)={r['chiEr_over_psiE_at_1']:.3g}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open(os.path.join(C.CACHE, 'window_map.json'), 'w'), indent=1)
    with open(os.path.join(C.CACHE, 'window_map.txt'), 'w') as fh:
        fh.write(__doc__ + "\n")
        fh.write("J  m/mgap  beta  bhat  F_LR  hel  eps_ds_B eps_ds_10 eps_ds_20  eps_raw_B eps_raw_10 eps_raw_20  Gtol_B Gtol_10 Gtol_20  nuB  Er/E(nu=1)\n")
        for r in out:
            fh.write(f"{r['J']} {r['mfac']} {r['beta']:.5f} {r['bhat']:.2f} {r['F_LR']:.6e} {r['helicity']} "
                     f"{r['eps_ds_B']:.4e} {r['eps_ds_10']:.4e} {r['eps_ds_20']:.4e} {r['eps_raw_B']:.4e} {r['eps_raw_10']:.4e} "
                     f"{r['eps_raw_20']:.4e} {r['Gtol_B']:.4g} {r['Gtol_10']:.4g} {r['Gtol_20']:.4g} {r['nuB']:.3f} {r['Er_over_E_at_1']:.4g}\n")


if __name__ == '__main__':
    main()
