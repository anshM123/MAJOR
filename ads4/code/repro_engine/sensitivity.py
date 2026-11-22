"""repro_engine/sensitivity.py -- how large an error in rho(nu) moves s_d of Q's binding state to the delta = 0.2 fail point.

Quantities (engine K = 5 rho as the baseline; H0 gauge; Dg = 300 and 1000):
  (a) worst-case shape error: rho -> rho (1 + eps(nu)), |eps(nu) - eps(0.02)| <= e0, sign chosen against the margin; linear in e0
      (exact: s_d is linear in Gamma_K, Gamma_K linear in rho at fixed c; c re-fixed at nu = 0.02 -- the shift eps(0.02) is a constant
      and drops out).  e0_fail = (s_fail - s_d)/A,  A = d s_d/d e0 (worst sign) -- also restricted to nu < 2 and to nu < 4;
  (b) uniform-in-nu relative error profile (eps = e0 const for nu >= 0.5, 0 at 0.02): s_d(e0) is linear; e0_fail;
  (c) scaling of the beyond-leading-Regge content: rho = rho_LR (1 + lambda (rho/rho_LR - 1)), rho_LR = T_st/T_ref; lambda_fail;
  (d) the same expressed as the FF/N shift per unit error.
Usage: python3 sensitivity.py -> cache/sensitivity.txt
"""
import os
import numpy as np
import sdchain as C

HERE = os.path.dirname(os.path.abspath(__file__))


def sd_of_rho(nus, rho, D):
    g = C.gammaK(nus, rho, D)
    return C.sd(nus, g['GK'])[0], g


def main():
    out = []
    pr = lambda s='': (print(s, flush=True), out.append(s))
    for Dg in (300, 1000):
        H = C.hyb_numbers(Dg)
        sfail = H['FF']/(2*H['a'])
        nus, rho = C.engine_ratio(Dg, 5)
        s0, g = sd_of_rho(nus, rho, Dg)
        slope = -2*H['a']/H['N']
        pr(f'# Dg = {Dg}: s_d(K=5 engine) = {s0:+.5f}, s_fail = {sfail:.4f}, margin = {sfail - s0:.4f}; dFF/N/ds_d = {slope:+.4f}')
        # (a) linear response: d s_d / d eps(nu_i) for a relative perturbation at one nu_i (i >= 1)
        m, _ = C.mass(Dg)
        nuw, cv, cd = C.chi('H0')
        k = min(len(nus), len(cd))
        # Gamma_exact = m^2 (c Tref rho - Tst), c = Tst0/(Tref0 rho0): d Gamma/d eps_i = m^2 c Tref_i rho_i (delta at i), the change of
        # c through eps(0.02) is a constant offset; GK subtracts (GK0/NH0) NH with GK0 = 0 unaffected.
        c = g['Tst'][0]/(g['Tref'][0]*rho[0])
        dG = m*m*c*g['Tref']*rho
        resp = 0.04*cd[:k]*dG[:k].real            # d s_d / d eps_i  (eps real)
        respi = -0.04*cd[:k]*dG[:k].imag          # d s_d / d (i eps_i) (phase error)
        for lo, hi in ((0, 41), (0, 2), (0, 4), (1, 2)):
            sel = (nus[:k] >= lo) & (nus[:k] < hi) & (np.arange(k) > 0)
            A = np.sum(np.abs(resp[sel])); Ai = np.sum(np.abs(respi[sel]))
            pr(f'  (a) worst-case |rel. error of rho| on nu in [{lo},{hi}): A = {A:.4g} per unit -> e0_fail = {(sfail - s0)/A:.2e} '
               f'(modulus);  phase error: {(sfail - s0)/Ai:.2e} rad')
        # (b) uniform relative error for nu >= 0.5 (smooth onset irrelevant at this level)
        for lo in (0.5, 0.06):
            sel = (nus[:k] >= lo)
            Ab = np.sum(resp[sel])
            pr(f'  (b) uniform eps = e0 on nu >= {lo}: d s_d/d e0 = {Ab:+.4g} -> e0_fail = {(sfail - s0)/Ab:+.2e}')
        # (c) scaling of the beyond-LR content
        rlr = g['Tst']/g['Tref']
        sLR, _ = sd_of_rho(nus, rlr, Dg)
        s2, _ = sd_of_rho(nus, rlr*(1 + 2.0*(rho/rlr - 1)), Dg)
        lam = 1 + (sfail - s0)/(s2 - s0)
        dev = np.abs(rho/rlr/(rho[0]/rlr[0]) - 1)
        pr(f'  (c) pure leading-Regge rho: s_d = {sLR:+.4f}; scaling of (rho/rho_LR - 1) by lambda: s_d(2) = {s2:+.4f}, '
           f'lambda_fail = {lam:.3f};  max|rho/rho_LR - 1| (normalised at 0.02) on nu<4: {np.max(dev[nus < 4]):.2e}, '
           f'in units of 1/m^2: {np.max(dev[nus < 4])*m*m:.1f}')
        pr(f'  (d) FF/N per unit s_d: {slope:+.4f};  FF/N = 0 at e0 = (a) values above')
    open(os.path.join(HERE, 'cache', 'sensitivity.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
