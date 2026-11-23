r"""
ee5_sd_report (EESTEP_CLOSE.md §5): s_d(D) for J = 0 states of dimension D under the FIXED de-shifted derivative functional (Stage 6 H0
chi_d, Delta_g-independent), from the direct 7-d quadrature (repro_engine, read-only; cached D = 300, 1000, reference D_ref = 2e4) plus the new
D = 3000, 6000 runs of ee5_sd_run.py (eestep/cache), all on the same grid n8-8-18-18-6-20-20, nu < 12 (direct part of s_d).

Reading of the numbers.  The chain measures Gamma = m^2 (c T_ref rho - T_st), rho = I(D)/I(D_ref).  If the true image is
I(D) = I_inf (1 + g(nu)/m^2 + O(m^-4)), then Gamma -> g (1 - m^2/m_ref^2): the c2-unit correction s_d is m-independent and the finite-reference
factor (1 - m^2/m_ref^2) is the only D-dependence.  F_NORM §2.6: exact cancellation of the g5 light term (P) needs the ++ J <= 3 correction
in c2 units to VANISH as m -> infinity at fixed functional; the double-counting reading predicts s_d ~ 0.51 unchanged.
Output: cache/ee5_sd_report.txt.
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
RE = os.path.join(HERE, '..', 'repro_engine')
sys.path.insert(0, RE)
import numpy as np
import direct_rho as DR
import sdchain as C

OUT = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


NS = '8-8-18-18-6-20-20'
DREF = 2e4


def path(D):
    for d in (os.path.join(RE, 'cache'), os.path.join(HERE, 'cache')):
        fn = os.path.join(d, 'direct_D%g_n%s_f64.npz' % (D, NS))
        if os.path.exists(fn):
            return fn
    return None


if __name__ == '__main__':
    pr('# ee5_sd_report: m-dependence of the engine\'s ++ J = 0 derivative-piece correction at fixed functional (F_NORM §2.6 deciding test)')
    nus, Iref, _ = DR.load(path(DREF))
    mref, _ = C.mass(DREF)
    rows = []
    pr(f'  grid {NS}, nu < 12, reference D_ref = {DREF:g} (m_ref = {mref:.1f}); Lambda_300 = 298.50')
    pr('     D      m     m/Lambda_300   s_d(nu<12)   s_v(nu<12)   s_d/(1 - m^2/m_ref^2)   windows [0,1) [1,2) [2,3) [3,4) [4,12)')
    for D in (300.0, 1000.0, 3000.0, 6000.0):
        fn = path(D)
        if fn is None:
            pr(f'  D = {D:g}: no run on this grid'); continue
        nu2, Ist, _ = DR.load(fn)
        assert np.allclose(nu2, nus)
        g = C.gammaK(nus, Ist/Iref, D)
        sd, sv = C.sd(nus, g['GK'], 'H0', 0.0, 12.0)
        win = [C.sd(nus, g['GK'], 'H0', lo, hi)[0] for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 12))]
        m = g['m']; corr = 1 - m*m/(mref*mref)
        rows.append((D, m, sd, sd/corr))
        pr(f'  {D:6.0f} {m:7.1f}   {m/298.50:8.2f}      {sd:+.5f}     {sv:+.4f}        {sd/corr:+.5f}            '
           + ' '.join(f'{w:+.4f}' for w in win))
    pr('  (the column s_d/(1 - m^2/m_ref^2) is NOT the right correction: the chain subtracts the m-independent 0.75 sum R_j (nu-ij)^2 and')
    pr('   const N^4 after the m^2 scaling, so a finite reference shifts s_d by -x (s_inf + S_sub), x = m^2/m_ref^2; done properly below)')
    pr('(b) infinitely heavy reference.  r_D(nu) = I(D; nu)/I(D; 0.02) removes D-dependent constants; r_inf = Richardson of r_6000 and r_20000')
    pr('    in 1/m^2 (error O(1/(m_a^2 m_b^2)) ~ 1e-16 relative); then the IDENTICAL chain with rho = r_D/r_inf and T_ref -> T_inf = sum_j R_j.')
    rr = {}
    for D in (300.0, 1000.0, 3000.0, 6000.0, DREF):
        fn = path(D)
        if fn is None:
            continue
        nu2, I, _ = DR.load(fn)
        rr[D] = I/I[0]
    if 6000.0 in rr:
        ma2 = C.mass(DREF)[0]**2; mb2 = C.mass(6000.0)[0]**2
        rinf = (ma2*rr[DREF] - mb2*rr[6000.0])/(ma2 - mb2)
        R = C.Rj(nus); Tinf = sum(R.values())
        NH = sum(R.values())
        res = []
        pr('     D      m/Lambda_300   s_d^inf(nu<12)   s_v^inf(nu<12)   windows [0,1) [1,2) [2,3) [3,4) [4,12)      (s_d with ref 2e4)')
        for D in (300.0, 1000.0, 3000.0, 6000.0):
            m, eta = C.mass(D)
            Tst = sum(R[j]*C.legP(nus - 1j*j, eta) for j in R)
            tot = Tinf*rr[D]/rinf; c = Tst[0]/tot[0]
            Gex = m*m*(c*tot - Tst)
            GK = Gex - 0.75*sum(R[j]*(nus - 1j*j)**2 for j in R)
            GK = GK - (GK[0]/NH[0])*NH
            sd, sv = C.sd(nus, GK, 'H0', 0.0, 12.0)
            win = [C.sd(nus, GK, 'H0', lo, hi)[0] for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 12))]
            res.append((D, sd))
            old = [r[2] for r in rows if r[0] == D][0]
            pr(f'  {D:6.0f}     {m/298.5:8.2f}        {sd:+.5f}         {sv:+.4f}       ' + ' '.join(f'{w:+.4f}' for w in win) + f'     ({old:+.5f})')
        Ds = np.array([r[0] for r in res]); ys = np.array([r[1] for r in res])
        A = np.vstack([np.ones_like(Ds), 1/Ds**2]).T
        coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
        pr(f'  fit s_d^inf(D) = {coef[0]:+.5f} + {coef[1]:.1f}/D^2  (max residual {np.max(np.abs(A @ coef - ys)):.1e});  '
           f'max - min over m/Lambda = 1..20: {np.ptp(ys):.1e}')
        # robustness: r_inf from (3000, 2e4) instead of (6000, 2e4)
        mb2b = C.mass(3000.0)[0]**2
        rinf2 = (ma2*rr[DREF] - mb2b*rr[3000.0])/(ma2 - mb2b)
        alt = []
        for D in (300.0, 1000.0, 6000.0):
            m, eta = C.mass(D)
            Tst = sum(R[j]*C.legP(nus - 1j*j, eta) for j in R)
            tot = Tinf*rr[D]/rinf2; c = Tst[0]/tot[0]
            GK = m*m*(c*tot - Tst) - 0.75*sum(R[j]*(nus - 1j*j)**2 for j in R)
            GK = GK - (GK[0]/NH[0])*NH
            alt.append(f'D={D:g}: {C.sd(nus, GK, "H0", 0.0, 12.0)[0]:+.5f}')
        pr('  robustness, r_inf from (3000, 2e4) instead: ' + ', '.join(alt) + f';  max |r_inf - r_inf\'| = {np.max(np.abs(rinf/rinf2 - 1)):.1e}')
        pr('  REPRO_ENGINE (full chain incl. engine K = 5 tail nu >= 12, +0.00159 / +0.00133): s_d = +0.50967 (D = 300), +0.50776 (D = 1000)')
    with open(os.path.join(HERE, 'cache', 'ee5_sd_report.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')
