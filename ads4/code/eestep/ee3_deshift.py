r"""
ee3_deshift (EESTEP_CLOSE.md §2): does the de-shift cancellation survive beyond leading Regge for the EE kernel itself?

By ee2 the exact EE per-nu kernel at fixed polarisations is  K^EE_nu = kappa(nu) r^2 E_nu,
    E_nu(rho, beta) = sum_{b=0}^{4} C_b(x) v^{b/2} Gh_b(rho, beta; nu),   C_0 = c_0, C_b = c_b + c_{-b},  (Psi^(b) = Psi^(-b))
with nu-INDEPENDENT C_b(x) (ee_quad.laurent_cb) and Gh_b the (0, b) continuous-spin blocks (ee_blocks; positive b only, no cancellation).
The leading-Regge kernel is L_nu(beta) = lim_{r -> 0, beta fixed} E_nu, extracted exactly as the r^0 Taylor coefficient by a Cauchy average
over the circle |r| = 0.3 e^{-beta} (E is analytic in r there; the c_b have poles up to r^-8 that cancel in the sum).

Test (as lrc2 for the scalar): S(rho, beta) = eta^n int_0^40 w(nu) E_nu dnu,  S_LR = eta^n int w L_nu,  S_raw = eta^n int |w E_nu|;
n = 1 (k = 2), n = 3 (k = 4); weights w2g, w2s, w4g, w4p of lrc2.  Local exponents d log|S - S_LR|/d beta at rho = 0.2, 0.5, 0.85.
Also: each b-block alone (E -> Gh_b, LR = P_l for every b), and the w4p residue formula for Gh_b:
    S_b - S_LR -> eta^3 Re[2 pi i Res_{nu = 3i/2} w c_l e^{l beta} F_b(rho; nu)].
Output: cache/ee3_deshift.txt (+ cache/ee3_deshift.npz).  ~10-20 min, one core.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import ee_blocks as EB
import ee_quad as Q

OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def W(name, nu):
    if name == 'w2g':
        return nu**2*np.exp(-nu**2/9)
    if name == 'w2s':
        return nu**2/np.cosh(np.pi*nu/2)
    if name == 'w4g':
        return nu**2*(nu**2 + 1)*(nu**2 + 4)*np.exp(-nu**2/9)
    if name == 'w4p':
        return nu**2*(nu**2 + 1)*(nu**2 + 4)*np.exp(-nu**2/9)/(nu**2 + 2.25)
    raise KeyError(name)


NAMES = ('w2g', 'w2s', 'w4g', 'w4p')
NPOW = {'w2g': 1, 'w2s': 1, 'w4g': 3, 'w4p': 3}
RHOS = np.array([0.2, 0.5, 0.85])
BETAS = np.arange(1.0, 17.0, 1.0)
BMAX = 4
NRAD = 260
MC = 32                                  # Cauchy points
rng = np.random.default_rng(23)
POLS = [tuple(rng.uniform(0, 2*np.pi, 4)) for _ in range(2)]


def Cb_of(r, beta, pol):
    cb, tail, odd, q = Q.laurent_cb(r, beta, pol)
    return np.array([cb[0]] + [cb[b] + cb[-b] for b in range(1, BMAX + 1)])


def main():
    xg, wg = np.polynomial.legendre.leggauss(700)
    numax = 40.0
    nus = 0.5*numax*(xg + 1); wq = 0.5*numax*wg
    Wv = {nm: W(nm, nus)*wq for nm in NAMES}
    nb, nr, npol = len(BETAS), len(RHOS), len(POLS)
    # nu-independent coefficients: real points and Cauchy circle
    C_real = np.zeros((npol, nb, nr, BMAX + 1)); v_real = np.zeros((nb, nr)); cancel = np.zeros((npol, nb, nr))
    phis = 2*np.pi*np.arange(MC)/MC
    C_circ = np.zeros((npol, nb, MC, BMAX + 1), dtype=complex); v_circ = np.zeros((nb, MC), dtype=complex); rcirc = np.zeros((nb, MC), dtype=complex)
    for ib, beta in enumerate(BETAS):
        for ir, rho in enumerate(RHOS):
            r = rho*np.exp(-beta); v_real[ib, ir] = EB.v_cross(rho, beta)
            for ip, pa in enumerate(POLS):
                C_real[ip, ib, ir] = Cb_of(r, beta, Q.pols(beta, *pa))
        R = 0.3*np.exp(-beta)
        for j, ph in enumerate(phis):
            rc = R*np.exp(1j*ph); rcirc[ib, j] = rc
            rho_c = rc*np.exp(beta); rb_c = rc*np.exp(-beta)
            v_circ[ib, j] = ((1 - rho_c)*(1 - rb_c)/((1 + rho_c)*(1 + rb_c)))**2
            for ip, pa in enumerate(POLS):
                C_circ[ip, ib, j] = Cb_of(rc, beta, Q.pols(beta, *pa))
    pr(f'  nu-independent coefficients done [{time.time() - T0:.0f}s]; sizes |C_b| at beta = 8: '
       + ', '.join(f'rho={rho}: {np.max(np.abs(C_real[0, 7, ir])):.1e}' for ir, rho in enumerate(RHOS)))
    Sx = {nm: np.zeros((npol, nb, nr)) for nm in NAMES}; Sl = {nm: np.zeros((npol, nb)) for nm in NAMES}
    Sa = {nm: np.zeros((npol, nb, nr)) for nm in NAMES}; Sfl = {nm: np.zeros((npol, nb, nr)) for nm in NAMES}
    Sb = {nm: np.zeros((BMAX + 1, nb, nr)) for nm in NAMES}; SP = {nm: np.zeros(nb) for nm in NAMES}
    imax = 0.0; lrim = 0.0
    for k, nu in enumerate(nus):
        coeffs = [EB.block_coeffs_b(nu, b, NRAD) for b in range(BMAX + 1)]
        for ib, beta in enumerate(BETAS):
            Pt = EB.legendre_scaled(nu, beta, NRAD + 1)
            Ys = [EB.radial_Y(B, m, Pt, beta) for (B, m) in coeffs]
            g = np.array([EB.Gh(Y, RHOS) for Y in Ys])                         # (BMAX+1, nr)
            imax = max(imax, np.max(np.abs(g.imag)/np.abs(g)))
            g = g.real
            gc = np.array([EB.Gh(Y, rcirc[ib]*np.exp(beta)) for Y in Ys])     # (BMAX+1, MC), complex r
            for nm in NAMES:
                Sb[nm][:, ib, :] += Wv[nm][k]*g
                SP[nm][ib] += Wv[nm][k]*Pt[0].real
            for ip in range(npol):
                vb = v_real[ib][None, :]**(0.5*np.arange(BMAX + 1)[:, None])    # (BMAX+1, nr)
                terms = C_real[ip, ib].T*vb*g                                    # (BMAX+1, nr)
                E = terms.sum(0); Eabs = np.abs(terms).sum(0)
                vcb = v_circ[ib][None, :]**(0.5*np.arange(BMAX + 1)[:, None])
                Lc = (C_circ[ip, ib].T*vcb*gc).sum(0)                            # E on the circle
                L = Lc.mean()
                lrim = max(lrim, abs(L.imag)/max(abs(L), 1e-300))
                for nm in NAMES:
                    Sx[nm][ip, ib] += Wv[nm][k]*E
                    Sl[nm][ip, ib] += Wv[nm][k]*L.real
                    Sa[nm][ip, ib] += np.abs(Wv[nm][k]*E)
                    Sfl[nm][ip, ib] += np.abs(Wv[nm][k])*Eabs
        if k % 100 == 99:
            pr(f'    nu node {k + 1}/700 [{time.time() - T0:.0f}s]')
    pr(f'  done: max |Im Gh_b|/|Gh_b| = {imax:.1e}; max |Im L|/|L| (LR extraction) = {lrim:.1e}')
    np.savez(os.path.join(HERE, 'cache', 'ee3_deshift.npz'), betas=BETAS, rhos=RHOS,
             **{f'Sx_{nm}': Sx[nm] for nm in NAMES}, **{f'Sl_{nm}': Sl[nm] for nm in NAMES}, **{f'Sa_{nm}': Sa[nm] for nm in NAMES},
             **{f'Sfl_{nm}': Sfl[nm] for nm in NAMES}, **{f'Sb_{nm}': Sb[nm] for nm in NAMES}, **{f'SP_{nm}': SP[nm] for nm in NAMES})
    return Sx, Sl, Sa, Sfl, Sb, SP


def slopes(y, betas):
    lg = np.log(np.abs(y))
    return np.diff(lg)/np.diff(betas)


def residue_pred(beta, rho, b, nm='w4p', nn=3):
    mp.mp.dps = 30
    nus = mp.mpc(0, 1.5)

    def g(nu):
        l = mp.mpf(-0.5) + 1j*nu
        w = nu**2*(nu**2 + 1)*(nu**2 + 4)*mp.exp(-nu**2/9)
        if nm == 'w4p':
            w = w/(nu**2 + mp.mpf(2.25))
        h = (mp.mpf(0.5) + 1j*nu)/2
        z = 4*mp.mpf(rho)/(1 + mp.mpf(rho))**2
        F = (z/(4*mp.mpf(rho)))**h*mp.hyp2f1(h, h + b, 2*h, z)
        return w*mp.gamma(l + mp.mpf(0.5))/(mp.sqrt(mp.pi)*mp.gamma(l + 1))*mp.exp(l*beta)*F

    rad = mp.mpf(0.2)
    res = mp.quad(lambda t: g(nus + rad*mp.exp(1j*t))*rad*mp.exp(1j*t), [0, mp.pi, 2*mp.pi])/(2*mp.pi)
    return float(mp.re(2j*mp.pi*res))*np.cosh(beta)**nn


def report(Sx, Sl, Sa, Sfl, Sb, SP):
    eta = np.cosh(BETAS)
    sel = slice(2, 13)
    pr('(1) each b-shifted block alone: local exponents of eta^n [int w Gh_b - int w P_l] (beta = 3.5..13.5); LR target -0.50, raw +0.5 / +2.5')
    for nm in NAMES:
        nn = NPOW[nm]
        for b in range(BMAX + 1):
            rows = []
            for ir, rho in enumerate(RHOS):
                d = (Sb[nm][b, :, ir] - SP[nm])*eta**nn
                rows.append(f'rho={rho}: ' + ' '.join(f'{x:+.2f}' for x in slopes(d, BETAS)[sel]))
            pr(f'    {nm}, b = {b}: ' + ' ; '.join(rows))
    pr('  residue check (w4p): (S_b - S_LR)/prediction eta^3 Re[2 pi i Res_{3i/2} w c_l e^{l beta} F_b] at beta = 6, 8, 10')
    for b in range(BMAX + 1):
        rows = []
        for ir, rho in enumerate(RHOS):
            vals = []
            for beta in (6.0, 8.0, 10.0):
                ib = list(BETAS).index(beta)
                d = (Sb['w4p'][b, ib, ir] - SP['w4p'][ib])*eta[ib]**3
                vals.append(d/residue_pred(beta, rho, b))
            rows.append(f'rho={rho}: ' + ', '.join(f'{x:.6f}' for x in vals))
        pr(f'    b = {b}: ' + ' ; '.join(rows))
    pr('  the same residue for the regular weight w4g (must vanish): ' + ', '.join(
        f'b={b}: {abs(residue_pred(8.0, 0.5, b, nm="w4g")):.0e}' for b in range(BMAX + 1)) + '  (rho = 0.5, beta = 8)')
    pr('(2) the EE kernel itself (two random polarisation sets): local exponents d log|S - S_LR|/d beta, beta = 3.5..13.5')
    for ip, pa in enumerate(POLS):
        pr(f'  polarisation set {ip}: angles (a1, a3, a2, a4) = {np.round(pa, 3).tolist()}')
        for nm in NAMES:
            nn = NPOW[nm]
            rows = []
            for ir, rho in enumerate(RHOS):
                d = (Sx[nm][ip, :, ir] - Sl[nm][ip])*eta**nn
                rows.append(f'rho={rho}: ' + ' '.join(f'{x:+.2f}' for x in slopes(d, BETAS)[sel]))
            pr(f'    {nm} (eta^{nn}): ' + ' ; '.join(rows))
            lr = slopes(Sl[nm][ip]*eta**nn, BETAS)[sel]
            raw = slopes(Sa[nm][ip, :, -1]*eta**nn, BETAS)[sel]
            pr(f'      LR target S_LR: ' + ' '.join(f'{x:+.2f}' for x in lr) + '   |  raw (rho = 0.85): ' + ' '.join(f'{x:+.2f}' for x in raw))
        # precision floor
        nm = 'w4g'
        fl = [1e-16*Sfl[nm][ip, ib, ir]*eta[ib]**3 for ib in range(len(BETAS)) for ir in range(len(RHOS))]
        dd = [abs(Sx[nm][ip, ib, ir] - Sl[nm][ip, ib])*eta[ib]**3 for ib in range(len(BETAS)) for ir in range(len(RHOS))]
        ok = [(BETAS[ib], RHOS[ir]) for ib in range(len(BETAS)) for ir in range(len(RHOS))
              if abs(Sx[nm][ip, ib, ir] - Sl[nm][ip, ib])*eta[ib]**3 < 30*1e-16*Sfl[nm][ip, ib, ir]*eta[ib]**3]
        pr(f'    double-precision floor (1e-16 x sum_b|terms| x eta^3, w4g) reached within x30 at (beta, rho): {ok[:12]}')
    pr('  table: (S - S_LR) eta^n and S_LR eta^n at rho = 0.5 (polarisation set 0)')
    for nm in NAMES:
        nn = NPOW[nm]
        pr(f'    {nm}: beta: ' + ' '.join(f'{b:5.0f}' for b in BETAS[2:14:2]))
        pr(f'          diff: ' + ' '.join(f'{x:+.2e}' for x in ((Sx[nm][0, :, 1] - Sl[nm][0])*eta**nn)[2:14:2]))
        pr(f'          LR  : ' + ' '.join(f'{x:+.2e}' for x in (Sl[nm][0]*eta**nn)[2:14:2]))
        pr(f'          raw : ' + ' '.join(f'{x:+.2e}' for x in (Sa[nm][0, :, 1]*eta**nn)[2:14:2]))


if __name__ == '__main__':
    pr('# ee3_deshift: de-shift cancellation beyond leading Regge for the EE kernel (EESTEP_CLOSE.md §2)')
    if len(sys.argv) > 1 and sys.argv[1] == 'report':
        z = np.load(os.path.join(HERE, 'cache', 'ee3_deshift.npz'))
        args = [{nm: z[f'{p}_{nm}'] for nm in NAMES} for p in ('Sx', 'Sl', 'Sa', 'Sfl', 'Sb', 'SP')]
    else:
        args = main()
    report(*args)
    pr(f'  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'ee3_deshift.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')
