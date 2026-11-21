"""
v4d_wu (RESULT_Q_V4DPRIME.md §1): the u-channel identity corner contact term with the ACTUAL de-shift weights chi(nu) (Stage 4 stored,
Stage 7 XP; value and derivative pieces), heavy-normalised (weight on the v4c position kernel w = chi N_com(0)/N_com(nu)).
Usage: python3 v4d_wu.py  -> cache/v4d_wu.txt   (~1 min, 1 core)
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import v4d_core as D

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
L = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)


def main():
    t0 = time.time()
    mp.mp.dps = 30
    pr('# v4d_wu: u-channel identity contact term at the light-cone-edge/Regge corner with the actual C4 de-shift weights')
    # ---------------------------------------------------------------- (0) checks
    pr('(0) checks')
    import v4c_wu as W
    W.T42 = W.t42()
    mx = 0
    for eps in (1e-3, 1e-5):
        for dph in (0.3*np.sqrt(eps), np.sqrt(eps), 5*np.sqrt(eps), 0.2, 2.0):
            u = W.U_eps(np.array([dph]), eps)[0]
            mx = max(mx, abs(u/D.U_exact(dph, eps) - 1))
    pr(f'  closed form U = (pi^2/64) T42 eps^4 (1+c)^2 (1-X^2)/X^6 vs v4c Cauchy-FFT residue: max rel. dev. {mx:.1e} (FFT precision)')
    mx = 0
    for nu in (0.5, 1.0, 3.0):
        a, b = D.ab(nu)
        for dph in (1e-2, 1e-3):
            K = W.Kfun(dph, nu)
            mx = max(mx, abs(K - complex(a + b*mp.mpf(dph)**(-2j*nu)))/abs(K))
    pr(f'  K_nu(dphi) = a + b dphi^(-2 i nu) + O(dphi^2) vs v4c Kfun (nu = .5, 1, 3; dphi = 1e-2, 1e-3): max rel. dev. {mx:.1e}')
    for nu, ref in ((1e-9, -1.516347e-2), (1.0, -(4.3518 + 3.9975j)*1e-3)):
        v = D.pair_exact(nu, lambda x: D.RtEE(np.cos(x)), N=2**14)
        pr(f'  I_t(nu = {nu:.0f}) by exact Fourier pairing {v:.6e}  vs v4c trapezoid {ref:.6e}  (rel {abs(v - ref)/abs(ref):.1e}, v4c grid error)')
    pr('  exact (Fourier) G_nu(eps) vs leading small-eps closed form, tube x1 -> x1 + i eps e:')
    for nu in (1.0, 3.0, 8.0, 14.0):
        row = []
        for eps in (1e-4, 1e-5, 1e-6):
            v = D.pair_exact(nu, lambda x: D.U_exact(x, eps), N=2**17 if eps < 3e-6 else 2**15)
            l = D.G_lead(nu, eps)
            row.append(f'eps={eps:.0e}: |G|={abs(v):.3e} rel.dev {abs(v - l)/abs(l):.1e}')
        pr(f'    nu = {nu:5.1f}: ' + '; '.join(row))
    # ---------------------------------------------------------------- (1) nu-structure
    pr('(1) per-nu structure: G_nu(eps) = (pi^2/16) T42 eps^{-1/4} [c+ eps^{i nu/2} + c- eps^{-i nu/2}];  heavy normalisation N_com(0)/N_com(nu)')
    Wt = D.weights()
    pr('    nu    |c+| (tube +)  |c-| (tube +)   N(0)/N(nu)    |c+|N0/N     |chi_xp_value|  |chi c+ N0/N| (xp value)')
    lN0 = D.logNcom(0.0)
    nuX, chiX, dX = Wt['xp_value']
    for nu in (0.5, 1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 60, 80, 100, 120):
        cp, cm = D.coef_lead(nu, +1)
        h = np.exp(lN0 - D.logNcom(nu))
        k = np.argmin(np.abs(nuX - nu)); sl = slice(max(0, k - 12), k + 13)
        ch = np.abs(chiX[sl]).max()
        pr(f'   {nu:5.1f}  {abs(cp):.3e}   {abs(cm):.3e}   {h:.3e}   {abs(cp)*h:.3e}   {ch:.3e}   {abs(cp)*h*ch:.3e}')
    pr('    (tube -: c+ <-> c- magnitudes exchanged; the dominant branch always decays only like nu^{9/2} e^{-pi nu/4}, the heavy'
       ' normalisation grows like nu^{-9} e^{pi nu}, the weights decay like e^{-(pi/2) nu} nu^{~18}: the product grows like e^{pi nu/4})')
    # local exponential rates
    for nm in ('stored_value', 'xp_value', 'xp_deriv'):
        nu, chi, dnu = Wt[nm]
        rr = []
        for a, b in ((20, 30), (30, 40), (50, 60), (80, 90), (100, 110)):
            if b > nu[-1]:
                continue
            ia = np.argmin(np.abs(nu - a)); ib = np.argmin(np.abs(nu - b))
            ma = np.abs(chi[max(0, ia - 12):ia + 13]).max(); mb = np.abs(chi[max(0, ib - 12):ib + 13]).max()
            rr.append(f'[{a},{b}]: {np.log(ma/mb)/(b - a):.3f}')
        pr(f'    decay rate -dlog|chi|/dnu ({nm}): ' + ', '.join(rr) + '   (pi/2 = 1.571, 3 pi/4 = 2.356 needed for convergence)')
    # chi near nu = 0
    for nm in Wt:
        nu, chi, dnu = Wt[nm]
        pr(f'    chi/nu^2 at the first nodes ({nm}): ' + ' '.join(f'{chi[k]/nu[k]**2:+.4e}' for k in range(4)) +
           f'   chi/nu: ' + ' '.join(f'{chi[k]/nu[k]:+.2e}' for k in range(3)))
    # ---------------------------------------------------------------- (2) partial sums
    pr('(2) de-shifted contact term  G_ds(eps; nu_c) = sum_{nu_k <= nu_c} dnu w_k G_{nu_k}(eps),  w = chi N_com(0)/N_com(nu)  [unit: v4c kernel at nu = 0]')
    pr('    reference sizes (same unit): raw per-nu |G_0(eps)| = 0.0143 / 0.062 / 0.125 at eps = 1e-2 / 1e-4 / 1e-5;  |I_t(0)| = 0.0151')
    coefs = {}
    for nm in Wt:
        nu, chi, dnu = Wt[nm]
        for sgn in (+1, -1):
            cc = np.array([D.coef_lead(v, sgn) for v in nu])
            coefs[(nm, sgn)] = cc
    EPS = (1e-2, 1e-4, 1e-6, 1e-8, 1e-12)
    for nm in Wt:
        nu, chi, dnu = Wt[nm]
        w = chi*np.exp(lN0 - D.logNcom(nu))
        for sgn in (+1, -1):
            cc = coefs[(nm, sgn)]
            pr(f'  {nm}, tube {"+" if sgn > 0 else "-"}:   |G_ds| at eps = ' + ' / '.join(f'{e:.0e}' for e in EPS))
            for nuc in (2, 4, 6, 10, 15, 20, 30, 41, 60, 90, 123):
                if nuc > nu[-1] + dnu:
                    continue
                m = nu <= nuc
                vals = []
                for e in EPS:
                    G = np.pi**2/16*D.T42*e**(-0.25)*(cc[m, 0]*e**(0.5j*nu[m]) + cc[m, 1]*e**(-0.5j*nu[m]))
                    vals.append(abs(dnu*np.sum(w[m]*G)))
                pr(f'      nu_c = {nuc:5.1f}: ' + '  '.join(f'{v:.3e}' for v in vals))
    # ---------------------------------------------------------------- (3) eps-scaling of the low-nu part
    pr('(3) eps-scaling of the truncated sums (tube +): local exponent p in |G_ds| ~ eps^p between successive eps (decades)')
    epsg = 10.0**-np.arange(2, 25, 1.0)
    for nm in ('stored_value', 'stored_deriv', 'xp_value', 'xp_deriv'):
        nu, chi, dnu = Wt[nm]
        w = chi*np.exp(lN0 - D.logNcom(nu)); cc = coefs[(nm, +1)]
        for nuc in (6, 10, 20):
            m = nu <= nuc
            v = np.array([abs(dnu*np.sum(w[m]*np.pi**2/16*D.T42*e**(-0.25)*(cc[m, 0]*e**(0.5j*nu[m]) + cc[m, 1]*e**(-0.5j*nu[m]))))
                          for e in epsg])
            p = np.diff(np.log(v))/np.diff(np.log(epsg))
            pr(f'    {nm:13s} nu_c = {nuc:3d}: |G_ds| = ' + ' '.join(f'{x:.1e}' for x in v[::3]) + '  (eps = 1e-2, 1e-5, ...);  p (running mean over'
               f' 4 decades) = ' + ' '.join(f'{np.mean(p[i:i + 4]):+.2f}' for i in range(0, len(p) - 3, 4)))
    # ---------------------------------------------------------------- (4) per-nu raw vs the claimed cure; analytic structure
    pr('(4) analytic structure of the dominant branch c+(nu) = a(nu) J(5/2 + i nu) in the lower strip (the contour shift of V4prime §3):')
    for z in (-0.5j, -1j, -1.5j, -2j, -2.5j, -3j):
        vals = []
        for h in (mp.mpf('1e-8'), mp.mpf('1e-9')):
            d = mp.mpf(1)/2 + 1j*(mp.mpc(z) + h)
            vals.append(abs(2*mp.pi*mp.gamma(1 - 2*d)*mp.rgamma(1 - d)**2))
        kind = 'POLE' if vals[1] > 5*vals[0] else ('zero' if vals[1] < 0.2*vals[0] else f'regular |a| = {float(vals[1]):.3e}')
        ex = -0.25 + (0.5j*complex(z)).real
        pr(f'    nu = {complex(z)}: a(nu) {kind};  a residue here would contribute eps^(-1/4) eps^(i nu/2) = eps^({ex:+.2f})')
    pr('    a(nu) = 2 pi Gamma(1-2delta)/Gamma(1-delta)^2: simple poles at nu = 0, -i, -2i, -3i (delta = 1/2, 3/2, 5/2, 7/2), zeros at'
       ' nu = -i/2, -3i/2, -5i/2.  A contour shift by 3 units therefore yields at best eps^{+1/4} (residue at nu = -i), not eps^{5/4},'
       ' unless w(-i) = 0; and it is not available at all, because w c+ grows like e^{pi nu/4} on the real axis.')
    # ---------------------------------------------------------------- (4b) smooth (Gaussian) taper: order of limits
    pr('(4b) Gaussian taper w -> w exp(-(nu/Lam)^2) (XP value, tube +): |G_ds| vs eps and Lam (non-uniform limits)')
    nu, chi, dnu = Wt['xp_value']
    w = chi*np.exp(lN0 - D.logNcom(nu)); cc = coefs[('xp_value', +1)]
    for Lam in (3, 6, 10, 20, 40):
        tap = np.exp(-(nu/Lam)**2)
        vals = []
        for e in (1e-2, 1e-4, 1e-6, 1e-10, 1e-20):
            G = np.pi**2/16*D.T42*e**(-0.25)*(cc[:, 0]*e**(0.5j*nu) + cc[:, 1]*e**(-0.5j*nu))
            vals.append(abs(dnu*np.sum(w*tap*G)))
        pr(f'    Lam = {Lam:3d}: ' + '  '.join(f'{v:.2e}' for v in vals) + '   (eps = 1e-2, 1e-4, 1e-6, 1e-10, 1e-20)')
    # ---------------------------------------------------------------- (4d) the nu = 0 junction (closed form)
    pr('(4d) nu = 0 junction.  b(nu) = a(-nu) (1 + i sinh(pi nu))/cosh(pi nu) exactly (duplication + reflection formulas); check:')
    for v in (0.3, 1.0, 4.0):
        a, b = D.ab(v); am, bm = D.ab(-v)
        pr(f'    nu = {v}: b/a(-nu) = {complex(b/am):.12f},  (1 + i sinh pi nu)/cosh pi nu = {complex((1 + 1j*mp.sinh(mp.pi*v))/mp.cosh(mp.pi*v)):.12f}')
    pr('    near nu = 0: a = i/nu + 4 log 2 + O(nu), b = -i/nu + 4 log 2 + pi + O(nu): a0 - b0 = -pi.  The two branches of the nu > 0 integral'
       ' do not continue into each other, so the endpoint nu = 0 survives any contour shift: with chi = chi2 nu^2 + ... (checked above),')
    pr('    G_ds(eps) -> (pi^2/16) T42 eps^{-1/4} * 2 pi i chi2 J(5/2) (-1)/L^3 * (1 + O(1/L)),  L = log(1/eps)/2  (divergent).')
    J52 = complex(D.Jfun(mp.mpf(5)/2, +1))
    nu, chi, dnu = Wt['xp_value']
    chi2 = chi[0]/nu[0]**2
    w = chi*np.exp(lN0 - D.logNcom(nu)); m = nu < 12
    cc = np.array([D.coef_lead(v, +1) for v in nu[m]]); tap = np.exp(-(nu[m]/3)**2)
    for Lv in (10, 20, 40, 60):
        sdir = dnu*np.sum(w[m]*tap*(cc[:, 0]*np.exp(-1j*nu[m]*Lv) + cc[:, 1]*np.exp(1j*nu[m]*Lv)))
        sas = 2j*np.pi*chi2*J52*(-1)/Lv**3
        pr(f'    L = {Lv:3d} (eps = e^-{2*Lv}): bracket (taper Lam = 3, xp value) = {abs(sdir):.3e};  junction asymptote |2 pi chi2 J(5/2)|/L^3 = {abs(sas):.3e}  (ratio {abs(sdir/sas):.2f})')
    # ---------------------------------------------------------------- (4c) ratio to an O(1) light family under the same truncated weight
    pr('(4c) contact term relative to an O(1) regular light family under the same (truncated) weight: R = |G_ds|/sum dnu |w||I_t|')
    for nm in ('stored_value', 'xp_value', 'xp_deriv'):
        nu, chi, dnu = Wt[nm]
        w = chi*np.exp(lN0 - D.logNcom(nu)); cc = coefs[(nm, +1)]
        for nuc in (6, 10, 20):
            m = nu <= nuc
            sel = np.where(m)[0][::max(1, int(0.2/dnu))]
            It = np.array([abs(D.pair_exact(float(nu[k]), lambda x: D.RtEE(np.cos(x)), N=2**12)) for k in sel])
            ref = np.sum(np.abs(w[sel])*It)*dnu*max(1, int(0.2/dnu))
            rr = []
            for e in (1e-2, 1e-6, 1e-12):
                G = np.pi**2/16*D.T42*e**(-0.25)*(cc[m, 0]*e**(0.5j*nu[m]) + cc[m, 1]*e**(-0.5j*nu[m]))
                rr.append(abs(dnu*np.sum(w[m]*G))/ref)
            pr(f'    {nm:13s} nu_c = {nuc:3d}: sum |w||I_t| = {ref:.2e};  R = ' + ' / '.join(f'{x:.2e}' for x in rr) + '  (eps = 1e-2 / 1e-6 / 1e-12)')
    # ---------------------------------------------------------------- (5) regular family I_t in heavy normalisation
    pr('(5) regular light family I_t(nu) (J-resummed [TT]_0, exact Fourier pairing): size and heavy-normalised size')
    for nu in (1, 3, 5, 8, 12, 16, 20, 25, 30):
        v = D.pair_exact(float(nu), lambda x: D.RtEE(np.cos(x)), N=2**14)
        h = np.exp(lN0 - D.logNcom(nu))
        pr(f'    nu = {nu:3d}: |I_t| = {abs(v):.3e}   |I_t| N0/N = {abs(v)*h:.3e}')
    pr('    (I_s(nu) = -I_t(nu) pointwise in the angles (V4prime §2), so the family SUM vanishes for every nu; each family alone is not'
       ' integrable against chi in heavy normalisation)')
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(CACHE, 'v4d_wu.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
