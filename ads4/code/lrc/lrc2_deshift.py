r"""
lrc2_deshift (LRC_CLOSE.md §2-§3): does the de-shift cancellation survive the O(rho) corrections of the EXACT kernel?

By lrc1 the exact scalar j'=0 per-nu kernel is N(nu) * (nu-independent prefactor) * Gh_{1, l}(rho, beta), l = -1/2 + i nu, the d = 3
continuous-spin Casimir eigenfunction (lrc_blocks).  Here, at fixed rho in the Regge box and growing beta (eta = cosh beta, r = rho e^-beta):

(a) Edge (light-cone) structure: Gh - 2 Re[c_l e^{l beta} F(rho; h)] = O(e^{-2 beta}) relative, with F = 2F1(1/2, h; h+1/2; rho^2), h = (Dp+l)/2.
    F is the exact all-orders-in-rho correction of the top shell; its poles in nu sit at nu = i(Dp + 1/2 + 2p) and coincide with zeros of c_l.
(b) Pole table of the radial coefficients B[n, m](nu): zeros of K(Dp+n, l+m) = n(2Dp+n-3) + m(m + 2 i nu).
(c) Smeared kernels  S(rho, beta) = eta^n int_0^inf w(nu) Gh dnu  (n = 1: k = 2;  n = 3: k = 4; the exact f-factors are nu-independent and
    bounded relative to eta^n on the box [LEMMAS §3.3 (c), P]), against the leading-Regge S_LR = eta^n int w P_l and the raw size
    S_raw = eta^n int |w Gh|.  Weights (position normalisation, de-shift = strip analyticity of w c_{l+j}):
      w2g = nu^2 e^{-nu^2/9}                       k = 2, entire
      w2s = nu^2 sech(pi nu/2)                     k = 2, poles at +-i (like the real D weights: target ~ beta eta^-1/2)
      w4g = nu^2 (nu^2+1)(nu^2+4) e^{-nu^2/9}      k = 4, regular at nu = 3i/2
      w4p = w4g/(nu^2 + 9/4)                       k = 4, POLE at nu = +-3i/2 (generic weight: f_pos(3i/2) != 0)
(d) Prediction for w4p: S - S_LR -> eta^3 Re[2 pi i Res_{nu = 3i/2} w c_l e^{l beta} F(rho; h)]  ~ rho^2 eta^{+1}.
Output: cache/lrc2_deshift.txt.  ~ 1-2 min.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import lrc_blocks as L

OUT = []
T0 = time.time()
DP = 1.0


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def W(name, nu):
    nu = np.asarray(nu, dtype=complex) if np.iscomplexobj(nu) else np.asarray(nu, dtype=float)
    if name == 'w2g':
        return nu**2*np.exp(-nu**2/9)
    if name == 'w2s':
        return nu**2/np.cosh(np.pi*nu/2)
    if name == 'w4g':
        return nu**2*(nu**2 + 1)*(nu**2 + 4)*np.exp(-nu**2/9)
    if name == 'w4p':
        return nu**2*(nu**2 + 1)*(nu**2 + 4)*np.exp(-nu**2/9)/(nu**2 + 2.25)
    raise KeyError(name)


def part_a():
    pr('(a) edge structure of the exact kernel: Gh - 2Re[c_l e^{l beta} F(rho;h)] relative to |Gh|, at fixed rho and growing beta')
    N = 260
    for nu in (0.7, 3.1):
        B, m = L.block_coeffs(DP, nu, N)
        l = -0.5 + 1j*nu; cl = L.c_coef(mp.mpc(l))
        for rho in (0.3, 0.85):
            row = []
            Fv = L.edge_F(DP, nu, rho)
            for beta in (2.0, 4.0, 6.0, 8.0):
                Pt = L.legendre_scaled(nu, beta, N + 1)
                g = L.Gh(L.radial_Y(B, m, Pt, beta), rho)[0]
                ed = 2*(cl*np.exp(l*beta)*Fv).real
                row.append(f'b={beta:.0f}: {abs(g - ed)/abs(g):.2e} (e^-2b {np.exp(-2*beta):.1e})')
            pr(f'    nu = {nu}, rho = {rho}: |F - 1| = {abs(Fv - 1):.3f};  ' + '; '.join(row))


def part_b():
    pr('(b) poles in nu of the radial coefficients B[n, m] (Dp = 1): K(1+n, l+m) = n(n-1) + m(m + 2i nu) = 0  =>  nu = i[n(n-1)+m^2]/(2m)')
    ims = sorted(set(round((n*(n - 1) + mm*mm)/(2*abs(mm)), 4) for n in range(2, 41, 2) for mm in range(-n, n + 1, 2) if mm != 0))
    pr(f'    |Im nu| of all such poles for n <= 40: smallest {ims[:6]};  top shell (m = n) at nu = i(n - 1/2) = 1.5i, 3.5i, 5.5i, ...')
    pr('    top-shell resummation: F = 2F1(1/2, h; h+1/2; rho^2), poles at h + 1/2 = -p  <=>  nu = i(3/2 + 2p); there l = -2 - 2p and c_l = 0')
    pr('    (Legendre re-expansion denominators (2j+1)(2j+-1) give coefficient poles at nu = i k, k integer: they cancel between m and -m; the')
    pr('     edge function c_l e^{l beta} F has, besides the poles of c_l at nu = i, 2i, ... [these are the leading-Regge ones], none below 3.5i)')
    for p in range(3):
        nus = 1j*(1.5 + 2*p)
        l = -0.5 + 1j*nus
        pr(f'    check p = {p}: nu* = {nus}, l = {l.real:+.1f}, 1/Gamma(l+1) = {complex(mp.rgamma(l + 1)):.1e}')


def part_c():
    pr('(c) smeared kernels: exact (all orders in rho) vs leading Regge, k = 2 (eta^1) and k = 4 (eta^3)')
    N = 260
    xg, wg = np.polynomial.legendre.leggauss(700)
    numax = 40.0
    nus = 0.5*numax*(xg + 1); wq = 0.5*numax*wg
    betas = np.arange(1.0, 17.0, 1.0)
    rhos = np.array([0.2, 0.5, 0.85])
    names = ('w2g', 'w2s', 'w4g', 'w4p')
    Wv = {nm: W(nm, nus)*wq for nm in names}
    Sx = {nm: np.zeros((len(betas), len(rhos))) for nm in names}
    Sl = {nm: np.zeros(len(betas)) for nm in names}
    Sa = {nm: np.zeros((len(betas), len(rhos))) for nm in names}
    imax = 0.0
    for k, nu in enumerate(nus):
        B, m = L.block_coeffs(DP, nu, N)
        for ib, beta in enumerate(betas):
            Pt = L.legendre_scaled(nu, beta, N + 1)
            g = L.Gh(L.radial_Y(B, m, Pt, beta), rhos)
            imax = max(imax, np.max(np.abs(g.imag)/np.abs(g)))
            for nm in names:
                Sx[nm][ib] += Wv[nm][k]*g.real
                Sl[nm][ib] += Wv[nm][k]*Pt[0].real
                Sa[nm][ib] += np.abs(Wv[nm][k]*g)
    pr(f'    nu-grid: 700 Gauss-Legendre nodes on [0, 40]; radial order N = {N}; max |Im Gh|/|Gh| = {imax:.1e} (Gh real, as it must be)')
    out = {}
    for nm in names:
        nn = 1 if nm.startswith('w2') else 3
        eta = np.cosh(betas)
        pr(f'  weight {nm} (prefactor eta^{nn}):')
        pr('    beta |  S_LR        | ' + ' | '.join(f'S_ex - S_LR (rho={r})' for r in rhos) + ' | S_raw (rho=0.85) | eta^-1/2   | floor 1e-16*S_raw')
        for ib, beta in enumerate(betas):
            d = (Sx[nm][ib] - Sl[nm][ib])*eta[ib]**nn
            pr(f'    {beta:4.0f} | {Sl[nm][ib]*eta[ib]**nn:+.4e} | ' + ' | '.join(f'{x:+.4e}' for x in d)
               + f' | {Sa[nm][ib][-1]*eta[ib]**nn:.3e}       | {eta[ib]**-0.5:.3e}  | {1e-16*Sa[nm][ib][-1]*eta[ib]**nn:.1e}')
        out[nm] = (betas, (Sx[nm] - Sl[nm][:, None])*(np.cosh(betas)**nn)[:, None], Sl[nm]*np.cosh(betas)**nn, Sa[nm]*(np.cosh(betas)**nn)[:, None])
    return out, rhos


def residue_pred(beta, rho, nm='w4p', nn=3):
    """eta^n Re[2 pi i Res_{nu = 3i/2} w c_l e^{l beta} F]."""
    mp.mp.dps = 30
    nus = mp.mpc(0, 1.5)

    def g(nu):
        l = mp.mpf(-0.5) + 1j*nu
        w = nu**2*(nu**2 + 1)*(nu**2 + 4)*mp.exp(-nu**2/9)
        if nm == 'w4p':
            w = w/(nu**2 + mp.mpf(2.25))
        h = (DP - mp.mpf(0.5) + 1j*nu)/2
        F = mp.hyp2f1(mp.mpf(0.5), h, h + mp.mpf(0.5), mp.mpf(rho)**2)
        return w*mp.gamma(l + mp.mpf(0.5))/(mp.sqrt(mp.pi)*mp.gamma(l + 1))*mp.exp(l*beta)*F

    rad = mp.mpf(0.2)
    res = mp.quad(lambda t: g(nus + rad*mp.exp(1j*t))*rad*mp.exp(1j*t), [0, mp.pi, 2*mp.pi])/(2*mp.pi)
    return float(mp.re(2j*mp.pi*res))*np.cosh(beta)**nn


def part_d(out, rhos):
    pr('(d) local growth exponents d log|S_ex - S_LR| / d beta  (raw rate: k = 2: +0.5, k = 4: +2.5; de-shifted target: -0.5)')
    for nm in ('w2g', 'w2s', 'w4g', 'w4p'):
        betas, D, SL, SA = out[nm]
        rows = []
        for ir, rho in enumerate(rhos):
            lg = np.log(np.abs(D[:, ir]))
            sl = np.diff(lg)/np.diff(betas)
            rows.append(f'rho={rho}: ' + ' '.join(f'{x:+.2f}' for x in sl[2:13]))
        pr(f'    {nm} (beta = 3.5 ... 13.5): ' + ' ;  '.join(rows))
    pr('    LR target S_LR: local exponents ' + ' '.join(
        f'{nm}: ' + ' '.join(f'{x:+.2f}' for x in (np.diff(np.log(np.abs(out[nm][2])))/np.diff(out[nm][0]))[2:10]) for nm in ('w2g', 'w4g')))
    pr('  prediction for w4p: residue at nu = 3i/2 (pole of F times pole of w times zero of c_l)')
    betas, D, SL, SA = out['w4p']
    for ir, rho in enumerate(rhos):
        rows = []
        for beta in (4.0, 6.0, 8.0, 10.0, 12.0):
            ib = list(betas).index(beta)
            pred = residue_pred(beta, rho)
            rows.append(f'b={beta:.0f}: meas/pred = {D[ib, ir]/pred:.6f}')
        pr(f'    rho = {rho}: ' + '; '.join(rows))
    pr('  same residue for the regular weight w4g (must be ~0):  ' + ', '.join(
        f'rho={rho}: {residue_pred(8.0, rho, nm="w4g"):.1e}' for rho in rhos) + '  (at beta = 8)')


if __name__ == '__main__':
    pr('# lrc2_deshift: survival of the de-shift cancellation beyond leading Regge (LRC_CLOSE.md §2-§3), Dp = 1 (lrc1)')
    part_a(); pr(f'  [{time.time() - T0:.0f}s]')
    part_b(); pr(f'  [{time.time() - T0:.0f}s]')
    out, rhos = part_c(); pr(f'  [{time.time() - T0:.0f}s]')
    part_d(out, rhos); pr(f'  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'lrc2_deshift.txt'), 'w') as f:
        f.write('\n'.join(OUT) + '\n')
