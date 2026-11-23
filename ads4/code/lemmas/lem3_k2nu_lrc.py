r"""
lem3_k2nu_lrc (LEMMAS_CLOSE.md §3): K2nu and LRC for the spinning de-shifted kernels (Ĉ'_EE, Ĉ4^ds, dĈ4^ds).  Read-only imports:
resultq/v4c_wu.py (Kfun), resultq/stage4_hj.py (exact R_j), h7i14 via stage4_hj.

(a) K2nu, celestial factor.  C_nu(phi1, phi3) = A13^{(2+delta)/2} K_nu(dphi), K_nu = int dphib [A10 A30]^{-delta/2}, delta = 1/2 + i nu,
    A_ij = 2 - 2 cos(phi_ij) >= 0.  |A^{-delta/2}| = A^{-1/4} for every real nu  =>  |K_nu(dphi)| <= Kbar(dphi) := K_0-modulus,
    nu-independent (gamma = 0).  Kbar ~ 2 log(1/dphi) as dphi -> 0.  K_nu evaluated by a log-substitution trapezoid (any nu), checked against
    v4c_wu.Kfun (mpmath) at small nu and against V4DPRIME's corner expansion a + b dphi^{-2 i nu} (|a|, |b| ~ sqrt(pi/nu)).
(b) K2nu / LRC geometry in the H7 §18.1 Regge frame (exact).  y = |y| (1,0,0), x = |x| (cosh B, sinh B, 0), x2 = x/x^2 = -x4,
    x0 = |y| n_theta, n_theta = (0, cos th, sin th), z0 = (1, cos ph, sin ph); r = |x||y|, eta = cosh B, rho = r e^B.
    Every nu-dependent base keeps a fixed sign on the box rho <= c < 1, and its x0-correction relative to the x0 = 0 (leading-Regge) value
    is bounded by rho (NOT by r): X20 |x|^2 = -(1 - 2 r sinh B cos th - r^2), and |y z0.n|/(-z0.x/|x|) <= r e^B.
(c) f-factor (Stage 3 Lemma B): R = 16 r^2 (eta^2 (1+r^2)^2 - 4 r^2)/(s sb)^2 = [4 (rho+rhob)^2 (1 + rho rhob)^2 - 64 rho^2 rhob^2]/((1-rho^2)(1-rhob^2))^2
    in radial variables: nu-INDEPENDENT, corrections O(rho^2) relative to 16 r^2 eta^2.
(d) LRC structure against the stored exact R_j (stage4_hj): the heavy image is a DIAGONAL image of the position kernel eta^3 P_l:
    R_j(nu) = const * a_j(l) E_9(nu - i j) * s(nu - i j) / N_com(nu),  with a_j the Legendre-recurrence coefficients of eta^3 P_l,
    E_9 the plane-wave eigenvalue of v3_regge_kernel/stage12_rawkernel, and s = 1 (scalar) or the EE harmonic factor.  Tested: sc3, EE3, and
    k = 2 (sc1, EE1).  Also the growth in nu of the heavy-normalised coefficients |R_j| (polynomial: gamma_heavy = 0 once N_com is in chi).
Output: cache/lem3_k2nu_lrc.txt.  ~ 1-2 min.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import numpy as np
import mpmath as mp
import sympy as sp

OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


# ------------------------------------------------------------------ (a)
def K_trap(dphi, nu, power=None, nper=40, smin=-40.0):
    """K = int_0^{2pi} dphib [A(phib - h) A(phib + h)]^{-delta/2}, h = dphi/2 (power = -delta/2 by default; real power gives Kbar).
    Split at the singular points +-h; near each, u = +-e^s (log variable) removes the |u|^{-1/2} and resolves the |u|^{-i nu} oscillation."""
    d = 0.5 + 1j*nu
    p = -d/2 if power is None else power
    h = dphi/2
    A = lambda x: 4*np.sin(x/2)**2
    # singular points at pb = +-h (mod 2 pi); parametrise each half-interval by its distance u to the singular end
    total = 0
    for (lo, hi, dlo, dhi) in ((-h, h, (-2*h, 0.0), (0.0, 2*h)), (h, 2*np.pi - h, (0.0, 2*h), (2*np.pi - 2*h, 0.0))):
        half = 0.5*(hi - lo)
        smax = np.log(half)
        ds = min(0.05, 0.25/max(nu, 1e-9)) if nu > 0 else 0.05
        ds = min(ds, 1.0/nper)
        sgrid = np.arange(smin, smax, ds); sgrid = np.append(sgrid, smax)
        u = np.exp(sgrid)
        for (base, sg) in ((dlo, +1), (dhi, -1)):
            d1 = base[0] + sg*u; d2 = base[1] + sg*u          # pb - h, pb + h
            g = ((A(d1)*A(d2)).astype(complex))**p*u
            total = total + np.trapezoid(g, sgrid)
    return total


def part_a():
    pr('(a) K2nu, celestial factor: |K_nu(dphi)| <= Kbar(dphi) (nu-independent), and the nu -> inf decay')
    import v4c_wu as W
    for dphi in (0.3, 1.5):
        for nu in (0.5, 3.0):
            kt = K_trap(dphi, nu); km = W.Kfun(dphi, nu)
            pr(f'    check vs v4c_wu.Kfun: dphi = {dphi}, nu = {nu}: trapezoid {kt:.8f}, Kfun {km:.8f}, rel {abs(kt - km)/abs(km):.1e}')
    pr('    dphi     Kbar(dphi)   2log(1/dphi)   max_nu |K_nu|/Kbar over nu in {0.5,1,2,5,10,20,40,80}   |K_nu| at nu = 10 / 40 / 80')
    for dphi in (1e-3, 1e-2, 0.1, 0.5, 1.5, 3.0):
        kb = K_trap(dphi, 0.0, power=-0.25).real
        vals = {nu: abs(K_trap(dphi, nu)) for nu in (0.5, 1, 2, 5, 10, 20, 40, 80)}
        pr(f'    {dphi:7.0e}  {kb:10.5f}   {2*np.log(1/dphi):10.5f}   {max(vals.values())/kb:.4f}   '
           f'{vals[10]:.4e} / {vals[40]:.4e} / {vals[80]:.4e}')
    # corner expansion of V4DPRIME: a = 2 pi Gamma(1 - 2 delta)/Gamma(1 - delta)^2,  b = 2^{1+2 i nu}[B(1/2, 1-delta) + B(delta-1/2, 1-delta)]/2
    pr('    V4DPRIME corner expansion K ~ a + b dphi^{-2 i nu}: |a| sqrt(nu/pi), |b| sqrt(nu/pi), and trapezoid vs expansion at dphi = 1e-3')
    for nu in (2.0, 10.0, 40.0):
        dl = mp.mpf(0.5) + 1j*mp.mpf(nu)
        a = 2*mp.pi*mp.gamma(1 - 2*dl)/mp.gamma(1 - dl)**2
        b = 2**(1 + 2j*nu)*(mp.beta(0.5, 1 - dl) + mp.beta(dl - 0.5, 1 - dl))/2
        ex = complex(a + b*mp.mpf(1e-3)**(-2j*nu))
        kt = K_trap(1e-3, nu)
        pr(f'      nu = {nu}: |a| sqrt(nu/pi) = {float(abs(a))*np.sqrt(nu/np.pi):.4f}, |b| sqrt(nu/pi) = {float(abs(b))*np.sqrt(nu/np.pi):.4f}; '
           f'|K_trap - (a + b dphi^-2inu)| / |K| = {abs(kt - ex)/abs(kt):.1e}')


# ------------------------------------------------------------------ (b)
def mdot(u, v):
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]


def part_b():
    pr('(b) Regge-frame geometry: sign of the nu-dependent bases and size of their x0-corrections on the box rho = r e^B <= c')
    th = np.linspace(0, 2*np.pi, 721)[:, None]; ph = np.linspace(0, 2*np.pi, 721)[None, :]
    for (B, rho) in ((1.0, 0.3), (5.0, 0.3), (12.0, 0.3), (12.0, 0.9), (16.0, 0.5)):
        r = rho*np.exp(-B); xa = 1.0; ya = r/xa
        xh = np.array([np.cosh(B), np.sinh(B), 0.0]); yh = np.array([1.0, 0.0, 0.0])
        x2 = -xh/xa
        n = np.stack([0*th, np.cos(th), np.sin(th)], -1)[:, 0, :]          # (721, 3)
        x0 = ya*n
        x20 = x2[None, :] - x0; x40 = -x2[None, :] - x0
        X20 = mdot(x20, x20); X40 = mdot(x40, x40); X20_0 = mdot(x2, x2)
        z0 = np.stack([np.ones_like(ph), np.cos(ph), np.sin(ph)], -1)[0]     # (721, 3)
        zx20 = (z0[None, :, :]*np.array([-1, 1, 1])*x20[:, None, :]).sum(-1)  # z0 . x20
        zx20_0 = (z0*np.array([-1, 1, 1])*x2[None, :]).sum(-1)
        rel_X = np.max(np.abs(X20/X20_0 - 1)); rel_z = np.max(np.abs(zx20/zx20_0[None, :] - 1))
        pr(f'    B = {B:5.1f} (eta = {np.cosh(B):.2e}), rho = {rho}, r = {r:.2e}: sign(X20) constant: {bool(np.all(X20 < 0))}, '
           f'sign(z0.x20) constant: {bool(np.all(zx20 > 0) or np.all(zx20 < 0))}; max rel corr X20: {rel_X:.3f} (2r sinhB = {2*r*np.sinh(B):.3f}), '
           f'z0.x20: {rel_z:.3f} (r e^B = {rho:.3f}); r alone = {r:.1e}')


# ------------------------------------------------------------------ (c)
def part_c():
    pr('(c) f-factor (Stage 3 Lemma B) in radial variables: nu-independent, corrections O(rho^2)')
    r, eta, rho, rhob = sp.symbols('r eta rho rhob', positive=True)
    s = 1 - 2*r*eta + r**2; sb = 1 + 2*r*eta + r**2
    R = 16*r**2*(eta**2*(1 + r**2)**2 - 4*r**2)/(s*sb)**2
    Rrad = R.subs({r: sp.sqrt(rho*rhob), eta: (rho + rhob)/(2*sp.sqrt(rho*rhob))})
    target = (4*(rho + rhob)**2*(1 + rho*rhob)**2 - 64*rho**2*rhob**2)/((1 - rho**2)*(1 - rhob**2))**2
    pr('    R(r, eta) in (rho, rhob) minus closed form:', sp.simplify(Rrad - target))
    for (B, rho_) in ((3.0, 0.1), (10.0, 0.1), (10.0, 0.5), (20.0, 0.5)):
        rr = rho_*np.exp(-B); ee = np.cosh(B)
        Rv = 16*rr**2*(ee**2*(1 + rr**2)**2 - 4*rr**2)/((1 - 2*rr*ee + rr**2)*(1 + 2*rr*ee + rr**2))**2
        pr(f'      B = {B}, rho = {rho_}: R/(16 r^2 eta^2) - 1 = {Rv/(16*rr**2*ee**2) - 1:.4e};  rho^2 = {rho_**2:.4e};  r = {rr:.1e}')


# ------------------------------------------------------------------ (d)
def legendre_coeffs(n):
    mu = sp.Symbol('mu'); coeffs = {0: sp.Integer(1)}
    for _ in range(n):
        new = {}
        for j, c in coeffs.items():
            m = mu + j
            new[j + 1] = new.get(j + 1, 0) + c*(m + 1)/(2*m + 1)
            new[j - 1] = new.get(j - 1, 0) + c*m/(2*m + 1)
        coeffs = {j: sp.simplify(c) for j, c in new.items() if sp.simplify(c) != 0}
    return mu, coeffs


def part_d():
    pr('(d) LRC structure vs the stored exact R_j (resultq/stage4_hj): heavy image = diagonal image of the position kernel eta^n P_l')
    import stage4_hj as HJ
    D = HJ.load(); nu = HJ.nu_s
    mp.mp.dps = 30
    E = lambda b, v: 2**(2*b - 2)*mp.pi*mp.gamma((b - mp.mpf(1)/2 - 1j*v)/2)**2*mp.gamma((b - mp.mpf(1)/2 + 1j*v)/2)**2
    Ncom = lambda v: 2**11*mp.pi**2*abs(mp.gamma(mp.mpf(11)/4 + 1j*v/2))**4
    for key, n, b in (('sc1', 1, 7), ('EE1', 1, 7), ('sc3', 3, 9), ('EE3', 3, 9)):
        mu, a = legendre_coeffs(n)
        Rj = {j: sp.lambdify(nu, v, 'mpmath') for j, v in D[key].items()}
        rows = []
        for v in (mp.mpf('0.37'), mp.mpf('2.9'), mp.mpf('11.3'), mp.mpf('41.0')):
            l = -mp.mpf(1)/2 + 1j*v
            rat = {}
            for j in sorted(Rj):
                aj = sp.lambdify(mu, a[j], 'mpmath')(l)
                rat[j] = Rj[j](v)*Ncom(v)/(aj*E(b, v - 1j*j))
            rows.append((v, rat))
        # scalar: ratio constant in j and nu; EE: ratio / ratio_scalar = harmonic factor s(nu - i j)?
        if key.startswith('sc'):
            allr = [complex(r) for _, rat in rows for r in rat.values()]
            c0 = allr[0]
            pr(f'    {key}: R_j N_com/(a_j(l) E_{b}(nu - i j)) over j and nu in {{0.37, 2.9, 11.3, 41}}: const = {c0.real:.10e}, '
               f'max rel spread {max(abs(x/c0 - 1) for x in allr):.1e}')
            globals()['_c_' + key] = c0
        else:
            c0 = globals()['_c_sc' + key[-1]]
            dev = []
            for v, rat in rows:
                for j, r in rat.items():
                    sv = complex(r)/c0
                    w = v - 1j*j
                    cf = complex(32*(4*w**2 + 1)*(4*w**2 + 9)/((4*w**2 + 25)*(4*w**2 + 81)))
                    dev.append(abs(sv/cf - 1))
            pr(f'    {key}: test [R_j N_com/(a_j E_{b}(nu - i j))]/const_scalar == EE harmonic factor s(nu - i j), '
               f's(w) = 32(4w^2+1)(4w^2+9)/((4w^2+25)(4w^2+81)): max rel dev over j, nu = {max(dev):.1e}'
               + ('  -> diagonal (per-harmonic) structure HOLDS' if max(dev) < 1e-8 else '  -> NOT per-harmonic (rank-3 EE contraction mixes harmonics)'))
            if n == 3:
                for j in sorted(D[key]):
                    rr = sp.factor(sp.simplify(D[key][j]/D['sc3'][j]))
                    pr(f'      R_{j}^EE3 / R_{j}^sc3 = {rr}')
        # growth in nu of the heavy-normalised coefficients
        degs = {j: sp.degree(sp.numer(sp.together(v)), nu) - sp.degree(sp.denom(sp.together(v)), nu) for j, v in D[key].items()}
        pr(f'      heavy-normalised |R_j(nu)| ~ nu^p at large nu: p = {degs}  (polynomial: gamma_heavy = 0)')


if __name__ == '__main__':
    pr('# lem3_k2nu_lrc: K2nu and LRC for the spinning kernels, LEMMAS_CLOSE.md §3')
    part_a(); pr(f'  [{time.time() - T0:.0f}s]')
    part_b(); pr(f'  [{time.time() - T0:.0f}s]')
    part_c(); pr(f'  [{time.time() - T0:.0f}s]')
    part_d(); pr(f'  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'lem3_k2nu_lrc.txt'), 'w') as f:
        f.write('\n'.join(OUT) + '\n')
