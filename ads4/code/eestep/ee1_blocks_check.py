r"""
ee1_blocks_check (EESTEP_CLOSE.md §1.2): checks of the b-shifted continuous-spin recursion of ee_blocks.py.

(a) b = 0 reproduces lrc/lrc_blocks.block_coeffs (imported read-only) and its Gh.
(b) The Dolan-Osborn Casimir PDE with (a, b) = (0, b), on g = (4r)^Dp Gh_b, by 4th-order finite differences in (z, zb) at Lorentzian
    Regge-box points (0 < rhob < rho < 1): residual |(DO - C) g| / |C g|.
(c) Polynomiality in b: B^(b)[n, m] is a polynomial in b of degree <= n (finite differences of order n+1 over b = -6..6 vanish).
(d) Pole content in 0 < Im nu <= 3: max_{n <= 8, m} |B^(b)[n, m](nu)| at nu = i y + eps for y in {1/2, 3/2, 5/2}, eps -> 0 (a genuine pole
    grows like 1/eps).  Prediction (a = 0): only y = 3/2 (n = 2, m = 2 and descendants).
(e) Edge (light-cone) structure: Gh_b - 2 Re[c_l e^{l beta} F_b(rho)] = O(e^{-2 beta}) relative, F_b = (z/4rho)^h 2F1(h, h+b; 2h; z).
Output: cache/ee1_blocks_check.txt.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'lrc'))
import numpy as np
import mpmath as mp
import ee_blocks as EB

OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def part_a():
    import lrc_blocks as L
    pr('(a) b = 0 against lrc_blocks (read-only import)')
    for nu in (0.4, 2.7):
        N = 60
        B, m = EB.block_coeffs_b(nu, 0, N)
        B0, m0 = L.block_coeffs(1.0, nu, N)
        dev = 0.0
        for lev in range(N//2 + 1):
            n = 2*lev
            for i, mm in enumerate(m0):
                if abs(mm) <= n:
                    dev = max(dev, abs(B[n, mm + N] - B0[lev, i])/max(1e-300, np.max(np.abs(B0[lev]))))
        odd = np.max(np.abs(B[1::2]))
        g1 = EB.Gh_point(nu, 0, 0.6, 1.7, N=160)
        Pt = L.legendre_scaled(nu, 1.7, 161); Y = L.radial_Y(*L.block_coeffs(1.0, nu, 160), Pt, 1.7); g0 = L.Gh(Y, 0.6)[0]
        pr(f'    nu = {nu}: max rel diff of B (even levels) = {dev:.1e}; max |B| on odd levels = {odd:.1e}; Gh(0.6, 1.7): {g1:.12f} vs {g0:.12f}')


def g_of_zzb(nu, b, z, zb, cache):
    rho = z/(1 + np.sqrt(1 - z))**2; rhob = zb/(1 + np.sqrt(1 - zb))**2
    r = np.sqrt(rho*rhob); beta = 0.5*np.log(rho/rhob)
    key = (nu, b)
    if key not in cache:
        cache[key] = EB.block_coeffs_b(nu, b, 120)
    B, m = cache[key]
    Pt = EB.legendre_scaled(nu, beta, 121)
    return (4*r)*EB.Gh(EB.radial_Y(B, m, Pt, beta), rho)[0]


def part_b():
    pr('(b) Dolan-Osborn Casimir PDE, (a, b) = (0, b), Dp = 1: residual |(DO - C) g|/|C g| (finite differences, h = 2e-3 z)')
    cache = {}
    for nu in (0.6, 2.3):
        for b in (-2, -1, 1, 2, 3):
            res = []
            for (z, zb) in ((0.45, 0.08), (0.7, 0.2), (0.3, 0.02)):
                hz, hzb = 2e-3*z, 2e-3*zb
                G = {}
                for i in range(-2, 3):
                    for j in range(-2, 3):
                        if i*j != 0 and not (abs(i) == 1 and abs(j) == 1):
                            continue
                        G[(i, j)] = g_of_zzb(nu, b, z + i*hz, zb + j*hzb, cache)
                c4 = np.array([1, -8, 0, 8, -1])/12.0; c2 = np.array([-1, 16, -30, 16, -1])/12.0
                gz = sum(c4[i + 2]*G[(i, 0)] for i in range(-2, 3))/hz
                gzz = sum(c2[i + 2]*G[(i, 0)] for i in range(-2, 3))/hz**2
                gb = sum(c4[j + 2]*G[(0, j)] for j in range(-2, 3))/hzb
                gbb = sum(c2[j + 2]*G[(0, j)] for j in range(-2, 3))/hzb**2
                g = G[(0, 0)]
                Dz = z**2*(1 - z)*gzz - (b + 1)*z**2*gz
                Dzb = zb**2*(1 - zb)*gbb - (b + 1)*zb**2*gb
                lhs = 2*Dz + 2*Dzb + 2*z*zb/(z - zb)*((1 - z)*gz - (1 - zb)*gb)
                l = -0.5 + 1j*nu
                C = 1*(1 - 3) + l*(l + 1)
                # control: same g with the b = 0 operator
                lhs0 = lhs + 2*b*(z**2*gz + zb**2*gb)
                res.append((abs(lhs - C*g)/abs(C*g), abs(lhs0 - C*g)/abs(C*g)))
            pr(f'    nu = {nu}, b = {b:+d}: residual ' + ', '.join(f'{a:.1e}' for a, _ in res)
               + '   [control, b = 0 operator on the same g: ' + ', '.join(f'{c:.1e}' for _, c in res) + ']')


def part_c():
    pr('(c) polynomiality in b: max over m of |Delta_b^{n+1} B^(b)[n, m]| / max|B^(b)[n, m]|  (b = -6..6, unit steps)')
    for nu in (0.37, 3.1):
        bs = np.arange(-6, 7)
        N = 10
        Bs = np.array([EB.block_coeffs_b(nu, int(b), N)[0] for b in bs])   # (13, N+1, 2N+1)
        row = []
        for n in range(1, 9):
            d = np.diff(Bs[:, n, :], n=n + 1, axis=0)
            dlow = np.diff(Bs[:, n, :], n=n, axis=0)
            row.append(f'n={n}: {np.max(np.abs(d))/np.max(np.abs(Bs[:, n, :])):.0e} (order n: {np.max(np.abs(dlow))/np.max(np.abs(Bs[:, n, :])):.0e})')
        pr(f'    nu = {nu}: ' + '; '.join(row))


def part_d():
    pr('(d) poles of B^(b)[n, m](nu), n <= 8, in 0 < Im nu <= 3: max|B| at nu = i y + eps (a pole shows |B| ~ 1/eps)')
    for b in (0, 1, -2, 3):
        rows = []
        for y in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
            vals = []
            for eps in (1e-3, 1e-6):
                B, m = EB.block_coeffs_b(1j*y + eps, b, 8)
                vals.append(np.max(np.abs(B)))
            rows.append(f'y={y}: {vals[0]:.1e} -> {vals[1]:.1e}')
        pr(f'    b = {b:+d}: ' + '; '.join(rows))
    pr('    (each candidate y is a zero of some K(1+n, l+m) = n(n-1) + m(m + 2 i nu); the growth x1000 marks a surviving pole)')


def part_e():
    pr('(e) edge structure: |Gh_b - 2 Re[c_l e^{l beta} F_b(rho)]| / |2 c_l e^{l beta} F_b| (envelope) at fixed rho, growing beta  (expect ~ e^{-2 beta})')
    N = 220
    for nu in (0.7, 3.1):
        l = -0.5 + 1j*nu; cl = EB.c_coef(mp.mpc(l))
        for b in (-3, 1, 4):
            B, m = EB.block_coeffs_b(nu, b, N)
            for rho in (0.3, 0.85):
                Fv = EB.edge_Fb(nu, b, rho)
                row = []
                for beta in (2.0, 4.0, 6.0, 8.0, 10.0):
                    Pt = EB.legendre_scaled(nu, beta, N + 1)
                    g = EB.Gh(EB.radial_Y(B, m, Pt, beta), rho)[0]
                    ed = 2*(cl*np.exp(l*beta)*Fv).real
                    row.append(f'b={beta:.0f}: {abs(g - ed)/abs(2*cl*np.exp(l*beta)*Fv):.1e}')
                pr(f'    nu = {nu}, b = {b:+d}, rho = {rho}: |F_b - 1| = {abs(Fv - 1):.3f}, Im Gh/|Gh| = {abs(g.imag)/abs(g):.0e};  ' + '; '.join(row)
                   + f'   (e^-2b at 8, 10: {np.exp(-16):.0e}, {np.exp(-20):.0e})')


if __name__ == '__main__':
    pr('# ee1_blocks_check: b-shifted continuous-spin radial recursion (EESTEP_CLOSE.md §1.2)')
    for f in (part_a, part_b, part_c, part_d, part_e):
        f(); pr(f'  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'ee1_blocks_check.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')
