r"""
lem1_cprime (LEMMAS_CLOSE.md §1): Lemma C' -- the domain on which the swapped-form integrand is ONE analytic function, derived from the
s-channel (12)(34) OPE, and checked against EIN's physical-route continuation (read-only import of resultq/ein_dbar.py).

Gauge (og_free / ein_dbar): a = 1/t1, b = 1/t3, y = 2 sin^2(dphi/2) in (0, 2],
    U = z zb = (2-a)(2-b)/(8y),   V = (1-z)(1-zb) = (2+a)(2+b)/(8y).
Claims checked here:
 (A) [P, sympy] lambda is a root (z or zb) iff (a + 4 lambda - 2)(b + 4 lambda - 2) = 8 (2-y) lambda (lambda - 1).
     Hence for Re a >= 0, Re b >= 0 and lambda >= 1: |lhs| >= (4 lambda - 2)^2 = rhs_max + 4 > rhs  =>  z, zb never in [1, inf).
     Sampled margin: min over the closed right half-planes of dist({z, zb}, [1, inf)).  Control: in Re a < 0 the "wall" z in [1, inf)
     IS reached.
 (B) On Omega = {Re a > 0, a not in [2, inf)}^2 the s-channel-sheet function is fixed by L = Log(2-a) + Log(2-b) - log(8y)
     (principal logs) and principal Li2 / log(1-z): closed form Phi_C' for Dbar_1111.  Compared with EIN's tracked physical-route
     continuation at the C x C nodes (kappa = 1) and on hairpin-family ellipses (kappa = 0.5, 0.2), both continuation orders.
 (C) Path tests: polylines inside Omega (via Re a -> 0+, around t = 1/2 through a in (1, 2)) reproduce Phi_C'; crossing the slit
     a in (2, inf) (= t in (0, 1/2)) changes the value by exactly the s-channel monodromy L -> L +- 2 pi i.
 (D) Outside Omega (Re a < 0, still in (C \ ((-inf,0] u [2,inf)))^2): two paths from the Euclidean point to the same target, passing on
     the two sides of a discriminant point a* (z = zb), give different values once one of them has crossed the wall -- the
     non-single-valuedness EIN found, explained.
Output: cache/lem1_cprime.txt.   Runtime ~ 2 min, 1 core.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'resultq'))
import numpy as np
import sympy as sp
import ein_dbar as E          # read-only

OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def phi_cprime(a, b, y):
    """C' closed form of Phi on Omega: principal Li2, log(1-z); log U continued as Log(2-a) + Log(2-b) - log(8y)."""
    U, V = E.UV_of(a, b, y)
    z, zb = E.roots(U, V)
    L = np.log(2 - a + 0j) + np.log(2 - b + 0j) - np.log(8*y)
    num = 2*E.li2p(z) - 2*E.li2p(zb) + L*(np.log(1 - z) - np.log(1 - zb))
    return num/(z - zb), z, zb


def wall_dist(z):
    """distance of z from the ray [1, inf)."""
    return np.where(z.real >= 1, np.abs(z.imag), np.abs(z - 1))


def part_A():
    pr('(A) algebra of the wall')
    a, b, y, lam = sp.symbols('a b y lam')
    U = (2 - a)*(2 - b)/(8*y); V = (2 + a)*(2 + b)/(8*y)
    quad = lam**2 - (1 + U - V)*lam + U                  # roots: z, zb
    ident = sp.expand(8*y*quad - ((a + 4*lam - 2)*(b + 4*lam - 2) - 8*(2 - y)*lam*(lam - 1)))
    pr('    8y [lam^2 - (1+U-V) lam + U] - [(a+4lam-2)(b+4lam-2) - 8(2-y) lam(lam-1)] =', ident)
    pr('    (4lam-2)^2 - 16 lam(lam-1) =', sp.expand((4*lam - 2)**2 - 16*lam*(lam - 1)), ' (> 0: strict separation, also at y -> 0)')
    rng = np.random.default_rng(1)
    N = 400000
    for lab, re_lo, re_hi in (('Re a, Re b in [0, 3] (closed right half-planes; C is Re = 1)', 0.0, 3.0),
                              ('Re a, Re b in [0, 0.02] (near t = i R, the boundary of Omega)', 0.0, 0.02)):
        A = rng.uniform(re_lo, re_hi, N) + 1j*rng.normal(0, 3, N)*rng.choice([0.01, 1, 10], N)
        B = rng.uniform(re_lo, re_hi, N) + 1j*rng.normal(0, 3, N)*rng.choice([0.01, 1, 10], N)
        Y = 2*rng.uniform(1e-6, 1, N)**2
        U, V = E.UV_of(A, B, Y); z, zb = E.roots(U, V)
        d = np.minimum(wall_dist(z), wall_dist(zb))
        # scale-free margin: |(A+4l-2)(B+4l-2)| - 8(2-y) l (l-1) >= 4 on the wall; report min over samples of d/(1+|z|+|zb|)
        pr(f'    {lab}: {N} samples, min dist(z, zb; [1,inf)) = {d.min():.3e}, min of dist/(1+|z|) = {(d/(1 + np.maximum(abs(z), abs(zb)))).min():.3e}; '
           f'any on the wall: {bool((d < 1e-12).any())}')
    # control: left half-plane
    A = -rng.uniform(0, 3, N) + 1j*rng.normal(0, 1, N); B = -rng.uniform(0, 3, N) + 1j*rng.normal(0, 1, N)
    Y = 2*rng.uniform(1e-6, 1, N)**2
    U, V = E.UV_of(A, B, Y); z, zb = E.roots(U, V)
    d = np.minimum(wall_dist(z), wall_dist(zb))
    pr(f'    control Re a, Re b in [-3, 0): fraction of samples within 1e-2 of the wall {np.mean(d < 1e-2):.3e}; min dist {d.min():.2e}')
    # exact wall points in the left half plane: b = conj(a) real-negative family
    lamv = 3.0; yv = 0.05; c = 1.0
    # (a + 4 lam - 2)^2 = 8(2-y) lam(lam-1) with a = b real: a = -4 lam + 2 + sqrt(8(2-y) lam(lam-1))
    av = -4*lamv + 2 + np.sqrt(8*(2 - yv)*lamv*(lamv - 1))
    U, V = E.UV_of(np.array([av + 0j]), np.array([av + 0j]), yv); z, zb = E.roots(U, V)
    pr(f'    explicit wall point: a = b = {av:.6f} (t = {1/av:.4f}), y = {yv}: roots z, zb = {z[0]:.6f}, {zb[0]:.6f} (lam = {lamv})')


def part_B():
    pr('(B) Phi_C\' (closed form on Omega) vs EIN physical-route continuation (ein_dbar.continue_ellipse)')
    for dphi in (0.1, 1.0, 2.5):
        y = 2*np.sin(dphi/2)**2
        for kappa in (1.0, 0.5, 0.2):
            t, w, th = E.contour_ellipse(24, kappa)
            sel = np.arange(0, 24, 1)
            for order in ('ab', 'ba'):
                st, ms = E.continue_ellipse(th[sel][:, None], th[sel][None, :], y, kappa=kappa, order=order)
                phE = E.phi_of(st)
                A = 1/t[sel][:, None] + 0*t[sel][None, :]; B = 0*t[sel][:, None] + 1/t[sel][None, :]
                phC, z, zb = phi_cprime(A, B, y)
                ok = np.abs(z - zb) > 1e-6
                rel = (np.abs(phE - phC)/np.maximum(np.abs(phC), 1e-300))[ok]
                pr(f'    dphi = {dphi}, kappa = {kappa}, order {order}: {ok.sum()} node pairs, max rel |Phi_EIN - Phi_C\'| = {rel.max():.1e}; '
                   f'min wall distance on the grid {np.minimum(wall_dist(z), wall_dist(zb)).min():.2e}; min |z-zb| on route {ms.min():.2e}')


def part_C():
    pr('(C) path tests inside / across Omega (ein_dbar.continue_path from the Euclidean point (eps, eps))')
    y = 2*np.sin(0.7/2)**2
    tgt = (1 + 2.3j, 1 - 1.7j)
    ref = phi_cprime(np.array([tgt[0]]), np.array([tgt[1]]), y)[0][0]
    paths = {
        'via Re a = 0.03 (t near the imaginary axis)': [(0.03 + 2.0j, 1e-3), (0.03 + 2.0j, 0.03 - 4.0j), (1 + 2.3j, 0.03 - 4.0j), tgt],
        'a through (1,2) (t through (1/2,1)) twice': [(1.5 - 1.0j, 1e-3), (1.5 + 1.0j, 1e-3), (1.5 + 1.0j, 1.5 + 2.0j), (1.5 + 1.0j, 1.5 - 1.0j), tgt],
        'large loop around t = 1/2 inside Re a > 0, avoiding the slit': [(1 + 5j, 1e-3), (6 + 5j, 1e-3), (6 + 0.5j, 1e-3), (1 + 2.3j, 1e-3), tgt],
    }
    for lab, pts in paths.items():
        st = E.continue_path(pts, y, nsub=600)
        v = E.phi_of(st)[0]
        pr(f'    {lab}: |Phi - Phi_C\'|/|Phi_C\'| = {abs(v - ref)/abs(ref):.1e}')
    # crossing the slit a in (2, inf) (t in (0, 1/2)): predicted monodromy L -> L - 2 pi i (from above to below)
    pts = [(3 + 0.5j, 1e-3), (3 - 0.5j, 1e-3), (1 - 1.7j, 1e-3), (1 - 1.7j, tgt[1])]
    st = E.continue_path(pts, y, nsub=600)
    v = E.phi_of(st)[0]
    a0, b0 = np.array([1 - 1.7j]), np.array([tgt[1]])
    base, z, zb = phi_cprime(a0, b0, y)
    for m in (+1, -1):
        pred = base[0] + m*2j*np.pi*(np.log(1 - z[0]) - np.log(1 - zb[0]))/(z[0] - zb[0])
        pr(f'    crossing the slit a = 3 (t = 1/3) downward, then to (1-1.7i, 1-1.7i): |Phi - Phi_C\'| = {abs(v - base[0]):.3e}; '
           f'vs monodromy L -> L {"+" if m > 0 else "-"} 2 pi i: |Phi - pred| = {abs(v - pred):.1e}')


def part_D():
    pr('(D) outside Omega: (C \\ ((-inf,0] u [2,inf)))^2 is simply connected, but Phi is not single-valued on it')
    a = sp.symbols('a')

    def wallsteps(pts, y):
        cur = (1e-3 + 0j, 1e-3 + 0j); n = 0
        for (pa, pb) in pts:
            for s_ in np.linspace(0, 1, 801)[1:]:
                aa = cur[0] + s_*(pa - cur[0]); bb = cur[1] + s_*(pb - cur[1])
                Uu, Vv = E.UV_of(np.array([aa]), np.array([bb]), y); z, zb = E.roots(Uu, Vv)
                n += int(min(wall_dist(z)[0], wall_dist(zb)[0]) < 3e-3)
            cur = (pa, pb)
        return n

    for (dphi, b1) in ((0.4, 0.8 - 0.9j), (1.0, 0.5 - 1.2j), (0.4, 1.0 + 0.0j)):
        y = 2*np.sin(dphi/2)**2
        U = (2 - a)*(2 - b1)/(8*y); V = (2 + a)*(2 + b1)/(8*y)
        roots = [complex(r) for r in sp.Poly(sp.expand((1 + U - V)**2 - 4*U), a).nroots()]
        pr(f'    dphi = {dphi}, b = {b1}: discriminant points a* = ' + ', '.join(f'{r:.4f}' for r in roots))
        for r in roots:
            if abs(r.imag) < 0.2:
                continue
            sg = np.sign(r.imag); h = min(0.5, 0.8*abs(r.imag))
            base = [(1e-3 + 0j, b1), (1e-3 + 1j*sg*(abs(r.imag) + 1.0), b1)] if r.real < 0 else [(1e-3 + 0j, b1)]
            P = r + 0.6; T = r - 0.6
            if r.real < 0:
                pre = base + [(P + 1j*sg*1.0, b1), (P, b1)]
            else:
                pre = base + [(P, b1)]
            p_up = pre + [(P + 1j*h, b1), (T + 1j*h, b1), (T, b1)]
            p_dn = pre + [(P - 1j*h, b1), (T - 1j*h, b1), (T, b1)]
            allp = [p[0] for p in p_up + p_dn]
            touches = any((x.real <= 0 and abs(x.imag) < 1e-9) or (x.real >= 2 and abs(x.imag) < 1e-9) for x in allp)
            v1 = E.phi_of(E.continue_path(p_up, y, nsub=1500))[0]; v2 = E.phi_of(E.continue_path(p_dn, y, nsub=1500))[0]
            where = 'in Omega' if r.real > 0 else 'Re a* < 0 (outside Omega)'
            extra = ''
            if T.real > 0 and min(x.real for x in allp) > 0:
                ref = phi_cprime(np.array([T]), np.array([b1]), y)[0][0]; extra = f'; vs Phi_C\' {abs(v1 - ref)/abs(ref):.1e}'
            pr(f'      a* = {r:.4f} [{where}]: Phi(T) above / below a* = {v1:.6f} / {v2:.6f}, |diff|/|Phi| = {abs(v1 - v2)/abs(v1):.2e}; '
               f'near-wall steps above/below = {wallsteps(p_up, y)}/{wallsteps(p_dn, y)}; path on a cut: {touches}{extra}')


def part_E():
    pr('(E) an EIN-type "joint interior" route: (t1, t3) -> (0.8 - 0.4i, 0.8 + 0.4i), then straight lines in t to near-endpoint nodes')
    y = 2*np.sin(0.7/2)**2
    t, w, th = E.contour_ellipse(30, 1.0)
    for (i1, i3) in ((2, 27), (27, 2), (2, 2), (5, 24)):
        t1n, t3n = t[i1], t[i3]
        pts = [(1/(0.8 - 0.4j), 1e-3), (1/(0.8 - 0.4j), 1/(0.8 + 0.4j))]
        for s_ in np.linspace(0, 1, 401)[1:]:
            pts.append((1/((0.8 - 0.4j) + s_*(t1n - (0.8 - 0.4j))), 1/(0.8 + 0.4j)))
        for s_ in np.linspace(0, 1, 401)[1:]:
            pts.append((1/t1n, 1/((0.8 + 0.4j) + s_*(t3n - (0.8 + 0.4j)))))
        cross = []
        for (p, q) in (((0.8 - 0.4j), t1n), ((0.8 + 0.4j), t3n)):
            if (p.imag > 0) != (q.imag > 0):
                sx = p.imag/(p.imag - q.imag); cross.append((p + sx*(q - p)).real)
            else:
                cross.append(None)
        v = E.phi_of(E.continue_path(pts, y, nsub=4))[0]
        base, z, zb = phi_cprime(np.array([1/t1n]), np.array([1/t3n]), y)
        dl = 2j*np.pi*(np.log(1 - z[0]) - np.log(1 - zb[0]))/(z[0] - zb[0])
        m = np.round(((v - base[0])/dl).real)
        pr(f'    nodes t1 = {t1n:.4f}, t3 = {t3n:.4f}: real-axis crossings at Re t = {cross}; |Phi_route - Phi_C\'|/|Phi_C\'| = '
           f'{abs(v - base[0])/abs(base[0]):.2e}; = Phi_C\' + m * (2 pi i monodromy term) with m = {m:+.0f}, residual {abs(v - base[0] - m*dl):.1e}')


if __name__ == '__main__':
    pr('# lem1_cprime: Lemma C\' (domain of one analytic integrand), LEMMAS_CLOSE.md §1')
    part_A(); pr(f'  [{time.time() - T0:.0f}s]')
    part_B(); pr(f'  [{time.time() - T0:.0f}s]')
    part_C(); pr(f'  [{time.time() - T0:.0f}s]')
    part_D(); pr(f'  [{time.time() - T0:.0f}s]')
    part_E(); pr(f'  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'lem1_cprime.txt'), 'w') as f:
        f.write('\n'.join(OUT) + '\n')
