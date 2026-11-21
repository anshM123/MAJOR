"""
v4e_lt (RESULT_Q_V4FINAL.md §1): what the CLS-D / KS-D definitions say about the collinear point phi1 -> phi3 (H_LT).

(1) [sympy, exact] On the COMPLEXIFIED light-transform manifold
        x1 = x0 + a z1,  x3 = x0 + b z3,   x0 in C^3, a, b in C, z1, z3 complex null (z = (1, (1-m^2)/(1+m^2), 2m/(1+m^2)), m in C),
    with the light-ray polarisations z1, z3 (CLS-D: O_1, O_3 polarised along x_10, x_30), the numerator
        N = (z1.z3) x13^2 - 2 (z1.x13)(z3.x13)
    of <T3 T1> vanishes IDENTICALLY, and x13^2 = -2ab z1.z3.  Hence W_u = <T4T2><T3T1> times ANY function of the points that is
    regular on the manifold off z1.z3 = 0 (in particular f4, whose only (13)-singularity is 1/x13^2) is identically zero there.  Every
    regularisation that keeps 1 and 3 on their (complexified) light-transform contours from a COMMON base point -- KS-D's continuation
    v -> v + i Im(v) z of the light-ray parameter, the cone-preserving t -> t + i eps, complex x0, complex celestial angles -- gives 0.
    The timelike tube x1 -> x1 + i eps e leaves the manifold: x10^2 != 0 (violates CLS-D's delta(x_10^2)) and z1 is no longer along
    x10; N = eps^2 (...) != 0.
(2) [numerics] Prescription independence of the absolutely integrable O(G) corner pieces: for the twist-1 light-ray structures
    U_mn (v4d_og; the singular (13) factors of the u-channel T family) the tube limit eps -> 0 and the light-ray (pointwise) value
    coincide.  Pointwise value: soft cutoff (x^2/(x^2+d^2))^2 on the eps = 0+ function, d -> 0.  Both limits are extracted by the
    least-squares fit  G(s) = G0 + A s^(1/4 + i nu/2) + B s^(1/4 - i nu/2)  (s = eps, resp. s = d^2).
Usage: python3 v4e_lt.py  -> cache/v4e_lt.txt  (~1 min, < 1 GB)
"""
import os, sys, time
sys.dont_write_bytecode = True
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import v4d_core as D
from v4d_og import U_mn

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
L = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)


def mdot(a, b):
    return -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def part1():
    pr('(1) exact algebra on the complexified light-transform manifold')
    m1, m3, a, b, eps, t1, t3 = sp.symbols('m1 m3 a b eps t1 t3')
    X0 = sp.symbols('X0:3')

    def z(m):
        return sp.Matrix([1, (1 - m**2)/(1 + m**2), 2*m/(1 + m**2)])
    z1, z3 = z(m1), z(m3)
    assert sp.simplify(mdot(z1, z1)) == 0 and sp.simplify(mdot(z3, z3)) == 0
    x0 = sp.Matrix(X0)
    x1 = x0 + a*z1; x3 = x0 + b*z3; x13 = x1 - x3
    X = sp.simplify(mdot(x13, x13))
    N = sp.simplify(mdot(z1, z3)*X - 2*mdot(z1, x13)*mdot(z3, x13))
    pr(f'    x13^2 = {sp.factor(X)}   (= -2ab z1.z3, z1.z3 = {sp.factor(sp.simplify(mdot(z1, z3)))})')
    pr(f'    N = {N}   (identically zero for complex a, b, m1, m3, x0)  [P]')
    # tube: x1 -> x1 + i eps e, real gauge of v4c (x0 = 0, x1 = -t1 z1, x3 = t3 z3), real angles
    p1, p3 = sp.symbols('p1 p3', real=True)
    w1 = sp.Matrix([1, sp.cos(p1), sp.sin(p1)]); w3 = sp.Matrix([1, sp.cos(p3), sp.sin(p3)])
    e = sp.Matrix([1, 0, 0])
    x1t = -t1*w1 + sp.I*eps*e; x3t = t3*w3; x13t = x1t - x3t
    Xt = sp.expand(mdot(x13t, x13t))
    Nt = sp.simplify(sp.expand(mdot(w1, w3)*Xt - 2*mdot(w1, x13t)*mdot(w3, x13t)))
    pr(f'    tube x1 -> x1 + i eps e:  N = {sp.factor(sp.trigsimp(Nt))}   (non-zero, t-independent; V4DPRIME: eps^2 (1+c))')
    pr(f'    tube: x10^2 = {sp.factor(sp.expand(mdot(x1t, x1t)))}  (!= 0: point 1 leaves the null cone of x0, i.e. the support of CLS-D\'s'
       f' delta(x_10^2); the polarisation z1 is no longer along x10)')
    # cone-preserving: t -> t + i eps
    x1c = -(t1 + sp.I*eps)*w1; x3c = (t3 + sp.I*eps)*w3; x13c = x1c - x3c
    Nc = sp.simplify(sp.expand(mdot(w1, w3)*mdot(x13c, x13c) - 2*mdot(w1, x13c)*mdot(w3, x13c)))
    pr(f'    cone-preserving t -> t + i eps:  N = {Nc},  x10^2 = {sp.simplify(mdot(x1c, x1c))}  [P]')
    # f4's (13) pole is the only (13)-singularity of f4: f4 ~ A (A^2 - x13^4)/(4 B x13^2)
    A_, B_, Y = sp.symbols('A B Y')        # Y = x13^2, A = x14^2 x23^2 - x12^2 x34^2 (over x24^2), B = x12^2 x34^2 x14^2 x23^2/x24^4
    up = sp.Symbol('s12')*sp.Symbol('s34')/(Y*sp.Symbol('s24')); vp = sp.Symbol('s14')*sp.Symbol('s23')/(Y*sp.Symbol('s24'))
    f4 = sp.factor(sp.together(((vp - up)**3 - (vp - up))/(4*up*vp)))
    pr(f'    f4 as a function of Y = x13^2:  {f4}   -> simple pole at Y = 0 only; Lemma C assigns Wightman i-eps to the (12), (14), (23),'
       f' (34) poles; the (13) pole has none (x13^2 > 0 on the whole real domain).')


def fit_limit(s, G, nu):
    """least squares G(s) = G0 + A s^(1/4 + i nu/2) + B s^(1/4 - i nu/2); returns G0 and residual."""
    s = np.asarray(s, float)
    M = np.stack([np.ones_like(s) + 0j, s**(0.25 + 0.5j*nu), s**(0.25 - 0.5j*nu)], axis=1)
    c, *_ = np.linalg.lstsq(M, np.asarray(G), rcond=None)
    res = np.max(np.abs(M @ c - G))
    return c[0], res


def part2():
    pr('(2) prescription independence of the O(G) twist-1 corner pieces (tube + limit vs light-ray pointwise value)')
    pr('    U_mn: v4d_og structures (edge residue of t^3 t^3 f4 x (13)-leg D1^m D3^n X^(-1/2)); pairing: v4d_core.pair_exact (exact Fourier)')
    EPS = np.array([1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7])
    DD = np.sqrt(EPS)                                   # soft-cutoff radii on the same scale
    for (m, n) in ((1, 1), (2, 2), (0, 0)):
        for nu in (2.0, 6.0):
            Gt, Gc = [], []
            for e, d in zip(EPS, DD):
                N = 2**18 if e < 2e-6 else 2**17
                Gt.append(D.pair_exact(nu, lambda x: U_mn(x, e, m, n, +1), N=N))

                def gc(x, d=d):
                    with np.errstate(all='ignore'):
                        u = U_mn(x, 1e-300, m, n, +1)*(x**2/(x**2 + d**2))**2
                    return np.where(x == 0, 0.0, u)
                Gc.append(D.pair_exact(nu, gc, N=N))
            Gt, Gc = np.array(Gt), np.array(Gc)
            G0t, rt = fit_limit(EPS, Gt, nu)
            G0c, rc = fit_limit(DD**2, Gc, nu)
            pr(f'    (m,n)=({m},{n}) nu={nu}:  tube eps->0: G0 = {G0t:.6e} (fit res {rt:.1e});  light-ray (cutoff d->0): G0 = {G0c:.6e}'
               f' (fit res {rc:.1e});  |diff|/|G0| = {abs(G0t - G0c)/abs(G0c):.1e}')
            pr(f'        raw sequences |G_tube| = ' + ' '.join(f'{abs(x):.4e}' for x in Gt) + ';  |G_cut| = ' + ' '.join(f'{abs(x):.4e}' for x in Gc))
    pr('    identity W_u: pointwise U(dphi, eps = 0) = 0 for dphi != 0 (factor eps^4) -> light-ray value exactly 0; the tube value diverges'
       ' (V4DPRIME §1): the only prescription-dependent piece is the universal C_T^2 u-channel identity.')


def main():
    t0 = time.time()
    pr('# v4e_lt: H_LT from the definitions (exact algebra) + prescription independence of the integrable O(G) corner pieces')
    part1()
    part2()
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(CACHE, 'v4e_lt.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
