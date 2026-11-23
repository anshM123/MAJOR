r"""
hk1_heavy (HK_RESIDUAL.md §1): the heavy image of the k = 4 EE position kernel of EESTEP, and the nu-normalisation behind f_pos(+-3i/2) = 0.

Two links, then the exact algebra at +-3i/2.
(a) POSITION SPACE.  EESTEP's leading-Regge EE kernel L_nu(beta; pols) (ee6_lrprofile.L: the r^0 Cauchy coefficient of the b-shifted block
    sum sum_b C_b v^{b/2} Gh_b, scalar normalisation E = K/(N(nu) r^2), scalar L = P_l) is compared with the split representation
    whose plane-wave transform stage4_hj computed:
        K0^split[pols] = int dphi A(xh, z0) B(yh, z0),   A_EE = 16 (z2'.z4')^2 A_sc,  A_sc = (-z0.xh)^{-(1-delta)},
        B_EE = v3_regge_kernel.B_monomials(spin=True) (theta-average of tau13 with |sin|^{-delta}, divided by the common <1>),
    i.e. the prediction   L_nu/P_l  =  16 (z2.z4)^2 int A_sc B_EE / int A_sc.     A nu-independent ratio (= 1) means that the
    kernel stage 4 transformed IS the EESTEP kernel, with the scalar nu-normalisation.
(b) HEAVY IMAGE.  The rank-3 plane-wave transform of that split kernel, L_4 = -16 <d^3 Ahat, d^3 Bhat> (stage4_hj §1.1), evaluated
    NUMERICALLY: a truncated 5-variable Taylor series (e.d, f.d, d_0, d_1, d_2; orders 2, 2, 3) of S^-lam R^kap on a z0 grid, full
    eta-contraction, trapezoid in phi, then a least-squares fit L_4(beta) = sum_{j = +-1, +-3} h_j P_{l+j}(cosh beta) over 14 betas.
    No step of stage4_hj's symbolic route (equivariant decomposition, Gram contraction, to_harmonics) is used.  Compared with the
    stored exact R_j (stage4_hj.pkl): h_j^EE / h_{-1}^sc  vs  R_j^EE3 / R_{-1}^sc3, and h_j^sc / h_{-1}^sc vs R_j^sc3 / R_{-1}^sc3.
(c) EXACT ALGEBRA at nu = +-3i/2: R_{-3}^EE3, R_{-1}^EE3, R^EE3/R^sc3 there; the position/heavy normalisation N(nu)/N_com of the scalar
    (= R_j^sc3 / (a_j(l) E_9(nu - i j)) x const, LEMMAS §3.2) is regular and non-zero at +-3i/2.
Output: cache/hk1_heavy.txt.   ~ 3-5 min, one core.
"""
import os, sys, time, math, itertools
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..')
for sub in ('eestep', 'h7i14', 'spinfunc', 'resultq'):
    sys.path.insert(0, os.path.join(CODE, sub))
import numpy as np
import mpmath as mp
import sympy as sp

OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def mdot(u, v):
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]


# ------------------------------------------------------------------------------------------------------------------ (a)
def part_a():
    import ee6_lrprofile as E6          # read-only import (EESTEP)
    import ee_quad as Q
    import ee3_deshift as D3
    import v3_regge_kernel as V
    Bm = V.B_monomials(True)
    nphi = 2**17                         # (-z0.xh)^(delta-1) peaks with width ~2e^-beta: resolved to beta = 6
    phi = 2*np.pi*np.arange(nphi)/nphi
    z0 = np.stack([np.ones(nphi), np.cos(phi), np.sin(phi)], -1)
    yh = np.array([1.0, 0.0, 0.0])
    rng = np.random.default_rng(7)
    POLS = list(D3.POLS) + [tuple(rng.uniform(0, 2*np.pi, 4)) for _ in range(2)]
    pr('(a) EESTEP leading-Regge EE kernel vs the split representation of stage 4 (v3 B_monomials, A_EE = 16 e24^2 A_sc)')
    pr('    ratio := [L_nu/P_l (EESTEP, b-shifted blocks + Cauchy r^0)] / [16 (z2.z4)^2 int A_sc B_EE / int A_sc]')
    allr = []
    for ip, pa in enumerate(POLS):
        for beta in (1.0, 3.0, 6.0):
            w1, w3, z2, z4 = Q.pols(beta, *pa)
            xh = np.array([np.cosh(beta), np.sinh(beta), 0.0])
            row = []
            for nu in (0.5, 1.5, 3.0, 6.0):
                delta = 0.5 + 1j*nu
                L, Pl = E6.L(nu, beta, pa)
                Asc = (-mdot(z0, xh))**(delta - 1)
                A1 = mdot(z0, w1); A3 = mdot(z0, w3); E13 = mdot(w1, w3); U1 = mdot(w1, yh); U3 = mdot(w3, yh)
                Bv = 0
                for (cf, kf, mf, i, j) in Bm:          # Y = |yh| = 1, t = -z0.yh = 1
                    Bv = Bv + cf(A1, A3, E13, delta)*U1**i*U3**j
                pred = 16*mdot(z2, z4)**2*np.sum(Asc*Bv)/np.sum(Asc)
                rat = (L/Pl)/pred
                row.append(rat); allr.append(rat)
            pr(f'    pols {ip} {np.round(pa, 2).tolist()}, beta = {beta}: nu = 0.5, 1.5, 3, 6: '
               + ' '.join(f'{x.real:+.8f}{x.imag:+.1e}i' for x in row))
    allr = np.array(allr)
    pr(f'    max |ratio - 1| over 4 pol sets x 3 betas x 4 nus: {np.max(np.abs(allr - 1)):.2e}   [{time.time() - T0:.0f}s]')
    return allr


# ------------------------------------------------------------------------------------------------------------------ (b)
SH = (3, 3, 4, 4, 4)
IDX = [k for k in itertools.product(range(3), range(3), range(4), range(4), range(4)) if k[2] + k[3] + k[4] <= 3]
MASK = np.zeros(SH, bool)
for k in IDX:
    MASK[k] = True
FACT = np.zeros(SH)
for k in IDX:
    FACT[k] = np.prod([math.factorial(x) for x in k])
S3 = [(a, b, c) for a in range(4) for b in range(4) for c in range(4) if a + b + c == 3]


def smul(A, B):
    out = np.zeros_like(A)
    for k in IDX:
        a = A[k]
        if not np.any(a):
            continue
        i, j, s0, s1, s2 = k
        out[i:, j:, s0:, s1:, s2:] += a[None, None, None, None, None]*B[:3 - i, :3 - j, :4 - s0, :4 - s1, :4 - s2]
    out[~MASK] = 0
    return out


def gbinom(al, n):
    v = 1.0 + 0j
    for k in range(n):
        v *= (al - k)/(k + 1)
    return v


class Side:
    """the truncated Taylor powers of S = -z0.(p+h), R = -(p+h)^2 around p, h = te e + tf f + s_mu E_mu; vectorised over z0."""
    def __init__(self, p, z, e, f):
        n = z.shape[0]
        E = np.eye(3)
        dS = np.zeros(SH + (n,), complex); dR = np.zeros(SH + (n,), complex)
        S0 = -mdot(z, p[None, :]); R0 = -mdot(p, p)
        dirs = {(1, 0, 0, 0, 0): e, (0, 1, 0, 0, 0): f, (0, 0, 1, 0, 0): E[0], (0, 0, 0, 1, 0): E[1], (0, 0, 0, 0, 1): E[2]}
        for k, v in dirs.items():
            dS[k] = -mdot(z, v[None, :])
            dR[k] = -2*mdot(p, v)
        for (k1, v1), (k2, v2) in itertools.combinations_with_replacement(list(dirs.items()), 2):
            k = tuple(x + y for x, y in zip(k1, k2))
            c = -mdot(v1, v2)*(1 if k1 == k2 else 2)
            dR[k] = dR[k] + c
        self.S0, self.R0 = S0, R0
        X = dS/S0[None, None, None, None, None, :]; Y = dR/R0
        one = np.zeros_like(X); one[0, 0, 0, 0, 0] = 1
        self.Xn = [one]; self.Yn = [one]
        for _ in range(7):
            self.Xn.append(smul(self.Xn[-1], X)); self.Yn.append(smul(self.Yn[-1], Y))
        self.cache = {}

    def deriv(self, lam, kap):
        key = (complex(lam), complex(kap))
        if key not in self.cache:
            SX = sum(gbinom(-lam, n)*self.Xn[n] for n in range(8))
            SY = sum(gbinom(kap, n)*self.Yn[n] for n in range(8))
            T = smul(SX, SY)*(self.S0**(-lam)*self.R0**kap)[None, None, None, None, None, :]
            self.cache[key] = T*FACT[..., None]
        return self.cache[key]


def side_rank3(side, monos, coef_args, delta, b=6, lshift=0.0):
    """dict (a,b,c) -> array over z0: d0^a d1^b d2^c of the plane-wave transform of |x|^{b-3} x (monomial sum)."""
    import v3_regge_kernel as V
    T = {s: 0 for s in S3}
    for (cf, kf, mf, i, j) in monos:
        k = complex(kf(delta)); m = complex(mf(delta))
        c = (b - 3) + k; lam = -m; kap = (2*lam - c - 3)/2
        Kv = np.exp(V.logKc(c, lam) - lshift)
        D = side.deriv(lam, kap)
        cv = cf(*coef_args, delta)
        for s in S3:
            T[s] = T[s] + cv*Kv*D[(i, j) + s]
    return T


def contract3(TA, TB):
    tot = 0
    for (a, b, c) in S3:
        mult = math.factorial(3)/(math.factorial(a)*math.factorial(b)*math.factorial(c))
        tot = tot + mult*(-1)**a*TA[(a, b, c)]*TB[(a, b, c)]
    return tot


def heavy_rank3(nu, beta, kind, nphi=192):
    import v3_regge_kernel as V
    import s1_offdiag_kernel as S1
    delta = 0.5 + 1j*nu
    ph, pbh, ep_p, em_p, ep_pb, em_pb = V.frames(beta)
    e2, e4, e1, e3 = ep_p, em_p, ep_pb, em_pb                 # h1 = h2 = +1 (helicity-blind, stage4 (H))
    phi = 2*np.pi*np.arange(nphi)/nphi
    z = np.stack([np.ones(nphi), np.cos(phi), np.sin(phi)], -1)
    lsA = V.logKc(6 - 2 - delta, 1 - delta).real; lsB = V.logKc(6 - 3 + delta, delta).real   # as v3.heavy_image (common nu-only)
    sA = Side(ph, z, e2, e4); sB = Side(pbh, z, e1, e3)
    cA = (mdot(z, e2[None, :]), mdot(z, e4[None, :]), np.full(nphi, mdot(e2, e4)))
    cB = (mdot(z, e1[None, :]), mdot(z, e3[None, :]), np.full(nphi, mdot(e1, e3)))
    if kind == 'sc':
        TA = side_rank3(sA, V.A_monomials(False), cA, delta, lshift=lsA)
        TB = side_rank3(sB, V.B_monomials(False), (0*phi, 0*phi, 0*phi), delta, lshift=lsB)
    else:
        TA = side_rank3(sA, S1.A_monos_struct('EE'), cA, delta, lshift=lsA)
        TB = side_rank3(sB, V.B_monomials(True), cB, delta, lshift=lsB)
    return -16*2*np.pi*np.mean(contract3(TA, TB))


def Pcon(nu, beta):
    return complex(mp.legenp(-0.5 + 1j*nu, 0, mp.cosh(beta), type=3))


def part_b():
    import stage4_hj as S4
    Rex = S4.load()
    fn = {key: {j: sp.lambdify(sp.Symbol('nu'), v, 'mpmath') for j, v in Rex[key].items()} for key in ('sc3', 'EE3')}
    betas = np.linspace(0.15, 3.0, 14)
    JS = (-3, -1, 1, 3)
    pr('\n(b) numerical rank-3 plane-wave transform of the split kernel (independent of stage4_hj\'s symbolic route), fitted to'
       ' sum_j h_j P_{l+j}(cosh beta), j = -3, -1, 1, 3')
    worst = 0.0
    for nu in (0.37, 1.3, 2.9, 5.5):
        M = np.array([[Pcon(nu - 1j*j, b) for j in JS] for b in betas])
        h = {}
        for kind in ('sc', 'EE'):
            v = np.array([heavy_rank3(nu, b, kind) for b in betas])
            sol, *_ = np.linalg.lstsq(M, v, rcond=None)
            res = np.linalg.norm(M @ sol - v)/np.linalg.norm(v)
            h[kind] = dict(zip(JS, sol)); h[kind + '_res'] = res
        line = []
        for kind, key in (('sc', 'sc3'), ('EE', 'EE3')):
            for j in JS:
                num = h[kind][j]/h['sc'][-1]
                ex = complex(fn[key][j](mp.mpf(nu))/fn['sc3'][-1](mp.mpf(nu)))
                dev = abs(num/ex - 1); worst = max(worst, dev)
                line.append(f'{kind} j={j:+d}: {dev:.1e}')
        pr(f'    nu = {nu}: fit resid sc {h["sc_res"]:.1e}, EE {h["EE_res"]:.1e};  |(h_j/h_-1^sc)/(R_j/R_-1^sc3) - 1|: ' + ', '.join(line)
           + f'   [{time.time() - T0:.0f}s]')
        # the k = 2 factorisation "EE = scalar x F(nu - i j)" fails: print the numerical h^EE_j/h^sc_j next to the exact ratio
        pr('        h_j^EE/h_j^sc (num | exact): ' + ', '.join(
            f'j={j:+d}: {complex(h["EE"][j]/h["sc"][j]):.6g} | {complex(fn["EE3"][j](mp.mpf(nu))/fn["sc3"][j](mp.mpf(nu))):.6g}' for j in JS))
    pr(f'    worst relative deviation from the stored exact R_j (EE3 and sc3, both relative to R_-1^sc3): {worst:.1e}')
    return worst


# ------------------------------------------------------------------------------------------------------------------ (c)
def part_c():
    import stage4_hj as S4
    Rex = S4.load()
    nu = sp.Symbol('nu')
    pr('\n(c) exact algebra at nu = +-3i/2 (stored R_j)')
    for nv in (-sp.Rational(3, 2)*sp.I, sp.Rational(3, 2)*sp.I):
        vals = {j: sp.nsimplify(sp.simplify(Rex['EE3'][j].subs(nu, nv))) for j in sorted(Rex['EE3'])}
        pr(f'    R_j^EE3({nv}): ' + ', '.join(f'j={j:+d}: {v}' for j, v in vals.items()))
        vals = {j: sp.simplify(Rex['sc3'][j].subs(nu, nv)) for j in sorted(Rex['sc3'])}
        pr(f'    R_j^sc3({nv}): ' + ', '.join(f'j={j:+d}: {v}' for j, v in vals.items()))
    # N(nu)/N_com of the scalar, from R_j^sc3 = const N(nu)/N_com a_j(l) E_9(nu - i j)
    l = -sp.Rational(1, 2) + sp.I*nu

    def a3(j):
        co = {0: sp.Integer(1)}
        for _ in range(3):
            new = {}
            for jj, c in co.items():
                m = l + jj
                new[jj + 1] = new.get(jj + 1, 0) + c*(m + 1)/(2*m + 1)
                new[jj - 1] = new.get(jj - 1, 0) + c*m/(2*m + 1)
            co = new
        return sp.factor(sp.simplify(co[j]))
    E9 = lambda w: sp.gamma((9 - sp.Rational(1, 2) - sp.I*w)/2)**2*sp.gamma((9 - sp.Rational(1, 2) + sp.I*w)/2)**2
    rat = {}
    for j in (-3, -1, 1, 3):
        rat[j] = sp.simplify(sp.gammasimp(Rex['sc3'][j]/(a3(j)*E9(nu - sp.I*j))))
    # all four must be the same function of nu (up to one constant): print ratios rat_j/rat_-1 (exactly 1 expected)
    pr('    scalar normalisation n(nu) := R_j^sc3/(a_j(l) E_9(nu - i j)): ratios n_j/n_-1 at nu = 0.7, 2.3+0.4i: ' + ', '.join(
        f'j={j:+d}: {complex(sp.N((rat[j]/rat[-1]).subs(nu, 0.7))):.3g}, {complex(sp.N((rat[j]/rat[-1]).subs(nu, 2.3 + 0.4*sp.I))):.3g}'
        for j in (-3, 1, 3)))
    nfun = sp.simplify(sp.gammasimp(rat[-1]))
    pr(f'    n(nu) (from j = -1) = {nfun}')
    for nv in (-sp.Rational(3, 2)*sp.I, sp.Rational(3, 2)*sp.I):
        pr(f'    n({nv}) = {sp.N(nfun.subs(nu, nv), 12)}   (finite and non-zero => N(nu)/N_com regular, non-zero there)')
    # EE3/sc3 ratio at the shifted points: R_-3^EE3 vs R_-3^sc3
    r3 = sp.factor(sp.simplify(Rex['EE3'][-3]/Rex['sc3'][-3]))
    pr(f'    R_-3^EE3/R_-3^sc3 = {r3}')


if __name__ == '__main__':
    pr('# hk1_heavy')
    part_a()
    part_b()
    part_c()
    pr(f'\n[{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'hk1_heavy.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')
