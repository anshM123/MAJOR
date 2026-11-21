"""Result Q, Stage 0 (RESULT_Q_STAGE0.md): flat D = 4 graviton sum rules with forward terms, for the quartic contact g4.
Floating point only: every LP number produced with this module is [NE].

Sum rules (arXiv:2201.06602 (2.26), (2.32), heavy densities (2.33)-(2.35a); conventions of FLAT_GRAVITON_BOUND.md section 1).
Units M = 1, u = 1/m^2, JJ = J(J+1), sigma_{++} = +1, sigma_{+-} = -1.  The functional family used here:

  F = int_0^1 [psi2 B2 + psi3 B3] dp  - a (d_{p^2}B4|_0 - B5|_0)  - c B4|_0          (a >= 0, c free)

 * low side (-F|_low = F|_high, tree level, no light matter):
       16 pi G L(eps) - pi G k |g3|^2 - 2 c g4 ,   k = -2 mu2_6 + 2 mu3_4 + 4 a ,  L(eps) = int_eps^1 psi2/p^2 = x0 log(1/eps) + Lfin
   (the pair -a dB4 + a B5 removes g5;  -B4|_0 low = 2 g4,  -dB4|_0 low = 4 pi G|g3|^2 + g5,  -B5|_0 low = g5)
 * heavy action per |c^h_{J,m^2}|^2 in normal form FF = m^4 F:
       FF = int psi2 (2 - u p^2) D_h(1 - 2 u p^2) + sigma u int psi3 D_h + a (2JJ - c4_h + sigma) u^3 - 2 c u^2 ,
   c4_{++} = 1, c4_{+-} = 41 (rcommon.phi_coef), D_{++} = P_J (J even), D_{+-} = P^{(0,8)}_{J-4} (J >= 4).
With c = 1/2 and heavy positivity:  g4/(8 pi G) <= 2 L(eps) - (k/8) |g3|^2 , i.e.  g4 M^6/(8 pi G) <= A log(M/m_IR) + B - kappa |g3|^2 M^8,
A = 2 x0, B = 2 Lfin, kappa = k/8.  AdS IR dictionary (Result D): m_IR -> 3/(2 R_AdS), log(M/m_IR) -> log(2 Lam/3), Lam = M R_AdS.

Only read-only imports from ads4/code/robust/rcommon.py (partial-wave recursion, quadrature, polynomial helpers).
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'robust'))
from rcommon import (pmul, peval, mom, flat_integrals, nodes, nquad, jmax_for, scaling_B2, SIGMA, C4, phi_coef)  # noqa: E402

CACHE = os.path.join(HERE, 'cache')
os.makedirs(CACHE, exist_ok=True)


def omp(r):
    c = np.array([1.])
    for _ in range(r):
        c = pmul(c, [1., -1.])
    return c


def basis2(d2, r):
    """psi2: (1-p)^r p^{i+1}, i = 0..d2 (psi2(0) = 0 is required by the graviton pole)."""
    o = omp(r)
    return [pmul(o, np.eye(i + 2)[i + 1]) for i in range(d2 + 1)]


def basis3(d3, r, j0=0):
    """psi3: (1-p)^r p^j, j = j0..d3 (j0 = 1 imposes psi3(0) = 0)."""
    o = omp(r)
    return [pmul(o, np.eye(j + 1)[j]) for j in range(j0, d3 + 1)]


def u_grid(nu=70):
    return np.unique(np.concatenate([np.logspace(np.log10(2.5e-5), 0, nu), np.linspace(0.5, 1, 26), np.linspace(0.8, 1, 21)]))


def parts(h, u, Jmax, b2, b3):
    n = nquad(Jmax)
    p, w = nodes(n)
    cols = [w * peval(c, p) * (2 - u * p * p) for c in b2] + [w * SIGMA[h] * u * peval(c, p) for c in b3]
    Js, R = flat_integrals(h, u, Jmax, np.array(cols).T, p)
    return Js, R[:, :len(b2)], R[:, len(b2):]


class Grid:
    """rows over (u, h, J): B2 basis columns, B3 basis columns, FWc (per unit a), CB4 = 2u^2 (per unit c, enters with minus)."""

    def __init__(self, d2, d3, r=3, j0=0, bmax=200.0, Jcap=2500, nu=70, verbose=False):
        self.d2, self.d3, self.r, self.j0 = d2, d3, r, j0
        self.b2 = basis2(d2, r); self.b3 = basis3(d3, r, j0)
        self.n2, self.n3 = len(self.b2), len(self.b3)
        fn = os.path.join(CACHE, 'grid_d%d_%d_r%d_j%d_b%d_J%d_n%d.npz' % (d2, d3, r, j0, bmax, Jcap, nu))
        if os.path.exists(fn):
            z = np.load(fn)
            for k in ['B2', 'B3', 'u', 'J', 'h']:
                setattr(self, k, z[k])
        else:
            L = {k: [] for k in ['B2', 'B3', 'u', 'J', 'h']}
            for u in u_grid(nu):
                Jm = jmax_for(u, bmax, Jcap)
                for h in SIGMA:
                    Js, B2, B3 = parts(h, u, Jm, self.b2, self.b3)
                    if h == '++':
                        s = Js % 2 == 0; Js, B2, B3 = Js[s], B2[s], B3[s]
                    L['B2'].append(B2); L['B3'].append(B3); L['u'].append(np.full(len(Js), u)); L['J'].append(Js)
                    L['h'].append(np.full(len(Js), SIGMA[h]))
                if verbose:
                    print('  grid u=%.3e' % u, flush=True)
            for k in L:
                setattr(self, k, np.concatenate(L[k]))
            np.savez(fn, **{k: getattr(self, k) for k in L})
        JJ = self.J * (self.J + 1.0)
        c4 = np.where(self.h > 0, 1.0, 41.0)
        self.FWc = (2 * JJ - c4 + self.h) * self.u ** 3
        self.CB4 = 2 * self.u ** 2
        # scaling limit u -> 0 at fixed bhat: 2 int psi2 J0(bhat p) (forward terms vanish like u^2)
        self.bsc = np.unique(np.concatenate([np.linspace(0, 60, 1201), np.linspace(60, 400, 1700)]))
        self.Bsc = scaling_B2(self.b2, self.bsc)
        # large-bhat asymptotics at fixed u (Lemma-M form of rcommon/lp.py):  B2 ~ -psi2''(0)/bhat^3,
        # B3 ~ sigma u (psi3(0)/bhat - psi3''(0)/(2 bhat^3)),  FW ~ a (bhat^2 u^2/2 - u^3/2 - (c4 - sigma) u^3),  CB4 = 2u^2
        bb = np.logspace(np.log10(0.75 * bmax), 4, 120); uu = np.logspace(-9, np.log10(0.25), 160)
        Bb, Uu = np.meshgrid(bb, uu); self.ab, self.au = Bb.ravel(), Uu.ravel()
        d2pp = np.array([2 * c[2] for c in self.b2])
        p30 = np.array([c[0] for c in self.b3]); p3pp = np.array([2 * c[2] if len(c) > 2 else 0. for c in self.b3])
        self.asym = {}
        for h, s in SIGMA.items():
            A2 = np.outer(1 / self.ab ** 3, -d2pp)
            A3 = s * self.au[:, None] * (np.outer(1 / self.ab, p30) - np.outer(1 / (2 * self.ab ** 3), p3pp))
            Fc = self.ab ** 2 * self.au ** 2 / 2 - self.au ** 3 / 2 - (C4[h] - s) * self.au ** 3
            self.asym[h] = (A2, A3, Fc, 2 * self.au ** 2)
        self.mu2_6 = np.array([mom(c, 6) for c in self.b2]); self.mu3_4 = np.array([mom(c, 4) for c in self.b3])
        self.x0 = np.array([c[1] for c in self.b2])
        self.Lfin = np.array([sum(c[i] / (i - 1) for i in range(2, len(c))) for c in self.b2])
        self.pp0 = np.array([2 * c[2] for c in self.b2])

    def LAdS(self, Lam):
        """int_0^1 phi/(p^2 + (3/(2Lam))^2) per psi2 basis function (exact AdS-regularised log, as in Result D)."""
        a = 1.5 / Lam
        edges = np.unique(np.concatenate([[0.0], a * np.logspace(-3, np.log10(1 / a), 60), [1.0]]))
        edges = edges[edges <= 1.0]
        x, w = np.polynomial.legendre.leggauss(40)
        tot = np.zeros(self.n2)
        for lo, hi in zip(edges[:-1], edges[1:]):
            p = lo + (x + 1) * (hi - lo) / 2; ww = w * (hi - lo) / 2
            for i, c in enumerate(self.b2):
                tot[i] += np.sum(ww * peval(c, p) / (p * p + a * a))
        return tot

    def matter_rows(self, nm=201):
        """light scalar / spin-2 matter (FLAT_GRAVITON_BOUND.md 1.3, B.1) coefficients of a0, a2 in -F|_low must be <= 0... in
        F|_high = -F|_low they appear with a plus sign on the RIGHT of the bound; requiring them <= 0 makes the bound matter-robust.
        Columns (x, y, a, c).  B4|_0 matter: -B4 low = 12 a2 p^2/m^4 -> contributes 0 at p = 0;  dB4|_0: -12 a2/m^4."""
        rows = []
        for m in np.linspace(1e-3, 1, nm):
            m2 = m * m
            r = np.zeros(self.n2 + self.n3 + 2)
            r[:self.n2] = [mom(c, 2) - 2 * m2 * mom(c, 0) for c in self.b2]
            r[self.n2:self.n2 + self.n3] = [-mom(c, 0) for c in self.b3]
            rows.append(r)
            r = np.zeros(self.n2 + self.n3 + 2)
            r[:self.n2] = [mom(c, 2) * m2 * m2 - 6 * m2 * mom(c, 4) + 6 * mom(c, 6) - 2 * m2 * (m2 * m2 * mom(c, 0) - 6 * m2 * mom(c, 2) + 6 * mom(c, 4)) for c in self.b2]
            r[self.n2:self.n2 + self.n3] = [-(m2 * m2 * mom(c, 0) - 6 * m2 * mom(c, 2) + 6 * mom(c, 4)) for c in self.b3]
            r[-2] = -12
            rows.append(r)
        return np.array(rows)


def heavy_action(psi2, psi3, a, c, h, u, Jmax):
    """FF for explicit polynomials (monomial coefficients) for J = J0..Jmax at fixed u (allowed J only)."""
    Js, B2, B3 = parts(h, u, Jmax, [np.asarray(psi2, float)], [np.asarray(psi3, float)])
    if h == '++':
        s = Js % 2 == 0; Js, B2, B3 = Js[s], B2[s], B3[s]
    FW = a * phi_coef(h, Js, u)
    return Js, B2[:, 0] + B3[:, 0] + FW - 2 * c * u * u, np.abs(B2[:, 0]) + np.abs(B3[:, 0]) + np.abs(FW) + 2 * abs(c) * u * u
