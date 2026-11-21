"""
bg_ee (BELOW_FOLD2.md §2.3): the EE graviton integrand (w0_ee_common.make_build_even_EE: C^_EE = Psi_2[EE], j' = 0, f = v' - u',
even J) on the a-series algebra of bg_ser / bg_geom, i.e. the EE polarisation prefactor as Laurent polynomials in the radial
variables (exact coefficients at fixed Delta), plus the Lemma-B factor R = ((v'-u')^2 - 1)/(4u'v') for k = 4.

Mirrors l5_spin (Vfun, Hfun, Xfun, tto_samples, null_circle, kappa, fourier) and w0_ee_common (e_factor, _gs, _contract) term by
term; the conserved <TTO> coefficient vectors are l5_spin.conserved (read-only import).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import comb
from bg_ser import Poly, ASer
import bg_geom as GG
from bg_geom import dot, sub, msq

M = Poly.mono


def null_circle(l):
    n = 2*l + 1; ph = 2*np.pi*np.arange(n)/n
    return np.stack([1 - np.cos(ph), 1 + np.cos(ph), np.sin(ph)], axis=1)


def kappa(l):
    n = np.arange(-l, l + 1)
    return (-1.0)**(l + n)*2.0**(-l)*comb(2*l, l + n)


def dft(l):
    """matrix F (2l+1 x 2l+1): coefficients n = -l..l (rows) from samples on null_circle(l) (cols) (= l5_spin.fourier)."""
    nn = 2*l + 1
    k = np.arange(nn)
    full = np.exp(-2j*np.pi*np.outer(k, k)/nn)/nn            # numpy fft convention
    idx = np.arange(-l, l + 1) % nn
    return full[idx]


def lin(coefs, objs):
    """sum_i coefs[i] objs[i] (objs ASer)."""
    out = None
    for c, o in zip(coefs, objs):
        if c == 0:
            continue
        t = o.scale(complex(c)) if abs(np.imag(c)) > 0 else o.scale(float(np.real(c)))
        out = t if out is None else out + t
    return out


class EEGeom:
    def __init__(self, K, G=None):
        self.K = K
        G = GG.geometry(K) if G is None else G
        self.G = G
        self.x = {0: G['x0'], 1: G['x1'], 2: G['x2'], 3: G['x3'], 4: G['x4'], 5: G['e']}
        self.z0 = G['z0']
        self._sq = {}
        self._inv = {}
        self.basis = [[ASer.const(1.0 if i == m else 0.0, K) for i in range(3)] for m in range(3)]

    def sq(self, i, j):
        """(x_i - x_j)^2 (= Xfun)."""
        key = (min(i, j), max(i, j))
        if key not in self._sq:
            r = sub(self.x[i], self.x[j]); self._sq[key] = dot(r, r)
        return self._sq[key]

    def inv_sq(self, i, j):
        key = (min(i, j), max(i, j))
        if key not in self._inv:
            self._inv[key] = self.sq(i, j).recip(f'x{key[0]}{key[1]}^2')
        return self._inv[key]

    def zx(self, z, i, j):
        """Z_i . X_j = z . (x_j - x_i)."""
        return dot(z, sub(self.x[j], self.x[i]))

    def V(self, i, z, j, k):
        """V_{i,jk} = (Z_i.X_j X_i.X_k - Z_i.X_k X_i.X_j)/(X_j.X_k), X_a.X_b = -(x_a - x_b)^2/2."""
        num = (self.zx(z, i, j)*self.sq(i, k) - self.zx(z, i, k)*self.sq(i, j))
        return num*self.inv_sq(j, k)             # (-1/2)/(-1/2) = 1

    def H(self, i, zi, j, zj):
        """H_ij = -2(Z_i.Z_j X_i.X_j - Z_i.X_j Z_j.X_i)."""
        return (dot(zi, zj)*self.sq(i, j)).scale(1.0) + (self.zx(zi, i, j)*self.zx(zj, j, i)).scale(2.0)


def tto(Gm, ia, za, ib, zbs, i5, z5s, J, mons, sol):
    """l5_spin.tto_samples: array [struct][nb][n5] of ASer (numerators only)."""
    B = Gm.basis
    V1 = Gm.V(ia, za, ib, i5)
    V2b = [Gm.V(ib, B[m], i5, ia) for m in range(3)]
    V5b = [Gm.V(i5, B[m], ia, ib) for m in range(3)]
    H12b = [Gm.H(ia, za, ib, B[m]) for m in range(3)]
    H15b = [Gm.H(ia, za, i5, B[m]) for m in range(3)]
    H25b = [[Gm.H(ib, B[m], i5, B[k]) for k in range(3)] for m in range(3)]
    nb, n5 = zbs.shape[0], z5s.shape[0]
    V2 = [lin(zbs[i], V2b) for i in range(nb)]
    V5 = [lin(z5s[i], V5b) for i in range(n5)]
    H12 = [lin(zbs[i], H12b) for i in range(nb)]
    H15 = [lin(z5s[i], H15b) for i in range(n5)]
    H25 = [[lin([zbs[i, m]*z5s[j, k] for m in range(3) for k in range(3)], [H25b[m][k] for m in range(3) for k in range(3)])
            for j in range(n5)] for i in range(nb)]
    K = Gm.K
    one = ASer.const(1.0, K)

    def pw(x, kmax):
        out = [one]
        for _ in range(kmax):
            out.append(out[-1]*x)
        return out
    V1p = pw(V1, 2)
    V2p = [pw(v, 2) for v in V2]; V5p = [pw(v, J) for v in V5]; H12p = [pw(v, 2) for v in H12]; H15p = [pw(v, 2) for v in H15]
    H25p = [[pw(v, 2) for v in row] for row in H25]
    ns = sol.shape[1]
    res = [[[None]*n5 for _ in range(nb)] for _ in range(ns)]
    used = [i for i in range(len(mons)) if not np.all(np.abs(sol[i]) < 1e-14)]
    for ib_ in range(nb):
        for i5_ in range(n5):
            acc = [None]*ns
            for i in used:
                a1, a2, a5, b12, b15, b25 = mons[i]
                t = V1p[a1]*V2p[ib_][a2]*H12p[ib_][b12]*V5p[i5_][a5]*H15p[i5_][b15]*H25p[ib_][i5_][b25]
                for s_ in range(ns):
                    c = float(sol[i, s_])
                    if c != 0:
                        acc[s_] = t.scale(c) if acc[s_] is None else acc[s_] + t.scale(c)
            for s_ in range(ns):
                res[s_][ib_][i5_] = acc[s_]
    return res


def fourier2(arr, l1, l2):
    """arr[n1][n2] of ASer (samples) -> Fourier coefficients along both axes (n = -l..l)."""
    F1, F2 = dft(l1), dft(l2)
    n1, n2 = len(arr), len(arr[0])
    tmp = [[lin(F1[a], [arr[i][j] for i in range(n1)]) for j in range(n2)] for a in range(n1)]
    return [[lin(F2[b], [tmp[a][j] for j in range(n2)]) for b in range(n2)] for a in range(n1)]


def e_factor(Gm, W2):
    """E(z2, z4) samples (5 x 5): V042^2 (H24 - 2 V_{2,04} V_{4,02})/(X20 X40)."""
    B = Gm.basis
    v042 = Gm.V(0, Gm.z0, 4, 2)
    H24b = [[Gm.H(2, B[m], 4, B[k]) for k in range(3)] for m in range(3)]
    V204b = [Gm.V(2, B[m], 0, 4) for m in range(3)]
    V402b = [Gm.V(4, B[m], 0, 2) for m in range(3)]
    inv = Gm.inv_sq(2, 0)*Gm.inv_sq(4, 0)
    vv = v042*v042*inv
    n = W2.shape[0]
    V204 = [lin(W2[i], V204b) for i in range(n)]
    V402 = [lin(W2[i], V402b) for i in range(n)]
    E = [[vv*(lin([W2[i, m]*W2[j, k] for m in range(3) for k in range(3)], [H24b[m][k] for m in range(3) for k in range(3)])
               - (V204[i]*V402[j]).scale(2.0)) for j in range(n)] for i in range(n)]
    return E, v042


def ee_pieces(K, Delta, J, kind='EE', verbose=False):
    """EE k = 2 integrand on bg_ser: logs (as bg_geom.scalar_pieces, w0 _gs exponents), polys (v'-u', pref[0,0]) and R (k = 4)."""
    import l5_spin as S5
    mons, sol, _ = S5.conserved(float(Delta), int(J))
    _, solt, _ = S5.conserved(3.0 - float(Delta), int(J))
    logs_sc, polys_sc, G, nm = GG.scalar_pieces(K)
    Gm = EEGeom(K, G)
    W2, WJ = null_circle(2), null_circle(J)
    x10, x30 = sub(Gm.x[1], Gm.x[0]), sub(Gm.x[3], Gm.x[0])
    F = tto(Gm, 1, x10, 2, W2, 5, WJ, J, mons, sol)
    if verbose:
        print('   F done', flush=True)
    Gt = tto(Gm, 3, x30, 4, W2, 5, WJ, J, mons, solt)
    if verbose:
        print('   Gt done', flush=True)
    F = [fourier2(f, 2, J) for f in F]; Gt = [fourier2(g, 2, J) for g in Gt]
    E, v042 = e_factor(Gm, W2)
    D = [[E[i][j]*E[i][j] for j in range(5)] for i in range(5)]
    D = fourier2(D, 2, 2)
    k2, kJ = kappa(2), kappa(J)
    n5 = 2*J + 1
    na, nc = len(F), len(Gt)
    pref = [[None]*nc for _ in range(na)]
    for a in range(na):
        Fa = F[a]
        Em = [[lin([1.0/k2[i] for i in range(5)], [D[i][j]*Fa[4 - i][k] for i in range(5)]) for k in range(n5)] for j in range(5)]
        for c in range(nc):
            Gc = Gt[c]
            terms = []; co = []
            for j in range(5):
                for k in range(n5):
                    terms.append(Em[j][k]*Gc[4 - j][n5 - 1 - k]); co.append(1.0/(k2[j]*kJ[k]))
            pref[a][c] = lin(co, terms)
    if verbose:
        print('   pref done', [len(p) for p in pref[0][0].c], flush=True)
    # slow part (w0 _gs): scalar logs without 'sub' and with the EE 3-pt exponents; f = v' - u' (power 1); -4 log V_{0,42}
    keep = ('A', 'xsq', 'FP1', 'FP2', 'FP3', 'dual', 'cel', 'meas')
    logs = [l for l in logs_sc if l[0] in keep]
    logs += [('3pt', nm['m12'], -(10 - J)/2, 0.0), ('3pt', nm['m15'], -J/2, 0.0), ('3pt', nm['m25'], -J/2, 0.0),
             ('3pt', nm['m34'], -(7 - J)/2, 0.0), ('3pt', nm['m35'], -(3 + J)/2, 0.0), ('3pt', nm['m45'], -(3 + J)/2, 0.0),
             ('v042', v042, -4.0, 0.0)]
    vmu = nm['vp'] - nm['up']
    Rk4 = (vmu*vmu - 1.0)*(nm['up']*nm['vp']).recip('u\'v\'').scale(0.25)
    return dict(logs=logs, vmu=vmu, pref=pref, R=Rk4, G=G, nm=nm)
