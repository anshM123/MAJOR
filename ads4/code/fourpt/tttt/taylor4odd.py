"""
Taylor-mode application of PARITY-ODD differential words to the scalar seed (extension of taylor4.Engine4; PARITY_ODD_HEAVY.md sect. 4).

State: dict key -> series, key = zexp (18 exponents of the z-symbols, as taylor4) + (eps,), eps = canonical tuple of epsilon factors.
An epsilon factor is a tuple of 5 slots; a slot is a vector name ('P1'..'Z4') or a contraction tag ('t', n) shared by exactly two factors
(meaning the sum eps(.., e_A, ..) eta^{AB} eps(.., e_B, ..)).  A pending factor (one derivative slot still open) is kept apart while the
second derivative of a two-derivative operator acts.

Odd operators at a vertex (i, j) (letters as odd3; the right vertex uses (3, 4)):
   'O' = eps(Zi, Zj, Pi, Pj, d/dPi)    'Q' = eps(Zi, Zj, Pi, Pj, d/dPj)          (first order, applied first: seed is Z-independent)
   'R' = eps(Zi, Pi, Pj, d/dPi, d/dPj) 'S' = eps(Zj, Pi, Pj, d/dPi, d/dPj)     (second order, applied first)
   'U' = eps(Zi, Pi, Pj, Zj, d/dZj)    'V' = eps(Zj, Pj, Pi, Zi, d/dZi)        (first order in Z, gauge covariant)
Even operators (A, B, C, E, H) act on epsilon slots by replacement (V.d/dW: slot W -> V; a repeated vector kills the factor).
d/dP of a symbol sigma(P, X) puts the partner X (factor 2 for s = 2P.P) into the derivative slot; d/dP of an epsilon factor holding P
in a slot creates a contraction tag.

Evaluation (Gram form, valid for the Euclidean-defined structures continued to any kinematics):
   eps(a1..a_{5-k}, T1..Tk) eps(b1..b_{5-k}, T1..Tk) = k! det[a_i.b_j] / det(eta_E),   det(eta_E) = -1/4
(eta_E: embedding metric of the Euclidean frame used for the three-point matching in odd3, P.P' = -(x-x')^2/2).  A single factor
(waves W[E,O]) is evaluated as the determinant of the Lorentzian embedding components (P = (1, x^2, x), Z = (0, 2x.z, z)).
"""
import itertools, math
import numpy as np
from taylor4 import Engine4, SN, ZN, SIDX, ZIDX, PAIR, dot, st_add, sderiv, smulvar, tables

DET_ETA_E = -0.25
VEC = ['P1', 'P2', 'P3', 'P4', 'Z1', 'Z2', 'Z3', 'Z4']

def _sym_partner(name, W):
    a, b = PAIR[name]
    if W == a: return b
    if W == b: return a
    return None

def _canon(eps):
    """canonical form: relabel tags in order of appearance; sort factors (tuple order is kept inside a factor)"""
    fs = sorted(eps, key=lambda f: tuple(str(s) for s in f))
    ren = {}; out = []
    for f in fs:
        g = []
        for s in f:
            if isinstance(s, tuple):
                if s not in ren: ren[s] = ('t', len(ren))
                g.append(ren[s])
            else:
                g.append(s)
        out.append(tuple(g))
    return tuple(out)

class Engine4Odd(Engine4):
    def __init__(self, M, s0, dtype=float):
        super().__init__(M, s0, dtype)
        self.ntag = 0

    def _newtag(self):
        self.ntag += 1
        return ('t', 1000 + self.ntag)

    # --- keys: zexp (18) + (eps,)
    @staticmethod
    def zx(key): return key[:18]
    @staticmethod
    def ep(key): return key[18] if len(key) > 18 else ()

    def mul_sym(self, A, name, fac):
        out = {}
        if name in SIDX:
            k = SIDX[name]
            for key, v in A.items():
                out[key] = smulvar(self.T, v, k, self.s0[k])*fac
        else:
            m = ZIDX[name]
            for key, v in A.items():
                kk = list(key); kk[m] += 1; kk = tuple(kk)
                out[kk] = out[kk] + v*fac if kk in out else v*fac
        return out

    def _eps_replace(self, A, V, W):
        """V.d/dW acting on the epsilon factors only"""
        out = {}
        for key, v in A.items():
            eps = self.ep(key)
            for fi, f in enumerate(eps):
                for si, s in enumerate(f):
                    if s != W: continue
                    if V in f and V != W: continue
                    g = list(f); g[si] = V
                    new = list(eps); new[fi] = tuple(g)
                    nk = self.zx(key) + (_canon(new),)
                    out[nk] = out[nk] + v if nk in out else v
        return out

    def dirder(self, A, V, W):
        out = super().dirder(A, V, W)        # z- and s-symbols (keys keep their eps part: base code copies the key tuple)
        if any(self.ep(k) for k in A):
            out = st_add(out, self._eps_replace(A, V, W))
        return out

    def _deriv_into_slot(self, A, W, base_eps_fn):
        """d/dW with the derivative index placed into a new (pending) factor: returns list of (state, slotvalue) pieces as a state whose
        keys carry ('pend', slotvalue) as an extra element.  base_eps_fn not used (kept for symmetry)."""
        out = {}
        # (a) symbols
        for name in list(SIDX) + list(ZIDX):
            X = _sym_partner(name, W)
            if X is None: continue
            fac = 2.0 if name in SIDX else 1.0
            d = self.d_sym(A, name)
            for key, v in d.items():
                nk = key + (('pend', X),) if len(key) == 19 else key[:18] + ((),) + (('pend', X),)
                if len(key) == 19: nk = key + (('pend', X),)
                out[nk] = out[nk] + fac*v if nk in out else fac*v
        # (b) existing epsilon factors holding W
        for key, v in A.items():
            eps = self.ep(key)
            for fi, f in enumerate(eps):
                for si, s in enumerate(f):
                    if s != W: continue
                    t = self._newtag()
                    g = list(f); g[si] = t
                    new = list(eps); new[fi] = tuple(g)
                    nk = self.zx(key) + (tuple(new),) + (('pend', t),)
                    out[nk] = out[nk] + v if nk in out else v
        return out

    def odd_op(self, A, kind, i, j):
        Pi, Pj, Zi, Zj = 'P%d' % i, 'P%d' % j, 'Z%d' % i, 'Z%d' % j
        A = {(k if len(k) == 19 else k + ((),)): v for k, v in A.items()}
        if kind in ('O', 'Q', 'U', 'V'):
            fixed = {'O': [Zi, Zj, Pi, Pj], 'Q': [Zi, Zj, Pi, Pj], 'U': [Zi, Pi, Pj, Zj], 'V': [Zj, Pj, Pi, Zi]}[kind]
            W = {'O': Pi, 'Q': Pj, 'U': Zj, 'V': Zi}[kind]
            D = self._deriv_into_slot(A, W, None)
            out = {}
            for key, v in D.items():
                slot = key[19][1]
                f = tuple(fixed + [slot])
                if not isinstance(slot, tuple) and slot in fixed: continue
                eps = tuple(key[18]) + (f,)
                nk = key[:18] + (_canon(eps),)
                out[nk] = out[nk] + v if nk in out else v
            return out
        if kind in ('R', 'S'):
            Zs = Zi if kind == 'R' else Zj
            fixed = [Zs, Pi, Pj]
            D1 = self._deriv_into_slot(A, Pi, None)                 # keys: zexp + (eps,) + (('pend', X),)
            out = {}
            # second derivative d/dPj acting on everything except the pending slot
            groups = {}
            for key, v in D1.items():
                groups.setdefault(key[19], {})[key[:19]] = v
            for pend, St in groups.items():
                X = pend[1]
                if not isinstance(X, tuple) and X in fixed: continue
                D2 = self._deriv_into_slot(St, Pj, None)
                for key, v in D2.items():
                    Y = key[19][1]
                    if not isinstance(Y, tuple) and (Y in fixed or Y == X): continue
                    f = tuple(fixed + [X, Y])
                    eps = tuple(key[18]) + (f,)
                    nk = key[:18] + (_canon(eps),)
                    out[nk] = out[nk] + v if nk in out else v
            return out
        raise ValueError(kind)

    def op(self, A, kind, i, j):
        if kind in 'OQRSUV':
            return self.odd_op(A, kind, i, j)
        A = {(k if len(k) == 19 else k + ((),)): v for k, v in A.items()}
        return super().op(A, kind, i, j)

def word_budget_odd(wL, wR):
    n = 0
    for c in wL + wR:
        if c in 'HUV': continue
        n += 2 if c in 'RS' else 1
    return n

def seed_shifts(word):
    s1 = sum(word.count(c) for c in 'AEHQUV'); s2 = sum(word.count(c) for c in 'BCHOUV')
    return s1, s2

# ------------------------------------------------------------------ evaluation of the epsilon factors on the polarisation grid
def gram_value(eps, dots, embed=None):
    """value of the epsilon product eps (tuple of factors) given dots(u, v) -> array; embed(name) -> component array [..., 5] for a single
    factor (Lorentzian determinant)."""
    if len(eps) == 0:
        return 1.0
    if len(eps) == 1:
        f = eps[0]
        assert not any(isinstance(s, tuple) for s in f)
        M = np.stack([embed(s) for s in f], axis=-2)
        return np.linalg.det(M)
    assert len(eps) == 2, eps
    f, g = eps
    tags = [s for s in f if isinstance(s, tuple)]
    assert sorted(tags) == sorted(s for s in g if isinstance(s, tuple))
    # move tags to the end, in a common order
    def reorder(h):
        free = [s for s in h if not isinstance(s, tuple)]
        perm = [h.index(s) for s in free] + [h.index(t) for t in tags]
        sign = 1
        for a in range(len(perm)):
            for b in range(a + 1, len(perm)):
                if perm[a] > perm[b]: sign = -sign
        return free, sign
    fa, sa = reorder(f); fb, sb = reorder(g)
    k = len(tags)
    n = 5 - k
    G = np.stack([np.stack([np.broadcast_to(dots(a, b), np.broadcast_shapes(*(np.shape(dots(x, y)) for x in fa for y in fb)))
                            for b in fb], axis=-1) for a in fa], axis=-2)
    return sa*sb*math.factorial(k)*np.linalg.det(G)/DET_ETA_E
