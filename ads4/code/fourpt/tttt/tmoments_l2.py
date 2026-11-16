"""
<TTTT> exact partial waves by Taylor mode (taylor4) on the moment grid of e2e/moments.py, and their tensor Regge moments
(same grid, polarisation grid, Fourier modes and kernel as fourpt/tmoments4 for <JJJJ>).

Waves:  W[X_hL, X_hR] for hL, hR in {same, opp}, X_h = CHL's helicity structures (chl_xspace) with exact word coefficients at (Delta, J)
(t8_kappa_exact.decompose; the d = 3 relation coefficient set to 0, irrelevant since the relation vanishes identically).
Per point:  s13^3 s24^3 W = exp(logscale) * sum_m C[m] * monomial_m(z),  logscale = Delta log(4r) + 3 (log s13 + log s24),
C = (-1)^J D_wL D_wR [Kt gt/(4r)^Delta] (Taylor mode),  gt jets from cblock12 (double radial sums) and taylor4.pq_jets_vec.
Usage: python3 tmoments_l2.py beta m [dtype=float64|longdouble] [chunk]     (cache: cache/modes_l2_*.npz, log to stdout)
"""
import os, sys, time, ctypes, itertools, functools
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
import blocks
from blocks import CBlock
from moments import make_grid, frame_data, frame_positions, series_converged, series_X
from cases import cases_fixed_beta
from moments4 import pol_grid4, fourier4, kernel_np
from spinblock2 import invariants4
from embed2 import shifts
from structures2 import words
from taylor4 import Engine4, seed, word_budget, ZN, SN, tables, _pmul2
CACHE = os.path.join(HERE, 'cache')
L, DX = 2, 3

_L12 = None
def lib12():
    global _L12
    if _L12 is None:
        L_ = ctypes.CDLL(os.path.join(HERE, 'libcblock12_d.so')); L0 = blocks._lib('d')
        L_.block_eval.argtypes = L0.block_eval.argtypes
        _L12 = L_
    return _L12

def series12(blk, r, chi, kmax, lmax, cut=30.0, pwcut=10.0):
    """CBlock(prec='d').series with the order-11 library (float path)"""
    r = np.ascontiguousarray(r, dtype=float); chi = np.ascontiguousarray(chi, dtype=float)
    r2 = np.zeros((r.size, 2)); r2[:, 0] = r; c2 = np.zeros((r.size, 2)); c2[:, 0] = chi
    out = np.zeros((r.size, kmax + 1, lmax + 1, 2))
    lib12().block_eval(blk.ptr, blk.off.ctypes.data, blk.D2.ctypes.data, blk.J, blk.N, r.size, r2.ctypes.data, c2.ctypes.data,
                       1, int(kmax), int(lmax), float(cut), float(pwcut), out.ctypes.data)
    return out[..., 0] + out[..., 1]

@functools.lru_cache(None)
def _acosh_np(M):
    c = sp.Symbol('c')
    return [sp.lambdify(c, sp.simplify(sp.diff(sp.acosh(c), c, k)), 'numpy') for k in range(1, M + 1)]

def pq_jets_vec(S, r, A0, E0, M, vpow=0.0):
    """vectorised taylor4.pq_jets_M (float64): S[npts, k, l] -> G[npts, i, j]"""
    import math
    npts = S.shape[0]
    Z = lambda: np.zeros((npts, M + 1, M + 1))
    dA = Z(); dE = Z()
    for i in range(M + 1):
        for j in range(M + 1 - i):
            if i + j > 0:
                dA[:, i, j] = A0*(-0.5)**i*0.5**j/(math.factorial(i)*math.factorial(j))
        if i > 0:
            dE[:, 0, i] = E0/math.factorial(i)
    dc = dA + dE; dx = dE - dA
    c0 = (r + 1/r)/2
    fk = [f(c0)*np.ones(npts) for f in _acosh_np(M)]
    dL = Z(); p = Z(); p[:, 0, 0] = 1
    for k in range(1, M + 1):
        p = _pmul2(p, dc, M)
        dL = dL - p*(fk[k - 1]/math.factorial(k))[:, None, None]
    one = Z(); one[:, 0, 0] = 1
    powL = [one]; powX = [one]
    for _ in range(M):
        powL.append(_pmul2(powL[-1], dL, M)); powX.append(_pmul2(powX[-1], dx, M))
    G = Z()
    for k in range(M + 1):
        for l in range(M + 1 - k):
            G = G + _pmul2(powL[k], powX[l], M)*(S[:, k, l]/(math.factorial(k)*math.factorial(l)))[:, None, None]
    if vpow != 0:
        ex = Z()
        for i in range(M + 1):
            for j in range(M + 1 - i):
                ex[:, i, j] = vpow**(i + j)/(math.factorial(i)*math.factorial(j))
        G = _pmul2(G, ex, M)*((E0/A0)**(2*vpow))[:, None, None]
    for i in range(M + 1):
        for j in range(M + 1 - i):
            G[:, i, j] *= math.factorial(i)*math.factorial(j)
    return G

def helicity_words(Dv, Jv):
    fn = os.path.join(CACHE, 'xwords_%.10f_%d.npz' % (Dv, Jv))
    if os.path.exists(fn):
        d = np.load(fn); return {'same': dict(zip(words(L), d['same'])), 'opp': dict(zip(words(L), d['opp']))}
    from t8_kappa_exact import decompose
    res = {}
    for kind in ('same', 'opp'):
        ws, sol, params = decompose(L, sp.Rational(Dv), sp.Integer(Jv), kind)
        sol = sol.subs({p: 0 for p in params})
        res[kind] = np.array([complex(sp.N(sol[i], 30)) for i in range(len(ws))])
    np.savez(fn, **res)
    return {k: dict(zip(words(L), v)) for k, v in res.items()}

COMBOS = [('same', 'same'), ('same', 'opp'), ('opp', 'same'), ('opp', 'opp')]

class TaylorWaves:
    def __init__(self, Delta, J, N, dt=np.float64):
        self.Delta, self.J, self.N, self.dt = Delta, J, N, dt
        self.ws = words(L)
        self.hw = helicity_words(Delta, J)
        self.blk = {}; self.pairs = {}
        for wL in self.ws:
            for wR in self.ws:
                s1, s2 = shifts(wL); s3, s4 = shifts(wR)
                d = (DX + s1, DX + s2, DX + s3, DX + s4)
                a = -(d[0] - d[1])/2; b = (d[2] - d[3])/2
                self.pairs.setdefault((a, b), []).append((wL, wR, d))
                if (a, b) not in self.blk:
                    c1, c2 = -a, b
                    self.blk[(a, b)] = (CBlock(Delta, J, c1, c2, N, prec='d'), 0.0) if c1 + c2 >= 0 else (CBlock(Delta, J, -c1, -c2, N, prec='d'), c1 + c2)

    def chunk(self, s, rs, chi, A0, E0):
        """returns dict combo -> dict monomial -> array[npts] of C (sign (-1)^J included)"""
        dt = self.dt
        acc = {c: {} for c in COMBOS}
        sv = {n: s[n].astype(dt) for n in SN}
        for ab, plist in self.pairs.items():
            bl, vpow = self.blk[ab]
            Mmax = max(word_budget(wL, wR) for wL, wR, d in plist)
            S = series12(bl, rs, chi, Mmax, Mmax)
            G = pq_jets_vec(S, rs, A0, E0, Mmax, vpow)
            d0 = plist[0][2]
            assert all(d == d0 for (_, _, d) in plist)
            Gl = [[G[:, i, j].astype(dt) for j in range(Mmax + 1)] for i in range(Mmax + 1)]
            sd = seed(Mmax, sv, tuple(dt(x) for x in d0), Gl, dtype=dt, one=dt(1))[(0,)*18]
            nb = tables(Mmax)['nb']
            for (wL, wR, d) in plist:
                M = word_budget(wL, wR)
                E = Engine4(M, sv)
                res = E.apply_words({(0,)*18: sd[:nb[M]]}, wL, wR)
                for c in COMBOS:
                    w = self.hw[c[0]][wL]*self.hw[c[1]][wR]
                    if w == 0: continue
                    A = acc[c]
                    for mono, v in res.items():
                        if mono in A:
                            A[mono][0] += w.real*v; A[mono][1] += w.imag*v
                        else:
                            A[mono] = [w.real*v, w.imag*v]
        sg = (-1)**self.J
        return {c: {mono: sg*(v[0].astype(float) + 1j*v[1].astype(float)) for mono, v in A.items()} for c, A in acc.items()}

def contract(Cd, z, npts):
    """sum_m C[m] prod z^e  on the polarisation grid: returns [npts, 5,5,5,5] complex"""
    pw = {}
    out = 0
    for mono, c in Cd.items():
        v = c.reshape((npts,) + (1,)*4).astype(complex)
        for n, e in zip(ZN, mono):
            if e:
                if (n, e) not in pw: pw[(n, e)] = z[n]**e
                v = v*pw[(n, e)]
        out = out + v
    return out

def compute(beta, mm, dt=np.float64, chunk=200, N=None, maxpts=None):
    (Dv, Jv), = cases_fixed_beta(beta, [mm])
    N = N or int(min(8000, max(4000, 400*Dv)))
    tag = 'modes_l2_D%.6f_J%d_%s.npz' % (Dv, Jv, 'ld' if dt == np.longdouble else 'd')
    fn = os.path.join(CACHE, tag)
    g = make_grid(Dv, Jv); ok = g['ok']
    fd = frame_data(g['U'][ok], g['S'][ok])
    conv = series_converged(fd['rs'], fd['chi'], Jv, N, cut=25.0)
    convfull = np.zeros(g['S'].shape, dtype=bool); convfull[ok] = conv
    fd = {k: v[conv] for k, v in fd.items()}
    npts = fd['rs'].size if maxpts is None else min(maxpts, fd['rs'].size)
    print(f"# <TTTT> Taylor-mode waves: beta={beta} m~{mm}: Delta={Dv:.6f} J={Jv} N={N} dtype={np.dtype(dt).name}; {npts} grid points", flush=True)
    t0 = time.time()
    TW = TaylorWaves(Dv, Jv, N, dt)
    print(f"  setup (helicity words, 25 seed blocks): {time.time() - t0:.0f} s", flush=True)
    modes = {c: {} for c in COMBOS}
    logscale = np.zeros(npts)
    for c0 in range(0, npts, chunk):
        t = time.time()
        sl = slice(c0, min(npts, c0 + chunk)); n = sl.stop - sl.start
        x1, x2, x3, x4 = frame_positions(fd['ru'][sl], fd['xi'][sl])
        z1, z3, z2, z4 = pol_grid4(fd['xi'][sl], L)
        s, z = invariants4({1: x1, 2: x2, 3: x3, 4: x4}, {1: z1, 2: z2, 3: z3, 4: z4})
        Cs = TW.chunk(s, fd['rs'][sl], fd['chi'][sl], fd['A0'][sl], fd['E0'][sl])
        logscale[sl] = Dv*np.log(4*fd['rs'][sl]) + DX*(np.log(s['s13']) + np.log(s['s24']))
        for c in COMBOS:
            Mv = contract(Cs[c], z, n).reshape(n, -1).reshape((n,) + (5,)*4)
            for q, v in fourier4(Mv, L).items():
                modes[c].setdefault(q, []).append(v)
        print(f"  points {sl.start}-{sl.stop}: {time.time() - t:.0f} s (total {time.time() - t0:.0f} s)", flush=True)
    keys, arrs = [], []
    for ci, c in enumerate(COMBOS):
        for q, lst in modes[c].items():
            arr = np.zeros(g['S'].shape, dtype=complex)
            full = np.concatenate(lst)
            idx = np.flatnonzero(convfull.ravel())[:npts]
            arr.ravel()[idx] = full
            keys.append((ci,) + q); arrs.append(arr)
    ls = np.full(g['S'].shape, -np.inf)
    ls.ravel()[np.flatnonzero(convfull.ravel())[:npts]] = logscale
    conv2 = np.zeros(g['S'].shape, dtype=bool); conv2.ravel()[np.flatnonzero(convfull.ravel())[:npts]] = True
    np.savez_compressed(fn, keys=np.array(keys), modes=np.array(arrs), us=g['us'], ts=g['ts'], U=g['U'], S=g['S'], W=g['W'], ok=g['ok'],
                        conv=conv2, logscale=ls, Delta=Dv, J=Jv, N=N, npts=npts, time=time.time() - t0)
    print(f"  saved {fn} ({time.time() - t0:.0f} s)", flush=True)
    return fn

if __name__ == "__main__":
    beta = float(sys.argv[1]); mm = float(sys.argv[2])
    dt = np.longdouble if (len(sys.argv) > 3 and sys.argv[3] == 'longdouble') else np.float64
    chunk = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    maxpts = int(sys.argv[5]) if len(sys.argv) > 5 else None
    compute(beta, mm, dt, chunk, maxpts=maxpts)
