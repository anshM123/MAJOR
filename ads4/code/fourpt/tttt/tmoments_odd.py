"""
<TTTT> waves with a PARITY-ODD exchanged structure (PARITY_ODD_HEAVY.md sect. 4), Taylor mode, on the moment grid of e2e/moments.py.

Structures at each vertex: E = CHL's parity-even opposite-helicity structure X^E_opp (even words, coefficients of t8_kappa_exact.decompose
as in tmoments_l2), O = CHL's parity-odd opposite-helicity structure X^O_opp = i[T^{+-} - T^{-+}] (odd words of odd3, real coefficients
c_r with X^O = i sum c_r F_w; the factor i is dropped: the wave uses the real-coefficient structure sum c_r F_w).
Waves:  W[E,E], W[E,O], W[O,E], W[O,O]  (left vertex (12), right vertex (34)).  Same polarisation grid, Fourier modes, cache format as
tmoments_l2 (combo labels ('same','same') -> E,E etc. are NOT used: keys (ci,) + q with ci in 0..3 = EE, EO, OE, OO).
Resumable: points are processed in chunks and every chunk is saved (cache/odd_<tag>_chunks/), the final npz is assembled at the end.
Usage: python3 tmoments_odd.py J [chunk] [maxpts]      (beta = 1, (Delta, J) = e2e/cases.case(1, J); J odd for the identical-T case)
"""
import os, sys, time, pickle
import numpy as np, sympy as sp, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from blocks import CBlock
from moments import make_grid, frame_data, frame_positions, series_converged
from cases import case
from moments4 import pol_grid4, fourier4
from spinblock2 import invariants4
from embed2 import shifts as shifts_even
from structures2 import words as even_words
from taylor4 import seed, tables, ZN, SN
from taylor4odd import Engine4Odd, word_budget_odd, seed_shifts, gram_value
from tmoments_l2 import series12, pq_jets_vec, helicity_words
import odd3
CACHE = os.path.join(HERE, 'cache')
L, DX = 2, 3
COMBOS = [('E', 'E'), ('E', 'O'), ('O', 'E'), ('O', 'O')]

# basis of the 12-dim space of parity-odd (2,2,J) structures reached by the odd words, chosen greedily by Taylor budget (odd3; rank 12)
ODD_BASIS = ['HO', 'HQ', 'ACO', 'ACQ', 'AEO', 'AEQ', 'BCO', 'BCQ', 'AACS', 'AAES', 'ABCS', 'BBCS']

def odd_coeffs(Dv, Jv, kind='opp'):
    fn = os.path.join(CACHE, 'oddwords_%.10f_%d_%s.pkl' % (Dv, Jv, kind))
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    ws, c, res, res2, rk = odd3.decompose_odd(Dv, Jv, kind, npts=50, dps=60, words=ODD_BASIS)
    out = dict(coef={w: float(mp.im(x)) for w, x in zip(ws, c) if abs(x) > 1e-25}, fit=float(res), fresh=float(res2), rank=rk)
    pickle.dump(out, open(fn, 'wb'))
    return out

class OddWaves:
    def __init__(self, Delta, J, N, kind='opp', only=None):
        self.Delta, self.J, self.N = Delta, J, N
        hwE = helicity_words(Delta, J)[kind]
        self.co = {'E': {w: complex(v).real for w, v in hwE.items() if abs(v) > 0},
                   'O': odd_coeffs(Delta, J, kind)['coef']}
        assert all(abs(complex(v).imag) < 1e-12*max(1, abs(v)) for v in hwE.values())
        self.combos = only or COMBOS
        self.pairs = {}; self.blk = {}
        for cL, cR in self.combos:
            for wL in self.co[cL]:
                for wR in self.co[cR]:
                    sL = shifts_even(wL) if cL == 'E' else seed_shifts(wL)
                    sR = shifts_even(wR) if cR == 'E' else seed_shifts(wR)
                    d = (DX + sL[0], DX + sL[1], DX + sR[0], DX + sR[1])
                    a = -(d[0] - d[1])/2; b = (d[2] - d[3])/2
                    self.pairs.setdefault((a, b), []).append((cL, cR, wL, wR, d))
                    if (a, b) not in self.blk:
                        c1, c2 = -a, b
                        self.blk[(a, b)] = ((c1, c2), 0.0) if c1 + c2 >= 0 else ((-c1, -c2), c1 + c2)
        self.S = {}

    def budget(self, wL, wR, cL, cR):
        return word_budget_odd(wL, wR)

    def precompute(self, rs, chi):
        for ab, ((c1, c2), vpow) in self.blk.items():
            Mmax = max(self.budget(wL, wR, cL, cR) for cL, cR, wL, wR, d in self.pairs[ab])
            bl = CBlock(self.Delta, self.J, c1, c2, self.N, prec='d')
            self.S[ab] = series12(bl, rs, chi, Mmax, Mmax)
            del bl

    def chunk(self, s, rs, chi, A0, E0, sl):
        """factorised: within a class of equal seed dimensions (d1,d2,d3,d4), the left words are applied once each and summed with their
        coefficients, then every right word is applied to the sum (cost |L| + |R| instead of |L| x |R| per class)."""
        from taylor4 import st_add
        acc = {c: {} for c in self.combos}
        sv = {n: s[n] for n in SN}
        for ab, plist in self.pairs.items():
            _, vpow = self.blk[ab]
            classes = {}
            for (cL, cR, wL, wR, d) in plist:
                classes.setdefault((cL, cR, d), (set(), set()))
                classes[(cL, cR, d)][0].add(wL); classes[(cL, cR, d)][1].add(wR)
            Mab = max(max(word_budget_odd(w, '') for w in Ls) + max(word_budget_odd(w, '') for w in Rs) for (Ls, Rs) in classes.values())
            G = pq_jets_vec(self.S[ab][sl], rs, A0, E0, Mab, vpow)
            Gl = [[G[:, i, j] for j in range(Mab + 1)] for i in range(Mab + 1)]
            for (cL, cR, d), (Ls, Rs) in classes.items():
                M = max(word_budget_odd(w, '') for w in Ls) + max(word_budget_odd(w, '') for w in Rs)
                sd = seed(M, sv, d, [row[:M + 1 - i] for i, row in enumerate(Gl[:M + 1])])[(0,)*18]
                nb = tables(M)['nb']
                E = Engine4Odd(M, sv)
                AL = {}
                for wL in sorted(Ls):
                    St = {(0,)*18 + ((),): sd[:nb[M]]}
                    for kind in reversed(wL): St = E.op(St, kind, 1, 2)
                    AL = st_add(AL, St, self.co[cL][wL]) if AL else {k: v*self.co[cL][wL] for k, v in St.items()}
                A = acc[(cL, cR)]
                for wR in sorted(Rs):
                    St = AL
                    for kind in reversed(wR): St = E.op(St, kind, 3, 4)
                    for key, v in St.items():
                        if not v.shape[0]: continue
                        val = self.co[cR][wR]*v[0]
                        A[key] = A[key] + val if key in A else val
        sg = (-1)**self.J
        return {c: {k: sg*v for k, v in A.items()} for c, A in acc.items()}

def contract(Cd, s, z, x, zpol, npts):
    """sum over keys C * prod z^e * eps-value on the polarisation grid -> [npts, 5,5,5,5]"""
    shp = (npts,) + (5,)*4
    def dots(u, v):
        if u[0] == 'P' and v[0] == 'P':
            i, j = sorted((int(u[1]), int(v[1])))
            return 0.0 if i == j else s['s%d%d' % (i, j)].reshape((npts,) + (1,)*4)/2
        if u[0] == 'P': u, v = v, u
        if v[0] == 'P':
            i, j = int(u[1]), int(v[1])
            return np.zeros(shp) if i == j else z['z%d%d' % (i, j)]
        i, j = sorted((int(u[1]), int(v[1])))
        return np.zeros(shp) if i == j else z['z%dz%d' % (i, j)]
    def embed(name):
        i = int(name[1]); xi = x[i]
        xx = (-xi[:, 0]**2 + xi[:, 1]**2 + xi[:, 2]**2).reshape((npts,) + (1,)*4)
        if name[0] == 'P':
            comp = [np.ones(shp), np.broadcast_to(xx, shp)] + [np.broadcast_to(xi[:, k].reshape((npts,) + (1,)*4), shp) for k in range(3)]
        else:
            zz = zpol[i]
            xz = (-xi[:, 0].reshape((npts,) + (1,)*4)*zz[..., 0] + xi[:, 1].reshape((npts,) + (1,)*4)*zz[..., 1]
                  + xi[:, 2].reshape((npts,) + (1,)*4)*zz[..., 2])
            comp = [np.zeros(shp), 2*xz] + [zz[..., k] for k in range(3)]
        return np.stack(comp, axis=-1)
    pw = {}; ev = {}
    out = 0
    for key, c in Cd.items():
        mono, eps = key[:18], key[18]
        v = c.reshape((npts,) + (1,)*4).astype(complex)
        for n, e in zip(ZN, mono):
            if e:
                if (n, e) not in pw: pw[(n, e)] = z[n]**e
                v = v*pw[(n, e)]
        if eps:
            if eps not in ev: ev[eps] = gram_value(eps, dots, embed)
            v = v*ev[eps]
        out = out + v
    return out

GRIDS = {'full': {}, 'half': dict(upanels=((0.0, 0.01, 4), (0.01, 0.03, 5), (0.03, 0.1, 7), (0.1, 0.3, 8), (0.3, 1.0, 12)),
                                   tpanels=((0.0, 6.0, 8), (6.0, 28.0, 13), (28.0, 85.0, 10)))}

def compute(J, chunk=100, maxpts=None, only=None, grid='full', worker=0, nworkers=1):
    Dv, Jv = case(1.0, J)
    N = int(min(8000, max(4000, 400*Dv)))
    tag = 'odd_D%.6f_J%d' % (Dv, Jv) + ('' if only is None else '_' + '_'.join(a + b for a, b in only)) + ('' if maxpts is None else '_n%d' % maxpts) + ('' if grid == 'full' else '_' + grid)
    cdir = os.path.join(CACHE, tag + '_chunks'); os.makedirs(cdir, exist_ok=True)
    g = make_grid(Dv, Jv, **GRIDS[grid]); ok = g['ok']
    fd = frame_data(g['U'][ok], g['S'][ok])
    conv = series_converged(fd['rs'], fd['chi'], Jv, N, cut=25.0)
    convfull = np.zeros(g['S'].shape, dtype=bool); convfull[ok] = conv
    fd = {k: v[conv] for k, v in fd.items()}
    npts = fd['rs'].size if maxpts is None else min(maxpts, fd['rs'].size)
    print(f"# <TTTT> odd-exchange waves: beta=1 Delta={Dv:.6f} J={Jv} N={N}; {npts} grid points", flush=True)
    t0 = time.time()
    OW = OddWaves(Dv, Jv, N, only=only)
    oc = odd_coeffs(Dv, Jv)
    print(f"  odd structure: {len(oc['coef'])} odd words, 3-pt fit residual {oc['fit']:.1e}, fresh {oc['fresh']:.1e}, rank {oc['rank']}", flush=True)
    print(f"  word pairs: " + ", ".join(f"{c}: {sum(1 for p in OW.pairs.values() for q in p if (q[0], q[1]) == c)}" for c in OW.combos), flush=True)
    OW.precompute(fd['rs'][:npts], fd['chi'][:npts])
    print(f"  setup: {time.time() - t0:.0f} s ({len(OW.blk)} seed blocks)", flush=True)
    for c0 in range(0, npts, chunk):
        if (c0//chunk) % nworkers != worker: continue
        fnc = os.path.join(cdir, 'c%06d.pkl' % c0)
        if os.path.exists(fnc): continue
        t = time.time()
        sl = slice(c0, min(npts, c0 + chunk)); n = sl.stop - sl.start
        x1, x2, x3, x4 = frame_positions(fd['ru'][sl], fd['xi'][sl])
        z1, z3, z2, z4 = pol_grid4(fd['xi'][sl], L)
        s, z = invariants4({1: x1, 2: x2, 3: x3, 4: x4}, {1: z1, 2: z2, 3: z3, 4: z4})
        Cs = OW.chunk(s, fd['rs'][sl], fd['chi'][sl], fd['A0'][sl], fd['E0'][sl], sl)
        modes = {}
        for ci, c in enumerate(OW.combos):
            Mv = contract(Cs[c], s, z, {1: x1, 2: x2, 3: x3, 4: x4}, {1: z1, 2: z2, 3: z3, 4: z4}, n)
            for q, v in fourier4(Mv.reshape((n,) + (5,)*4), L).items():
                modes[(ci,) + q] = v
        logscale = Dv*np.log(4*fd['rs'][sl]) + DX*(np.log(s['s13']) + np.log(s['s24']))
        pickle.dump(dict(modes=modes, logscale=logscale), open(fnc, 'wb'))
        print(f"  points {sl.start}-{sl.stop}: {time.time() - t:.0f} s (total {time.time() - t0:.0f} s)", flush=True)
    # assemble (only when every chunk exists)
    if not all(os.path.exists(os.path.join(cdir, 'c%06d.pkl' % c0)) for c0 in range(0, npts, chunk)):
        print("  chunks missing (other worker): not assembled", flush=True); return None
    keys = None; parts = {}; ls = []
    for c0 in range(0, npts, chunk):
        d = pickle.load(open(os.path.join(cdir, 'c%06d.pkl' % c0), 'rb'))
        ls.append(d['logscale'])
        for k, v in d['modes'].items(): parts.setdefault(k, []).append(v)
    idx = np.flatnonzero(convfull.ravel())[:npts]
    keys, arrs = [], []
    for k, lst in parts.items():
        arr = np.zeros(g['S'].shape, dtype=complex); arr.ravel()[idx] = np.concatenate(lst)
        keys.append(k); arrs.append(arr)
    L_ = np.full(g['S'].shape, -np.inf); L_.ravel()[idx] = np.concatenate(ls)
    conv2 = np.zeros(g['S'].shape, dtype=bool); conv2.ravel()[idx] = True
    fn = os.path.join(CACHE, 'modes_' + tag + '.npz')
    np.savez_compressed(fn, keys=np.array(keys), modes=np.array(arrs), us=g['us'], ts=g['ts'], U=g['U'], S=g['S'], W=g['W'], ok=g['ok'],
                        conv=conv2, logscale=L_, Delta=Dv, J=Jv, N=N, npts=npts, time=time.time() - t0)
    print(f"  saved {fn} ({time.time() - t0:.0f} s)", flush=True)
    return fn

if __name__ == "__main__":
    J = int(sys.argv[1]); chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    maxpts = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != 'all' else None
    only = [tuple(c) for c in sys.argv[4].split(',')] if len(sys.argv) > 4 and sys.argv[4] != 'all' else None
    grid = sys.argv[5] if len(sys.argv) > 5 else 'full'
    worker, nworkers = (int(sys.argv[6]), int(sys.argv[7])) if len(sys.argv) > 7 else (0, 1)
    compute(J, chunk, maxpts, only, grid, worker, nworkers)
