"""Explicit spin-2 contact diagrams near the collinear limit y -> 0 (test of Lemma Y's conclusion), and the
y-orders of epsilon structures.  Independent of the r4heft code (no Mellin representation is used).

Configuration (3d metric (+,-,-); embedding P=(1,x^2,x), Z=(0,2x.z,z), P.Q = -(P+Q- + P-Q+)/2 + p.q):
  x1 = -t1 z(a), x3 = t3 z(b), x2 = T e0, x4 = -T e0, z(phi) = (1, cos phi, sin phi),
  a = phibar - dphi/2, b = phibar + dphi/2, y = 1 - cos(dphi);  polarisations z1 = z(a), z3 = z(b), z2, z4 fixed null.
  X_ij = -2 P_i.P_j = x_ij^2 > 0 for t1, t3 < T/2;  X13 = 2 t1 t3 y.

Bulk contact of four spin-2 Delta=3 fields: K_i = (V_i.W)^2 / W_i^5, W_i = -2P_i.X, U_i = Z_i.X,
V_i = W_i Z_i + 2 U_i P_i (tangent and null).  Vertices: index contractions
   cyc1234 = (V1.V2)(V2.V3)(V3.V4)(V4.V1),  cyc1324 = (V1.V3)(V3.V2)(V2.V4)(V4.V1),
   pair12  = (V1.V2)^2 (V3.V4)^2,  pair13 = (V1.V3)^2 (V2.V4)^2,  pair14 = (V1.V4)^2 (V2.V3)^2.
Reduction (Schwinger): int dX prod U_i^n_i W_i^-b_i = C * G(sum(b-n)/2 - h)/prod G(b_i)
      * int prod dt t_i^(b_i-1) prod_i (sum_{j!=i} t_j Z_i.P_j)^n_i exp(-sum t_i t_j X_ij)
so the contact = sum_k (monomial in invariants)_k * I_{Delta'_k}(X) with the Symanzik I of dfun.py.
"""
import itertools, math, os, sys, time
import numpy as np
import mpmath as mp
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from dfun import I_fast as I_sym

OUT = os.path.join(os.path.dirname(__file__), 'cache', 'B2_contact.txt')
lines = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

G3 = np.diag([1.0, -1.0, -1.0])
def dot3(u, v): return u @ G3 @ v
def zvec(phi): return np.array([1.0, math.cos(phi), math.sin(phi)])
def P_of(x): return np.concatenate([[1.0, dot3(x, x)], x])
def Z_of(x, z): return np.concatenate([[0.0, 2 * dot3(x, z)], z])
GE = np.zeros((5, 5)); GE[0, 1] = GE[1, 0] = -0.5; GE[2:, 2:] = G3
def dotE(A, B): return A @ GE @ B

def config(dphi, t1=0.3, t3=0.2, T=1.0, phibar=0.4, a2=1.1, a4=-2.3):
    a, b = phibar - dphi / 2, phibar + dphi / 2
    e0 = np.array([1.0, 0, 0])
    xs = {1: -t1 * zvec(a), 3: t3 * zvec(b), 2: T * e0, 4: -T * e0}
    zs = {1: zvec(a), 3: zvec(b), 2: zvec(a2), 4: zvec(a4)}
    P = {i: P_of(xs[i]) for i in xs}
    Z = {i: Z_of(xs[i], zs[i]) for i in xs}
    return P, Z

def invariants(P, Z):
    ZZ = {(i, j): dotE(Z[i], Z[j]) for i in P for j in P}
    ZP = {(i, j): dotE(Z[i], P[j]) for i in P for j in P}
    PP = {(i, j): dotE(P[i], P[j]) for i in P for j in P}
    X = {f'{i}{j}': -2 * PP[(i, j)] for i in P for j in P if i < j}
    return ZZ, ZP, PP, X

# ---- polynomial algebra: dict exponent-tuple -> coeff ; variables W1..W4, U1..U4 (8 exps)
def pmul(A, B):
    C = defaultdict(float)
    for ka, va in A.items():
        for kb, vb in B.items():
            C[tuple(x + y for x, y in zip(ka, kb))] += va * vb
    return C

def VV(i, j, ZZ, ZP, PP):
    # V_i.V_j = W_i W_j Zi.Zj + 2 W_i U_j Zi.Pj + 2 U_i W_j Pi.Zj + 4 U_i U_j Pi.Pj
    def e(Wi=(), Ui=()):
        v = [0] * 8
        for k in Wi: v[k - 1] += 1
        for k in Ui: v[4 + k - 1] += 1
        return tuple(v)
    return {e(Wi=(i, j)): ZZ[(i, j)], e(Wi=(i,), Ui=(j,)): 2 * ZP[(i, j)],
            e(Wi=(j,), Ui=(i,)): 2 * ZP[(j, i)], e(Ui=(i, j)): 4 * PP[(i, j)]}

VERTS = {'cyc1234': [(1, 2), (2, 3), (3, 4), (4, 1)], 'cyc1324': [(1, 3), (3, 2), (2, 4), (4, 1)],
         'pair12': [(1, 2), (1, 2), (3, 4), (3, 4)], 'pair13': [(1, 3), (1, 3), (2, 4), (2, 4)],
         'pair14': [(1, 4), (1, 4), (2, 3), (2, 3)]}

def reduce_contact(vert, ZZ, ZP, PP):
    """returns dict Delta' (4-tuple) -> coefficient (float); common constant C*G(6-3/2) dropped."""
    poly = {(0,) * 8: 1.0}
    for (i, j) in VERTS[vert]:
        poly = pmul(poly, VV(i, j, ZZ, ZP, PP))
    out = defaultdict(float)
    for ex, cf in poly.items():
        if abs(cf) == 0: continue
        p, n = ex[:4], ex[4:]
        assert all(p[k] + n[k] == 2 for k in range(4))
        pref = cf / math.prod(math.gamma(3 + n[k]) for k in range(4))
        # expand prod_i (sum_{j != i} t_j ZP(i,j))^{n_i}
        tpoly = {(0, 0, 0, 0): pref}
        for i in range(1, 5):
            for _ in range(n[i - 1]):
                new = defaultdict(float)
                for kt, vt in tpoly.items():
                    for j in range(1, 5):
                        if j == i: continue
                        kk = list(kt); kk[j - 1] += 1
                        new[tuple(kk)] += vt * ZP[(i, j)]
                tpoly = new
        for kt, vt in tpoly.items():
            D = tuple(3 + n[k] + kt[k] for k in range(4))
            out[D] += vt
    return out

def sigma(D): return (D[0] + D[2] - D[1] - D[3]) / 2

def run_contacts(dphis, verts):
    Icache = {}
    res = {v: [] for v in verts}
    for dphi in dphis:
        P, Z = config(dphi)
        ZZ, ZP, PP, X = invariants(P, Z)
        y = 1 - math.cos(dphi)
        for v in verts:
            terms = reduce_contact(v, ZZ, ZP, PP)
            tot = mp.mpf(0); maxterm = 0; worst = None
            for D, cf in terms.items():
                if cf == 0: continue
                key = (D, dphi)
                if key not in Icache:
                    Icache[key] = I_sym(D, X)
                val = cf * Icache[key]
                tot += val
                if abs(val) > maxterm: maxterm = abs(val); worst = D
            res[v].append((dphi, y, float(tot), float(maxterm), worst, len(terms)))
    return res, Icache

def main():
    t0 = time.time()
    mp.mp.dps = 30
    # ---- sanity of the configuration
    P, Z = config(0.3)
    ZZ, ZP, PP, X = invariants(P, Z)
    y = 1 - math.cos(0.3)
    say('configuration check at dphi = 0.3: P_i^2, Z_i^2, Z_i.P_i =',
        max(abs(dotE(P[i], P[i])) for i in P), max(abs(dotE(Z[i], Z[i])) for i in P), max(abs(dotE(Z[i], P[i])) for i in P))
    say('  X13/(2 t1 t3 y) =', X['13'] / (2 * 0.3 * 0.2 * y), '; all X_ij > 0:', all(v > 0 for v in X.values()))
    say('  Z1.Z3/y, Z1.P3/(t3 y), Z3.P1/(t1 y), P1.P3/(t1 t3 y):', ZZ[(1, 3)] / y, ZP[(1, 3)] / (0.2 * y), ZP[(3, 1)] / (0.3 * y), PP[(1, 3)] / (0.06 * y))
    # ---- y-orders of the invariants and of epsilon structures
    dps = [0.2 / 2**k for k in range(0, 8)]
    def order(fun):
        vals = [abs(fun(d)) for d in dps]
        ys = [1 - math.cos(d) for d in dps]
        return [math.log(vals[k + 1] / vals[k]) / math.log(ys[k + 1] / ys[k]) if vals[k] > 0 and vals[k+1] > 0 else float('nan')
                for k in range(len(dps) - 1)]
    say('\n[eps] local y-orders of eps(V1..V5) (det of the 5 embedding vectors), dphi = 0.2 -> 0.2/128:')
    vecnames = [('P', i) for i in range(1, 5)] + [('Z', i) for i in range(1, 5)]
    stats = defaultdict(list)
    for combo in itertools.combinations(vecnames, 5):
        k13 = sum(1 for t, i in combo if i in (1, 3))
        def f(d, combo=combo):
            P, Z = config(d)
            M = np.array([(P if t == 'P' else Z)[i] for t, i in combo])
            return np.linalg.det(M)
        o = order(f)
        stats[k13].append(o[-1])
    for k in sorted(stats):
        arr = np.array(stats[k])
        fin = arr[np.isfinite(arr)]
        say(f'   k13 = {k}: {len(arr)} placements; final local order min {fin.min() if len(fin) else float("nan"):.4f} '
            f'max {fin.max() if len(fin) else float("nan"):.4f}  (non-finite/zero: {len(arr) - len(fin)})')
    # k13 = 4 : absolute size
    P, Z = config(0.2)
    m4 = max(abs(np.linalg.det(np.array([(P if t == 'P' else Z)[i] for t, i in c])))
             for c in itertools.combinations(vecnames, 5) if sum(1 for t, i in c if i in (1, 3)) == 4)
    say(f'   k13 = 4: max |eps| at dphi = 0.2: {m4:.2e} (identically zero up to roundoff)')

    # ---- explicit contacts
    verts = list(VERTS)
    dphis = [0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625]
    say('\n[contacts] reduction to sum monomial x I_Delta\'(X); value of the full contact and of its largest single term')
    P, Z = config(0.1); ZZ, ZP, PP, X = invariants(P, Z)
    for v in verts:
        terms = reduce_contact(v, ZZ, ZP, PP)
        smax = max(sigma(D) for D in terms)
        say(f'   {v}: {len(terms)} distinct Delta\' ; max sigma = {smax}; Delta\' with max sigma: '
            f'{sorted([D for D in terms if sigma(D) == smax])[:4]}')
    res, Icache = run_contacts(dphis, verts)
    for v in verts:
        say(f'\n   vertex {v}:')
        say('     dphi        y           G(y)                 slope dlog|G|/dlog y   G/log(y)    max|term|      ')
        rows = res[v]
        for k, (d, y, g, mt, w, nt) in enumerate(rows):
            sl = '' if k == 0 else f'{math.log(abs(g / rows[k-1][2])) / math.log(y / rows[k-1][1]):+.4f}'
            say(f'     {d:<10.6g} {y:<11.4e} {g:+.12e}  {sl:>10}   {g / math.log(y):+.6e}  {mt:.3e} (D\'={w})')
        # fit G = A + B log y on the last 4 points
        ys = np.array([r[1] for r in rows[-4:]]); gs = np.array([r[2] for r in rows[-4:]])
        Mx = np.vstack([np.ones_like(ys), np.log(ys)]).T
        coef, *_ = np.linalg.lstsq(Mx, gs, rcond=None)
        resid = np.max(np.abs(Mx @ coef - gs)) / np.max(np.abs(gs))
        Mx2 = np.vstack([np.ones_like(ys), np.log(ys), ys * np.log(ys), ys]).T
        coef2, *_ = np.linalg.lstsq(Mx2, gs, rcond=None)
        say(f'     fit last 4 pts G = A + B log y: A = {coef[0]:+.6e}, B = {coef[1]:+.6e}, rel. residual {resid:.1e}; '
            f'with y log y, y terms: A = {coef2[0]:+.6e}, B = {coef2[1]:+.6e}')
    # ---- the D-functions alone blow up: show I_Delta' for the max-sigma Delta' of pair13
    say('\n[control] single I_Delta\'(X) with sigma > 0 alone (no monomial): local y-order should be -sigma')
    for D in [(7, 3, 7, 3), (5, 3, 5, 3), (4, 3, 4, 3), (6, 4, 5, 3)]:
        vals = []
        for d in dphis[2:]:
            P, Z = config(d); _, _, _, X = invariants(P, Z)
            vals.append((1 - math.cos(d), float(I_sym(D, X))))
        sl = [math.log(vals[k+1][1] / vals[k][1]) / math.log(vals[k+1][0] / vals[k][0]) for k in range(len(vals) - 1)]
        say(f'   D\'={D} sigma={sigma(D)}: local orders ' + ' '.join(f'{s:+.4f}' for s in sl))
    say(f'\n(elapsed {time.time() - t0:.0f} s)')
    open(OUT, 'w').write('\n'.join(lines) + '\n')

if __name__ == '__main__':
    main()
