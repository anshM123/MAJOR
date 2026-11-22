"""
rs2_uniform (HRS_CLOSE.md §3, check B6): uniformity of the analytic-radial 1/Delta expansion near every divisor of the Regge-scaled
integrand.  The contiguity argument (HRS_CLOSE §2) needs every factor that the engine logs, raises to a non-integer power or inverts
to be  c * atoms^e * (1 + eps)  with eps = sum_k a^k eps_k BOUNDED as any divisor atom -> 0 (no boundary layer: no zero of 1 + eps
approaching a divisor as a -> 0).  Exponent tables cannot decide this (the atoms u, v, ca, sa and cf, sf are redundant), so it is
tested numerically: for each factor, each order k <= K and each divisor (p1, p2, q1, q2 -> 0; w -> 0 [s -> 1]; cf -> 0 [phi -> pi/2];
u -> 0 [alpha -> 0]; v -> 0 [alpha -> pi]) the growth exponent  g = d log max|eps_k| / d log t  at t = 1e-2 .. 1e-6 is fitted
(g >= 0: uniform; g < 0: boundary layer).
Factors: all log factors of the scalar (l1_crystal via bg_geom.scalar_pieces) and of the EE k = 2/4 integrand (bg_ee.ee_pieces:
3-pt powers, V_{0,42}), every inverted distance of the EE polarisation prefactor (inv_sq), u'v' (inverted in the scalar f and in
Lemma-B R), m13 m24 and -x^2.
Usage: python3 rs2_uniform.py K   -> cache/rs2_uniform_K{K}.txt
"""
import os, sys, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_ser as S, bg_geom as GG, bg_ee as BE

C = os.path.join(HERE, 'cache')
TS = 10.0**-np.arange(2, 7)


def approach(div, t, n=6, seed=21):
    P = S.atom_points(n, seed).astype(float)
    if div in ('p1', 'p2', 'q1', 'q2'):
        P[:, S.NAMES.index(div)] = t
    elif div == 'w':
        s = np.sqrt(1 - t); P[:, 4] = s; P[:, 5] = t
    elif div == 'cf':
        P[:, 6] = t; P[:, 7] = np.sqrt(1 - t*t)
    elif div in ('u', 'v'):
        ca = 1 - t if div == 'u' else -1 + t
        P[:, 8] = ca; P[:, 9] = np.sqrt(1 - ca*ca); P[:, 10] = 1 - ca; P[:, 11] = 1 + ca
    return P


DIVS = ('p1', 'p2', 'q1', 'q2', 'w', 'cf', 'u', 'v')


def growth(eps_k):
    out = {}
    for d in DIVS:
        vals = np.array([np.max(np.abs(eps_k.ev(approach(d, t)))) for t in TS])
        if np.all(vals < 1e-300):
            out[d] = np.inf
            continue
        g = np.polyfit(np.log(TS), np.log(np.maximum(vals, 1e-300)), 1)[0]
        out[d] = g
    return out


def main(K):
    t0 = time.time()
    lines = [f'rs2_uniform K={K}: growth exponent g of max|eps_k| as the divisor atom t -> 0 (g >= 0 uniform; g < 0 boundary layer)']
    logs, polys, G, nm = GG.scalar_pieces(K)
    P = BE.ee_pieces(K + 1, 300.0, 0)
    facs = []
    for (n, f, p0, p1) in logs:
        facs.append((f'scalar log {n} (p0={p0:+.2f}, p1={p1:+.1f})', f))
    for (n, f, p0, p1) in P['logs']:
        if n in ('3pt', 'v042'):
            facs.append((f'EE log {n} (p0={p0:+.2f})', f))
    Gm = BE.EEGeom(K + 1, P['G'])
    for (i, j) in ((1, 2), (1, 5), (2, 5), (3, 4), (3, 5), (4, 5), (2, 4), (2, 0), (4, 0), (1, 3), (1, 4), (2, 3)):
        facs.append((f'inverted x{i}{j}^2', Gm.sq(i, j)))
    facs.append(("inverted u'v'", nm['up']*nm['vp']))
    facs.append(('inverted m13 m24', nm['m13']*nm['m24']))
    facs.append(('inverted -x^2', -G['xsq']))
    worst = np.inf
    for name, f in facs:
        c0, e0, eps = f.trunc(K).rel(name)
        row = []
        for k, p in enumerate(eps.c):
            if k == 0 or len(p) == 0:
                continue
            g = growth(p)
            gmin = min(g.values()); worst = min(worst, gmin)
            row.append(f'k={k}: ' + ' '.join(f'{d}:{v:+.2f}' if np.isfinite(v) else f'{d}:0' for d, v in g.items()))
        lines.append(f'{name}: lead ' + str({S.NAMES[i]: int(e0[i]) for i in range(12) if e0[i]}))
        lines += ['    ' + r for r in row] if row else ['    (exact monomial: eps = 0)']
    lines.append(f'WORST growth exponent over all factors, orders, divisors: {worst:+.3f}  ({time.time() - t0:.0f}s)')
    txt = '\n'.join(lines)
    print(txt, flush=True)
    open(os.path.join(C, f'rs2_uniform_K{K}.txt'), 'w').write(txt + '\n')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
