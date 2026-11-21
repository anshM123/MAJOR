"""Result Q Stage 0: LP for the upper bound on g4 (and the |g3|^2-dependence) from the minimal flat sum-rule families.  [NE]

Family (common_q.py):  F = int psi2 B2 + int psi3 B3 - a (dB4|_0 - B5|_0) - c B4|_0, with c = 1/2 fixed (normalisation), so
    g4 M^6/(8 pi G) <= A log(M/m_IR) + B - kappa |g3|^2 M^8,   A = 2 x0, B = 2 Lfin, kappa = k/8.
Switches: --no-b3 (psi3 = 0), --no-fwd (a = 0), --psi3zero (psi3(0) = 0), --matter (light-matter robust), delta (margin, relative
to N = |B2| + |B3| + |FW| + |CB4|), --G3 g (objective A log(LREF) + B - kappa g at fixed |g3|^2 M^8 = g), --lref (log reference).
Objective default: minimise A (tie-break 1e-3 B); with --lref L: minimise A log(L) + B (- kappa G3).
Usage: python3 lp_g4.py d2 d3 r [options]    -> prints JSON line; appends to cache/lp_g4.jsonl
"""
import sys, os, json, time, argparse
import numpy as np
from scipy.optimize import linprog
from common_q import Grid, CACHE, SIGMA


def solve(G, delta=0.0, use_b3=True, use_fwd=True, matter=False, lref=None, G3=0.0, maxit=60, vb=1e5, verbose=False,
          kappa_min=None):
    n2, n3 = G.n2, G.n3
    nv = n2 + n3 + 2                                   # x, y, a, c
    cfix = 0.5
    # objective
    obj = np.zeros(nv)
    if lref is None:
        obj[:n2] = 2 * G.x0 + 2e-3 * G.Lfin
    else:
        obj[:n2] = 2 * G.x0 * np.log(lref) + 2 * G.Lfin
    # -kappa G3 = -(k/8) G3, k = -2 mu2_6 + 2 mu3_4 + 4a
    obj[:n2] += -G3 / 8 * (-2 * G.mu2_6); obj[n2:n2 + n3] += -G3 / 8 * (2 * G.mu3_4); obj[n2 + n3] += -G3 / 8 * 4
    bounds = [(-vb, vb)] * n2 + [((-vb, vb) if use_b3 else (0, 0))] * n3 + [((0, vb) if use_fwd else (0, 0)), (cfix, cfix)]
    fixed, bfix = [], []
    # scaling limit: 2 int psi2 J0 >= 0
    R = np.zeros((len(G.bsc), nv)); R[:, :n2] = -G.Bsc; fixed.append(R); bfix.append(np.zeros(len(G.bsc)))
    # asymptotic rows (large bhat, fixed u)
    for h in SIGMA:
        A2, A3, Fc, Cb = G.asym[h]
        R = np.zeros((len(G.ab), nv)); R[:, :n2] = -A2; R[:, n2:n2 + n3] = -A3; R[:, n2 + n3] = -Fc; R[:, -1] = Cb
        fixed.append(R); bfix.append(np.zeros(len(G.ab)))
    if matter:
        M = G.matter_rows(); fixed.append(M); bfix.append(np.zeros(len(M)))
    if kappa_min is not None:                          # require kappa >= kappa_min  (i.e. -k/8 <= -kappa_min)
        r = np.zeros(nv); r[:n2] = 2 * G.mu2_6 / 8; r[n2:n2 + n3] = -2 * G.mu3_4 / 8; r[n2 + n3] = -4 / 8
        fixed.append(r[None]); bfix.append(np.array([-kappa_min]))
    Afix = np.concatenate(fixed); bfix = np.concatenate(bfix)
    sc = np.maximum(np.abs(Afix).max(axis=1), 1e-300); Afix = Afix / sc[:, None]; bfix = bfix / sc

    def rows(idx, s2, s3):
        sF = np.sign(G.FWc[idx])
        Rr = np.zeros((len(idx), nv))
        Rr[:, :n2] = (1 - delta * s2)[:, None] * G.B2[idx]
        Rr[:, n2:n2 + n3] = (1 - delta * s3)[:, None] * G.B3[idx]
        Rr[:, n2 + n3] = (1 - delta * sF) * G.FWc[idx]
        Rr[:, -1] = -(1 + delta) * G.CB4[idx]
        return -Rr

    init = np.unique(np.concatenate([np.arange(0, len(G.u), max(1, len(G.u) // 3000)), np.where(G.J <= 12)[0]]))
    cuts = []
    for s2 in (1, -1):
        for s3 in (1, -1):
            Rr = rows(init, np.full(len(init), s2), np.full(len(init), s3))
            cuts.append(Rr / np.maximum(np.abs(Rr).max(axis=1), 1e-300)[:, None])
    for it in range(maxit):
        Aub = np.concatenate([Afix] + cuts); bub = np.concatenate([bfix, np.zeros(sum(len(x) for x in cuts))])
        res = linprog(obj, A_ub=Aub, b_ub=bub, bounds=bounds, method='highs')
        if res.status != 0:
            return {'status': res.status, 'msg': res.message, 'it': it}
        v = res.x
        x, y, a, c = v[:n2], v[n2:n2 + n3], v[n2 + n3], v[-1]
        B2 = G.B2 @ x; B3 = G.B3 @ y; FW = a * G.FWc; CB = c * G.CB4
        FF = B2 + B3 + FW - CB; N = np.abs(B2) + np.abs(B3) + np.abs(FW) + CB
        g = FF - delta * N
        viol = np.where(g < -1e-7 * np.maximum(N, 1e-300))[0]
        if verbose:
            print('  it %d obj %.5f viol %d min g/N %.3e' % (it, res.fun, len(viol), np.min(g / np.maximum(N, 1e-300))), flush=True)
        if len(viol) == 0:
            break
        order = viol[np.argsort(g[viol] / N[viol])][:4000]
        Rr = rows(order, np.sign(B2[order]), np.sign(B3[order]))
        cuts.append(Rr / N[order][:, None])
    else:
        return {'status': 9, 'msg': 'no convergence'}
    k = -2 * G.mu2_6 @ x + 2 * G.mu3_4 @ y + 4 * a
    A = 2 * G.x0 @ x; B = 2 * G.Lfin @ x
    # minimum of FF/N on the grid, and where
    rat = FF / np.maximum(N, 1e-300); i = int(np.argmin(rat))
    psi2 = sum(xi * cc for xi, cc in zip(x, [np.pad(cc, (0, 40 - len(cc))) for cc in G.b2]))
    psi3 = sum(yi * cc for yi, cc in zip(y, [np.pad(cc, (0, 40 - len(cc))) for cc in G.b3])) if n3 else np.zeros(40)
    return {'status': 0, 'A': A, 'B': B, 'kappa': k / 8, 'a': a, 'c': c, 'it': it,
            'psi2pp0': float(G.pp0 @ x), 'psi3_0': float(psi3[0]),
            'minFF_N': float(rat[i]), 'at': [float(1 / np.sqrt(G.u[i])), int(G.J[i]), int(G.h[i])],
            'x': x.tolist(), 'y': y.tolist(), 'psi2': np.trim_zeros(psi2, 'b').tolist(), 'psi3': np.trim_zeros(psi3, 'b').tolist()}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('d2', type=int); ap.add_argument('d3', type=int); ap.add_argument('r', type=int)
    ap.add_argument('--delta', type=float, default=0.0)
    ap.add_argument('--no-b3', action='store_true'); ap.add_argument('--no-fwd', action='store_true')
    ap.add_argument('--psi3zero', action='store_true'); ap.add_argument('--matter', action='store_true')
    ap.add_argument('--lref', type=float, default=None); ap.add_argument('--G3', type=float, default=0.0)
    ap.add_argument('--kmin', type=float, default=None)
    ap.add_argument('-v', action='store_true')
    o = ap.parse_args()
    t0 = time.time()
    G = Grid(o.d2, o.d3, o.r, j0=1 if o.psi3zero else 0)
    out = solve(G, o.delta, not o.no_b3, not o.no_fwd, o.matter, o.lref, o.G3, verbose=o.v, kappa_min=o.kmin)
    out.update(dict(d2=o.d2, d3=o.d3, r=o.r, delta=o.delta, b3=not o.no_b3, fwd=not o.no_fwd, psi3zero=o.psi3zero, matter=o.matter,
                    lref=o.lref, G3=o.G3, kmin=o.kmin, secs=round(time.time() - t0, 1)))
    short = {k: (round(v, 5) if isinstance(v, float) else v) for k, v in out.items() if k not in ('x', 'y', 'psi2', 'psi3')}
    print(json.dumps(short))
    with open(os.path.join(CACHE, 'lp_g4.jsonl'), 'a') as f:
        f.write(json.dumps(out) + '\n')
