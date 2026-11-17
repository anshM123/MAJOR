"""
Estimate of the moment contribution of the grid points dropped by the radial-series convergence mask (tiny t, i.e. r -> 1), for the
cached <JJJJ> mode grids of fourpt/tmoments4: per u-row, power-law extrapolation in t of the moment density from the two smallest
accepted t-nodes, summed over the dropped nodes of that row (rows without two accepted nodes are skipped and reported).
Returns (moment, power-law estimate of the dropped part, monotonic upper bound |d(first accepted)| x weights, log scale, #skipped)
in the units of TensorMoments4.moment (the dropped region is where the density decreases steeply towards t -> 0).
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from moments4 import kernel_np

def density(mom, k, nu, G):
    g = mom.g; U, S = g['U'], g['S']
    xi = -np.log(g['us'])
    rho, rhob = S**2, (S*U)**2
    with np.errstate(divide='ignore', invalid='ignore'):
        logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k - 4)*np.log(rho*rhob) + mom.logscale
    return logw, G

def moment_and_corner(mom, k, nu, ab, q, cL=None):
    g = mom.g
    q1, q3, q2, q4 = q
    F = kernel_np(nu, [(q1 + q3, q2 + q4)], -np.log(g['us']))[(q1 + q3, q2 + q4)]
    G = mom.modes[ab][q] if cL is None else sum(w*mom.modes[k_][q] for k_, w in cL.items())
    logw, _ = density(mom, k, nu, G)
    logw = np.where(mom.conv, logw, -np.inf)
    sh = np.max(logw)
    d = np.where(mom.conv, np.exp(logw - sh), 0.0)*F[:, None]*G          # density without quadrature weights
    val = np.sum(d*g['W'])
    est = 0j; bound = 0.0; skipped = 0
    ts = g['ts']; us = g['us']
    rowsum = np.zeros(d.shape[0], dtype=complex); rowok = np.zeros(d.shape[0], dtype=bool); skiprows = []
    for iu in range(d.shape[0]):
        rowsum[iu] = np.sum(d[iu]*g['W'][iu])
        acc = np.where(mom.conv[iu])[0]
        drop = np.where(g['ok'][iu] & ~mom.conv[iu])[0]
        if len(drop) == 0: continue
        if len(acc) < 2:
            skipped += len(drop); skiprows.append(iu); continue
        i1, i2 = acc[0], acc[1]
        d1, d2 = d[iu, i1], d[iu, i2]
        if d1 == 0 or d2 == 0:
            skipped += len(drop); continue
        p = np.log(abs(d2)/abs(d1))/np.log(ts[i2]/ts[i1])
        for it in drop:
            if it < i1:
                # bound: density no larger than at the first accepted node (monotonic decrease towards t -> 0, as observed);
                # estimate: power law with exponent p if p > 0, else the bound
                bound += abs(d1)*g['W'][iu, it]
                e_ = d1*(ts[it]/ts[i1])**max(p, 0.0)*g['W'][iu, it]
                est += e_; rowsum[iu] += e_
        rowok[iu] = True
    # rows with fewer than two accepted nodes (smallest u): power law in u of the (completed) row integrals of the two nearest complete
    # rows (u-weights included via the row sums; rows are Gauss nodes of the same panel); counted fully in the bound
    if skiprows:
        good = [i for i in range(d.shape[0]) if rowok[i] and i > max(skiprows)][:2]
        if len(good) == 2:
            r1, r2 = rowsum[good[0]], rowsum[good[1]]
            # per unit u-weight: rowsum includes W = uw*tw/scale; divide by the u-quadrature weight ratio implicitly via u-power fit
            pu = np.log(abs(r2)/abs(r1))/np.log(us[good[1]]/us[good[0]]) if r1 != 0 and r2 != 0 else 0.0
            for iu in skiprows:
                uw_ratio = g['W'][iu].max()/g['W'][good[0]].max() if g['W'][good[0]].max() > 0 else 1.0
                e_ = r1*(us[iu]/us[good[0]])**max(pu, 0.0)*uw_ratio
                est += e_ - rowsum[iu]*0; bound += abs(r1)*uw_ratio
    return val, est, bound, sh, skipped
