"""
<TTTT> beyond one point (H1 of RESULT_D.md): covariant-kernel moments of the exact Taylor-mode waves (tmoments_l2 caches) against the
factorised prediction (t13: T^{hLhR} = 8 sqrt(lambda_hL lambda_hR)/m^8, a = k + 5), with a ROBUST extraction: a fixed kernel contraction
of the modes (no fit), as recommended by REFEREE_SPINNING_ODDK.md M1.

Kernels at each Regge pair (frames of SPINNING_4PT.md sect. 4.5: on the polarisation circles D_1 = cos(th_out - th_in), D_2 = -1,
D_5 = sin(th_out - th_in)):
   E  = D_1^2               (project even kernel "D1 x D1" for l = 2)      -> modes (q_in, q_out) = (0,0): 1/2, (+-2,-+2): 1/4
   FF = D_1^2 + D_1 D_2 = D_1 (z_in.z_out)   (flip-free even kernel S_ff of SUPERCONVERGENT_UPLIFT.md)  -> E and (+-1,-+1): -1/2
   O  = D_5 D_1             (odd kernel A)                                   -> (2,-2): 1/(4i), (-2,2): -1/(4i)
(the average of K(th) e^{i q.th} over the circles picks K's coefficient at -q).  Two-pair kernels E(x)E, FF(x)FF, O(x)O.
Waves: W[X_s,X_s], W[X_o,X_o], W[X_s,X_o].  Ratio R = measured/predicted kernel moment; corner (dropped tiny-t points) estimated with
oddk/corner_est (power law; second estimate held at the first accepted node; NOT a bound).
Usage: python3 t16_kernels_ext.py [files...]   (default: all cache/modes_l2_*_d.npz)    Output: out_t16_kernels_ext.txt
"""
import os, sys, glob, itertools
import numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', 'oddk')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from t13_analysis_l2 import Modes, design, Tpred, lam, CLASSES, QS
from corner_est import moment_and_corner
from predict4 import weyl
import t5_analysis_l1 as T5
mp.mp.dps = 20
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

KER = {
    'E': {(0, 0): 0.5, (2, -2): 0.25, (-2, 2): 0.25},
    'FF': {(0, 0): 0.5, (2, -2): 0.25, (-2, 2): 0.25, (1, -1): -0.5, (-1, 1): -0.5},
    'O': {(2, -2): 1/(4j), (-2, 2): -1/(4j)},
}
KS = (2, 3, 4, 5, 6, 7); NUS = (0.0, 1.0, 3.0)
WAVES = [('same', 'same'), ('opp', 'opp'), ('same', 'opp')]

LEAD = {('same', 'same'): (2, -2, -2, 2), ('opp', 'opp'): (2, -2, 2, -2), ('same', 'opp'): (2, -2, 2, 2)}

def kernel_moment(mo, k, nu, wave, kern, A, Tp):
    """raw kernel moment: (measured, predicted) in a common scale (measured includes (-1)^J); no corner estimate (the per-row
    extrapolation of corner_est fails for modes that are small in the bulk, e.g. (0,0,0,0) of W[o,o], where it exceeds the mode 30x)"""
    K = KER[kern]
    terms = []; pred = 0j
    for (q1, q3), w13 in K.items():
        for (q2, q4), w24 in K.items():
            q = (q1, q3, q2, q4)
            v, sh = mo.moment(k, nu, wave, q)
            terms.append((w13*w24, v, sh))
            pred += w13*w24*sum(A[q][c]*Tp[c] for c in CLASSES)
    shm = max(t[2] for t in terms)
    meas = sum(w*v*np.exp(sh - shm) for w, v, sh in terms)*(-1)**mo.J
    return complex(meas), complex(pred), shm

def lead_mode(mo, k, nu, wave, A, Tp):
    q = LEAD[wave]
    v, e, bd, sh, sk = moment_and_corner(mo, k, nu, wave, q)
    p = sum(A[q][c]*Tp[c] for c in CLASSES)
    f = (-1)**mo.J*np.exp(sh - np.log(abs(p)))/(p/abs(p))
    return complex(v*f), complex((v + e)*f)

def rel(meas, pred, shm):
    """measured/predicted with the scales restored"""
    return complex(meas*np.exp(shm - np.log(abs(pred)))/(pred/abs(pred))) if pred != 0 else None

def run(fn):
    mo = Modes(fn)
    D, J = mo.Delta, mo.J
    Ap, Am, m, beta = weyl(D, J)
    res = {}
    for k in KS:
        for nu in NUS:
            A = design(D, J, k, nu)
            for wave in WAVES:
                Tp = Tpred(wave[0], wave[1], D, J, 8.0)
                for kern in KER:
                    res[(k, nu, wave, kern)] = kernel_moment(mo, k, nu, wave, kern, A, Tp)
                if wave[0] == wave[1]:
                    res[(k, nu, wave, 'lead')] = lead_mode(mo, k, nu, wave, A, Tp)
    return float(m), float(beta), D, J, res

def R_of(res, k, nu, wave, kern):
    if kern == 'lead':
        return res[(k, nu, wave, 'lead')][1].real
    meas, pred, shm = res[(k, nu, wave, kern)]
    return rel(meas, pred, shm).real

def mixed_size(res, k, nu, kern):
    """|K(W[s,o])| / sqrt(|K(W[s,s])| |K(W[o,o])|), measured and predicted (log-scale safe)"""
    ms, ps, ss_ = res[(k, nu, ('same', 'same'), kern)]; mo_, po, so_ = res[(k, nu, ('opp', 'opp'), kern)]
    mx, px, sx = res[(k, nu, ('same', 'opp'), kern)]
    lm = np.log(abs(mx) + 1e-300) + sx - 0.5*(np.log(abs(ms)) + ss_ + np.log(abs(mo_)) + so_)
    lp = np.log(abs(px) + 1e-300) - 0.5*(np.log(abs(ps)) + np.log(abs(po)))
    return np.exp(lm), np.exp(lp)

if __name__ == "__main__":
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(HERE, 'cache', 'modes_l2_*_d.npz')))
    pr("# <TTTT> kernel moments (robust: fixed kernel contraction) / factorised prediction (fourpt/tttt/t16_kernels_ext.py)")
    pr("# R = corner-corrected ratio [raw ratio] (corner shift in units of 1e-4); flip-free test: |FF(W[s,o])| relative to FF(W[s,s]) pred")
    data = []
    for fn in files:
        m, beta, D, J, res = run(fn)
        data.append((beta, m, D, J, res))
        pr(f"\n## beta={beta:.3f} m={m:.2f} (Delta,J)=({D:.4f},{J})  [{os.path.basename(fn)}]")
        pr("#  per wave: leading mode R (corner-corrected) [raw] | kernel ratios E (D1^2 x D1^2), FF (S_ff x S_ff), O (A x A), raw;  mixed wave: |K(so)|/sqrt|K(ss)K(oo)| meas/pred for E, FF")
        for k in KS:
            for nu in NUS:
                cells = []
                for wave in (('same', 'same'), ('opp', 'opp')):
                    lr, lc = res[(k, nu, wave, 'lead')]
                    cells.append(f"{wave[0][0]}{wave[1][0]}: lead {lc.real:.4f} [{lr.real:.4f}] E {R_of(res, k, nu, wave, 'E'):.4f} FF {R_of(res, k, nu, wave, 'FF'):.4f} O {R_of(res, k, nu, wave, 'O'):.4f}")
                mE, pE = mixed_size(res, k, nu, 'E'); mF, pF = mixed_size(res, k, nu, 'FF')
                cells.append(f"so: E {mE:.1e}/{pE:.1e} FF {mF:.1e}/{pF:.1e}")
                pr(f"  k={k} nu={nu}: " + " | ".join(cells))
    # convergence in m per (beta, k, nu, wave, kernel)
    pr("\n# convergence: m^2 (R - 1) per mass (corrected), and c1 of R - 1 = c1/m + c2/m^2 + c4/m^4 (three masses)")
    for beta in sorted(set(round(d[0], 1) for d in data)):
        rows = sorted([d for d in data if round(d[0], 1) == beta], key=lambda d: d[1])
        if len(rows) < 2: continue
        for k in KS:
            for nu in NUS:
                for wave in (('same', 'same'), ('opp', 'opp'), ('same', 'opp')):
                    if wave[0] != wave[1]: continue
                    for kern in ('lead', 'E', 'FF', 'O'):
                        ms = [r[1] for r in rows]; Rs = [R_of(r[4], k, nu, wave, kern) for r in rows]
                        line = f"  beta~{beta} k={k} nu={nu} {wave[0][0]}{wave[1][0]}-{kern:2s}: R = " + " ".join(f"{x:.4f}" for x in Rs) + \
                               " | m^2(R-1) = " + ", ".join(f"{mm*mm*(x-1):+.1f}" for mm, x in zip(ms, Rs))
                        if len(ms) >= 3:
                            f = T5.fits(ms, Rs)
                            line += f" | c1 = {f['A'][0]:+.3f}"
                        pr(line)
    open(os.path.join(HERE, 'out_t16_kernels_ext.txt'), 'w').write("\n".join(out) + "\n")
