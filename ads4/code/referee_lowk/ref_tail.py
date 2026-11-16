"""
Referee test of the Darboux-tail corner completion (ads4/code/lowk/corner.py): rebuild one l = 1 case from scratch (settings A, B
as in the note -> reproduction of the table entries), and re-evaluate the same corner points with further tail settings:
changed fit window, changed number of asymptotic terms, shorter exact sum, log terms in the alternating family, and a
deliberately wrong (mutated) Casimir exponent.  For every setting the *same* acceptance mask (from A/B) is used and the
k = 2, 3 helicity ratios are recomputed with the note's error machinery.
Usage: python3 ref_tail.py Delta J [l]      -> cache/reftail_l{l}_D{Delta}.pkl and stdout
"""
import os, sys, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
LOWK = os.path.join(HERE, '..', 'lowk'); E2E = os.path.join(HERE, '..', 'e2e')
sys.path.insert(0, E2E); sys.path.insert(0, LOWK)
import numpy as np, mpmath as mp
import corner, lcblock
from corner import TensorMomentsLC, refined_base, moments_with_err, _Pre
from moments import frame_data, frame_positions, pol_grid, fourier_modes
from spinblock import SpinBlock, invariants
from lcblock import LCBlock, poch_over_fact
import t1_lowk_e2e as T1

class LCLog(LCBlock):
    """alternating family with q^e and q^e log q terms for BOTH exponents (resonant case c1 - c2 integer -> logs at zbar = oo)."""
    def _basis(self, q, l, t):
        cols, kinds = [], []
        for m in range(self.M):
            cols.append(poch_over_fact(2*self.sig + l + t - m, q)); kinds.append(('main', 2*self.sig + l + t - m))
        alt = (-1.0)**q
        for c in sorted(set([self.c1, self.c2])):
            for m in range(self.Malt):
                a = -2*c + l + t - m
                cols.append(alt*poch_over_fact(a, q)); kinds.append(('alt', a))
                e = a - 1
                cols.append(alt*q**e*np.log(q)); kinds.append(('altlog', e))
        return np.stack(cols, 1), kinds

class LCMut(LCBlock):
    """MUTATION: main Darboux exponent 2 sigma + l + t - m replaced by 2 sigma + l + t - m + 0.25 (wrong Casimir exponent)."""
    def _basis(self, q, l, t):
        A, kinds = LCBlock._basis(self, q, l, t)
        cols = []; nk = []
        for j, kd in enumerate(kinds):
            if kd[0] == 'main':
                a = kd[1] + 0.25
                cols.append(poch_over_fact(a, q)); nk.append(('main', a))
            else:
                cols.append(A[:, j]); nk.append(kd)
        return np.stack(cols, 1), nk

A = dict(N=10000, I=1500, Q=8000, M=8)
SETTINGS = [
    ('A', LCBlock, A),
    ('B', LCBlock, dict(N=8000, I=1200, Q=6000, M=6)),
    ('win_hi', LCBlock, dict(A, qlo=6000)),            # fit window [6000, 8000]
    ('win_lo', LCBlock, dict(A, qlo=1000)),            # fit window [1000, 8000]
    ('M4', LCBlock, dict(A, M=4, Malt=2)),             # fewer asymptotic terms
    ('M12', LCBlock, dict(A, M=12, Malt=5)),           # more asymptotic terms
    ('Q3000', LCBlock, dict(A, Q=3000)),               # exact sum only to 3000, tail from there
    ('altlog', LCLog, A),                              # log terms in the alternating family
    ('MUT', LCMut, A),                                 # mutation: wrong main exponent (+0.25)
]

STASH = {}

def evaluate_all(l, Delta, J, Dphi, U, S):
    fd = frame_data(U, S)
    z1, z3 = pol_grid(len(U))
    x1, x2, x3, x4 = frame_positions(fd['ru'], fd['xi'])
    s, z = invariants(x1, x2, x3, x4, z1, z3)
    sb = SpinBlock(l, Delta, J, Dphi, J + 60)
    keys = {ab: (tuple(sorted((blk.a, blk.b))), vpow) for ab, (blk, vpow) in sb.blocks.items()}
    print("seed (c1,c2) pairs:", sorted(set(k for k, v in keys.values())), flush=True)
    Ssets = {name: {} for name, _, _ in SETTINGS}; itails = {name: [] for name, _, _ in SETTINGS}; fitres = {name: {} for name, _, _ in SETTINGS}
    for key in sorted(set(k for k, v in keys.values())):
        t0 = time.time()
        big = LCBlock(Delta, J, key[0], key[1], **A)
        for name, cls, par in SETTINGS:
            lb = big if (name == 'A') else cls(Delta, J, key[0], key[1], coeffs=big.C, **par)
            Ssets[name][key] = lb.series(fd['rs'], fd['chi'])
            itails[name].append(lb.i_tail)
            fitres[name][key] = float(lb.fitres[:301].max())
        print(f"  seed {key}: {time.time()-t0:.0f} s; fit residuals (i<=300): " + ", ".join(f"{n}={fitres[n][key]:.1e}" for n, _, _ in SETTINGS), flush=True)
        del big
    out = {}
    for name, _, _ in SETTINGS:
        for ab, (key, vpow) in keys.items():
            sb.blocks[ab] = (_Pre(Ssets[name][key], *key), vpow)
        ls, M = sb.eval(s, z, fd['rs'], fd['chi'], fd['A0'], fd['E0'])
        out[name] = (ls, fourier_modes(M, lmax=l), np.max(itails[name], axis=0), fitres[name])
    STASH['out'] = out
    return [out['A'], out['B']]

if __name__ == "__main__":
    Delta, J = float(sys.argv[1]), int(sys.argv[2]); l = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    corner.evaluate_points = lambda l_, D_, J_, Dphi_, U_, S_, pars=None: evaluate_all(l_, D_, J_, Dphi_, U_, S_)
    t0 = time.time()
    base = refined_base(l, Delta, J)
    fn = os.path.join(HERE, 'cache', f'rebuild_l{l}_D{Delta:.6f}_J{J}.npz')
    if os.path.exists(fn): os.remove(fn)
    mom = TensorMomentsLC.build(l, Delta, J, base=base, fn=fn, cache=True)
    print(f"build done {time.time()-t0:.0f} s; info {mom.info}", flush=True)
    # re-derive the point selections exactly as in build
    g = base.g; todo = g['ok'] & ~base.conv
    out = STASH['out']
    keys = sorted(mom.modes.keys())
    # positions of the todo points among the evaluated ones (pts_all = todo | vmask): recover vmask from build's rng
    cand = np.argwhere(g['ok'] & base.conv & (g['U'] < 0.25)); rng = np.random.default_rng(1)
    sel = cand[rng.choice(len(cand), size=min(300, len(cand)), replace=False)]
    vmask = np.zeros_like(todo); vmask[tuple(sel.T)] = True
    pts_all = todo | vmask; tsel = todo[pts_all]
    good = mom.lcpts[todo]
    gi = tuple(np.argwhere(todo)[good].T)
    lsA = out['A'][0][tsel]
    results = {}
    for name, _, _ in SETTINGS:
        ls, fm, it, fr = out[name]
        m2 = TensorMomentsLC.__new__(TensorMomentsLC)
        m2.__dict__.update({k: v for k, v in mom.__dict__.items() if k not in ('modes',)})
        m2.modes = {k: mom.modes[k].copy() for k in keys}
        for k in keys:
            V = fm[k][tsel]*np.exp(ls[tsel] - lsA)
            m2.modes[k][gi] = V[good]
        res, errs = moments_with_err(m2, (2, 3), (0.0, 1.0, 3.0))
        for k in (2, 3):
            for nu in (0.0, 1.0, 3.0):
                lc1, dr1, sh1, es1 = errs[(k, nu, l, -l)]; lc2, dr2, sh2, es2 = errs[(k, nu, l, l)]
                e1 = lc1 + dr1 + T1.ERR_RAD; e2 = lc2 + dr2 + T1.ERR_RAD
                Rp, Rf, dRp, dRf, worst = T1.components_with_error(l, Delta, J, k, nu, res, e1, e2)
                results[(name, k, nu)] = (Rp, Rf, dRp, dRf, lc1, lc2, dr1, dr2)
                print(f"{name:7s} k={k} nu={nu}: R_pres {Rp:.9f} +- {dRp:.1e}  R_flip {Rf:.9f} +- {dRf:.1e}   errlc {lc1:.1e}/{lc2:.1e} drop {dr1:.1e}/{dr2:.1e}", flush=True)
    pickle.dump(dict(results=results, info=mom.info, fitres={n: out[n][3] for n in out}), open(os.path.join(HERE, 'cache', f'reftail_l{l}_D{Delta:.6f}.pkl'), 'wb'))
    print(f"total {time.time()-t0:.0f} s")
