"""
Completion of the tensor-mode data of ads4/code/e2e (data_modes/) in the corner u = sqrt(rhobar_u/rho_u) -> 0.

The e2e runs evaluate the stripped spinning partial wave on a (u, t) Gauss grid with the radial series (N <= 10^4 levels) and
drop the grid points where that series has not converged (moments.series_converged): the region u <~ 30 alpha/(N t) and
very small t.  Here the dropped points are evaluated with the lightcone + asymptotic-tail representation of the scalar
seed blocks (lcblock.LCBlock), fed into the *same* spinning machinery (spinblock.SpinBlock.eval: differential operators,
(P,Q) jets, generated tensor code) - nothing else is changed.

Error control (per grid point, then per moment):
  * every LC point is evaluated with two independent settings A = (N, I, Q, M) = (10000, 1500, 8000, 8) and
    B = (8000, 1200, 6000, 6) (different radial depth, fit window, number of asymptotic terms, lightcone depth);
    A is used, 3|A - B| is the error estimate (summed over points, i.e. no cancellation assumed); points with |A - B| > tol_pt * max|A| (relative to the largest mode at that
    point) are rejected: these are the deepest corner points, where the conserved spinning combination cancels by many orders
    of magnitude between the seed blocks and double precision is insufficient;
  * rejected points and points where the lightcone sum over i has not converged (tiny t) are replaced by a power-law
    extrapolation in u, column by column in t, from the two smallest accepted u-nodes; the whole replaced contribution is
    counted as error (conservative);
  * a validation set of already converged points at u < 0.25 is re-evaluated with A and compared with the cached
    radial-series values (TensorMomentsLC.valid).
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E2E = os.path.join(HERE, '..', 'e2e')
if E2E not in sys.path:
    sys.path.insert(0, E2E)
from moments import TensorMoments, frame_data, frame_positions, pol_grid, fourier_modes
from spinblock import SpinBlock, invariants
from formulas import Fnu_np
from lcblock import LCBlock
CACHE = os.path.join(HERE, 'cache')
E2E_DATA = os.path.join(E2E, 'data_modes')
PAR_A = dict(N=10000, I=1500, Q=8000, M=8)
PAR_B = dict(N=8000, I=1200, Q=6000, M=6)
SAFETY = 3.0     # err_lc = SAFETY * |A - B| (block-level tests: |A - B| tracks the deviation of A from a Q = 14000 reference within ~2x)
_COEF = {}

def e2e_file(l, Delta, J):
    return os.path.join(E2E_DATA, f"modes_l{l}_D{Delta:.6f}_J{J}.npz")

def lc_file(l, Delta, J):
    return os.path.join(CACHE, f"lcmodes_l{l}_D{Delta:.6f}_J{J}.npz")

class _Pre:
    """stand-in for a scalar seed block: returns precomputed derivative sums S[p, k, l]."""
    def __init__(self, S, a, b):
        self.S, self.a, self.b = S, a, b
    def series(self, r, chi, sig=+1, kmax=4, lmax=4, cut=None):
        assert len(r) == self.S.shape[0]
        return self.S

def evaluate_points(l, Delta, J, Dphi, U, S, pars=(PAR_A, PAR_B)):
    """stripped tensor (log scale, Fourier modes) at grid points (U, S arrays) with LC seeds, for each LC setting in pars.
    Seeds g^{c1,c2} = g^{c2,c1} (the Casimir operator depends on c1 + c2 and c1 c2 only) are computed once; the coefficients
    C^{(l)}_{i,q} (which depend only on B_{n,j}, n <= I + Q - J + l, not on the truncation N) are computed once for the largest
    setting and re-used (sub-blocks) for the others.  Returns list of (logscale, modes, i_tail, info)."""
    fd = frame_data(U, S)
    x1, x2, x3, x4 = frame_positions(fd['ru'], fd['xi'])
    z1, z3 = pol_grid(len(U))
    s, z = invariants(x1, x2, x3, x4, z1, z3)
    sb = SpinBlock(l, Delta, J, Dphi, 60)              # small radial N: its CBlocks are replaced below
    keys = {}
    for ab, (blk, vpow) in sb.blocks.items():
        keys[ab] = (tuple(sorted((blk.a, blk.b))), vpow)
    Ssets = [dict() for _ in pars]; itails = [[] for _ in pars]; info = [dict() for _ in pars]
    big = max(range(len(pars)), key=lambda j: (pars[j]['I'], pars[j]['Q']))
    for key in sorted(set(k for k, v in keys.values())):
        t0 = time.time()
        lbig = LCBlock(Delta, J, key[0], key[1], **pars[big])
        for j, par in enumerate(pars):
            lb = lbig if j == big else LCBlock(Delta, J, key[0], key[1], coeffs=lbig.C, **par)
            Ssets[j][key] = lb.series(fd['rs'], fd['chi'])
            itails[j].append(lb.i_tail)
            info[j][key] = dict(fitres_i300=float(lb.fitres[:min(par['I'], 300) + 1].max()))
        info[big][key]['time'] = round(time.time() - t0, 1)
        del lbig
    out = []
    for j in range(len(pars)):
        for ab, (key, vpow) in keys.items():
            sb.blocks[ab] = (_Pre(Ssets[j][key], *key), vpow)
        ls, M = sb.eval(s, z, fd['rs'], fd['chi'], fd['A0'], fd['E0'])
        out.append((ls, fourier_modes(M, lmax=l), np.max(itails[j], axis=0), info[j]))
    return out

class TensorMomentsLC(TensorMoments):
    @classmethod
    def build(cls, l, Delta, J, Dphi=1.75, nvalid=300, itail_max=1e-13, tol_pt=1e-2, cache=True, verbose=True, base=None, fn=None):
        """base: a TensorMoments object (radial series on a grid); default: the cached e2e data_modes file."""
        fn = lc_file(l, Delta, J) if fn is None else fn
        if cache and os.path.exists(fn):
            return cls.load(fn)
        base = TensorMoments.load(e2e_file(l, Delta, J)) if base is None else base
        g = base.g
        todo = g['ok'] & ~base.conv
        cand = np.argwhere(g['ok'] & base.conv & (g['U'] < 0.25))
        rng = np.random.default_rng(1)
        sel = cand[rng.choice(len(cand), size=min(nvalid, len(cand)), replace=False)] if len(cand) else np.zeros((0, 2), int)
        vmask = np.zeros_like(todo); vmask[tuple(sel.T)] = True
        t0 = time.time()
        pts_all = todo | vmask
        (lsA, fmA, itA, infoA), (lsB_all, fmB_all, itB_all, infoB) = evaluate_points(l, Delta, J, Dphi, g['U'][pts_all], g['S'][pts_all])
        tA = time.time() - t0; tB = 0.0
        tsel_all = todo[pts_all]
        lsB = lsB_all[tsel_all]; fmB = {k: v[tsel_all] for k, v in fmB_all.items()}; itB = itB_all[tsel_all]
        pts = todo | vmask
        vsel = vmask[pts]; tsel = todo[pts]
        keys = sorted(base.modes.keys())
        # validation against the cached radial-series values
        dev = 0.0
        if vsel.any():
            Aval = np.array([base.modes[k][vmask]*np.exp(base.logscale[vmask]) for k in keys])
            Bval = np.array([fmA[k][vsel]*np.exp(lsA[vsel]) for k in keys])
            dev = float(np.max(np.abs(Aval - Bval)/np.max(np.abs(Aval), axis=0)))
        # point-wise comparison of the two LC settings
        ref = lsA[tsel]
        VA = np.array([fmA[k][tsel] for k in keys])
        VB = np.array([fmB[k]*np.exp(lsB - ref) for k in keys])
        with np.errstate(invalid='ignore', divide='ignore'):
            ptdev = np.max(np.abs(VA - VB), axis=0)/np.max(np.abs(VA), axis=0)
        ptdev = np.where(np.isfinite(ptdev), ptdev, np.inf)
        good0 = (itA[tsel] < itail_max) & (itB < itail_max) & (ptdev < tol_pt)
        # contiguity: in every t-column accept LC points only downwards from the radial-series region, stopping at the first
        # point that fails the A/B test or whose size jumps relative to its upper neighbour by more than a power u^{+-8}
        # (isolated points inside a rejected region are where A and B are both corrupted by the cancellation of the seeds)
        idx = np.argwhere(todo)
        okpt = np.zeros_like(todo); okpt[tuple(idx[good0].T)] = True
        mag = np.zeros(g['S'].shape)
        mag[base.conv] = np.max(np.abs([base.modes[k][base.conv] for k in keys]), axis=0)*np.exp(base.logscale[base.conv] - np.max(base.logscale[base.conv]))
        magA = np.max(np.abs(VA), axis=0)*np.exp(ref - np.max(base.logscale[base.conv]))
        tmp = np.zeros(g['S'].shape); tmp[tuple(idx.T)] = magA; mag = np.where(todo, tmp, mag)
        acc = base.conv.copy()
        us = g['us']
        for j in range(g['S'].shape[1]):
            col_ok = g['ok'][:, j]
            for i in range(g['S'].shape[0] - 2, -1, -1):
                if not (todo[i, j] and okpt[i, j]):
                    continue
                if not acc[i + 1, j]:
                    continue
                if mag[i + 1, j] > 0 and mag[i, j] > 0:
                    p = np.log(mag[i + 1, j]/mag[i, j])/np.log(us[i + 1]/us[i])
                    if not (-8.0 < p < 8.0):
                        continue
                acc[i, j] = True
        good = acc[todo] & good0
        self = cls.__new__(cls)
        self.l, self.Delta, self.J, self.Dphi = l, float(Delta), int(J), float(Dphi)
        self.g = g; self.X = base.X
        idx = np.argwhere(todo)
        gi = tuple(idx[good].T)
        self.conv = base.conv.copy(); self.conv[gi] = True
        self.lcpts = np.zeros_like(self.conv); self.lcpts[gi] = True
        self.modes = {k: base.modes[k].copy() for k in keys}
        self.dmodes = {k: np.zeros_like(base.modes[k]) for k in keys}        # |A - B| at LC points (same log scale)
        for j, k in enumerate(keys):
            self.modes[k][gi] = VA[j][good]
            self.dmodes[k][gi] = np.abs(VA[j] - VB[j])[good]
        self.logscale = base.logscale.copy(); self.logscale[gi] = ref[good]
        self.valid = dev
        self.info = dict(tA=round(tA), tB=round(tB), npts=int(todo.sum()), nfill=int(good.sum()),
                         n_rejected=int(((ptdev >= tol_pt) & (itA[tsel] < itail_max)).sum()),
                         n_itail=int((itA[tsel] >= itail_max).sum()), valid=dev,
                         seedsA=list(infoA.values()), seedsB=list(infoB.values()))
        if verbose:
            print(f"   [LC l={l} Delta={Delta:.4f} J={J}] {self.info['npts']} corner points (A {tA:.0f} s, B {tB:.0f} s): filled {self.info['nfill']}, "
                  f"rejected {self.info['n_rejected']} (|A-B| > {tol_pt}), {self.info['n_itail']} with unconverged i-sum; "
                  f"validation vs radial at {int(vmask.sum())} converged points: max dev {dev:.1e}", flush=True)
        if cache:
            os.makedirs(CACHE, exist_ok=True)
            self.save(fn)
        return self

    def save(self, fn):
        g = self.g
        keys = sorted(self.modes.keys())
        np.savez_compressed(fn, l=self.l, Delta=self.Delta, J=self.J, Dphi=self.Dphi, us=g['us'], ts=g['ts'], U=g['U'], S=g['S'], W=g['W'],
                            ok=g['ok'], conv=self.conv, lcpts=self.lcpts, X=self.X, logscale=self.logscale, keys=np.array(keys),
                            modes=np.array([self.modes[k] for k in keys]), dmodes=np.array([self.dmodes[k] for k in keys]),
                            valid=self.valid, info=np.array(repr(self.info)))

    @classmethod
    def load(cls, fn):
        d = np.load(fn, allow_pickle=False)
        self = cls.__new__(cls)
        self.l, self.Delta, self.J, self.Dphi = int(d['l']), float(d['Delta']), int(d['J']), float(d['Dphi'])
        self.g = dict(us=d['us'], ts=d['ts'], U=d['U'], S=d['S'], W=d['W'], ok=d['ok'])
        self.conv = d['conv']; self.lcpts = d['lcpts']; self.logscale = d['logscale']; self.X = d['X']
        self.modes = {tuple(int(x) for x in k): m for k, m in zip(d['keys'], d['modes'])}
        self.dmodes = {tuple(int(x) for x in k): m for k, m in zip(d['keys'], d['dmodes'])}
        self.valid = float(d['valid']); self.info = str(d['info'])
        return self

    # ------------------------------------------------------------------ moments with error budget
    def _density(self, k, nu, q1, q3, modes=None):
        """integrand * quadrature weight at every grid point, common scale exp(sh)."""
        g = self.g; U, S = g['U'], g['S']
        xi = -np.log(g['us'])
        if not hasattr(self, '_Fc'): self._Fc = {}
        if nu not in self._Fc:
            self._Fc[nu] = Fnu_np(nu, list(range(0, 2*self.l + 1)), xi)
        F = self._Fc[nu][abs(q1 + q3)]
        rho, rhob = S**2, (S*U)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k - 4)*np.log(rho*rhob) + self.logscale
        logw[~g['ok']] = -np.inf
        sh = np.max(logw[np.isfinite(logw)])
        md = self.modes[(q1, q3)] if modes is None else modes[(q1, q3)]
        return np.where(np.isfinite(logw), np.exp(logw - sh), 0.0)*F[:, None]*md*g['W'], sh

    DROP_REL_ERR = 0.5       # relative error assigned to the extrapolated estimate of the missing points

    def moment_err(self, k, nu, q1, q3):
        """(value, err_lc, err_drop, lc_share).  value = Pihat (complex mantissa, log scale) including the extrapolated estimate of
        the grid points that are still missing; err_lc = SAFETY * sum over LC points of |A - B|-contributions; err_drop =
        DROP_REL_ERR * |estimate of the missing points| (+ the full size of any part bounded only by the t^{a-1} law);
        errors relative to |value|.
        Missing points are of two kinds: (i) tiny t (lightcone sum not converged, radial series not converged): estimated along
        their u-row by the power law in t through the two smallest accepted t-nodes of that row; (ii) deepest corner (A/B test
        failed): estimated along their t-column by the power law in u through the two smallest accepted u-nodes.  Exponents are
        fitted (clipped to [0, 20]) and the estimate carries the phase of the nearest accepted node."""
        dens, sh = self._density(k, nu, q1, q3)
        derr, _ = self._density(k, nu, q1, q3, modes=self.dmodes)
        g = self.g; W = g['W']; us = g['us']; ts = g['ts']
        acc = self.conv & g['ok']
        miss = g['ok'] & ~self.conv
        est = np.zeros_like(dens); bound = 0.0
        f = np.where(W != 0, dens/np.where(W != 0, W, 1), 0)
        def plaw(x1, x2, f1, f2, x):
            a1, a2 = abs(f1), abs(f2)
            if a1 <= 0 or a2 <= 0:
                return f1*0
            p = min(max(np.log(a2/a1)/np.log(x2/x1), 0.0), 20.0)
            return f1*(x/x1)**p
        for (i, j) in np.argwhere(miss):
            ra = np.where(acc[i, :])[0]; ra = ra[ra > j]
            ca = np.where(acc[:, j])[0]; ca = ca[ca > i]
            if len(ra) >= 2 and (len(ca) < 2 or ts[ra[0]] - ts[j] <= 0):
                est[i, j] = plaw(ts[ra[0]], ts[ra[1]], f[i, ra[0]], f[i, ra[1]], ts[j])*W[i, j]
            elif len(ra) >= 2 and len(ca) >= 2:
                # both directions available: tiny-t points (row neighbours at small t) vs deep corner (column neighbours at small u)
                if ts[j] < 2.0:
                    est[i, j] = plaw(ts[ra[0]], ts[ra[1]], f[i, ra[0]], f[i, ra[1]], ts[j])*W[i, j]
                else:
                    est[i, j] = plaw(us[ca[0]], us[ca[1]], f[ca[0], j], f[ca[1], j], us[i])*W[i, j]
            elif len(ca) >= 2:
                est[i, j] = plaw(us[ca[0]], us[ca[1]], f[ca[0], j], f[ca[1], j], us[i])*W[i, j]
            else:
                # no usable neighbours: bound by the t^{a-1} law of the leading integrand (a = a_x + a_y) from the nearest column
                a = 2*self.Dphi + 2*k + 2*self.l
                full = [jj for jj in range(j + 1, len(ts)) if acc[:, jj].sum() >= 2]
                if not full:
                    bound = np.inf; continue
                j1 = full[0]
                bound += np.abs(f[:, j1]).max()*(ts[j]/ts[j1])**(a - 1)*W[i, j]
        val = dens.sum() + est.sum()
        err_lc = SAFETY*float(np.abs(derr).sum()/abs(val))
        lc_share = float(abs(dens[self.lcpts].sum())/abs(val))
        err_drop = float((self.DROP_REL_ERR*np.abs(est).sum() + bound)/abs(val))
        self.last_est_share = float(abs(est.sum())/abs(val))
        return (val, sh), err_lc, err_drop, lc_share

# ---------------------------------------------------------------------- corner-refined grid (quadrature check)
UP_E2E = ((0.0, 0.01, 8), (0.01, 0.03, 10), (0.03, 0.1, 14), (0.1, 0.3, 16), (0.3, 1.0, 24))
UP_REF = ((0.0, 1e-4, 6), (1e-4, 1e-3, 8), (1e-3, 0.004, 8), (0.004, 0.01, 8), (0.01, 0.02, 8), (0.02, 0.03, 8),
          (0.03, 0.1, 14), (0.1, 0.3, 16), (0.3, 1.0, 24))

UP_REF2 = ((0.0, 1e-4, 10), (1e-4, 1e-3, 14), (1e-3, 0.004, 14), (0.004, 0.01, 14), (0.01, 0.02, 14), (0.02, 0.03, 14),
           (0.03, 0.1, 14), (0.1, 0.3, 16), (0.3, 1.0, 24))

def refined_base(l, Delta, J, panels=UP_REF):
    """the e2e grid with the two u-panels below u = 0.03 (18 nodes) replaced by 46 nodes on six panels; the radial-series data for
    u >= 0.03 are copied from the e2e cache (identical nodes), all points with u < 0.03 are left to the LC evaluation."""
    from moments import make_grid
    old = TensorMoments.load(e2e_file(l, Delta, J))
    g = make_grid(Delta, J, upanels=panels)
    n_old = 18; n_new = sum(p[2] for p in panels[:6])
    assert np.allclose(g['us'][n_new:], old.g['us'][n_old:])
    new = TensorMoments.__new__(TensorMoments)
    new.l, new.Delta, new.J, new.Dphi = old.l, old.Delta, old.J, old.Dphi
    new.g = g
    sh = g['S'].shape
    new.conv = np.zeros(sh, bool); new.conv[n_new:] = old.conv[n_old:]
    new.logscale = np.full(sh, -np.inf); new.logscale[n_new:] = old.logscale[n_old:]
    new.X = np.zeros(sh); new.X[n_new:] = old.X[n_old:]
    new.modes = {}
    for k, v in old.modes.items():
        a = np.zeros(sh, dtype=v.dtype); a[n_new:] = v[n_old:]; new.modes[k] = a
    return new

def moments_with_err(mom, ks, nus):
    """res[(k, nu, q1, q3)] = (-1)^J Pihat^{(q1,q3)}_{k,nu} (mpmath) for all modes, as e2e/t4_spinning_e2e.moments_of, but with the two
    transverse modes (l,-l), (l,l) from TensorMomentsLC.moment_err (missing points estimated); errs[(k, nu, q1, q3)] =
    (err_lc, err_drop, lc_share, estimated share) for those two."""
    import mpmath as mp
    l, J = mom.l, mom.J
    res, errs = {}, {}
    for k in ks:
        for nu in nus:
            for q1 in range(-l, l + 1):
                for q3 in range(-l, l + 1):
                    if (q1, q3) in ((l, -l), (l, l)):
                        (val, sh), el, ed, share = mom.moment_err(k, nu, q1, q3)
                        errs[(k, nu, q1, q3)] = (el, ed, share, mom.last_est_share)
                    else:
                        val, sh = mom.moment(k, nu, q1, q3)
                    res[(k, nu, q1, q3)] = (-1)**J*mp.mpc(complex(val))*mp.exp(mp.mpf(sh))
    return res, errs
