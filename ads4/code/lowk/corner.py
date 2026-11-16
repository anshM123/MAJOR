"""
Completion of the tensor-mode data of ads4/code/e2e (data_modes/) in the corner u = sqrt(rhobar_u/rho_u) -> 0.

The e2e runs evaluate the stripped spinning partial wave on a (u, t) Gauss grid with the radial series (N <= 10^4 levels) and
drop the grid points where that series has not converged (moments.series_converged): the region u <~ 30 alpha/(N t) and
very small t.  Here the dropped points are evaluated with the lightcone + asymptotic-tail representation of the scalar
seed blocks (lcblock.LCBlock), fed into the *same* spinning machinery (spinblock.SpinBlock.eval: differential operators,
(P,Q) jets, generated tensor code) - nothing else is changed.  A validation set of already-converged points at small u is
re-evaluated the same way and compared with the cached radial-series values.

Points where the lightcone sum over i (powers of rho_s) has not converged with I levels (tiny t only) remain dropped; their
weight is reported (TensorMomentsLC.dropped_weight_bound).
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E2E = os.path.join(HERE, '..', 'e2e')
if E2E not in sys.path:
    sys.path.insert(0, E2E)
from moments import TensorMoments, frame_data, frame_positions, pol_grid, fourier_modes
from spinblock import SpinBlock, invariants
from lcblock import LCBlock
CACHE = os.path.join(HERE, 'cache')
E2E_DATA = os.path.join(E2E, 'data_modes')

def e2e_file(l, Delta, J):
    return os.path.join(E2E_DATA, f"modes_l{l}_D{Delta:.6f}_J{J}.npz")

def lc_file(l, Delta, J, tag):
    return os.path.join(CACHE, f"lcmodes_l{l}_D{Delta:.6f}_J{J}_{tag}.npz")

class LCSpinBlock(SpinBlock):
    """SpinBlock whose scalar seeds are LCBlocks (same differential operators, jets and tensor code)."""
    def __init__(self, l, Delta, J, Dphi, N, I, Q, M=8, Malt=3):
        super().__init__(l, Delta, J, Dphi, 4)          # tiny radial N: the CBlocks built here are replaced below
        self.lcinfo = {}
        for ab, (blk, vpow) in list(self.blocks.items()):
            c1, c2 = blk.a, blk.b
            t0 = time.time()
            lb = LCBlock(Delta, J, c1, c2, N, I, Q, M=M, Malt=Malt)
            self.blocks[ab] = (lb, vpow)
            self.lcinfo[ab] = dict(c=(c1, c2), build=time.time() - t0, fitres=float(lb.fitres[:min(I, 300) + 1].max()))

    def i_tail(self):
        return np.max([blk.i_tail for (blk, vp) in self.blocks.values()], axis=0)

def evaluate_points(l, Delta, J, Dphi, U, S, lcpar):
    """stripped tensor (log scale, Fourier modes) at grid points (U, S arrays) with LC seeds."""
    fd = frame_data(U, S)
    x1, x2, x3, x4 = frame_positions(fd['ru'], fd['xi'])
    z1, z3 = pol_grid(len(U))
    s, z = invariants(x1, x2, x3, x4, z1, z3)
    sb = LCSpinBlock(l, Delta, J, Dphi, **lcpar)
    ls, M = sb.eval(s, z, fd['rs'], fd['chi'], fd['A0'], fd['E0'])
    fm = fourier_modes(M, lmax=l)
    return ls, fm, sb.i_tail(), sb.lcinfo

class TensorMomentsLC(TensorMoments):
    """TensorMoments with the corner completed.  Attributes: conv (points used), lcpts (points from the LC evaluation),
    valid (validation: max relative deviation LC vs radial at re-evaluated converged points)."""
    @classmethod
    def build(cls, l, Delta, J, Dphi=1.75, lcpar=None, tag='A', nvalid=400, itail_max=1e-13, cache=True, verbose=True):
        lcpar = dict(N=10000, I=1500, Q=8000, M=10) if lcpar is None else lcpar
        fn = lc_file(l, Delta, J, tag)
        if cache and os.path.exists(fn):
            return cls.load(fn)
        base = TensorMoments.load(e2e_file(l, Delta, J))
        g = base.g
        todo = g['ok'] & ~base.conv
        # validation set: converged points with u < 0.25 (the LC representation is not used near u ~ 1, where chi -> 0)
        cand = np.argwhere(g['ok'] & base.conv & (g['U'] < 0.25))
        rng = np.random.default_rng(1)
        sel = cand[rng.choice(len(cand), size=min(nvalid, len(cand)), replace=False)] if len(cand) else np.zeros((0, 2), int)
        vmask = np.zeros_like(todo); vmask[tuple(sel.T)] = True
        pts = todo | vmask
        idx = np.argwhere(pts)
        t0 = time.time()
        ls, fm, itail, info = evaluate_points(l, Delta, J, Dphi, g['U'][pts], g['S'][pts], lcpar)
        dt = time.time() - t0
        self = cls.__new__(cls)
        self.l, self.Delta, self.J, self.Dphi = l, float(Delta), int(J), float(Dphi)
        self.g = g; self.X = base.X
        self.modes = {k: v.copy() for k, v in base.modes.items()}
        self.logscale = base.logscale.copy()
        # validation
        vsel = vmask[pts]
        dev = 0.0
        if vsel.any():
            A = {key: base.modes[key][vmask]*np.exp(base.logscale[vmask]) for key in base.modes}
            B = {key: fm[key][vsel]*np.exp(ls[vsel]) for key in base.modes}
            norm = np.max([np.abs(A[key]) for key in A], axis=0)
            dev = float(max(np.max(np.abs(A[key] - B[key])/norm) for key in A))
        self.valid = dev
        # fill the dropped points; keep those where the lightcone i-sum converged
        tsel = todo[pts]
        good = tsel & (itail < itail_max)
        conv = base.conv.copy()
        lcpts = np.zeros_like(conv)
        gi = idx[good]
        conv[tuple(gi.T)] = True; lcpts[tuple(gi.T)] = True
        for key in self.modes:
            self.modes[key][tuple(gi.T)] = fm[key][good]
        self.logscale[tuple(gi.T)] = ls[good]
        self.conv = conv; self.lcpts = lcpts
        self.info = dict(time=dt, npts=int(pts.sum()), nfill=int(good.sum()), nleft=int((tsel & ~good).sum()), valid=dev,
                         lcpar=lcpar, seeds=info)
        if verbose:
            print(f"   [LC l={l} Delta={Delta:.4f} J={J} tag={tag}] {self.info['npts']} points in {dt:.0f} s; filled {self.info['nfill']}, "
                  f"left {self.info['nleft']} (i-sum not converged); validation max dev {dev:.1e}", flush=True)
        if cache:
            os.makedirs(CACHE, exist_ok=True)
            self.save(fn)
        return self

    def save(self, fn):
        g = self.g
        keys = sorted(self.modes.keys())
        np.savez_compressed(fn, l=self.l, Delta=self.Delta, J=self.J, Dphi=self.Dphi, us=g['us'], ts=g['ts'], U=g['U'], S=g['S'], W=g['W'],
                            ok=g['ok'], conv=self.conv, lcpts=self.lcpts, X=self.X, logscale=self.logscale, keys=np.array(keys),
                            modes=np.array([self.modes[k] for k in keys]), valid=self.valid, info=np.array(repr(self.info)))

    @classmethod
    def load(cls, fn):
        d = np.load(fn, allow_pickle=False)
        self = cls.__new__(cls)
        self.l, self.Delta, self.J, self.Dphi = int(d['l']), float(d['Delta']), int(d['J']), float(d['Dphi'])
        self.g = dict(us=d['us'], ts=d['ts'], U=d['U'], S=d['S'], W=d['W'], ok=d['ok'])
        self.conv = d['conv']; self.lcpts = d['lcpts']; self.logscale = d['logscale']; self.X = d['X']
        self.modes = {tuple(int(x) for x in k): m for k, m in zip(d['keys'], d['modes'])}
        self.valid = float(d['valid']); self.info = str(d['info'])
        return self

    def moment_parts(self, k, nu, q1, q3):
        """moment split into (radial-series points, LC points), same normalisation as TensorMoments.moment."""
        v_all, sh = self.moment(k, nu, q1, q3)
        save = self.conv
        self.conv = save & ~self.lcpts
        lg = self.logscale
        self.logscale = np.where(self.lcpts, -np.inf, lg)
        v_rad, sh2 = self.moment(k, nu, q1, q3)
        self.conv, self.logscale = save, lg
        return v_all, v_rad*np.exp(sh2 - sh), sh

    def dropped_weight_bound(self, k):
        """crude bound on the relative weight of the points still dropped (all at tiny t): the fraction of
        int t^{a-1} e^{-t} dt (a = a_x + a_y, the t-dependence of the leading integrand) below the largest dropped t."""
        g = self.g
        left = g['ok'] & ~self.conv
        if not left.any():
            return 0.0
        tmax = float(np.max(g['S'][left]*2*(self.Delta - self.J - 2 + (self.Delta + self.J - 1)*g['U'][left])))
        from scipy.special import gammainc
        a = 2*self.Dphi + k - 1 + 2*self.l + k + 1
        return float(gammainc(a, tmax))
