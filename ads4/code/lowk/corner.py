"""
Completion of the tensor-mode data of ads4/code/e2e (data_modes/) in the corner u = sqrt(rhobar_u/rho_u) -> 0.

The e2e runs evaluate the stripped spinning partial wave on a (u, t) Gauss grid with the radial series (N <= 10^4 levels) and
drop the grid points where that series has not converged (moments.series_converged): the region u <~ 30 alpha/(N t) and
very small t.  Here the dropped points are evaluated with the lightcone + asymptotic-tail representation of the scalar
seed blocks (lcblock.LCBlock), fed into the *same* spinning machinery (spinblock.SpinBlock.eval: differential operators,
(P,Q) jets, generated tensor code) - nothing else is changed.

Error control (per grid point, then per moment):
  * every LC point is evaluated with two independent settings A = (N, I, Q, M) = (10000, 1500, 8000, 10) and
    B = (8000, 1200, 6000, 8) (different radial depth, fit window, number of asymptotic terms, lightcone depth);
    A is used, |A - B| is the error estimate; points with |A - B| > tol_pt * max|A| (relative to the largest mode at that
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
PAR_A = dict(N=10000, I=1500, Q=8000, M=10)
PAR_B = dict(N=8000, I=1200, Q=6000, M=8)

def e2e_file(l, Delta, J):
    return os.path.join(E2E_DATA, f"modes_l{l}_D{Delta:.6f}_J{J}.npz")

def lc_file(l, Delta, J):
    return os.path.join(CACHE, f"lcmodes_l{l}_D{Delta:.6f}_J{J}.npz")

class LCSpinBlock(SpinBlock):
    """SpinBlock whose scalar seeds are LCBlocks (same differential operators, jets and tensor code)."""
    def __init__(self, l, Delta, J, Dphi, N, I, Q, M=8, Malt=3):
        super().__init__(l, Delta, J, Dphi, 60)         # small radial N: the CBlocks built here are replaced below
        self.lcinfo = {}
        for ab, (blk, vpow) in list(self.blocks.items()):
            c1, c2 = blk.a, blk.b
            t0 = time.time()
            lb = LCBlock(Delta, J, c1, c2, N, I, Q, M=M, Malt=Malt)
            self.blocks[ab] = (lb, vpow)
            self.lcinfo[ab] = dict(c=(c1, c2), build=round(time.time() - t0, 1), fitres_i300=float(lb.fitres[:min(I, 300) + 1].max()))

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
    @classmethod
    def build(cls, l, Delta, J, Dphi=1.75, nvalid=300, itail_max=1e-13, tol_pt=1e-3, cache=True, verbose=True, base=None, fn=None):
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
        lsA, fmA, itA, infoA = evaluate_points(l, Delta, J, Dphi, g['U'][todo | vmask], g['S'][todo | vmask], PAR_A)
        tA = time.time() - t0
        lsB, fmB, itB, infoB = evaluate_points(l, Delta, J, Dphi, g['U'][todo], g['S'][todo], PAR_B)
        tB = time.time() - t0 - tA
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
        good = (itA[tsel] < itail_max) & (itB < itail_max) & (ptdev < tol_pt)
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

    def moment_err(self, k, nu, q1, q3):
        """(value, err_lc, err_drop, lc_share): value = Pihat (complex mantissa * e^sh as in TensorMoments.moment, returned as
        (mant, sh)); err_lc = sum over LC points of |A - B|-contributions; err_drop = |power-law extrapolation of the still
        missing points| (counted entirely as error); all errors relative to |value|."""
        dens, sh = self._density(k, nu, q1, q3)
        derr, _ = self._density(k, nu, q1, q3, modes=self.dmodes)
        val = dens.sum()
        err_lc = float(np.abs(derr).sum()/abs(val))
        lc_share = float(abs(dens[self.lcpts].sum())/abs(val))
        # extrapolation of the missing points (ok & ~conv): per t column, density/W ~ c u^p from the two smallest accepted u
        g = self.g; W = g['W']; us = g['us']
        miss = g['ok'] & ~self.conv
        extra = 0.0
        for j in np.where(miss.any(axis=0))[0]:
            acc = np.where(self.conv[:, j] & g['ok'][:, j])[0]
            mi = np.where(miss[:, j])[0]
            if len(acc) >= 2:
                i1, i2 = acc[0], acc[1]
                f1, f2 = abs(dens[i1, j])/W[i1, j], abs(dens[i2, j])/W[i2, j]
                p = np.log(max(f2, 1e-300)/max(f1, 1e-300))/np.log(us[i2]/us[i1]) if f1 > 0 and f2 > 0 else 0.0
                p = min(max(p, 0.0), 20.0)
                extra += sum(f1*(us[i]/us[i1])**p*W[i, j] for i in mi if us[i] < us[i1])
                # missing points above the first accepted one (should not occur): bound by the neighbour value
                extra += sum(f1*W[i, j] for i in mi if us[i] > us[i1])
            else:
                # (tiny t) no usable neighbours in u: bound the whole column by the t^{a-1} law of the leading integrand
                # (a = a_x + a_y), scaled from the nearest larger t column without missing points
                a = 2*self.Dphi + 2*k + 2*self.l
                ts = g['ts']
                full = [jj for jj in range(j + 1, len(ts)) if not miss[:, jj].any() and g['ok'][:, jj].any()]
                if not full:
                    extra = np.inf; continue
                j1 = full[0]
                c1 = np.abs(dens[:, j1]).sum()
                extra += c1*(ts[j]/ts[j1])**(a - 1)*(np.sum(W[:, j])/np.sum(W[:, j1]))
        err_drop = float(extra/abs(val))
        return (val, sh), err_lc, err_drop, lc_share
