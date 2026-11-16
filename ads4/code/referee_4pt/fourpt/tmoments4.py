"""
Fourier modes (two Regge pairs) of the exact stripped doubly-spinning partial waves on the (u, t) integration grid of e2e/moments.py,
and their nu-harmonic tensor Regge moments (kernel moments4.kernel_np).

For l = 1 the basis combinations  G^{ab} := s13^2 s24^2 W[v_a(Delta,J), v_b(Delta,J)]  (a, b = 1, 2; t1_structures_l1) are stored;
any coupling combination follows by bilinearity.
"""
import os, sys, pickle, numpy as np, mpmath as mp, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from moments import make_grid, frame_data, frame_positions, series_converged, series_X      # noqa: E402
from spinblock2 import SpinBlock4, invariants4, monomial_values
from moments4 import pol_grid4, fourier4, kernel_np
from embed2 import Delta as sD, J as sJ, CACHE

DATADIR = os.path.join(CACHE, 'modes4')

def basis_coeffs(l, Delta, J):
    ws, vecs = pickle.load(open(os.path.join(CACHE, 'consbasis_l%d.pkl' % l), 'rb'))
    sub = {sD: sp.nsimplify(Delta), sJ: J}
    return ws, [{w: float(sp.N(v[i].subs(sub), 30)) for i, w in enumerate(ws)} for v in vecs]

class TensorMoments4:
    def __init__(self, l, Delta, J, N, cut_mask=25.0, chunk=600, **gkw):
        self.l, self.Delta, self.J = l, Delta, J
        g = make_grid(Delta, J, **gkw); self.g = g
        ok = g['ok']
        fd = frame_data(g['U'][ok], g['S'][ok])
        conv = series_converged(fd['rs'], fd['chi'], J, N, cut=cut_mask)
        self.conv = np.zeros(g['S'].shape, dtype=bool); self.conv[ok] = conv
        self.X = np.zeros(g['S'].shape); self.X[ok] = series_X(fd['rs'], fd['chi'], J, N)
        fd = {key: v[conv] for key, v in fd.items()}
        npts = fd['rs'].size
        ws, vecs = basis_coeffs(l, Delta, J)
        nb = len(vecs)
        combos = {(a, b): (vecs[a], vecs[b]) for a in range(nb) for b in range(nb)}
        sb = SpinBlock4(l, Delta, J, N)
        fs, monos = sb.phi
        modes = {ab: {} for ab in combos}
        logscale = np.zeros(npts)
        for c0 in range(0, npts, chunk):
            sl = slice(c0, min(npts, c0 + chunk))
            x1, x2, x3, x4 = frame_positions(fd['ru'][sl], fd['xi'][sl])
            Xd = {1: x1, 2: x2, 3: x3, 4: x4}
            z1, z3, z2, z4 = pol_grid4(fd['xi'][sl], l)
            s, z = invariants4(Xd, {1: z1, 2: z2, 3: z3, 4: z4})
            ls, C = sb.eval_coeffs(s, fd['rs'][sl], fd['chi'][sl], fd['A0'][sl], fd['E0'][sl], combos)
            logscale[sl] = ls
            mv = monomial_values(monos, z)                      # [nmono, n, nth^4]
            for ab in combos:
                M = np.einsum('mp,mp...->p...', C[ab], mv)
                fm = fourier4(M, l)
                for q, v in fm.items():
                    modes[ab].setdefault(q, []).append(v)
        self.logscale = np.full(g['S'].shape, -np.inf); self.logscale[self.conv] = logscale
        self.modes = {}
        for ab in combos:
            self.modes[ab] = {}
            for q, lst in modes[ab].items():
                arr = np.zeros(g['S'].shape, dtype=complex); arr[self.conv] = np.concatenate(lst)
                self.modes[ab][q] = arr

    def moment(self, k, nu, ab, q, cL=None, xfac=None):
        """(complex mantissa, log scale) of Pihat^q_{k,nu}[G^{ab}]; if cL is given (dict (a,b)->weight), of the combination."""
        g = self.g; U, S = g['U'], g['S']
        xi = -np.log(g['us'])
        q1, q3, q2, q4 = q
        key = (nu, q1 + q3, q2 + q4)
        if not hasattr(self, '_K'): self._K = {}
        if key not in self._K:
            self._K[key] = kernel_np(nu, [(q1 + q3, q2 + q4)], xi)[(q1 + q3, q2 + q4)]
        F = self._K[key]
        rho, rhob = S**2, (S*U)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k-4)*np.log(rho*rhob) + self.logscale
        logw[~g['ok']] = -np.inf
        if xfac is not None:
            Xs = self.X/xfac
            with np.errstate(divide='ignore', invalid='ignore'):
                keep = (Xs > 8) & (Xs - 10*np.log(np.maximum(Xs, 1.0)) > 25.0)
            logw = np.where(keep, logw, -np.inf)
        sh = np.max(logw)
        G = self.modes[ab][q] if cL is None else sum(w*self.modes[k_][q] for k_, w in cL.items())
        val = np.sum(np.exp(logw - sh)*F[:, None]*G*g['W'])
        return val, sh

    def tail_estimate(self, k, nu, ab, q, xfac=1.5):
        v1, s1 = self.moment(k, nu, ab, q); v2, s2 = self.moment(k, nu, ab, q, xfac=xfac)
        diff = abs(v1 - v2*np.exp(s2 - s1))/abs(v1)
        return diff/(xfac**(k + 0.5) - 1)

    def save(self, fn):
        g = self.g
        d = dict(l=self.l, Delta=self.Delta, J=self.J, us=g['us'], ts=g['ts'], U=g['U'], S=g['S'], W=g['W'], ok=g['ok'], conv=self.conv,
                 X=self.X, logscale=self.logscale)
        keys = []
        arrs = []
        for ab, md in self.modes.items():
            for q, v in md.items():
                keys.append(ab + q); arrs.append(v)
        np.savez_compressed(fn, keys=np.array(keys), modes=np.array(arrs), **d)

    @classmethod
    def load(cls, fn):
        d = np.load(fn)
        self = cls.__new__(cls)
        self.l, self.Delta, self.J = int(d['l']), float(d['Delta']), int(d['J'])
        self.g = dict(us=d['us'], ts=d['ts'], U=d['U'], S=d['S'], W=d['W'], ok=d['ok'])
        self.conv, self.logscale, self.X = d['conv'], d['logscale'], d['X']
        self.modes = {}
        for key, v in zip(d['keys'], d['modes']):
            key = tuple(int(x) for x in key)
            self.modes.setdefault(key[:2], {})[key[2:]] = v
        return self

def default_N(l, Delta, c=400):
    return int(min(10000 if l == 1 else 8000, max(4000, c*Delta)))

def get_moments(l, Delta, J, N=None, cache=True):
    fn = os.path.join(DATADIR, f"modes4_l{l}_D{Delta:.6f}_J{J}.npz")
    if cache and os.path.exists(fn):
        return TensorMoments4.load(fn)
    N = default_N(l, Delta) if N is None else N
    mom = TensorMoments4(l, Delta, J, N)
    if cache:
        os.makedirs(DATADIR, exist_ok=True)
        mom.save(fn)
    return mom
