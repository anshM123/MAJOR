"""
Regge moments of exact blocks on the real section w, wbar < 0, in the symmetric frame of CHMRS (C.2)
   (x_1, x_2, x_3, x_4) = (-y, x/x^2, y, -x/x^2),  x, y future timelike,  r = |x||y|, eta = -x.y/(|x||y|) = cosh xi.

Scalar:  Pihat_{k,nu}[G] = int_1^oo [d eta] P_{-1/2+i nu}(eta) int_0^{r_max} dr r^{k-2} G(r, eta)        ([d eta] = 2 d eta, CHMRS (2.50))
Tensor:  Pihat^{(q1,q3)}_{k,nu}[G] = int [d eta] F^nu_{q1+q3}(xi) int dr r^{k-2} G_{q1 q3}(r, eta),
   G_{q1 q3} = Fourier modes of the stripped frame tensor  G(z_1(th_1), z_3(th_3)),  z(th) = yhat + cos th e_par + sin th e_perp
   (e_par: unit tangent at yhat pointing to xhat).  These are the nu-harmonic tensor Regge moments: see SPINNING_E2E.md sect. 3.

Integration variables (u, t):  sqrt(rho) = s = t/(2 (alpha + alphabar u)),  sqrt(rhobar) = s u   (rho >= rhobar are the
u-channel radial coordinates, r = sqrt(rho rhobar), e^{xi} = 1/u), so that a heavy block decays like e^{-t} uniformly in u:
   Pihat = (1/2) int du ds 4 s^3 u (rho - rhobar) (rho rhobar)^{(k-4)/2} F G.
"""
import numpy as np, mpmath as mp
from formulas import alpha_pair, Fnu_np

def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b-a)*x + 0.5*(b+a), 0.5*(b-a)*w

def make_grid(Delta, J, upanels=((0.0, 0.02, 10), (0.02, 0.25, 16), (0.25, 1.0, 28)), tpanels=((0.0, 6.0, 20), (6.0, 28.0, 28), (28.0, 85.0, 24))):
    us, uw = zip(*[gl(a, b, n) for (a, b, n) in upanels]); us, uw = np.concatenate(us), np.concatenate(uw)
    ts, tw = zip(*[gl(a, b, n) for (a, b, n) in tpanels]); ts, tw = np.concatenate(ts), np.concatenate(tw)
    al, alb = alpha_pair(Delta, J)
    U, T = np.meshgrid(us, ts, indexing='ij')
    scale = 2*(al + alb*U)
    S = T/scale
    W = np.outer(uw, tw)/scale                     # du ds
    ok = S < 0.985
    return dict(us=us, ts=ts, U=U, S=S, W=W, ok=ok)

def frame_data(U, S):
    """flattened arrays of everything needed at the grid points (only those with S < 1)."""
    s1 = S; s2 = S*U
    rs = (1-s1)*(1-s2)/((1+s1)*(1+s2))
    chi = np.log((1-s2)*(1+s1)/((1+s2)*(1-s1)))
    e1 = 4*s1/(1+s1)**2; e2 = 4*s2/(1+s2)**2
    A0 = e1*e2/(4*rs); E0 = (2-e1)*(2-e2)/(4*rs)
    ru = S*S*U; xi = -np.log(U)
    return dict(rs=rs, chi=chi, A0=A0, E0=E0, ru=ru, xi=xi)

def frame_positions(ru, xi):
    q = np.sqrt(ru)
    z0 = np.zeros_like(q)
    y = np.stack([q, z0, z0], axis=-1)
    x = np.stack([q*np.cosh(xi), q*np.sinh(xi), z0], axis=-1)
    x2 = x/(-q*q)[:, None]            # x/x^2, x^2 = -q^2
    return -y, x2, y, -x2

NTH = 5
THETAS = 2*np.pi*np.arange(NTH)/NTH
def pol_grid(npts, nth=NTH):
    th = 2*np.pi*np.arange(nth)/nth
    zz = np.stack([np.ones(nth), np.cos(th), np.sin(th)], axis=-1)
    z1 = np.broadcast_to(zz[None, :, None, :], (npts, nth, nth, 3)).copy()
    z3 = np.broadcast_to(zz[None, None, :, :], (npts, nth, nth, 3)).copy()
    return z1, z3

def fourier_modes(M, lmax=2):
    """M[npts, nth, nth] -> dict (q1, q3) -> array[npts] with  M = sum G_{q1q3} e^{i(q1 th1 + q3 th3)}"""
    nth = M.shape[1]
    th = 2*np.pi*np.arange(nth)/nth
    out = {}
    for q1 in range(-lmax, lmax+1):
        for q3 in range(-lmax, lmax+1):
            ph = np.exp(-1j*(q1*th[:, None] + q3*th[None, :]))
            out[(q1, q3)] = np.tensordot(M, ph, axes=([1, 2], [0, 1]))/nth**2
    return out

def series_converged(rs, chi, J, N, cut=25.0):
    """True where the radial series with N levels (and the derivative sums) has converged.  With X = -N log(rho_big) - J chi
    the observed truncation error of the conserved tensor is ~ e^{-X} X^5 (t4b_convergence.py: 3e-9 at X = 37, 1e-4 at X = 26,
    < 1e-16 at X > 50); points with  X - 6 log X < cut  are dropped."""
    X = -N*(np.log(rs) + chi) - J*chi
    with np.errstate(divide='ignore', invalid='ignore'):
        crit = X - 6*np.log(np.maximum(X, 1.0))
    return (X > 8) & (crit > cut)

class ScalarMoments:
    """Regge moments of the scalar partial wave of <phi_1 phi_2 phi_1 phi_2> (dimensions D1, D2) stripped as
       G = (s13)^{D1} (s24)^{D2} K g = u'^{-(D1+D2)/2} g^{a,b},  a = -(D1-D2)/2, b = (D1-D2)/2,  on w, wbar < 0."""
    def __init__(self, Delta, J, D1, D2, N, **gkw):
        from blocks import CBlock
        self.Delta, self.J, self.D1, self.D2 = Delta, J, D1, D2
        g = make_grid(Delta, J, **gkw); self.g = g
        ok = g['ok']
        fd = frame_data(g['U'][ok], g['S'][ok])
        b = (D1 - D2)/2
        blk = CBlock(Delta, J, b, b, N)           # reflected block g^{-a,b} = g^{b,b}
        Sr = blk.series(fd['rs'], fd['chi'], sig=+1)[:, 0, 0]
        # g^{a,b}(u,v) = (-1)^J v^{-b} gt;  u' = 1/A0^2, v' = (E0/A0)^2
        logG = (Delta*np.log(4*fd['rs']) + np.log(np.abs(Sr)) - 2*b*np.log(fd['E0']/fd['A0']) + (D1 + D2)*np.log(fd['A0']))
        self.logG = np.full(g['S'].shape, -np.inf); self.logG[ok] = logG      # times (-1)^J
        self.sgn = np.ones(g['S'].shape); self.sgn[ok] = np.sign(Sr)

    def moment(self, k, nu):
        """returns (sign, log|.|) of (-1)^J Pihat_{k,nu}[G]"""
        g = self.g; U, S = g['U'], g['S']
        F = Fnu_np(nu, [0], -np.log(g['us']))[0].real
        rho, rhob = S**2, (S*U)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k-4)*np.log(rho*rhob) + self.logG
        logw[~g['ok']] = -np.inf
        sh = np.max(logw)
        val = np.sum(np.exp(logw - sh)*self.sgn*F[:, None]*g['W'])
        return np.sign(val), np.log(abs(val)) + sh

class TensorMoments:
    """Fourier modes of the stripped frame tensor of the spinning partial wave on the (u, t) grid."""
    def __init__(self, l, Delta, J, Dphi, N, prec='d', cut_mask=25.0, **gkw):
        from spinblock import SpinBlock, invariants
        self.l, self.Delta, self.J, self.Dphi = l, Delta, J, Dphi
        g = make_grid(Delta, J, **gkw); self.g = g
        ok = g['ok']
        fd = frame_data(g['U'][ok], g['S'][ok])
        npts = fd['rs'].size
        x1, x2, x3, x4 = frame_positions(fd['ru'], fd['xi'])
        z1, z3 = pol_grid(npts)
        s, z = invariants(x1, x2, x3, x4, z1, z3)
        sb = SpinBlock(l, Delta, J, Dphi, N, prec=prec); self.sb = sb
        ls, M = sb.eval(s, z, fd['rs'], fd['chi'], fd['A0'], fd['E0'])
        # points where the radial series (and its 4th derivatives) is not converged with N levels are dropped:
        # the derivative sums there are garbage (not small), cf. the truncation rule in cblock.c
        conv = series_converged(fd['rs'], fd['chi'], J, N, cut=cut_mask)
        M = np.where(conv[:, None, None], M, 0.0)
        self.conv = np.zeros(g['S'].shape, dtype=bool); self.conv[ok] = conv
        self.modes = {}
        fm = fourier_modes(M, lmax=l)
        self.logscale = np.full(g['S'].shape, -np.inf); self.logscale[ok] = ls
        for key, v in fm.items():
            arr = np.zeros(g['S'].shape, dtype=complex); arr[ok] = v
            self.modes[key] = arr
        self.Mraw = M; self.fd = fd

    def moment(self, k, nu, q1, q3, fun=None):
        """Pihat^{(q1,q3)}_{k,nu} of the stripped tensor: returns (complex mantissa, log scale):  value = mant * exp(logscale).
        fun: optional replacement of the harmonic function (callable of xi array)."""
        g = self.g; U, S = g['U'], g['S']
        xi = -np.log(g['us'])
        F = Fnu_np(nu, [q1 + q3], xi)[q1 + q3] if fun is None else fun(xi)
        rho, rhob = S**2, (S*U)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k-4)*np.log(rho*rhob) + self.logscale
        logw[~g['ok']] = -np.inf
        sh = np.max(logw)
        val = np.sum(np.exp(logw - sh)*F[:, None]*self.modes[(q1, q3)]*g['W'])
        return val, sh

    def coverage(self):
        """smallest u and t at which all grid points are converged (diagnostic for the dropped region)"""
        g = self.g
        bad = g['ok'] & ~self.conv
        return dict(frac_dropped=float(bad.sum())/float(g['ok'].sum()),
                    u_all_ok=float(min([u for i, u in enumerate(g['us']) if not bad[i, len(g['ts'])//4:].any()] or [np.nan])))
