"""
lfcore.py -- shared light-front machinery for the free-field null-plane tests (d = 3).

One-body light-front QM (np1_kernels.py):  on L^2(dK dy) per particle,
    a(f) = K f(y),   h = p^2/(4K)  (p = -i d_y),   n(g) = -(1/2){g, p}.
A state is stored slice-wise: for each K-node (weight wK) a complex array chi(y) on a uniform periodic grid.
Covariant <-> LF:  p^0 = h + K,  p^x = h - K,  p^y = q,  h = q^2/(4K).
chi(K,q) = psi_cov(p)/sqrt(2K)  (covariant measure d^2p/(2p^0) = dK dq/(2K)); overall constants are irrelevant
because every quantity is divided by <psi|psi>.
"""
import numpy as np

# ------------------------------------------------------------------ profiles f_m and cutoffs
def fm(y, m):
    """f_m(y) = (1+y^2) e^{2 i m arctan y}, and its derivative f_m' = 2 (y + i m) e^{2 i m arctan y}."""
    ph = np.exp(2j*m*np.arctan(y))
    return (1 + y**2)*ph, 2*(y + 1j*m)*ph

def smooth_cut(y, R, p=8):
    """C^infinity super-Gaussian cutoff c(y) = exp(-(y/R)^p) and its derivative."""
    c = np.exp(-(y/R)**p)
    return c, -p*(y/R)**(p - 1)/R*c

def profile(y, m, part, R=None):
    """Real profile F = Re/Im f_m (optionally times the smooth cutoff at scale R); returns F, F'."""
    f, fp = fm(y, m)
    if R is not None:
        c, cp = smooth_cut(y, R)
        f, fp = f*c, fp*c + f*cp
    if part == 're':
        return f.real, fp.real
    return f.imag, fp.imag

# ------------------------------------------------------------------ slice grids
class Slice:
    def __init__(self, K, w, L, dy):
        N = int(2**np.ceil(np.log2(2*L/dy)))
        self.K, self.w, self.N = K, w, N
        self.dy = 2*L/N
        self.y = -L + self.dy*np.arange(N)
        self.q = 2*np.pi*np.fft.fftfreq(N, self.dy)
        self.h = self.q**2/(4*K)

    def from_q(self, chiq):
        """chi(y_j) = sum_k chiq(q_k) e^{i q_k y_j} (up to a constant)."""
        return np.fft.ifft(chiq*np.exp(1j*self.q*self.y[0]))

    def H(self, v):
        return np.fft.ifft(self.h*np.fft.fft(v))

    def ip(self, u, v):           # <u, v> on the slice (no K-weight)
        return np.vdot(u, v)*self.dy

    def hform(self, v):           # <v, h v>
        F = np.fft.fft(v)
        return float(np.sum(self.h*np.abs(F)**2))*self.dy/self.N

    def evolve(self, v, ts):      # e^{-i h t} v for an array of times -> (len(ts), N)
        F = np.fft.fft(v)
        return np.fft.ifft(np.exp(-1j*np.outer(ts, self.h))*F[None, :], axis=1)

def gauss_legendre(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b - a)*x + 0.5*(b + a), 0.5*(b - a)*w

# ------------------------------------------------------------------ regulator
def gh_nodes(n=80):
    return np.polynomial.hermite.hermgauss(n)

def I_eta_from_G(Gfun, eta, n=80):
    """I_eta = -int rho_eta'(t) G(t) dt,  G(t) = <psi,[A,A_t] psi>/i,  rho_eta Gaussian (hat rho >= 0, even).
    = (1/sqrt(pi)) sum_j w_j (sqrt2 x_j/eta) G(sqrt2 eta x_j)."""
    x, w = gh_nodes(n)
    t = np.sqrt(2)*eta*x
    return float(np.sum(w*np.sqrt(2)*x/eta*Gfun(t))/np.sqrt(np.pi))
