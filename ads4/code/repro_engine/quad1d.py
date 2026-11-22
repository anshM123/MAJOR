"""repro_engine/quad1d.py -- one-dimensional PRODUCT-INTEGRATION rules (own code) for the seven variables of the section D.

For a variable x with real weight W(x) and nu-dependent factor phase(x)^{i nu} (phase = atom^{E1}), the rule is
    int W(x) atom(x)^{i nu E1} f(x) dx  ~  sum_k omega_k(nu) f(x_k),     omega_k(nu) = int W atom^{i nu E1} l_k(x) dx,
with nodes x_k = Gauss nodes of the real weight (scipy roots_genlaguerre / roots_jacobi, mapped) and l_k the Lagrange basis on them
(barycentric form).  omega_k(nu) is computed with a fine double-exponential / exp-substituted trapezoid rule, NOT with closed-form
Gamma/Beta functions (so this route shares nothing with the engine's integrate_fast).  Interpolatory => exact for polynomial f of
degree < n at every nu; Gauss (degree < 2n) at nu = 0.

Kinds (variable, real weight, nu-atom):
  'lag'  : X in (0, inf),   X^c e^{-X},                    X^{i nu e}
  'tjac' : t in (0, 1),     (1 - t)^g t^{-1/2},            (1 - t)^{i nu e}          (t = s^2; integrand even in s)
  'cosw' : phi in (-pi/2, pi/2), cos(phi)^h,              cos(phi)^{i nu e}
  'sinw' : alpha in (0, pi), sin(alpha)^k,                sin(alpha)^{i nu e}
"""
import numpy as np
from scipy.special import roots_genlaguerre, roots_jacobi


def bary_weights(x):
    n = len(x)
    w = np.ones(n)
    for j in range(n):
        d = x[j] - np.delete(x, j)
        w[j] = 1.0/np.prod(d)
    return w


def lagrange_matrix(xn, t):
    """L[k, i] = l_k(t_i) for nodes xn (barycentric second form; exact at nodes)."""
    wb = bary_weights(xn)
    D = t[None, :] - xn[:, None]
    hit = np.abs(D) < 1e-300
    D = np.where(hit, 1.0, D)
    T = wb[:, None]/D
    L = T/np.sum(T, axis=0, keepdims=True)
    cols = np.where(hit.any(axis=0))[0]
    for i in cols:
        L[:, i] = hit[:, i].astype(float)
    return L


def _ts01(h=1/512., tmax=4.2):
    """tanh-sinh points on (0,1): t, 1 - t (accurate), weights."""
    tau = np.arange(-tmax, tmax + h/2, h)
    z = np.pi*np.sinh(tau)
    t = 1.0/(1.0 + np.exp(-z)); u = 1.0/(1.0 + np.exp(z))          # t, 1 - t
    w = h*np.pi*np.cosh(tau)*t*u
    keep = (t > 0) & (u > 0) & (w > 0)
    return t[keep], u[keep], w[keep]


class Rule:
    def __init__(self, kind, n, par, e):
        self.kind, self.n, self.par, self.e = kind, n, par, e
        if kind == 'lag':
            x, _ = roots_genlaguerre(n, par)
            tt = np.arange(-60.0, 5.6, 0.01)
            X = np.exp(tt)
            fine_x = X; lw = (par + 1)*tt - X                         # log(W(x) dx/dt)
            dw = 0.01*np.exp(lw); la = tt                             # log atom = log X
        elif kind == 'tjac':
            y, _ = roots_jacobi(n, par, -0.5)                         # weight (1-y)^g (1+y)^{-1/2} on [-1,1]
            x = (y + 1)/2
            t, u, w = _ts01()
            fine_x = t; dw = w*u**par/np.sqrt(t); la = np.log(u)       # atom = 1 - t
        elif kind == 'cosw':
            y, _ = roots_jacobi(n, par, par)
            x = y*np.pi/2
            t, u, w = _ts01()
            fine_x = np.pi*(t - 0.5); c = np.sin(np.pi*np.minimum(t, u))   # cos(phi)
            dw = np.pi*w*c**par; la = np.log(c)
        elif kind == 'sinw':
            y, _ = roots_jacobi(n, par, par)
            x = (y + 1)*np.pi/2
            t, u, w = _ts01()
            fine_x = np.pi*t; c = np.sin(np.pi*np.minimum(t, u))
            dw = np.pi*w*c**par; la = np.log(c)
        else:
            raise ValueError(kind)
        self.x = np.asarray(x, float)
        self._L = lagrange_matrix(self.x, fine_x)*dw[None, :]       # (n, nfine)
        self._la = la

    def weights(self, nus):
        """omega[k, nu] (complex)."""
        ph = np.exp(1j*np.outer(self._la, self.e*np.asarray(nus, float)))     # (nfine, nnu)
        return self._L @ ph


def selftest():
    from scipy.special import gamma, beta
    out = []
    nus = np.array([0.0, 0.5, 3.0, 12.0])
    # 'lag': int X^{c + j + i nu e} e^{-X} = Gamma(c + j + 1 + i nu e)
    import mpmath as mp
    for kind, par, e, exact in (
            ('lag', 2.25, -0.5, lambda j, nu: complex(mp.gamma(2.25 + j + 1 + 1j*nu*(-0.5)))),
            ('lag', -0.25, -0.5, lambda j, nu: complex(mp.gamma(-0.25 + j + 1 - 0.5j*nu))),
            ('tjac', 1.75, -0.5, lambda j, nu: complex(mp.beta(j + 0.5, 1.75 - 0.5j*nu + 1))),
            ('cosw', 0.5, 1.0, lambda j, nu: complex(mp.quad(lambda p: mp.cos(p)**(0.5 + 1j*nu)*p**j, [-mp.pi/2, 0, mp.pi/2]))),
            ('sinw', -0.5, -1.0, lambda j, nu: complex(mp.quad(lambda p: mp.sin(p)**(-0.5 - 1j*nu)*p**j, [0, mp.pi/2, mp.pi])))):
        R = Rule(kind, 8, par, e)
        W = R.weights(nus)
        err = 0
        for j in (0, 1, 3, 7):
            for i, nu in enumerate(nus):
                ex = exact(j, nu)
                err = max(err, abs(np.sum(W[:, i]*R.x**j) - ex)/abs(ex))
        out.append(f'  rule {kind} par={par} e={e}: max rel err on x^j (j<8), nu in {list(nus)}: {err:.1e}')
    return out


if __name__ == '__main__':
    print('\n'.join(selftest()))
