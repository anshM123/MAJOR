r"""
m1_remainder/m1_weight.py : the published endpoint-repaired k = 2 continuum weight of Omega^EE (read-only data
prl/data/k2/k2_continuum_D_<Dg>_rep.json), evaluated in mpmath (any precision) and in numpy (float64, vectorised).

  chi(nu)   = Q(nu) sum_signed x_k exp(-(nu - c)^2/(2 s^2)) + sum_j y_j nu^{2j} (nu^2+1) log(nu^2+1) e^{-nu^2/2},
              Q = nu^2/(nu^2+1/4)^2, signed centres = {c_k} u {-c_k : c_k > 0}   (exactly as fp_moments.py)
  theta(nu) = (1/2)[erf((nu+nu_c)/a) - erf((nu-nu_c)/a)],  nu_c = Lam + 100, a = 10
  R_-(nu)   from prl/data/kernels.json (D_k2_EE)
  Leading-Regge kernel identity (I2, fp_norm.py):  sum_k P_k z^k = F(z),
      F(z) (1-z)^3 = (2/Lam) int_0^oo chi theta 2 Re[ R_-(nu) P_{-3/2+i nu}(X) ] dnu,   X = (1+z)/(1-z).
"""
import os, json, math
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..', 'prl', 'data'))


class Weight:
    def __init__(self, Dg=1000):
        d = json.load(open(os.path.join(DATA, 'k2', f'k2_continuum_D_{Dg}_rep.json')))
        self.Dg = Dg
        self.lam = (Dg - 2)*(Dg - 1)
        self.Lam = math.sqrt(self.lam)
        cs = d['centres_sigma']
        self.c = np.array([float(c) for c, s in cs]); self.s = np.array([float(s) for c, s in cs])
        self.xs = d['x_re']; self.x = np.array([float(v) for v in self.xs])
        self.ys = d['y_re']; self.y = [float(v) for v in self.ys]
        self.nu_c = self.Lam + 100; self.a = 10.0
        kj = json.load(open(os.path.join(DATA, 'kernels.json')))['D_k2_EE']
        self.Rm_expr = kj['R_minus']

    # ---------------- numpy (float64)
    def chi_np(self, nu):
        nu = np.asarray(nu, dtype=float)
        out = np.zeros_like(nu)
        Q = nu**2/(nu**2 + 0.25)**2
        for i in range(len(nu)):
            v = nu[i]
            m = np.abs(self.c - v) < 12*self.s
            B = np.sum(self.x[m]*np.exp(-(v - self.c[m])**2/(2*self.s[m]**2)))
            m2 = (np.abs(self.c + v) < 12*self.s) & (self.c > 0)
            B += np.sum(self.x[m2]*np.exp(-(v + self.c[m2])**2/(2*self.s[m2]**2)))
            out[i] = Q[i]*B
        lg = (nu**2 + 1)*np.log1p(nu**2)*np.exp(-nu**2/2)
        out += (self.y[0] + self.y[1]*nu**2 + self.y[2]*nu**4)*lg
        return out

    def theta_np(self, nu):
        from scipy.special import erf
        return 0.5*(erf((nu + self.nu_c)/self.a) - erf((nu - self.nu_c)/self.a))

    def Rm_np(self, nu):
        I = 1j
        n = np.asarray(nu, dtype=complex)
        return -2*(2*n - I)*(2*n + I)**2*(2*n + 3*I)*(2*n + 5*I)*(2*n + 11*I)/(n*(2*n - 7*I)*(2*n - 3*I)*(2*n + 7*I))

    # ---------------- mpmath
    def chi_mp(self, nu):
        nu = mp.mpf(nu)
        B = mp.mpf(0)
        for (c, s, x) in zip(self.c, self.s, self.xs):
            c = mp.mpf(float(c)); s = mp.mpf(float(s)); x = mp.mpf(float(x))
            if abs(nu - c) < 40*s:
                B += x*mp.exp(-(nu - c)**2/(2*s*s))
            if c > 0 and abs(nu + c) < 40*s:
                B += x*mp.exp(-(nu + c)**2/(2*s*s))
        Q = nu**2/(nu**2 + mp.mpf(1)/4)**2
        y = [mp.mpf(v) for v in self.ys]
        return Q*B + (y[0] + y[1]*nu**2 + y[2]*nu**4)*(nu**2 + 1)*mp.log(1 + nu**2)*mp.exp(-nu**2/2)

    def theta_mp(self, nu):
        nu = mp.mpf(nu); nc = mp.sqrt(self.lam) + 100
        return (mp.erf((nu + nc)/10) - mp.erf((nu - nc)/10))/2

    def Rm_mp(self, nu):
        nu = mp.mpf(nu); I = mp.mpc(0, 1)
        return -2*(2*nu - I)*(2*nu + I)**2*(2*nu + 3*I)*(2*nu + 5*I)*(2*nu + 11*I)/(nu*(2*nu - 7*I)*(2*nu - 3*I)*(2*nu + 7*I))
