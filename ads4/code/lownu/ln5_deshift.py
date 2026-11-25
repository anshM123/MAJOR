"""Reproduction of the reviewer's analytic low-nu de-shift weight for the k = 2 EE piece (floating point, [NE]).
chi(z) = chi_sing(z) + P(z) sum_j c_j G_j(z),  G_j = Gaussian pair at +-mu_j (width s),  P = z^2(z^2+1)/(z^2+1/4)^q  (q = 2 or 1),
E(nu) = 2 Re[R_-(nu-i) chi(nu-i)]  must reproduce the target T(nu) = psi2_NF02(nu/L)  (dilated NF02, support [0, L]).
Reports: integrated relative error, max|c_j|, S-transform relative errors S[E-T](beta)/S[T](beta), the slope check, real-axis
concentration of chi, and the tail integral int_{nu_t}^inf |chi| (|R_-|+|R_+|) e^{pi nu} (log10).
Usage: python3 ln5_deshift.py L q [spacing] [s] [ridge]"""
import sys, json, numpy as np
from ln_common import nf02_float, guard
L = float(sys.argv[1]); q = int(sys.argv[2])
dmu = float(sys.argv[3]) if len(sys.argv) > 3 else 0.27
s = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
lam = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-10
p2, p3 = nf02_float(); k0 = 27.48509817
I = 1j
def Rm(z):
    return -2*(2*z-I)*(2*z+I)**2*(2*z+3*I)*(2*z+5*I)*(2*z+11*I)/(z*(2*z-7*I)*(2*z-3*I)*(2*z+7*I))
def T(nu):
    p = np.clip(nu / L, 0, 1); return np.where(nu < L, np.polyval(p2[::-1], p) / k0, 0.0)
kappa = p2[1] / k0 / L
def chi_sing(z):
    w = z*z + 1; c = -25*kappa/(36*np.pi*np.sqrt(np.e))
    return c*w*np.log(w)*(1 + 37/36*w)*np.exp(-z*z/2)
def P(z): return z*z*(z*z+1)/(z*z+0.25)**q
mus = np.arange(0, L + 6 + 1e-9, dmu)
def Gm(z):  # (len(z), nmu)
    return np.exp(-(z[:, None]-mus[None, :])**2/(2*s*s)) + np.exp(-(z[:, None]+mus[None, :])**2/(2*s*s))
nu = np.linspace(0, L + 25, int((L + 25) / 0.01) + 1); z = nu - I
Es = 2*np.real(Rm(z[1:]) * chi_sing(z[1:])); Es = np.concatenate([[0.0], Es])
Bm = 2*np.real((Rm(z)*P(z))[:, None] * Gm(z)); Bm[0] = 0.0
wq = np.full(len(nu), 0.01); wq[0] = wq[-1] = 0.005
# weighted ridge LS: minimise int (E - T)^2 + lam |c|^2
Aw = Bm * np.sqrt(wq)[:, None]; bw = (T(nu) - Es) * np.sqrt(wq)
cc = np.linalg.solve(Aw.T @ Aw + lam * np.eye(len(mus)), Aw.T @ bw)
E = Es + Bm @ cc; R = E - T(nu)
relint = np.sum(wq * np.abs(R)) / np.sum(wq * np.abs(T(nu)))
# slope at 0 (finite difference) vs kappa
sl = (E[1] - E[0]) / 0.01
# S transform via Mehler-Dirichlet: S[g](b) = (sqrt2/pi) int_0^b ghat(t)/sqrt(cosh b - cosh t) dt, t = b - y^2
def ghat(g, t):
    return (np.cos(np.outer(t, nu)) * (wq * g)[None, :]).sum(axis=1)
def S(g, b, ny=400):
    y, wy = np.polynomial.legendre.leggauss(ny); y = (y + 1) / 2 * np.sqrt(b); wy = wy / 2 * np.sqrt(b)
    t = b - y*y; jac = 2*y
    den = np.sqrt(np.maximum(np.cosh(b) - np.cosh(t), 1e-300))
    return np.sqrt(2)/np.pi * np.sum(wy * jac * ghat(g, t) / den)
rows = []
for b in (0.1, 0.5, 1, 2, 4, 8, 16, 32):
    sT = S(T(nu), b); sR = S(R, b); rows.append((b, sT, sR, abs(sR / sT)))
# real-axis concentration of chi and the e^{pi nu} tail
nr = np.linspace(0, 200, 200001); zr = nr.astype(complex)
chir = np.real(chi_sing(zr) + P(zr) * (Gm(zr) @ cc))
Rr = np.abs(Rm(zr + 1e-30)) + np.abs(Rm(-zr + 1e-30))
def tail_log10(nt):
    m = nr >= nt
    with np.errstate(over='ignore', divide='ignore'):
        lg = np.log(np.abs(chir[m]) * Rr[m] + 1e-320) + np.pi * nr[m]
    mx = lg.max(); return float((mx + np.log(np.sum(np.exp(lg - mx)) * (nr[1]-nr[0]))) / np.log(10))
res = dict(L=L, q=q, dmu=dmu, s=s, ridge=lam, n=len(mus), relint=float(relint), maxc=float(np.abs(cc).max()),
           slope=float(sl), kappa=float(kappa), Srel=[(b, float(r)) for b, _, _, r in rows],
           chi_max=float(np.abs(chir).max()), chi_at=[(x, float(np.abs(chir[int(x*1000)]))) for x in (L, L+5, L+10, L+20, L+30)],
           tail_log10={str(t): tail_log10(t) for t in (L+10, L+20, L+30, 1.6*L)})
print(json.dumps(res))
json.dump(res, open(guard('cache/ln5_L%.2f_q%d_s%.2f.json' % (L, q, s)), 'w'))
