"""rv5_build_chi.py -- double-precision tabulation of the NEW family chi = B q (note Sec.1), uncut, on Gauss-Legendre nodes,
for use as a discrete measure in the paper's moment formula (Eq. fp-moments).  Optional erf cutoff theta(L+2000,10) (note Sec.1).
Output: out/chi_<lam2>.npz with nu, wt (GL weights), chi, g = -9216/L (4nu^2+1)/((4nu^2+9)(4nu^2+49)) chi.
Usage: python3 rv5_build_chi.py LAM2 [EXTRA=150] [NGL=24]"""
import sys, math
import numpy as np
from scipy.special import loggamma, erfc
from rv_common import PSI2
import mpmath as mp

lam2 = int(sys.argv[1]); EXTRA = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0
NGL = int(sys.argv[3]) if len(sys.argv) > 3 else 24
L = math.sqrt(lam2)
from rv_common import X2, K0
cX = np.array([float(x) for x in X2]); k0 = float(K0)

def psi2(p):
    # factored form (1-p)^5 p X(p)/k0: no cancellation near p = 1
    return (1-p)**5*p*np.polynomial.polynomial.polyval(p, cX)/k0

def logF_real(v):
    """log F(v), v real (F>0): P(v) cosh(pi v/2) |Gamma(5/4+iv/2)|^2"""
    v = np.asarray(v, float)
    P = (v**2+0.25)**2*(v**2+2.25)*(v**2+6.25)*(v**2+20.25)
    a = np.pi*np.abs(v)/2
    lc = a + np.log1p(np.exp(-2*a)) - math.log(2)
    lg = 2*np.real(loggamma(1.25+0.5j*v))
    return np.log(P) + lc + lg

def B_real(nu):
    """B(nu) = F(nu+i)/R_-(nu) = [P(nu+i)/R_-(nu)] * i sinh(pi nu/2) * Gamma(3/4+i nu/2) Gamma(7/4-i nu/2)  (real, even)"""
    z = np.asarray(nu, float).astype(complex)
    I = 1j
    poly = -z*(2*z-7*I)**2*(2*z-3*I)**2*(2*z+3*I)*(2*z+7*I)**2/2048          # P(z+i)/R_-(z), the note's identity (checked in rv1)
    a = np.pi*np.abs(nu)/2
    lsh = a + np.log1p(-np.exp(-2*a)) - math.log(2)                         # log sinh(pi|nu|/2)
    lg = loggamma(0.75+0.5j*z) + loggamma(1.75-0.5j*z)
    val = poly*1j*np.sign(nu)*np.exp(lsh + lg)
    return val

# v-grid for the convolution: GL panels of width 0.25, 20 nodes, breakpoints 0 and +-L
tv, wv = np.polynomial.legendre.leggauss(20)
edges = np.concatenate([np.arange(0, L, 0.25), [L]])
V = np.concatenate([(a+b)/2+(b-a)/2*tv for a, b in zip(edges[:-1], edges[1:])])
WV = np.concatenate([(b-a)/2*wv for a, b in zip(edges[:-1], edges[1:])])
Wpos = psi2(V/L)*np.exp(-logF_real(V))            # W on (0,L); W even
# nu nodes
t, w = np.polynomial.legendre.leggauss(NGL)
nedges = list(np.arange(0, 1, 0.125)) + list(np.arange(1, 4, 0.25)) + list(np.arange(4, L+EXTRA, 0.5)) + [L+EXTRA]
if L not in nedges:
    nedges = sorted(set(nedges + [L]))
NU = np.concatenate([(a+b)/2+(b-a)/2*t for a, b in zip(nedges[:-1], nedges[1:])])
WT = np.concatenate([(b-a)/2*w for a, b in zip(nedges[:-1], nedges[1:])])

def q_of(nu):
    out = np.empty_like(nu)
    for i0 in range(0, len(nu), 256):
        x = nu[i0:i0+256][:, None]
        # (1/4) int_{-L}^{L} W(v) sech(pi(x-v)/2) dv = (1/4) int_0^L W(v)[sech(pi(x-v)/2)+sech(pi(x+v)/2)] dv
        k1 = 1/np.cosh(np.clip(np.pi*(x-V)/2, -700, 700)); k2 = 1/np.cosh(np.clip(np.pi*(x+V)/2, -700, 700))
        out[i0:i0+256] = 0.25*((k1+k2)*(Wpos*WV)).sum(1)
    return out

q = q_of(NU)
Bv = B_real(NU)
chi = (Bv*q)
print("max |Im B|/|B| =", float(np.max(np.abs(Bv.imag)/np.maximum(np.abs(Bv.real), 1e-300))))
chi = chi.real
theta = 1 - 0.5*erfc((NU + L + 2000)/10) - 0.5*erfc((L + 2000 - NU)/10)
g = -9216.0/L*(4*NU**2+1)/((4*NU**2+9)*(4*NU**2+49))*chi
np.savez(f'out/chi_{lam2}.npz', nu=NU, wt=WT, chi=chi, g=g, theta=theta, L=L)
print("sum wt*chi =", float((WT*chi).sum()))
