"""Referee r1 (SUPERCONVERGENT_UPLIFT.md): independent of laplace.py / s1 / s2.
Direct numerical Fourier-Laplace transform on H^2 (q=(1,0,0), l=(1,1,0), theta=0), as in the master formula of laplace.py docstring:
  M(h1,h3) = Gamma(a) int sinh(e) de dphi cosh(e)^{-a} (-yhat.l)^{-1/2-i nu} sum_q w_q t_q(yhat) / c(a,1/2+i nu).
Celestial weights from QUADRATURE of |sin psi|^{-delta} (j'=0) and sgn(sin psi)|sin psi|^{-delta} (j'=1), normalised to the 0-th coefficient.
Checks: (T1) N_5 at nu=1.3 (a=10); (T2) k=4 j'=0 pres/flip at nu=1 (a=9); (T3) D1^2 flip/pres at a=9, nu=10, 20."""
import numpy as np, mpmath as mp, sys
def dot(u, v): return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]
def perp(y, n): return np.cross(y, n)*np.array([-1.0, 1.0, 1.0])

def coef(M, nu, odd):
    dl = mp.mpf(1)/2 + 1j*mp.mpf(nu)
    def f(p):
        s = mp.sin(p)
        val = abs(s)**(-dl)
        if odd and s < 0: val = -val
        return val*mp.exp(-1j*M*p)
    return complex(mp.quad(f, [0, mp.pi/2, mp.pi, 3*mp.pi/2, 2*mp.pi])/(2*mp.pi))

def w_cel(jp, nu, L=2):
    c0 = coef(jp, nu, jp == 1)
    w = {}
    for q1 in range(-L, L+1):
        for q3 in range(-L, L+1):
            M = q1 + q3
            if M % 2 != jp: continue
            w[(q1, q3)] = (-1)**q3*coef(M, nu, jp == 1)/c0
    return w

def response(L, a, nu, w, h1, h3, Neta=200, Nu=256, etamax=18.0, npan=10):
    lv = np.array([1.0, 1.0, 0.0])
    eps = {+1: np.array([0, 1, 1j])/np.sqrt(2), -1: np.array([0, 1, -1j])/np.sqrt(2)}
    xg, wg = np.polynomial.legendre.leggauss(Neta)
    edges = np.linspace(0, etamax, npan+1)
    u = (np.arange(Nu) + 0.5)*2*np.pi/Nu - np.pi
    NT = 2*L + 4; thq = np.arange(NT)*2*np.pi/NT
    tot = 0j
    for A, B in zip(edges[:-1], edges[1:]):
        for e, ww in zip(0.5*(B-A)*xg + 0.5*(A+B), 0.5*(B-A)*wg):
            kap = np.exp(-e)
            phi = 2*np.arctan(kap*np.tan(u/2)); jac = kap*(1 + np.tan(u/2)**2)/(1 + kap**2*np.tan(u/2)**2)
            y = np.stack([np.cosh(e)*np.ones_like(phi), np.sinh(e)*np.cos(phi), np.sinh(e)*np.sin(phi)], -1)
            yl = 2*np.cosh(e)*np.sin(phi/2)**2 + np.exp(-e)*np.cos(phi)
            n = lv[None, :]/yl[:, None] - y; npp = perp(y, n)
            def cf(h):
                zt = y[:, None, :] + np.cos(thq)[None, :, None]*n[:, None, :] + np.sin(thq)[None, :, None]*npp[:, None, :]
                F = np.fft.fft(dot(zt, eps[h][None, None, :])**L, axis=1)/NT
                return {q: F[:, q % NT] for q in range(-L, L+1)}
            c1, c3 = cf(h1), cf(h3)
            Ft = sum(wq*c1[q[0]]*c3[q[1]] for q, wq in w.items())
            tot += ww*np.sinh(e)*np.cosh(e)**(-a)*np.sum(Ft*yl**(-0.5 - 1j*nu)*jac)*(2*np.pi/Nu)
    c = 2**(a-1)*np.sqrt(np.pi)*complex(mp.gamma((a-0.5-1j*nu)/2)*mp.gamma((a-0.5+1j*nu)/2))
    return complex(mp.gamma(a))*tot/c

def gam(aa, x): return float(abs(mp.gamma((aa - 0.5 + 1j*x)/2))**2)
WA = {(2, -2): 0.25j, (-2, 2): -0.25j}
WD = {(2, -2): 4.0, (-2, 2): 4.0, (0, 0): 2/9}
if __name__ == "__main__":
    mp.mp.dps = 20
    out = []
    def pr(s): print(s, flush=True); out.append(s)
    # T1
    nu, a = 1.3, 10
    w1 = w_cel(1, nu)
    sc = response(2, a, nu, w1, 1, -1); pj = response(2, a, nu, WA, 1, -1)
    scf = abs(response(2, a, nu, w1, 1, 1)); pjf = abs(response(2, a, nu, WA, 1, 1))
    Ncl = -(32/9)*(2*nu + 1j)*(4*nu**2 + 41)/(4*nu**2 + 393)
    g1 = gam(a-1, nu)/gam(a, nu)
    pr(f"T1 k=5 a=10 nu=1.3: proj A(+-)={pj:.8f} claimed {9j*g1*(4*nu**2+393)/(64*(4*nu**2+169)):.8f}; |flip| SC={scf:.1e} proj={pjf:.1e}")
    pr(f"   N5 direct = {sc/pj:.8f}  claimed {Ncl:.8f}  rel diff {abs(sc/pj/Ncl-1):.2e}")
    # T2
    nu, a = 1.0, 9
    w0 = w_cel(0, nu)
    pe = (response(2, a, nu, w0, 1, -1) + response(2, a, nu, w0, -1, 1))/2
    fe = (response(2, a, nu, w0, 1, 1) + response(2, a, nu, w0, -1, -1))/2
    pe_cl = 2*(16*nu**4 + 216*nu**2 + 263)/((4*nu**2 + 81)*(4*nu**2 + 169))
    x = nu
    fe_cl = 420*(2*x-5j)*(2*x-1j)/((2*x-13j)*(2*x-9j)*(2*x+3j)*(2*x+7j)*(2*x+9j)*(2*x+13j))
    pr(f"T2 k=4 a=9 nu=1 j'=0: pres {pe:.8f} claimed {pe_cl:.8f}; flip {fe:.8f} claimed {fe_cl:.8f}")
    # T3
    for nu, Nu in ((10.0, 1024), (20.0, 2048)):
        a = 9
        pD = (response(2, a, nu, WD, 1, -1, Nu=Nu, Neta=300) + response(2, a, nu, WD, -1, 1, Nu=Nu, Neta=300))/2
        fD = (response(2, a, nu, WD, 1, 1, Nu=Nu, Neta=300) + response(2, a, nu, WD, -1, -1, Nu=Nu, Neta=300))/2
        N1 = (4*nu**2 + 137)*(4*nu**2 + 561)
        cl = abs((2*nu-7j)*(2*nu-5j)*(2*nu-3j)*(2*nu-1j))/N1
        pr(f"T3 D1^2 a=9 nu={nu}: pres {pD:.6e} claimed {N1/(4*(4*nu**2+81)*(4*nu**2+169)):.6e}; |flip/pres| direct {abs(fD/pD):.5f} claimed {cl:.5f}")
    open(sys.argv[1], 'w').write("\n".join(out) + "\n")
