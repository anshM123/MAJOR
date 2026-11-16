"""Referee check r1 (OMEGA_PRIME.md section 6): N(nu) recomputed WITHOUT laplace.py / clsd_kernel.py.

Direct numerical Fourier-Laplace transform on H^2 (q = (1,0,0), l = (1,1,0), so -qhat.l = 1, theta = 0):
   M_K(h1,h3) = Gamma(a) * int sinh(eta) d eta d phi cosh(eta)^{-a} (-yhat.l)^{-1/2-i nu} F_K[t](yhat) / c(a, 1/2+i nu)
with t = (z1.eps_h1)^l (z3.eps_h3)^l, z(th) = yhat + cos th n_l + sin th n_l^perp (frame at yhat pointing to l),
F_K[t] = sum_q w_q t_q  (t_q = coefficient of e^{i(q1 th1 + q3 th3)}).
CLS-D odd kernel: w_q = (-1)^{q3} d_{(q1+q3-1)/2}/d_0, d_m Fourier coefficients of sgn(sin psi)|sin psi|^{-delta}
computed here by quadrature (not from the Pochhammer closed form).  Project kernels: D5 (l=1), A = D5 D1 (l=2).
"""
import numpy as np, mpmath as mp, sys

ETA = np.diag([-1.0, 1.0, 1.0])
def dot(u, v):  # u, v: (...,3)
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]

def levi_perp(y, n):
    # n_perp^mu = eta^{mu mu} eps_{mu a b} y^a n^b  (eps_{012} = 1)
    c = np.cross(y, n)               # Euclidean cross = eps_{mu a b} y^a n^b with lower mu
    return c*np.array([-1.0, 1.0, 1.0])

def dm_numeric(m, nu):
    dl = mp.mpf(1)/2 + 1j*mp.mpf(nu)
    f = lambda p: mp.sin(p)**(-dl)*mp.exp(-1j*(2*m+1)*p)
    # sgn(sin)|sin|^{-delta} over (0,2pi): (0,pi) part plus (pi,2pi) part
    I1 = mp.quad(f, [0, mp.pi/2, mp.pi])
    g = lambda p: -(mp.sin(p - mp.pi))**(-dl)*mp.exp(-1j*(2*m+1)*p)
    I2 = mp.quad(g, [mp.pi, 3*mp.pi/2, 2*mp.pi])
    return complex((I1 + I2)/(2*mp.pi))

def weights_clsd(l, nu):
    d0 = dm_numeric(0, nu)
    w = {}
    for q1 in range(-l, l+1):
        for q3 in range(-l, l+1):
            M = q1 + q3
            if M % 2 == 0: continue
            w[(q1, q3)] = (-1)**q3*dm_numeric((M-1)//2, nu)/d0
    return w

def weights_proj(l):
    if l == 1: return {(-1, 1): -0.5j, (1, -1): 0.5j}
    return {(-2, 2): -0.25j, (2, -2): 0.25j}

def response(l, a, nu, w, h1, h3, Neta=260, Nu=256, etamax=18.0):
    lv = np.array([1.0, 1.0, 0.0])
    eps = {+1: np.array([0, 1, 1j])/np.sqrt(2), -1: np.array([0, 1, -1j])/np.sqrt(2)}
    xg, wg = np.polynomial.legendre.leggauss(Neta)
    # split eta in panels for accuracy
    edges = np.linspace(0, etamax, 9)
    tot = 0j
    u = (np.arange(Nu) + 0.5)*2*np.pi/Nu - np.pi
    thq = np.arange(8)*2*np.pi/8
    for A, B in zip(edges[:-1], edges[1:]):
        et = 0.5*(B-A)*xg + 0.5*(A+B); we = 0.5*(B-A)*wg
        for e, ww in zip(et, we):
            kap = np.exp(-e)
            phi = 2*np.arctan(kap*np.tan(u/2))
            jac = kap*(1 + np.tan(u/2)**2)/(1 + kap**2*np.tan(u/2)**2)
            y = np.stack([np.cosh(e)*np.ones_like(phi), np.sinh(e)*np.cos(phi), np.sinh(e)*np.sin(phi)], -1)
            yl = 2*np.cosh(e)*np.sin(phi/2)**2 + np.exp(-e)*np.cos(phi)   # -yhat.l > 0, stable form
            n = lv[None, :]/yl[:, None] - y
            npp = levi_perp(y, n)
            # Fourier coefficients in theta of (z(th).eps)^l
            def coeffs(h):
                ev = eps[h]
                zt = y[:, None, :] + np.cos(thq)[None, :, None]*n[:, None, :] + np.sin(thq)[None, :, None]*npp[:, None, :]
                val = dot(zt, ev[None, None, :])**l       # (Nphi, 8)
                F = np.fft.fft(val, axis=1)/8              # coefficient of e^{-i k th}? fft gives sum val e^{-2pi i jk/8}
                # val(th) = sum_q c_q e^{i q th}  ->  c_q = (1/8) sum_j val_j e^{-i q th_j} = F[q mod 8]
                return {q: F[:, q % 8] for q in range(-l, l+1)}
            c1 = coeffs(h1); c3 = coeffs(h3)
            Ft = sum(wq*c1[q[0]]*c3[q[1]] for q, wq in w.items())
            integ = Ft*yl**(-0.5 - 1j*nu)
            tot += ww*np.sinh(e)*np.cosh(e)**(-a)*np.sum(integ*jac)*(2*np.pi/Nu)
    c = 2**(a-1)*np.sqrt(np.pi)*complex(mp.gamma((a-0.5-1j*nu)/2)*mp.gamma((a-0.5+1j*nu)/2))
    return complex(mp.gamma(a))*tot/c

def gam(aa, x): return float(abs(mp.gamma((aa - 0.5 + 1j*x)/2))**2)

if __name__ == "__main__":
    out = []
    def pr(*s):
        line = " ".join(str(x) for x in s); print(line, flush=True); out.append(line)
    for l, a in ((1, 6), (2, 8)):
        for nu in (0.7, 2.3):
            g1 = gam(a-1, nu)/gam(a, nu)
            wc = weights_clsd(l, nu); wp = weights_proj(l)
            R = {}
            for nm, w in (('CLSD', wc), ('proj', wp)):
                for hh in ((1, -1), (-1, 1), (1, 1), (-1, -1)):
                    R[(nm,)+hh] = response(l, a, nu, w, *hh)
            Nnum = R[('CLSD', 1, -1)]/R[('proj', 1, -1)]
            if l == 1:
                Nclaim = 0.4*(2*nu + 1j)
                pc = 5j*g1/8; cc = g1*(0.5j*nu - 0.25)
            else:
                Nclaim = -(32/7)*(2*nu + 1j)*(4*nu**2 + 9)/(4*nu**2 + 201)
                pc = 7j*g1*(4*nu**2 + 201)/(64*(4*nu**2 + 81))
                cc = g1*(1/648)*(-32j*nu**3 + 16*nu**2 - 72j*nu + 36)/(nu**2/20.25 + 1)
            pr(f"l={l} a={a} nu={nu}:")
            pr(f"   proj A(+-)={R[('proj',1,-1)]:.8f}  claimed {pc:.8f};  A(-+)={R[('proj',-1,1)]:.8f}  A(++)={abs(R[('proj',1,1)]):.2e} A(--)={abs(R[('proj',-1,-1)]):.2e}")
            pr(f"   CLSD A(+-)={R[('CLSD',1,-1)]:.8f}  claimed {cc:.8f};  A(-+)={R[('CLSD',-1,1)]:.8f}  A(++)={abs(R[('CLSD',1,1)]):.2e} A(--)={abs(R[('CLSD',-1,-1)]):.2e}")
            pr(f"   N numeric = {Nnum:.8f}   claimed N = {Nclaim:.8f}   rel.diff = {abs(Nnum/Nclaim-1):.2e}")
    open(sys.argv[1] if len(sys.argv) > 1 else 'cache/r1.log', 'w').write("\n".join(out) + "\n")
