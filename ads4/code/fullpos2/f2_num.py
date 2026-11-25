"""
fullpos2/f2_num.py (FULLPOS2.md): [DV] numerics for the finite-mass frequency kernels.

  g^{(p)}_k(nu) = [z^k] (1-z)^{-p} P_{-1/2+i nu}((1+z)/(1-z)),   recurrence (proved in FULLPOS2.md §1):
      (k+1)^2 g_{k+1} = (2k^2 + 2pk + p - s) g_k - (k+p-1)^2 g_{k-1},  s = nu^2 + 1/4,  g_0 = 1, g_{-1} = 0.
  Phi2(nu) = w^3 sum_{k=0}^{n+1} V_{b,k} g^{(3)}_k(nu)          (physical k = 2 de-shifted kernel; w = 1 - r, r = C/T)
      V_{b,k} = U_k [1 + ((k+3)^2/2 - u_{k-1})/T] + (k^2/(2T)) U_{k-1},  u_{-1} = (b-1) b   (ANGULAR_ROUTE §32)
  Phi3(nu) = sum_k U_k g^{(4)}_k(nu)                                 (bare k = 3 kernel; H3 = sigma_b Lam^2 w^3/T int psi3 Phi3(Lam p) dp)
  Hybrid (T0): Phi2_hyb(nu) = D_b(1 - 2nu^2/m^2) + P_{-1/2+i nu}(cosh beta) - J0((2J+1) nu/m),  sinh(beta/2) = (J+1/2)/m.
"""
import sys, os, json
sys.dont_write_bytecode = True
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
FP = os.path.join(HERE, '..', 'fullpos'); sys.path.insert(0, FP)


def state(J, b, T):
    T = mp.mpf(T); C = J*(J + 1); n = J - b
    U = [mp.mpf(1)]
    for k in range(n + 1):
        u = (b + k)*(b + k + 1); U.append(U[-1]*(C - u)/(T - u))   # U[n+1] = 0
    return dict(J=J, b=b, T=T, C=C, n=n, U=U, r=C/T, w=1 - C/T, m2=T - C)


def gcoef(p, nu, K):
    s = mp.mpf(nu)**2 + mp.mpf(1)/4
    g = [mp.mpf(1)]; gm = mp.mpf(0)
    for k in range(K):
        nx = ((2*k*k + 2*p*k + p - s)*g[k] - (k + p - 1)**2*(g[k - 1] if k else gm))/(k + 1)**2
        g.append(nx)
    return g


def Vb(st):
    b, T, U, n = st['b'], st['T'], st['U'], st['n']
    V = []
    for k in range(n + 2):
        um = (b - 1)*b if k == 0 else (b + k - 1)*(b + k)
        v = U[k]*(1 + (mp.mpf(k + 3)**2/2 - um)/T)
        if k >= 1:
            v += mp.mpf(k*k)/(2*T)*U[k - 1]
        V.append(v)
    return V


def Phi2(st, nu, V=None):
    V = V or Vb(st)
    g = gcoef(3, nu, st['n'] + 2)
    return st['w']**3*mp.fsum(V[k]*g[k] for k in range(st['n'] + 2))


def Phi3(st, nu):
    g = gcoef(4, nu, st['n'] + 1)
    return mp.fsum(st['U'][k]*g[k] for k in range(st['n'] + 1))


def Db(st, y):
    """2F1(-n, n+2b+1; 1; y)"""
    return mp.hyp2f1(-st['n'], st['n'] + 2*st['b'] + 1, 1, y)


def Phi2_hyb(st, nu):
    J = st['J']; m = mp.sqrt(st['m2']); nu = mp.mpf(nu)
    beta = 2*mp.asinh((J + mp.mpf(1)/2)/m)
    return Db(st, nu**2/m**2) + mp.legenp(-mp.mpf(1)/2 + 1j*nu, 0, mp.cosh(beta), type=3).real - mp.besselj(0, (2*J + 1)*nu/m)


_KJ = None
def Kshift(Phi, nu):
    """K[Phi](nu) = -(4x+9)(4x+49)/(4608(4x+1)) Re[R_-(nu) Phi(nu+i) + R_+(nu) Phi(nu-i)]"""
    global _KJ
    if _KJ is None:
        _KJ = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'kernels.json')))['D_k2_EE']
    nu = mp.mpf(nu); x = nu*nu
    ev = lambda e: mp.mpmathify(eval(e, {'nu': nu, 'I': mp.mpc(0, 1)}))
    return mp.re(-(4*x + 9)*(4*x + 49)/(4608*(4*x + 1))*(ev(_KJ['R_minus'])*Phi(nu + 1j) + ev(_KJ['R_plus'])*Phi(nu - 1j)))


def gcoef_c(p, nu, K):
    """complex nu"""
    s = mp.mpmathify(nu)**2 + mp.mpf(1)/4
    g = [mp.mpc(1)]; gm = 0
    for k in range(K):
        g.append(((2*k*k + 2*p*k + p - s)*g[k] - (k + p - 1)**2*(g[k - 1] if k else gm))/(k + 1)**2)
    return g


if __name__ == '__main__':
    import fp_kernel as FK
    for (J, Delta, b) in ((12, 40, 0), (30, 1030, 4), (100, 1100, 0)):
        T = (Delta - 1)*(Delta - 2)
        mp.mp.dps = 80
        st = state(J, b, T); V = Vb(st)
        for nu in (0.7, 5.0, 23.0):
            mp.mp.dps = int(0.75*nu + 60 + 0.1*J)
            st = state(J, b, T); V = Vb(st)
            ph = lambda zz: st['w']**3*mp.fsum(V[k]*g for k, g in enumerate(gcoef_c(3, zz, st["n"] + 1)))
            a = Kshift(ph, nu); kf = FK.kfin(nu, J, T, b, mp.mp.dps)[0]
            print(J, Delta, b, nu, mp.nstr(a, 20), mp.nstr(kf, 20), mp.nstr(a - kf, 3))
