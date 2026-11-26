"""
fullpos4/fp4_kern.py (FULLPOS4.md §1): per-frequency kernel K_p(nu) = L(F_{p,nu}) three ways [DV, mpmath].
  exact:  K_p = sum_{k<=n} U_k g_k,  U_k = (-n)_k(n+s)_k/((-N)_k(N+s)_k),  (k+1)^2 g_{k+1} = (2k^2+2pk+p-s_nu) g_k - (k-1+p)^2 g_{k-1}
  s59:    -(iC/2pi) int_{-pi/2-i eps}^{pi/2-i eps} e^{iB al}/sin(al) int_0^oo e^{-a v} xi^{s-1} F(xi^2) dv d al,
          xi = sinh(v/2)/sinh(v/2 + i al)   (ANGULAR_ROUTE §59; full conical function)
  angles: (al, ga) with gamma = al - i v;  xi = sin((ga-al)/2)/sin((ga+al)/2), e^beta = tan(ga/2)/tan(al/2)
Notation: n = J-b, N = Delta-b-2, a = N-n, s = 2b+1, B = 2J+1, S = B+a, C = a Gamma(n+1)Gamma(N+s)/(Gamma(n+s)Gamma(N+1)).
"""
import sys
sys.dont_write_bytecode = True
import mpmath as mp


def exact(p, nu, n, N, s):
    nu = mp.mpmathify(nu); sn = nu*nu + mp.mpf(1)/4
    g = [mp.mpf(1)]; gm = mp.mpf(0)
    tot = mp.mpf(0); U = mp.mpf(1)
    for k in range(n + 1):
        tot += U*g[k]
        g.append(((2*k*k + 2*p*k + p - sn)*g[k] - (k - 1 + p)**2*(g[k - 1] if k >= 1 else 0))/(k + 1)**2)
        U = U*(k - n)*(n + s + k)/((k - N)*(N + s + k))
    return tot


def Cconst(n, N, s):
    a = N - n
    return a*mp.exp(mp.loggamma(n + 1) + mp.loggamma(N + s) - mp.loggamma(n + s) - mp.loggamma(N + 1))


def Ffull(p, nu, xi):
    z = xi*xi
    X = (1 + z)/(1 - z)
    return (1 - z)**(-p)*mp.legenp(-mp.mpf(1)/2 + 1j*nu, 0, X, type=3)


def Pplus(nu, beta):
    """outgoing Legendre solution: P = P+ + P- (P-(beta) = conj P+(conj beta) for real nu)."""
    nu = mp.mpmathify(nu)
    pre = mp.gamma(1j*nu)/(mp.sqrt(mp.pi)*mp.gamma(mp.mpf(1)/2 + 1j*nu))
    return pre*mp.exp((-mp.mpf(1)/2 + 1j*nu)*beta)*mp.hyp2f1(mp.mpf(1)/2, mp.mpf(1)/2 - 1j*nu, 1 - 1j*nu, mp.exp(-2*beta))


def s59(p, nu, n, N, s, eps=None):
    a = N - n; J = n + (s - 1)//2; B = 2*J + 1
    S = B + a
    eps = eps or mp.mpf(1)/S
    C = Cconst(n, N, s)

    def inner(al):
        def f(v):
            xi = mp.sinh(v/2)/mp.sinh(v/2 + 1j*al)
            return mp.exp(-a*v)*xi**(s - 1)*Ffull(p, nu, xi)
        return mp.quad(f, [0, 0.1, 1, 5, mp.inf])

    def outer(x):
        al = x - 1j*eps
        return mp.exp(1j*B*al)/mp.sin(al)*inner(al)
    pts = [-mp.pi/2, -0.3, -10*eps, -eps, 0, eps, 10*eps, 0.3, mp.pi/2]
    I = mp.quad(outer, pts)
    return -1j*C/(2*mp.pi)*I


if __name__ == '__main__':
    mp.mp.dps = 20
    for (p, nu, J, b, Delta) in [(3, 1.2, 10, 0, 32.5), (3, 1.7, 3, 0, 9.3)]:
        n = J - b; N = mp.mpf(Delta) - b - 2; s = 2*b + 1
        print(p, nu, J, b, Delta, 'exact', exact(p, nu, n, N, s), 's59', s59(p, nu, n, N, s), flush=True)
