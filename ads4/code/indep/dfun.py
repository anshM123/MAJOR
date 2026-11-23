"""Symanzik integral I_Delta(X) = int_0^oo prod dt_i t_i^(Delta_i-1) exp(-sum_{i<j} t_i t_j X_ij), X_ij > 0,
evaluated WITHOUT Mellin-Barnes: t2, t4 integrated in closed form (Tricomi U), then the radial s-integral
in closed form (DLMF 13.10.7 -> 2F1), leaving a 1-d integral over c in (0,1).

  I = G2 G4 X24^-D4 * 1/2 * int_0^1 dc c^(D1-1)(1-c)^(D3-1) a^(D4-D2) lam^(-M/2) J(M/2, kappa/lam)
  a = c X12 + (1-c) X23, b = c X14 + (1-c) X34, lam = a b / X24, kappa = c(1-c) X13, M = D1+D3+D4-D2
  J(bb, z) = int_0^oo e^{-z w} w^{bb-1} U(D4, D4-D2+1, w) dw
           = G(bb) G(bb-cc+1)/G(D4+bb-cc+1) 2F1(bb, bb-cc+1; D4+bb-cc+1; 1-z)   (checked in selftest())
"""
import mpmath as mp

mp.mp.dps = 30

def J_closed(A, C, bb, z):
    # DLMF 13.10.7: int_0^oo e^{-z t} t^{b-1} U(a,c,t) dt = G(b)G(b-c+1)/G(a+b-c+1) 2F1(b, b-c+1; a+b-c+1; 1-z)
    return mp.gamma(bb) * mp.gamma(bb - C + 1) / mp.gamma(A + bb - C + 1) * mp.hyp2f1(bb, bb - C + 1, A + bb - C + 1, 1 - z)

def J_numeric(A, C, bb, z):
    return mp.quad(lambda w: mp.exp(-z * w) * w**(bb - 1) * mp.hyperu(A, C, w), [0, 1, 10, mp.inf])

def I_sym(D, X, dps=None):
    D1, D2, D3, D4 = [mp.mpf(d) for d in D]
    X12, X13, X14, X23, X24, X34 = [mp.mpf(X[k]) for k in ('12', '13', '14', '23', '24', '34')]
    M = D1 + D3 + D4 - D2
    A, C = D4, D4 - D2 + 1
    def f(c):
        if c < mp.mpf('1e-24') or 1 - c < mp.mpf('1e-24'):
            return mp.mpf(0)
        a = c * X12 + (1 - c) * X23
        b = c * X14 + (1 - c) * X34
        lam = a * b / X24
        kap = c * (1 - c) * X13
        return c**(D1 - 1) * (1 - c)**(D3 - 1) * a**(D4 - D2) * lam**(-M / 2) * J_closed(A, C, M / 2, kap / lam)
    val = mp.quad(f, [0, mp.mpf('1e-6'), mp.mpf('1e-3'), mp.mpf('0.5'), 1 - mp.mpf('1e-3'), 1 - mp.mpf('1e-6'), 1])
    return mp.gamma(D2) * mp.gamma(D4) * X24**(-D4) * val / 2

def I_simplex(D, X):
    """Feynman-parameter form: 1/2 G(S/2) int_simplex prod a^(D-1) Q(a)^(-S/2); 3-d, for checks only."""
    from scipy import integrate
    import math
    S = sum(D)
    Xs = {k: float(v) for k, v in X.items()}
    def f(a3, a2, a1):
        a4 = 1 - a1 - a2 - a3
        if a4 <= 0: return 0.0
        a = (a1, a2, a3, a4)
        Q = (a1*a2*Xs['12'] + a1*a3*Xs['13'] + a1*a4*Xs['14'] + a2*a3*Xs['23'] + a2*a4*Xs['24'] + a3*a4*Xs['34'])
        return math.prod(a[i]**(D[i]-1) for i in range(4)) * Q**(-S/2)
    v, e = integrate.nquad(f, [lambda a2, a1: [0, 1-a1-a2], lambda a1: [0, 1-a1], [0, 1]],
                           opts={'epsabs': 1e-13, 'epsrel': 1e-11, 'limit': 200})
    return 0.5 * math.gamma(S/2) * v

def selftest(log=print):
    # (1) DLMF form of J
    for (A, C, bb, z) in [(3, 1, 4.5, 0.3), (4, 2, 5, 1.7), (3, 0, 3.5, 0.05), (5, 3, 6, 2.5)]:
        a, b = J_closed(A, C, bb, z), J_numeric(A, C, bb, z)
        log(f'  J check A={A} C={C} b={bb} z={z}: closed {mp.nstr(a,15)} numeric {mp.nstr(b,15)} rel {mp.nstr(abs(a/b-1),3)}')
    # (2) full I vs Feynman-parameter simplex integral
    X = {'12': 1.3, '13': 0.7, '14': 0.9, '23': 1.1, '24': 2.0, '34': 1.6}
    for D in [(1, 1, 1, 1), (3, 3, 3, 3), (4, 3, 5, 3), (3, 4, 3, 5)]:
        a = I_sym(D, X); b = I_simplex(D, X)
        log(f'  I check D={D}: 1-d {mp.nstr(a,15)}  simplex {b:.15g}  rel {float(abs(a/b-1)):.2e}')
    # (3) box function: I_1111 X13 X24 = Phi(z, zbar)
    X = {'12': 1.0, '13': 1.0, '14': 1.0, '23': 1.0, '24': 1.0, '34': 1.0}
    u, v = mp.mpf(1), mp.mpf(1)
    # z, zbar from u = z zb, v = (1-z)(1-zb): z = e^{i pi/3}
    z = mp.exp(1j * mp.pi / 3); zb = mp.conj(z)
    Phi = (2 * mp.polylog(2, z) - 2 * mp.polylog(2, zb) + mp.log(z * zb) * mp.log((1 - z) / (1 - zb))) / (z - zb)
    a = I_sym((1, 1, 1, 1), X)
    log(f'  box check: I_1111 X13 X24 = {mp.nstr(a,15)}, Phi(e^(i pi/3)) = {mp.nstr(Phi.real,15)} (imag {mp.nstr(Phi.imag,2)})')

if __name__ == '__main__':
    selftest()

# ---- fast fixed-node version (double-exponential nodes on (0,1)), same formula
def _de_nodes(h, N):
    nodes = []
    for k in range(-N, N + 1):
        t = k * h
        u = mp.pi / 2 * mp.sinh(t)
        c = 1 / (1 + mp.exp(-2 * u))          # c = (1 + tanh u)/2
        w = h * mp.pi / 2 * mp.cosh(t) * 2 * c * (1 - c)   # dc/dt
        if c <= mp.mpf('1e-15') or 1 - c <= mp.mpf('1e-15'):
            continue
        nodes.append((c, w))
    return nodes
_NODES = {}

def I_fast(D, X, h=mp.mpf(1) / 12, N=48):
    key = (h, N, mp.mp.dps)
    if key not in _NODES:
        _NODES[key] = _de_nodes(h, N)
    D1, D2, D3, D4 = [mp.mpf(d) for d in D]
    X12, X13, X14, X23, X24, X34 = [mp.mpf(X[k]) for k in ('12', '13', '14', '23', '24', '34')]
    M = D1 + D3 + D4 - D2
    A, C = D4, D4 - D2 + 1
    bb = M / 2
    g = mp.gamma(bb) * mp.gamma(bb - C + 1) / mp.gamma(A + bb - C + 1)
    p1, p2, p3 = bb, bb - C + 1, A + bb - C + 1
    s = mp.mpf(0)
    for c, w in _NODES[key]:
        a = c * X12 + (1 - c) * X23
        b = c * X14 + (1 - c) * X34
        lam = a * b / X24
        kap = c * (1 - c) * X13
        s += w * c**(D1 - 1) * (1 - c)**(D3 - 1) * a**(D4 - D2) * lam**(-bb) * mp.hyp2f1(p1, p2, p3, 1 - kap / lam)
    return mp.gamma(D2) * mp.gamma(D4) * X24**(-D4) * g * s / 2
