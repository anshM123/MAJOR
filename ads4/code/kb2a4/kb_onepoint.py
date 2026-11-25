"""kb_onepoint.py -- one-point (Hofman-Maldacena) matrices of E(n) on the T family and the traces entering the k = 2 first-moment
terms (ads4/K_B2A4.md Sec. 4.2).  Mostly-minus metric, P in the forward cone, n-bar = (1, cos t, sin t).
  H_E(P,n) = omega Pi_TT(P),   H_W(P,n) = omega (|q><q| - Pi/4),  q = Pi_TT[N N],  N = unit rest-frame direction of n-bar,
  omega = (P^2)^{3/2}/(n-bar.P)^2;   H = H_E + t4 H_W  (flip entries t4/4: |flip|/diag = 1 at |t4| = 4, as C1_FINAL Sec. 2.3(e)).
Computes at the rest frame P = (1,0,0) the traces T31 = Tr(H3 H1), tau1 = Tr(H3 nbar1.d H1), tau3 = Tr((nbar3.d H3) H1) as trigonometric
polynomials in phi = phi3 - phi1 (DFT of exact-arithmetic evaluations, derivatives by central differences at 40 digits), and the
rest-frame helicity matrix of H_W (normalisation / HM phase check).
Usage: python3 -B kb_onepoint.py -> cache/kb_onepoint.txt"""
import mpmath as mp
mp.mp.dps = 40
G = [1, -1, -1]
def dot(a, b): return sum(G[i]*a[i]*b[i] for i in range(3))
def nbar(t): return [mp.mpf(1), mp.cos(t), mp.sin(t)]
def proj_data(P):
    P2 = dot(P, P)
    Pl = [G[i]*P[i] for i in range(3)]
    th_ud = [[(1 if a == b else 0) - P[a]*Pl[b]/P2 for b in range(3)] for a in range(3)]
    th_uu = [[(G[a] if a == b else 0) - P[a]*P[b]/P2 for b in range(3)] for a in range(3)]
    th_dd = [[(G[a] if a == b else 0) - Pl[a]*Pl[b]/P2 for b in range(3)] for a in range(3)]
    return P2, th_ud, th_uu, th_dd
def Pi(T, pd):
    P2, th_ud, th_uu, th_dd = pd
    A = [[sum(th_ud[a][c]*T[c][d]*th_ud[b][d] for c in range(3) for d in range(3)) for b in range(3)] for a in range(3)]
    tr = sum(th_dd[a][b]*T[a][b] for a in range(3) for b in range(3))
    return [[(A[a][b] + A[b][a])/2 - th_uu[a][b]*tr/2 for b in range(3)] for a in range(3)]
def qvec(P, t, pd):
    P2 = pd[0]
    n = nbar(t); nP = dot(n, P)
    N = [(n[i] - nP*P[i]/P2)*mp.sqrt(P2)/nP for i in range(3)]
    return Pi([[N[a]*N[b] for b in range(3)] for a in range(3)], pd)
def contract(A, B):  # A^{ab} B^{cd} g_ac g_bd
    return sum(G[a]*G[b]*A[a][b]*B[a][b] for a in range(3) for b in range(3))
BASIS = []
for a in range(3):
    for b in range(a, 3):
        E = [[mp.mpf(0)]*3 for _ in range(3)]; E[a][b] = mp.mpf(1); E[b][a] = mp.mpf(1)
        BASIS.append((a, b, E))
def Wmap(P, t, kind):
    """returns function e -> H e  (omega omitted), kind 'E' or 'W'"""
    pd = proj_data(P)
    if kind == 'E':
        return lambda e: Pi(e, pd)
    q = qvec(P, t, pd)
    def f(e):
        c = contract(q, e)
        pe = Pi(e, pd)
        return [[q[a][b]*c - pe[a][b]/4 for b in range(3)] for a in range(3)]
    return f
def omega(P, t):
    return dot(P, P)**mp.mpf(1.5)/dot(nbar(t), P)**2
def Hfull(P, t, kind):
    w = omega(P, t); f = Wmap(P, t, kind)
    return lambda e: [[w*x for x in row] for row in f(e)]
def trace(fmap):
    return sum(fmap(E)[a][b] for a, b, E in BASIS)
def compose(A, B): return lambda e: A(B(e))
REST = [mp.mpf(1), mp.mpf(0), mp.mpf(0)]
H_STEP = mp.mpf('1e-12')
def dmap(t, kind, direction):
    """directional derivative of H(P,t) along direction at rest frame, as map"""
    Pp = [REST[i] + H_STEP*direction[i] for i in range(3)]
    Pm = [REST[i] - H_STEP*direction[i] for i in range(3)]
    Hp = Hfull(Pp, t, kind); Hm = Hfull(Pm, t, kind)
    return lambda e: [[(x - y)/(2*H_STEP) for x, y in zip(rp, rm)] for rp, rm in zip(Hp(e), Hm(e))]
def data(phi, kind):
    t1 = mp.mpf(0); t3 = mp.mpf(phi)
    H1 = Hfull(REST, t1, kind); H3 = Hfull(REST, t3, kind)
    T31 = trace(compose(H3, H1))
    tau1 = trace(compose(H3, dmap(t1, kind, nbar(t1))))
    tau3 = trace(compose(dmap(t3, kind, nbar(t3)), H1))
    return T31, tau1, tau3
def dft(fun, K=10, M=48):
    vals = [fun(2*mp.pi*j/M) for j in range(M)]
    out = []
    for k in range(K+1):
        c = sum(v*mp.cos(2*mp.pi*j*k/M) for j, v in enumerate(vals))*2/M
        s = sum(v*mp.sin(2*mp.pi*j*k/M) for j, v in enumerate(vals))*2/M
        if k == 0: c /= 2
        out.append((k, c, s))
    return out
def helicity_matrix(t):
    """rest-frame matrix <eps_a|(Q - Pi/4)|eps_b>, eps_pm = e_pm e_pm, e_pm = (0,1,pm i)/sqrt2 (spatial, Euclidean pairing)."""
    e = {1: [0, 1/mp.sqrt(2), 1j/mp.sqrt(2)], -1: [0, 1/mp.sqrt(2), -1j/mp.sqrt(2)]}
    eps = {h: [[e[h][a]*e[h][b] for b in range(3)] for a in range(3)] for h in (1, -1)}
    f = Wmap(REST, mp.mpf(t), 'W')
    M = {}
    for a in (1, -1):
        for b in (1, -1):
            img = f(eps[b])
            M[(a, b)] = sum(mp.conj(eps[a][i][j])*img[i][j] for i in range(3) for j in range(3))
    return M


if __name__ == "__main__":
    import os
    out = []
    for kind in ('E', 'W'):
        for idx, name in ((0, 'T31'), (1, 'tau1'), (2, 'tau3')):
            co = dft(lambda ph: data(ph, kind)[idx], K=8, M=24)
            line = '%s %s: nonzero DFT coefficients (k, cos, sin) = %s' % (kind, name, [(k, mp.nstr(c, 15), mp.nstr(s, 3)) for k, c, s in co
                                                                                   if abs(c) > 1e-25 or abs(s) > 1e-25])
            print(line); out.append(line)
    for t in ('0', '0.37'):
        M = helicity_matrix(t)
        line = 'H_W helicity matrix at phi=%s: diag %s %s ; flip(+,-) %s (|.|=%s, arg=%s) ; flip(-,+) %s' % (
            t, mp.nstr(M[(1, 1)], 6), mp.nstr(M[(-1, -1)], 6), mp.nstr(M[(1, -1)], 6), mp.nstr(abs(M[(1, -1)]), 6),
            mp.nstr(mp.arg(M[(1, -1)]), 6), mp.nstr(M[(-1, 1)], 6))
        print(line); out.append(line)
    out.append('=> T31_W = cos(4phi)/8, tau1_W = tau3_W = 3cos(4phi)/8 (no extra first-moment term); T31_E = 2, tau_E = 6.')
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'kb_onepoint.txt'), 'w').write('\n'.join(out) + '\n')
