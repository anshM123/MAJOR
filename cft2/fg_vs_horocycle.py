"""Regulators in a Banados geometry (boundary-graviton dressing of Poincare AdS3).
Banados FG cutoff z = eps  vs  the 'frame horocyclic' cutoff (lift a_p = (f(x_p),1)/sqrt f'(x_p), same for fb).
Frame formula (exact for horocycles = exact for twist-field normalisation in the CFT):
    l_frame(1,2) = log[ (f1-f2)(fb1-fb2) / (eps^2 sqrt(f1' f2' fb1' fb2')) ]
Tests: (a) point-to-point FG length  d_pp = l_frame + c_pp(1,2) eps^2 + O(eps^4);
       (b) FG-truncated boundary-anchored geodesic d_tr = l_frame + c_tr(1,2) eps^2 + O(eps^4);
       (c) the O(eps^2) coefficients are NOT of Bell (endpoint-local) form a_i + a_j:
           X = c(1,3)+c(2,4)-c(1,4)-c(2,3) != 0.
Uses the map checked in banados_map_check.py."""
import mpmath as mp
mp.mp.dps = 40
f   = lambda x: mp.e**(2*x + mp.mpf('0.3')*mp.sin(x))
fb  = lambda y: mp.tan(y/2 + mp.mpf('0.2')*mp.sin(y)/3)
d1 = lambda F, x: mp.diff(F, x); d2 = lambda F, x: mp.diff(F, x, 2)
def to_poinc(z, x, y):
    F1, F2, G1, G2 = d1(f, x), d2(f, x), d1(fb, y), d2(fb, y)
    D = 4*F1*G1 + z**2*F2*G2
    return (4*z*(F1*G1)**mp.mpf(1.5)/D, f(x) - 2*z**2*F1**2*G2/D, fb(y) - 2*z**2*G1**2*F2/D)
def dist(P, Q):
    (u1, w1, b1), (u2, w2, b2) = P, Q
    s = 1 + ((u1-u2)**2 + (w1-w2)*(b1-b2))/(2*u1*u2)
    return mp.acosh(s)
def l_frame(p, q, eps):
    (x1, y1), (x2, y2) = p, q
    return mp.log((f(x1)-f(x2))*(fb(y1)-fb(y2))/(eps**2*mp.sqrt(d1(f, x1)*d1(f, x2)*d1(fb, y1)*d1(fb, y2))))
def z_of(P, guess):
    """invert the map: Poincare point -> Banados (z,x,y) by Newton"""
    v = mp.matrix(guess)
    for _ in range(60):
        F = mp.matrix(to_poinc(v[0], v[1], v[2])) - mp.matrix(P)
        if mp.norm(F) < mp.mpf(10)**(-34): break
        J = mp.matrix(3, 3); h = mp.mpf(10)**(-15)
        for k in range(3):
            e = mp.matrix(v); e[k] += h
            Fk = mp.matrix(to_poinc(e[0], e[1], e[2])) - mp.matrix(P)
            for i in range(3): J[i, k] = (Fk[i]-F[i])/h
        v = v - mp.lu_solve(J, F)
    return v
def d_trunc(p, q, eps):
    """Poincare geodesic between boundary points W_p, W_q, truncated where Banados z = eps."""
    (x1, y1), (x2, y2) = p, q
    W1, B1, W2, B2 = f(x1), fb(y1), f(x2), fb(y2)
    # geodesic: semicircle  u^2 + (w-wc)(b-bc) = R2 in the (w,b) plane through the two points, param s in (0,1):
    # w = W1 + s (W2-W1), b = B1 + s (B2-B1), u^2 = s(1-s)(W2-W1)(B2-B1)
    def P(s): return (mp.sqrt(s*(1-s)*(W2-W1)*(B2-B1)), W1 + s*(W2-W1), B1 + s*(B2-B1))
    def cut(end, sgn):
        # z(P(s)) = eps near the end; s measured from that end; first-order guess s0 = eps^2 f'fb'/(dW dB)
        s0 = eps**2*d1(f, end[0])*d1(fb, end[1])/((W2-W1)*(B2-B1))
        Pe = (lambda s: P(s)) if sgn > 0 else (lambda s: P(1-s))
        g = lambda s: z_of(Pe(s), (eps, end[0], end[1]))[0] - eps
        return mp.findroot(g, (s0*mp.mpf('0.5'), s0*2), solver='anderson')
    sa = cut(p, +1); sb = 1 - cut(q, -1)
    return dist(P(sa), P(sb))
pts = [(mp.mpf('-0.9'), mp.mpf('0.8')), (mp.mpf('-0.3'), mp.mpf('0.35')), (mp.mpf('0.25'), mp.mpf('-0.2')), (mp.mpf('0.85'), mp.mpf('-0.7'))]
# points (x, y) with x increasing and y decreasing: f increasing, fb increasing -> Delta f * Delta fb < 0 ...
# make spacelike: use y -> -y ordering so that (f1-f2)(fb1-fb2) > 0
pts = [(x, -y) for (x, y) in pts]
c_pp, c_tr = {}, {}
for i in range(4):
    for j in range(i+1, 4):
        p, q = pts[i], pts[j]
        r = []
        for eps in (mp.mpf('1e-3'), mp.mpf('5e-4')):
            Pp = to_poinc(eps, *p); Pq = to_poinc(eps, *q)
            r.append(((dist(Pp, Pq) - l_frame(p, q, eps))/eps**2, (d_trunc(p, q, eps) - l_frame(p, q, eps))/eps**2))
        c_pp[(i, j)] = r[1][0]; c_tr[(i, j)] = r[1][1]
        print(f"pair {i}{j}: (d_pp - l_frame)/eps^2 = {mp.nstr(r[0][0],8)}, {mp.nstr(r[1][0],8)} ;"
              f"  (d_trunc - l_frame)/eps^2 = {mp.nstr(r[0][1],8)}, {mp.nstr(r[1][1],8)}")
for name, c in (("point-to-point", c_pp), ("truncated", c_tr)):
    X = c[(0, 2)] + c[(1, 3)] - c[(0, 3)] - c[(1, 2)]
    print(f"{name}: non-Bell test X = c02+c13-c03-c12 = {mp.nstr(X, 8)}  (0 iff coefficients were endpoint-local)")
