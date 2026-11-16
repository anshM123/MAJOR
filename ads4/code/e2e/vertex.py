"""
Exact frame integrand of the Lorentzian shadow representation of the spinning partial wave and its heavy limit.

Frame (CHMRS (C.2)): (x_1,x_2,x_3,x_4;x_5) = (-y, x/x^2, y, -x/x^2; e), n = (1,1,0) [coordinates (t, z, perp)], x, y future timelike.
   T(z_1,z_3;x,y) := (-1)^J/(4 b) * s13^{l+1} s24^{Dphi} * L(z_1) * Fs(z_3),      b = b^{(1-Dphi)}_{Delta,J},
   L  = sum_w c_w(Delta,J)      D_w S3(1,2,5)          [left structure,  (Delta, J)],
   Fs = sum_w c_w(3-Delta,-1-J) D_w S3^S(3,4,5)        [the same analytic family at the Lorentzian-shadow point],
so that   Pihat^{K}[stripped W_diff] = 2 int d^3x d^3y (|x||y|)^{k-4} K^I(xhat,yhat) T_I(x,y)   (exact; SPINNING_E2E.md sect. 4).
Heavy limit (x = X/m, y = Y/m):
   T -> (-1)^J/(4b) * 2 * 4^{l+1+Dphi} |x|^{2 Dphi} |y|^{2(l+1)} e^{2 p.x + 2 pbar.y} V_l(z_1) V^R_l(z_3),
   p = (Delta, J, 0), pbar = (Delta, -J, 0);  n_par = (-J, Delta, 0)/m, e_perp = (0,0,1):
   V_1 = Delta J m (z.n_par), V^R_1 = -V_1;   V_2 = V^R_2 = (Delta^2 J^2 m^2/2) [ (z.n_par)^2 - (z.e_perp)^2 ].
"""
import mpmath as mp, sympy as sp
from structures import rel_structure, logMref, Dphi as symDphi, get as get_struct
from embed import Delta as symDelta, J as symJ
from formulas import log_b

def mdot(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def msub(u, v): return [u[i] - v[i] for i in range(3)]

_cache = {}
def struct_fun(l):
    """returns f(Delta, J, Dphi)(x1, x2, x5, z1, n) = (log of positive prefactor, polynomial part) of sum_w c_w D_w S3."""
    if l in _cache: return _cache[l]
    rs = rel_structure(l)
    sy, V5 = rs[0][2], rs[0][3]
    cs, se = get_struct(l)
    expr = sum(cs[w]*Q for (w, Q, _, _) in rs)
    names = ['s12', 's15', 's25', 'z12', 'z15', 'z51', 'z52', 'z1z5']
    args = [sy[n] for n in names] + [symDelta, symJ, symDphi]
    f = sp.lambdify(args, expr, modules='mpmath')
    # log M_ref with |V5|
    lm = logMref(l, sy, V5)
    fl = sp.lambdify(args, lm.subs(sp.log(V5), sp.log(sp.Abs(V5))), modules='mpmath')
    def F(D, Jv, Dp, x1, x2, x5, z1, n):
        s12 = -mdot(msub(x1, x2), msub(x1, x2)); s15 = -mdot(msub(x1, x5), msub(x1, x5)); s25 = -mdot(msub(x2, x5), msub(x2, x5))
        a = [s12, s15, s25, mdot(z1, msub(x2, x1)), mdot(z1, msub(x5, x1)), mdot(n, msub(x1, x5)), mdot(n, msub(x2, x5)), mdot(z1, n), D, Jv, Dp]
        return fl(*a), f(*a)
    _cache[l] = F
    return F

def integrand_exact(l, Delta, J, Dphi, x, y, z1, z3):
    """returns (log scale, value) with T = (-1)^J * exp(logscale) * value"""
    F = struct_fun(l)
    Delta = mp.mpf(Delta)
    e = [mp.mpf(1), mp.mpf(0), mp.mpf(0)]; n = [mp.mpf(1), mp.mpf(1), mp.mpf(0)]
    x2n = mdot(x, x)
    x1 = [-v for v in y]; x2 = [v/x2n for v in x]; x3 = list(y); x4 = [-v/x2n for v in x]
    lL, pL = F(Delta, J, Dphi, x1, x2, e, z1, n)
    lR, pR = F(3 - Delta, -1 - J, Dphi, x3, x4, e, z3, n)
    s13 = -mdot(msub(x1, x3), msub(x1, x3)); s24 = -mdot(msub(x2, x4), msub(x2, x4))
    logscale = lL + lR + (l+1)*mp.log(s13) + Dphi*mp.log(s24) - mp.log(4) - log_b(Delta, J, 1 - Dphi)
    return logscale, pL*pR

def integrand_leading(l, Delta, J, Dphi, x, y, z1, z3):
    Delta = mp.mpf(Delta)
    al, alb = Delta - J - 2, Delta + J - 1
    m = mp.sqrt(al*alb)
    p = [Delta, mp.mpf(J), mp.mpf(0)]; pb = [Delta, -mp.mpf(J), mp.mpf(0)]
    npar = [-mp.mpf(J)/m, Delta/m, mp.mpf(0)]; eperp = [mp.mpf(0), mp.mpf(0), mp.mpf(1)]
    ax = mp.sqrt(-mdot(x, x)); ay = mp.sqrt(-mdot(y, y))
    logscale = (mp.log(2) + (l + 1 + Dphi)*mp.log(4) + 2*Dphi*mp.log(ax) + 2*(l+1)*mp.log(ay) + 2*mdot(p, x) + 2*mdot(pb, y)
                - mp.log(4) - log_b(Delta, J, 1 - Dphi))
    if l == 1:
        V = lambda z: Delta*J*m*mdot(z, npar)
        return logscale, -V(z1)*V(z3)
    V = lambda z: Delta**2*J**2*m**2/2*(mdot(z, npar)**2 - mdot(z, eperp)**2)
    return logscale, V(z1)*V(z3)

def integrand_leading_weyl(l, Delta, J, Dphi, x, y, z1, z3):
    """the same leading form written in the Weyl-odd variables A_+ = Delta - 3/2 = (alpha+alphabar)/2, A_- = J + 1/2 = (alphabar-alpha)/2:
    P = (A_+, A_-, 0), Pbar = (A_+, -A_-, 0), |P| = m, n_par = (-A_-, A_+, 0)/m, prefactors (A_+ A_- m)^l."""
    Delta = mp.mpf(Delta)
    Ap, Am = Delta - mp.mpf(3)/2, J + mp.mpf(1)/2
    m = mp.sqrt(Ap*Ap - Am*Am)
    p = [Ap, Am, mp.mpf(0)]; pb = [Ap, -Am, mp.mpf(0)]
    npar = [-Am/m, Ap/m, mp.mpf(0)]; eperp = [mp.mpf(0), mp.mpf(0), mp.mpf(1)]
    ax = mp.sqrt(-mdot(x, x)); ay = mp.sqrt(-mdot(y, y))
    logscale = (mp.log(2) + (l + 1 + Dphi)*mp.log(4) + 2*Dphi*mp.log(ax) + 2*(l+1)*mp.log(ay) + 2*mdot(p, x) + 2*mdot(pb, y)
                - mp.log(4) - log_b(Delta, J, 1 - Dphi))
    if l == 1:
        V = lambda z: Ap*Am*m*mdot(z, npar)
        return logscale, -V(z1)*V(z3)
    V = lambda z: Ap**2*Am**2*m**2/2*(mdot(z, npar)**2 - mdot(z, eperp)**2)
    return logscale, V(z1)*V(z3)

if __name__ == "__main__":
    mp.mp.dps = 40
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    pr("# exact frame integrand / leading heavy form, x = X/m, y = Y/m, J/Delta fixed; several null polarisations (z_1, z_3)")
    X = [mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('-0.5')]; Y = [mp.mpf('0.9'), mp.mpf('-0.3'), mp.mpf('0.2')]
    def nullv(th): return [mp.mpf(1), mp.cos(th), mp.sin(th)]
    pols = [(0.3, 1.9), (2.2, 0.4), (4.0, 5.1)]
    for l in (1, 2):
        for Dp in (1.75,):
            for (D, Jv) in ((50.25, 20), (500.25, 200), (5000.25, 2000), (50000.25, 20000)):
                m = mp.sqrt((mp.mpf(D) - Jv - 2)*(mp.mpf(D) + Jv - 1))
                x = [v/m for v in X]; y = [v/m for v in Y]
                row = []
                for (t1, t3) in pols:
                    le, ve = integrand_exact(l, D, Jv, Dp, x, y, nullv(t1), nullv(t3))
                    ll, vl = integrand_leading(l, D, Jv, Dp, x, y, nullv(t1), nullv(t3))
                    row.append(mp.nstr(mp.exp(le - ll)*ve/vl, 8))
                pr(f"  l={l} Dphi={Dp} Delta={D} J={Jv} m={mp.nstr(m, 6)}: {row}")
    pr("# O(1/m) term of the integrand, pointwise, relative to the leading form in Weyl-odd variables (integrand_leading_weyl), times b-independent")
    pr("#   m*(exact/leading - 1)  for (z_1,z_3), for (z_3,z_1), and for the 1<->3 symmetrised tensor.")
    pr("#   Claim (Weyl reflection (x,y;alpha,alphabar) -> (-x,-y;-alpha,-alphabar) at fixed X = m x, Y = m y): the O(1/m) term is antisymmetric under")
    pr("#   z_1 <-> z_3, so the last column must decrease like 1/m.   [the factor 1/b is common and is tested separately in t3/t7]")
    for l in (1, 2):
        for (t1, t3) in ((0.3, 1.9), (4.0, 5.1)):
            for (D, Jv) in ((500.25, 200), (5000.25, 2000), (50000.25, 20000)):
                m = mp.sqrt((mp.mpf(D) - Jv - 2)*(mp.mpf(D) + Jv - 1))
                x = [v/m for v in X]; y = [v/m for v in Y]
                le, ve = integrand_exact(l, D, Jv, 1.75, x, y, nullv(t1), nullv(t3)); ll, vl = integrand_leading_weyl(l, D, Jv, 1.75, x, y, nullv(t1), nullv(t3))
                le2, ve2 = integrand_exact(l, D, Jv, 1.75, x, y, nullv(t3), nullv(t1)); ll2, vl2 = integrand_leading_weyl(l, D, Jv, 1.75, x, y, nullv(t3), nullv(t1))
                r1 = mp.exp(le - ll)*ve/vl; r2 = mp.exp(le2 - ll2)*ve2/vl2
                rs = (mp.exp(le - ll)*ve + mp.exp(le2 - ll2)*ve2)/(vl + vl2)
                pr(f"  l={l} pol=({t1},{t3}) m={mp.nstr(m, 6)}:  {mp.nstr(m*(r1-1), 6)}   {mp.nstr(m*(r2-1), 6)}   symmetrised: {mp.nstr(m*(rs-1), 6)}   m^2*(sym-1): {mp.nstr(m*m*(rs-1), 6)}")
    open('out_vertex.txt', 'w').write("\n".join(out) + "\n")
