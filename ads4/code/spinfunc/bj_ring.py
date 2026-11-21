"""
bj_ring (BELOW_FOLD_SPINJ.md §1): closed-form spin-J O-index pairing in the Fourier basis of the quartic remainder, on any coefficient
type (complex numbers for checks; bg_ser.ASer for the engine).

Every parity-even conserved <T T O_J> structure is f(z5) = (V.z)^{J-4} phi(z5), phi a quartic (a5 = J - b15 - b25 >= J - 4), and the
O-index contraction depends only on the restriction to the null cone, where phi = sum_{n=-4}^{4} phi_n e_n, e_n = (zeta_sgn(n).z)^{|n|} (u.z)^{4-|n|}
(zeta_+- = (-1, 1, +-i), u = (-1, -1, 0) in the lightcone components of bg_geom.dot: on null_circle, zeta_+-.z = e^{+-i phi}, u.z = 1).  So
   <f, g> = sum_{n,n'} phi_n gamma_n' B_{n n'},    B_{n n'} = <(V.z)^{J-4} e_n, (W.z)^{J-4} e_n'>
and by BELOW_FOLD3 §3 (Legendre identity + polarisation)
   B_{n n'} = N_J R^J ((J-4)!/J!)^2 |n|!(4-|n|)! |n'|!(4-|n'|)! [eps_s^|n| eps_u^(4-|n|) eta_s'^|n'| eta_u^(4-|n'|)] g ,   R = sqrt(V^2 W^2),
   g = (1+X)^{J/2} (1+Y)^{J/2} P_J((t+U)/sqrt((1+X)(1+Y))),  a = V/|V| + d_a, b = W/|W| + d_b, d_a = (eps_s zeta_s + eps_u u)/|V|, ...
In manifest form (q_a, q_b = O(1/Delta), everything else O(1)):
   s_a = v.d_a, sigma_a = d_a^2, q_a = w_perp.d_a (w_perp = (W - (V.W/V^2) V)/|W|), rho_a = (sigma_a - s_a^2)/(1+s_a)^2,
   (1+X)^{J/2} = (1+s_a)^J (1+rho_a)^{J/2},
   tau = t[(1+rho_a)^{-1/2}(1+rho_b)^{-1/2} - 1] + (q_a + q_b + c - t s_a s_b)(1+s_a)^{-1}(1+s_b)^{-1}(1+rho_a)^{-1/2}(1+rho_b)^{-1/2},
   g = (1+s_a)^J (1+rho_a)^{J/2} (1+s_b)^J (1+rho_b)^{J/2} sum_k P_J^(k)(t) tau^k / k!          (c = d_a.d_b).
=> B = N_J R^J sum_{k=0}^{8} P_J^(k)(t) Bk[k]; this module returns Bk (the factor N_J R^J is dropped: N_J is J-only, R^J cancels the J parts
of the 3-pt exponents exactly).  Bookkeeping for the a = 1/Delta expansion at fixed J/Delta: every J-polynomial coefficient ~ J^e is stored
as (value * a^e) * a^(-e) (Ctx.cst(value, -e)), so ASer truncation is uniform in J/Delta.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import factorial
try:
    import bj_ee  # noqa: F401  (installs the safe Poly.prune)
except ImportError:
    pass
from scipy.special import binom

ZETA = {+1: np.array([-1.0, 1.0, 1j]), -1: np.array([-1.0, 1.0, -1j])}
UVEC = np.array([-1.0, -1.0, 0.0])
DA = 4


def ldot_num(p, q):
    return -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]


class Ctx:
    """coefficient arithmetic.  mode 'num': python complex.  mode 'aser': bg_ser.ASer with depth K at Delta = Dl."""
    def __init__(self, mode='num', K=None, Dl=None):
        self.mode, self.K, self.Dl = mode, K, Dl
        if mode == 'aser':
            import bg_ser as S
            self.S = S

    def cst(self, x, shift=0):
        if self.mode == 'num':
            return complex(x)
        S = self.S
        return S.ASer(shift, [S.Poly.const(complex(x)*self.Dl**shift)] + [S.Poly() for _ in range(self.K)])

    def add(self, x, y):
        if x is None:
            return y
        if y is None:
            return x
        return x + y

    def mul(self, x, y):
        return x*y

    def sc(self, x, c):
        if self.mode == 'num':
            return x*c
        return x.scale(complex(c)) if np.iscomplexobj(c) and abs(np.imag(c)) > 0 else x.scale(float(np.real(c)))

    def dot(self, A, B):
        """lightcone dot product; A, B: 3-lists of coefficients or numbers."""
        terms = [(A[0], B[1], -0.5), (A[1], B[0], -0.5), (A[2], B[2], 1.0)]
        out = None
        for x, y, c in terms:
            if self.mode == 'aser':
                xn, yn = np.isscalar(x) or isinstance(x, np.number), np.isscalar(y) or isinstance(y, np.number)
                if xn and yn:
                    v = self.cst(x*y*c) if x*y != 0 else None
                elif xn:
                    v = self.sc(y, x*c) if x != 0 else None
                elif yn:
                    v = self.sc(x, y*c) if y != 0 else None
                else:
                    v = self.sc(x*y, c)
            else:
                v = x*y*c
            out = self.add(out, v)
        return out


# ------------------------------------------------------------------------------ truncated polynomial ring R_a (x) R_b
def _ok(m):
    return m[0] + m[1] <= DA and m[2] + m[3] <= DA


def rmul(ctx, X, Y):
    out = {}
    for m1, c1 in X.items():
        for m2, c2 in Y.items():
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2], m1[3] + m2[3])
            if not _ok(m):
                continue
            out[m] = ctx.add(out.get(m), ctx.mul(c1, c2))
    return out


def radd(ctx, X, Y, cy=1.0):
    out = dict(X)
    for m, c in Y.items():
        v = c if cy == 1.0 else ctx.sc(c, cy)
        out[m] = ctx.add(out.get(m), v)
    return out


def rsc(ctx, X, c):
    return {m: ctx.sc(v, c) for m, v in X.items()}


def rmulc(ctx, X, c):
    """ring element times a coefficient (ASer / number)."""
    return {m: ctx.mul(v, c) for m, v in X.items()}


ONE = (0, 0, 0, 0)


def rseries(ctx, x, coefs):
    """sum_m coefs[m] x^m, x nilpotent (no constant term); coefs[m] are coefficient objects (ctx.cst) or None."""
    out = {ONE: coefs[0]} if coefs[0] is not None else {}
    pw = None
    for mm in range(1, len(coefs)):
        pw = dict(x) if pw is None else rmul(ctx, pw, x)
        if not pw:
            break
        if coefs[mm] is not None:
            out = radd(ctx, out, rmulc(ctx, pw, coefs[mm]))
    return out


def binser(ctx, p, mmax, jpow):
    """coefficients C(p, m), m = 0..mmax, as ctx constants with shift -jpow*m (jpow = 1 if p ~ J, 0 if p = O(1))."""
    return [ctx.cst(binom(p, m), -jpow*m) for m in range(mmax + 1)]


def bk_matrices(ctx, J, V, W, verbose=False):
    """Bk[k] (k = 0..8): dict {(n, n'): coefficient} with B_{nn'} = N_J R^J sum_k P_J^(k)(t) Bk[k][(n, n')];  also returns t.
    V, W: 3-lists (lightcone components) of coefficient objects."""
    dot = ctx.dot
    VV, WW, VW = dot(V, V), dot(W, W), dot(V, W)
    if ctx.mode == 'num':
        R = np.sqrt(VV*WW)
        iV2, iW2, iR = 1/VV, 1/WW, 1/R
    else:
        import bg_ser as S
        iV2, iW2 = VV.recip('V^2'), WW.recip('W^2')
        c0, e0, eps = (VV*WW).rel('V^2W^2')
        assert np.all(e0 % 2 == 0), e0
        r = S.series_fn(eps, 'pow', -0.5)
        iR = S.ASer(-(VV*WW).d//2, [p.mul_mono(-e0//2, c0**-0.5) for p in r.c])
        assert (VV*WW).d % 2 == 0
    t = ctx.mul(VW, iR)
    # W_perp = W - (V.W/V^2) V, V_perp = V - (V.W/W^2) W  (manifestly O(a))
    fw, fv = ctx.mul(VW, iV2), ctx.mul(VW, iW2)
    Wp = [ctx.add(W[i], ctx.sc(ctx.mul(fw, V[i]), -1.0)) for i in range(3)]
    Vp = [ctx.add(V[i], ctx.sc(ctx.mul(fv, W[i]), -1.0)) for i in range(3)]
    if ctx.mode == 'aser':
        Wp = [x.norm() for x in Wp]; Vp = [x.norm() for x in Vp]
    pre = ctx.cst(1.0/float(J*(J - 1)*(J - 2)*(J - 3))**2, 8)                                    # ((J-4)!/J!)^2
    Bk = [dict() for _ in range(9)]
    for s in (+1, -1):
        for s2 in (+1, -1):
            A = [ZETA[s], UVEC]; Bv = [ZETA[s2], UVEC]
            Sa = [ctx.mul(dot(list(Ai), V), iV2) for Ai in A]
            Sb = [ctx.mul(dot(list(Bi), W), iW2) for Bi in Bv]
            Qa = [ctx.mul(dot(list(Ai), Wp), iR) for Ai in A]
            Qb = [ctx.mul(dot(list(Bi), Vp), iR) for Bi in Bv]
            # sigma_a = d_a^2: zeta.zeta = 0, zeta.u = 0, u.u = -1  -> sigma_a = -eps_u^2 / V^2
            sa = {(1, 0, 0, 0): Sa[0], (0, 1, 0, 0): Sa[1]}
            sb = {(0, 0, 1, 0): Sb[0], (0, 0, 0, 1): Sb[1]}
            qa = {(1, 0, 0, 0): Qa[0], (0, 1, 0, 0): Qa[1]}
            qb = {(0, 0, 1, 0): Qb[0], (0, 0, 0, 1): Qb[1]}
            siga = {(0, 2, 0, 0): ctx.sc(iV2, -1.0)}
            sigb = {(0, 0, 0, 2): ctx.sc(iW2, -1.0)}
            cc = {(0, 1, 0, 1): ctx.sc(iR, -1.0)}
            zz = ldot_num(ZETA[s], ZETA[s2])
            if abs(zz) > 0:
                cc[(1, 0, 1, 0)] = ctx.sc(iR, zz)
            side = {}
            for nm, sx, sig in (('a', sa, siga), ('b', sb, sigb)):
                inv1 = rseries(ctx, sx, [ctx.cst((-1.0)**m) for m in range(DA + 1)])            # (1+s)^-1
                inv2 = rmul(ctx, inv1, inv1)
                rho = rmul(ctx, radd(ctx, sig, rmul(ctx, sx, sx), -1.0), inv2)
                powJ = rseries(ctx, sx, binser(ctx, J, DA, 1))                                  # (1+s)^J
                rJ2 = rseries(ctx, rho, binser(ctx, J/2.0, DA//2, 1))                           # (1+rho)^{J/2}
                rmh = rseries(ctx, rho, binser(ctx, -0.5, DA//2, 0))                            # (1+rho)^{-1/2}
                side[nm] = dict(A=rmul(ctx, powJ, rJ2), beta=rmh, alpha=rmul(ctx, inv1, rmh))
            a_, b_ = side['a'], side['b']
            AB = rmul(ctx, a_['A'], b_['A'])
            bb = rmul(ctx, a_['beta'], b_['beta'])
            bb1 = {m: v for m, v in bb.items() if m != ONE}                                     # beta_a beta_b - 1
            # bb has constant term 1 exactly (both series start with cst(1)): removed above
            t1 = rmulc(ctx, bb1, t)
            lin = radd(ctx, radd(ctx, radd(ctx, qa, qb), cc), rmulc(ctx, rmul(ctx, sa, sb), t), -1.0)
            t2 = rmul(ctx, lin, rmul(ctx, a_['alpha'], b_['alpha']))
            tau = radd(ctx, t1, t2)
            keys = {}
            for n in range(0, DA + 1):
                for n2 in range(0, DA + 1):
                    keys[(n, n2)] = (n, DA - n, n2, DA - n2)
            tk = {ONE: ctx.cst(1.0)}
            for k in range(0, 9):
                if k:
                    tk = rmul(ctx, tk, tau)
                    tk = rsc(ctx, tk, 1.0/k)
                    if not tk:
                        break
                for (n, n2), key in keys.items():
                    nn, nn2 = s*n, s2*n2
                    if (n == 0 and s == -1) or (n2 == 0 and s2 == -1):
                        continue
                    acc = None
                    for mt, ct in tk.items():
                        ma = (key[0] - mt[0], key[1] - mt[1], key[2] - mt[2], key[3] - mt[3])
                        if min(ma) < 0 or ma not in AB:
                            continue
                        acc = ctx.add(acc, ctx.mul(AB[ma], ct))
                    if acc is None:
                        continue
                    fac = factorial(n)*factorial(DA - n)*factorial(n2)*factorial(DA - n2)
                    Bk[k][(nn, nn2)] = ctx.mul(ctx.sc(acc, fac), pre)
                if verbose:
                    print(f'   sign ({s:+d},{s2:+d}) k={k}: tau^k terms {len(tk)}', flush=True)
    return Bk, t


def legendre_derivs(J, t, kmax):
    import mpmath as mp
    return [complex(mp.fac2(2*k - 1)*mp.gegenbauer(J - k, k + mp.mpf(1)/2, mp.mpc(t))) if k <= J else 0.0 for k in range(kmax + 1)]


def _test():
    import mpmath as mp
    import bh_spinJ as BS
    rng = np.random.default_rng(5)
    ctx = Ctx('num')
    lines = []
    worst = 0
    for J in (4, 6, 8, 12, 24, 32):
        R = lambda: rng.normal(size=3) + 1j*rng.normal(size=3)
        V = R(); W = V + 0.05*R()
        Bk, t = bk_matrices(ctx, J, list(V), list(W))
        Pd = legendre_derivs(J, t, 8)
        NJ = float(mp.fac(J)/mp.fac2(2*J - 1))
        Rr = np.sqrt(ldot_num(V, V)*ldot_num(W, W))
        for (n, n2) in ((0, 0), (1, -1), (2, 1), (-3, 4), (4, -4), (-2, -2)):
            val = NJ*Rr**J*sum(Pd[k]*Bk[k].get((n, n2), 0) for k in range(9))
            As = [ZETA[1 if n >= 0 else -1]]*abs(n) + [UVEC]*(4 - abs(n))
            Bs = [ZETA[1 if n2 >= 0 else -1]]*abs(n2) + [UVEC]*(4 - abs(n2))
            ref = BS.pair_sampled_mp(J, V, As, W, Bs, dps=50)
            e = abs(val/ref - 1); worst = max(worst, e)
            lines.append(f'  J={J:3d} (n,n\')=({n:+d},{n2:+d}): closed {val:.6e}  sampled(mp) {ref:.6e}  rel {e:.1e}')
    lines.append(f'  worst {worst:.1e}')
    print('\n'.join(lines))
    return lines


if __name__ == '__main__':
    out = _test()
    open(os.path.join(HERE, 'cache', 'bj_ring_test.txt'), 'w').write('\n'.join(out) + '\n')


# ------------------------------------------------------------------------------ vectorised version (bj_xs), production
def bk_matrices_xs(J, V, W, K, Dl, verbose=False):
    """Bk as bk_matrices, with the ring arithmetic of bj_xs (one table per ring element).  Bookkeeping: all J-dependent numbers are
    kept EXACT (not expanded in 1/J); only geometric a-orders are truncated, relative to O(1) leading terms.  For that every geometric
    object entering the ring starts at a^0: the generators are rescaled by 1/J (s_a -> S/J, sigma_a, c -> /J^2) and the manifestly
    small q_a = (w_perp.A) = a qhat_a is entered as qhat_a (d = 0) times the number a/J; then every term of
    P^(k)(t) tau^k ~ J^{2k} (a/J)^f J^{-2s} = (J a)^f is O(1) at fixed J/Delta with a uniformly truncated geometric factor.
    Returns (Bk, t, extras) with Bk[k][(n, n')] ASer (k = 0..8)."""
    import bg_ser as S
    from bj_xs import XS
    import time
    t0 = time.time()
    ctx = Ctx('aser', K, float(Dl))
    dot = ctx.dot
    VV, WW, VW = dot(V, V), dot(W, W), dot(V, W)
    iV2, iW2 = VV.recip('V^2'), WW.recip('W^2')
    c0, e0, eps = (VV*WW).rel('V^2W^2')
    assert np.all(e0 % 2 == 0), e0
    r = S.series_fn(eps, 'pow', -0.5)
    assert (VV*WW).d % 2 == 0
    iR = S.ASer(-(VV*WW).d//2, [p.mul_mono(-e0//2, c0**-0.5) for p in r.c])
    t = VW*iR
    fw, fv = VW*iV2, VW*iW2
    Wp = [(W[i] - fw*V[i]).norm() for i in range(3)]
    Vp = [(V[i] - fv*W[i]).norm() for i in range(3)]
    BIG = 10**4
    Dl = float(Dl)

    def cst(x, shift=0):
        assert shift == 0
        return XS.const(complex(x), 0, BIG)

    def ser(x, coefs):
        out = coefs[0]
        pw = None
        for m in range(1, len(coefs)):
            pw = x if pw is None else pw*x
            if len(pw) == 0:
                break
            out = out + pw*coefs[m]
        return out

    def bser(p, mmax, jpow):
        return [cst(binom(p, m)) for m in range(mmax + 1)]

    def E(A, e, jp=0, small=False):
        """ASer A times eps-monomial e, times J^{-jp} (exact number); small: A = a Ahat (A.d >= 1), entered as Ahat (d - 1) times a."""
        if small:
            assert A.d >= 1, ('q not manifestly small', A.d)
            A = S.ASer(A.d - 1, A.c)
        X = XS.from_aser(A, e)
        fac = float(J)**(-jp)*((1.0/Dl) if small else 1.0)
        return X.scale(fac) if fac != 1.0 else X
    pre = float(J)**8/float(J*(J - 1)*(J - 2)*(J - 3))**2
    tX = E(t, (0, 0, 0, 0))
    Bk = [dict() for _ in range(9)]
    for s in (+1, -1):
        for s2 in (+1, -1):
            A = [ZETA[s], UVEC]; Bv = [ZETA[s2], UVEC]
            Sa = [dot(list(Ai), V)*iV2 for Ai in A]; Sb = [dot(list(Bi), W)*iW2 for Bi in Bv]
            Qa = [dot(list(Ai), Wp)*iR for Ai in A]; Qb = [dot(list(Bi), Vp)*iR for Bi in Bv]
            sa = E(Sa[0], (1, 0, 0, 0), 1) + E(Sa[1], (0, 1, 0, 0), 1)
            sb = E(Sb[0], (0, 0, 1, 0), 1) + E(Sb[1], (0, 0, 0, 1), 1)
            qa = E(Qa[0], (1, 0, 0, 0), 1, True) + E(Qa[1], (0, 1, 0, 0), 1, True)
            qb = E(Qb[0], (0, 0, 1, 0), 1, True) + E(Qb[1], (0, 0, 0, 1), 1, True)
            siga = E(iV2.scale(-1.0), (0, 2, 0, 0), 2); sigb = E(iW2.scale(-1.0), (0, 0, 0, 2), 2)
            cc = E(iR.scale(-1.0), (0, 1, 0, 1), 2)
            zz = ldot_num(ZETA[s], ZETA[s2])
            if abs(zz) > 0:
                cc = cc + E(iR.scale(complex(zz)), (1, 0, 1, 0), 2)
            side = {}
            for nm, sx, sig in (('a', sa, siga), ('b', sb, sigb)):
                inv1 = ser(sx, [cst((-1.0)**m) for m in range(DA + 1)])
                rho = (sig - sx*sx)*(inv1*inv1)
                powJ = ser(sx, bser(J, DA, 1))
                rJ2 = ser(rho, bser(J/2.0, DA//2, 1))
                rmh = ser(rho, bser(-0.5, DA//2, 0))
                side[nm] = dict(A=powJ*rJ2, beta=rmh, alpha=inv1*rmh)
            a_, b_ = side['a'], side['b']
            AB = a_['A']*b_['A']*cst(pre, 0)
            bb1 = a_['beta']*b_['beta'] - cst(1.0)
            t1 = bb1*tX
            lin = qa + qb + cc - sa*sb*tX
            tau = t1 + lin*(a_['alpha']*b_['alpha'])
            tk = cst(1.0)
            for k in range(0, 9):
                if k:
                    tk = (tk*tau).scale(1.0/k)
                    if len(tk) == 0:
                        break
                G = AB*tk
                for n in range(0, DA + 1):
                    for n2 in range(0, DA + 1):
                        if (n == 0 and s == -1) or (n2 == 0 and s2 == -1):
                            continue
                        key = (n, DA - n, n2, DA - n2)
                        A_ = G.eps_part(key)
                        if all(len(p) == 0 for p in A_.c):
                            continue
                        fac = factorial(n)*factorial(DA - n)*factorial(n2)*factorial(DA - n2)
                        Bk[k][(s*n, s2*n2)] = A_.scale(float(fac))
                if verbose:
                    print(f'   xs sign ({s:+d},{s2:+d}) k={k}: tau^k rows {len(tk)} d {tk.d} top {tk.top}  {time.time() - t0:.1f}s', flush=True)
    return Bk, t, dict(iR=iR, iV2=iV2, iW2=iW2)


def _test_xs(K=2, Dl=300.0, J=6, avals=(1.0,), flip=True, verbose=False):
    """XS/ASer pairing vs the exact numeric closed form at the same geometric point (V, W evaluated): sum_k P^(k)(t) Bk."""
    import bg_ser as S, bg_geom as GG, bg_ee as BE, bj_ee as BJ, time
    G = GG.geometry(K + 2); Gm = BE.EEGeom(K + 2, G); B = Gm.basis
    V = BJ.tv([Gm.V(5, B[m], 1, 2) for m in range(3)]); W = BJ.tv([Gm.V(5, B[m], 3, 4) for m in range(3)])
    if flip:
        W = [w.scale(-1.0) for w in W]
    VK = [S.ASer(v.d, v.c[:K + 2]) for v in V]; WK = [S.ASer(w.d, w.c[:K + 2]) for w in W]
    t0 = time.time(); Bx, tx, _ = bk_matrices_xs(J, VK, WK, K, Dl); t1 = time.time()
    pts = S.atom_points(4, 3)
    out = []
    for af in avals:
        a = af/Dl
        worst = 0
        for ip in range(len(pts)):
            P = pts[ip:ip + 1]
            Vn = [complex(v.ev(P, a)[0]) for v in V]; Wn = [complex(w.ev(P, a)[0]) for w in W]
            Bn, tn = bk_matrices(Ctx('num'), J, Vn, Wn)
            Pd = legendre_derivs(J, tn, 8)
            tt = complex(tx.ev(P, a)[0]); Pa = legendre_derivs(J, tt, 8)
            exs = {key: sum(Pd[k]*Bn[k].get(key, 0) for k in range(9)) for key in Bn[0]}
            aps = {key: sum(Pa[k]*complex(Bx[k][key].ev(P, a)[0]) for k in range(9) if key in Bx[k]) for key in Bn[0]}
            scale = max(abs(v) for v in exs.values())
            for key in exs:
                worst = max(worst, abs(aps[key] - exs[key])/scale)
                if verbose and ip == 0:
                    print(key, exs[key], aps[key])
        out.append(f'  J={J} K={K} Delta={Dl:g} a=1/{1/a:g}: max |diff|/max|B| (all n,n\') {worst:.2e}')
    out.append(f'  xs time {t1 - t0:.1f}s')
    print('\n'.join(out), flush=True)
    return out
