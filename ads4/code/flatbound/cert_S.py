"""Certificate S (scaling region): J >= JP, m >= 2, bhat=(2J+1)/m in [BS, BA]  (and the small-bhat part bhat <= BS).

Continuous extension.  With X = bhat^2/4, beta = J(J+1)u = X - u/4, lambda_i = 1/4 + e_h(i):
   FF_h(bhat,u) = G_h(bhat,u) + 10 X u^2 - (5/2 + c_h) u^3,
   G_h = sum_n (-1)^n/(n!)^2 prod_{i<=n} (X - lambda_i u) (T0_n + T1_n u) = sum_k a_k(X) u^k,
which equals FF_h at every physical point (J integer, u=1/m^2).  G_h is entire in bhat and analytic in |u|<1.
Physical points of this region satisfy u <= umax(bhat) := (bhat/(2JP+1))^2.

Per bhat-cell [bc-h, bc+h]: exact Taylor coefficients g_ik = d^i_b a_k(bc)/i! (i<I, k<K; exact rational series,
enclosed in arb, with a rigorous bound for the truncation n>NMAX), and a Cauchy bound on the polydisc
|bhat-bc|<=r, |u|<=rho from the Mehler-Dirichlet representation (Lemma S of the .md):
   |G_+| <= int_0^1 |Psi_+| cosh(Y(p)) theta/sin(theta) dp,
   |G_-| <= S_CG int_0^1 |Psi_-| (1-rho p^2)^-8 cosh(Y(p)+4 theta) theta/sin(theta) dp,
   Y(p) = r p s(rho p^2) + (bc+r) p (s(rho p^2)-1),  s(x)=arcsin(sqrt x)/sqrt x,  theta = 2 arcsin(p sqrt rho),
   S_CG >= sum_j |c_j(J)| for |J+1/2| >= (bc-r)/(2 sqrt rho).
Then FF >= Ptrunc(t,u) - M[(h/r)^I/((1-h/r)(1-u/rho)) + (u/rho)^K/((1-h/r)(1-u/rho))]  and the truncated
polynomial (plus the exact Phi part) is bounded below on sub-boxes by a mean-value form with interval gradients.
Usage: python3 cert_S.py channel BS BA JP
"""
import sys, time
import flint
from fractions import Fraction as Fr
from math import factorial
from common import PSI2, PSI3, peval, T0, T1, e_h, PHI_CONST, cg_abs_bound, set_prec

PREC = 768
set_prec(PREC)
A = flint.arb
NMAX = 420
K = 12          # u-order
I = 9           # bhat-order

def fq(x):
    x = Fr(x); return A(flint.fmpq(x.numerator, x.denominator))

def iv(lo, hi):
    return fq(lo).union(fq(hi))

# ------------------------------------------------------------------ exact a_k(X) coefficients
def build_ak(h):
    lam = [None] + [Fr(1, 4) + e_h(h, i) for i in range(1, NMAX + 1)]
    E = [Fr(1)] + [Fr(0)] * K           # E_k^{(n)}, k=0..K
    coef = [dict() for _ in range(K)]    # coef[k][j] : coefficient of X^j in a_k
    t0 = [T0(n) for n in range(NMAX + 1)]
    t1 = [T1(n, h) for n in range(NMAX + 1)]
    fact2 = Fr(1)
    for n in range(0, NMAX + 1):
        if n > 0:
            fact2 *= n * n
            for k in range(min(n, K), 0, -1):
                E[k] = E[k] + lam[n] * E[k - 1]
        sgn = -1 if n % 2 else 1
        for k in range(K):
            if k <= n:
                c = Fr(sgn * (-1) ** k) * E[k] * t0[n] / fact2
                coef[k][n - k] = coef[k].get(n - k, 0) + c
            if 1 <= k <= n + 1:
                c = Fr(sgn * (-1) ** (k - 1)) * E[k - 1] * t1[n] / fact2
                coef[k][n - k + 1] = coef[k].get(n - k + 1, 0) + c
    polys = []
    for k in range(K):
        deg = 2 * max(coef[k])
        cl = [A(0)] * (deg + 1)
        for j, c in coef[k].items():
            cl[2 * j] = fq(c / 4 ** j)
        polys.append(flint.arb_poly(cl))
    return polys

def tail_bound(h, Xbar):
    """Bound, valid for all k<K and |X|<=Xbar, of |sum_{n>NMAX} (n-term of a_k)| (geometric, ratio <= 1/2 checked)."""
    lamsum = lambda n: sum(Fr(1, 4) + e_h(h, i) for i in range(1, n + 1))
    Tmax = A(100)
    def term(n):
        Ln = fq(lamsum(n))
        return Tmax * (fq(Xbar) + 1) ** (n + 1) * 2 * (Ln + 1) ** K / A(factorial(n)) ** 2
    n0 = NMAX + 1
    t = term(n0)
    ratio = (fq(Xbar) + 1) * ((fq(lamsum(n0 + 1)) + 1) / (fq(lamsum(n0)) + 1)) ** K / A(n0 + 1) ** 2
    assert bool(ratio < A(1) / 2), 'tail ratio'
    return 2 * t

# ------------------------------------------------------------------ Cauchy bound
def s_fun(x):  # arcsin(sqrt x)/sqrt x for real x in [0,1), x ball >= 0
    if bool(x == 0): return A(1)
    sx = x.sqrt(); return sx.asin() / sx

def cauchy_M(h, bc, r, rho, NPC=160):
    sig = 1 if h == '++' else -1
    tot = A(0)
    for i in range(NPC):
        plo, phi = Fr(i, NPC), Fr(i + 1, NPC)
        p = iv(plo, phi); ph = fq(phi)
        psi2 = PSI2_ev(p); psi3 = PSI3_ev(p)
        Wabs = A(psi2.abs_upper()) * (2 + rho * ph * ph) + rho * A(psi3.abs_upper())
        x = rho * ph * ph
        sx = s_fun(x)
        Y = r * ph * sx + (bc + r) * ph * (sx - 1)
        th = 2 * (ph * rho.sqrt()).asin()
        fac = th / th.sin()
        if h == '+-':
            Wabs = Wabs / (1 - x) ** 8
            Y = Y + 4 * th
        tot += A((Wabs * Y.cosh() * fac).upper()) / NPC
    if h == '+-':
        Zmin = (bc - r) / (2 * rho.sqrt())
        Zlo = Fr(int((Zmin.lower() * 1000).floor().unique_fmpz()), 1000)
        assert Zlo > 5
        tot = tot * fq(cg_abs_bound(Zlo))
    return A(tot.upper())

def _mkev(poly):
    def f(x):
        rr = A(0)
        for c in reversed(poly): rr = rr * x + fq(c)
        return rr
    return f
PSI2_ev = _mkev(PSI2); PSI3_ev = _mkev(PSI3)

# ------------------------------------------------------------------ box certification
def eval2(g, t, u, Xc_poly):
    """Evaluate sum g[i][k] t^i u^k (interval Horner)."""
    r = A(0)
    for i in reversed(range(len(g))):
        ri = A(0)
        for k in reversed(range(len(g[i]))):
            ri = ri * u + g[i][k]
        r = r * t + ri
    return r

def deriv_t(g):
    return [[g[i][k] * i for k in range(len(g[i]))] for i in range(1, len(g))] or [[A(0)]]

def deriv_u(g):
    return [[g[i][k] * k for k in range(1, len(g[i]))] or [A(0)] for i in range(len(g))]

def certify_box(g, gt, gu, tail_fn, t0, t1, u0, u1, depth=0, stats=None):
    tm = (t0 + t1) / 2; um = (u0 + u1) / 2
    T = iv(t0, t1); U = iv(u0, u1)
    val = eval2(g, fq(tm), fq(um), None)
    dT = eval2(gt, T, U, None); dU = eval2(gu, T, U, None)
    low = val - A(dT.abs_upper()) * fq((t1 - t0) / 2) - A(dU.abs_upper()) * fq((u1 - u0) / 2)
    tl = tail_fn(u1)
    if bool(low - tl > 0):
        stats['boxes'] += 1
        return True
    if depth > 26:
        return False
    # split the dimension with the larger first-order loss
    lt = dT.abs_upper() * fq((t1 - t0) / 2); lu = dU.abs_upper() * fq((u1 - u0) / 2)
    if bool(A(lt) > A(lu)):
        return certify_box(g, gt, gu, tail_fn, t0, tm, u0, u1, depth + 1, stats) and \
               certify_box(g, gt, gu, tail_fn, tm, t1, u0, u1, depth + 1, stats)
    else:
        return certify_box(g, gt, gu, tail_fn, t0, t1, u0, um, depth + 1, stats) and \
               certify_box(g, gt, gu, tail_fn, t0, t1, um, u1, depth + 1, stats)

def cell(h, ak, bc, hh, JP, kappa=3):
    r = min(Fr(2), bc / 2)
    assert hh <= r / 10
    u0 = ((bc + hh) / (2 * JP + 1)) ** 2
    rho = fq(kappa * u0)
    rr = fq(r)
    M = cauchy_M(h, fq(bc), rr, rho)
    Xbar = ((bc + r) ** 2) / 4
    tb = tail_bound(h, Xbar)
    # Taylor coefficients g[i][k]
    g = [[None] * K for _ in range(I)]
    for k in range(K):
        pk = ak[k]
        fi = A(1)
        for i in range(I):
            v = pk(fq(bc))
            err = A((tb / rr ** i).upper())
            g[i][k] = v / fi + (-err).union(err)   # add truncation error (Cauchy on the n-tail)
            pk = pk.derivative(); fi = fi * (i + 1)
    # exact Phi part: 10 X u^2 - (5/2+c) u^3 with X = (bc+t)^2/4
    c = PHI_CONST[h]
    g[0][2] += fq(Fr(10) * bc * bc / 4); g[1][2] += fq(Fr(10) * 2 * bc / 4); g[2][2] += fq(Fr(10, 4))
    g[0][3] -= fq(Fr(5, 2) + c)
    q = fq(hh) / rr
    def tail_fn(ub):
        x = fq(ub) / rho
        return M * (q ** I + x ** K) / ((1 - q) * (1 - x))
    gt = deriv_t(g); gu = deriv_u(g)
    stats = {'boxes': 0, 'minlow': 0}
    ok = certify_box(g, gt, gu, tail_fn, -hh, hh, Fr(0), u0, 0, stats)
    return ok, stats['boxes'], M, u0

def small_b(h, BS, JP):
    """bhat <= BS, u <= (BS/(2JP+1))^2: crude absolute-series bound."""
    X = fq(Fr(BS) ** 2 / 4); u = fq((Fr(BS) / (2 * JP + 1)) ** 2)
    G0 = fq(T0(0))  # a lower bound uses T0_0 + T1_0 u >= T0_0 - |T1_0| u
    low = G0 - fq(abs(T1(0, h))) * u
    prod = A(1); tail = A(0)
    for n in range(1, 400):
        prod = prod * (X + fq(Fr(1, 4) + e_h(h, n)) * u) / (n * n)
        tail += prod * (fq(abs(T0(n))) + fq(abs(T1(n, h))) * u)
    # remaining n>=400 negligible: prod ratio < 1/2 already
    tail += 2 * prod * 100
    low = low - tail - (fq(Fr(5, 2)) + PHI_CONST[h]) * u ** 3
    return bool(low > 0), low

def main():
    h = sys.argv[1] if len(sys.argv) > 1 else '+-'
    BS = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(1)
    BA = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(120)
    JP = int(sys.argv[4]) if len(sys.argv) > 4 else 400
    lo = Fr(sys.argv[5]) if len(sys.argv) > 5 else BS
    hi = Fr(sys.argv[6]) if len(sys.argv) > 6 else BA
    t0 = time.time()
    ok_all = True
    if lo == BS:
        ok, low = small_b(h, BS, JP)
        print('S %s small-b (b<=%s): lower bound %s %s' % (h, BS, low.str(6), 'PASS' if ok else 'FAIL'), flush=True)
        ok_all &= ok
    ak = build_ak(h)
    print('S %s: a_k built (%.0fs)' % (h, time.time() - t0), flush=True)
    b = lo; ncell = 0
    while b < hi:
        r = min(Fr(2), (b + Fr(1, 10)) / 2)
        hh = min(r / 10, Fr(1, 10))
        # make the cell [b, b+2hh]; require hh <= r(bc)/10 with bc=b+hh
        bc = b + hh
        while hh > min(Fr(2), bc / 2) / 10:
            hh = hh / 2; bc = b + hh
        if bc + hh > hi:
            hh = (hi - b) / 2; bc = b + hh
        ok, nbox, M, u0 = cell(h, ak, bc, hh, JP)
        ncell += 1
        print('S %s cell b in [%s, %s] umax=%.3e M<=%s boxes=%d %s' % (
            h, float(b), float(b + 2 * hh), float(u0), M.str(4), nbox, 'PASS' if ok else 'FAIL'), flush=True)
        ok_all &= ok
        b = b + 2 * hh
    print('VERDICT cert_S %s b in [%s,%s] (%d cells, %.0fs): %s' % (h, lo, hi, ncell, time.time() - t0, 'PASS' if ok_all else 'FAIL'))

if __name__ == '__main__':
    main()
