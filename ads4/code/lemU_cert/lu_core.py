"""
lemU_cert/lu_core.py: shared machinery for the Lemma U certificate (ANGULAR_ROUTE §80, FULLPOS4 §8), python-flint balls.

Objects (notation of FULLPOS4 §8; fp4_lemU.py is the DV ancestor):
  * jets  y_j = d^p [t^j] F_{p,nu}(c + c d t),  F_{p,nu}(z) = (1-z)^{-p} P_{-1/2+i nu}((1+z)/(1-z)),  from the ODE recurrence
    (FULLPOS4 fp4_ref.jet) seeded by P_nu(cosh beta_c) and d_beta P_nu (Arb Legendre, type 3).
  * integrated jets  Y^A_j = 2 int_0^1 psi2(p) y^{(3)}_j(Lam p) dp,  Y^D_j = (Lam^2/d) int_0^1 psi3(p) y^{(4)}_j(Lam p) dp.
  * dressing contraction Y^B = D_b Y^A + sigma_b Y^D (exact, FULLPOS4 §8).
  * R recursion R_{m,j}(c, x) (x = x_* or any fixed x) and A_m = sum_j R_{m,j} Y_j,  P(y) = sum A_m y^m + d^2 x y sum (q_b A_m + B_m) y^m.
  * Chebyshev on y in [0,1]: Phat_l = (2 - delta_l0) sum_{m>=l} C(2m, m-l) 4^{-m} p_m.
Every function works in ball arithmetic; mid=True drops radii (DV mode only).
"""
import os, sys, json, math
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')
K0 = Fr('27.48509817')
LAM2 = Fr(998*999)


def A(q):
    q = Fr(q)
    return arb(fmpq(q.numerator, q.denominator))


def psi_coeffs(k):
    d = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    assert Fr(d['k0']) == K0
    return [Fr(c)/K0 for c in d['psi%d_exact' % k]]


PSI2 = psi_coeffs(2)
PSI3 = psi_coeffs(3)


def lam():
    return A(LAM2).sqrt()


def polyval(coef, x):
    """Horner with exact rational coefficients; x arb or acb."""
    v = A(coef[-1]) if not isinstance(x, acb) else acb(A(coef[-1]))
    for a in reversed(coef[:-1]):
        v = v*x + A(a)
    return v


def _m(x, mid):
    if not mid:
        return x
    return x.mid() if isinstance(x, arb) else acb(x.real.mid(), x.imag.mid())


# ------------------------------------------------------------------ per-frequency jets
def seeds(p, nu, c):
    """(y0, y1) = d^p (F, c d F') at z = c, for real or complex nu (arb/acb); real c (arb) -> real parts if nu real."""
    d = 1 - c
    X = (1 + c)/d
    mu = acb(-0.5) + acb(0, 1)*acb(nu)
    P0 = acb(X).legendre_p(mu, 0, 3); Pm = acb(X).legendre_p(mu - 1, 0, 3)
    dP = mu*(X*P0 - Pm)/(X*X - 1)          # dP/dX
    sh = 2*c.sqrt()/d                      # sinh beta_c
    y0 = P0; y1 = p*c*P0 + c.sqrt()*dP*sh
    if isinstance(nu, arb):
        return y0.real, y1.real            # true values are real: real-part enclosure is rigorous
    return y0, y1


def jets(p, nu, c, K, mid=False):
    """y_0..y_K (arb if nu arb, else acb)."""
    y0, y1 = seeds(p, nu, c)
    d = 1 - c
    A1, A2, A3 = d - 2*c, c*c - 2*c*d, c*c*d
    q = (1 + 2*p)*c
    B0, B1, B2 = 1 - q, -q*d - c*(1 - q), c*q*d
    nu2 = nu*nu
    C0 = c*(nu2 + arb(0.25) - p + p*p*c); C1 = c*p*p*c*d
    y = [_m(y0, mid), _m(y1, mid)]
    for k in range(0, K - 1):
        acc = A1*(k*(k + 1))*y[k + 1] + B0*(k + 1)*y[k + 1] + C0*y[k]
        if k >= 1:
            acc += A2*((k - 1)*k)*y[k] + B1*k*y[k] + C1*y[k - 1]
        if k >= 2:
            acc += A3*((k - 2)*(k - 1))*y[k - 1] + B2*(k - 1)*y[k - 1]
        y.append(_m(-acc/((k + 2)*(k + 1)), mid))
    return y[:K + 1]


# ------------------------------------------------------------------ integrated jets, DV (numpy GL nodes)
def integrated_jets_mid(c, K, Np, which=('A', 'D')):
    import numpy as np
    xg, wg = np.polynomial.legendre.leggauss(Np)
    LAM = lam(); d = 1 - c
    out = {k: [arb(0)]*(K + 1) for k in which}
    for xn, wn in zip(xg, wg):
        p = arb((xn + 1)/2); wp = arb(wn/2)
        nu = LAM*p
        if 'A' in which:
            f = 2*polyval(PSI2, p)*wp
            y = jets(3, nu, c, K, mid=True)
            o = out['A']
            for j in range(K + 1):
                o[j] = (o[j] + f*y[j]).mid()
        if 'D' in which:
            f = LAM*LAM*polyval(PSI3, p)*wp/d
            y = jets(4, nu, c, K, mid=True)
            o = out['D']
            for j in range(K + 1):
                o[j] = (o[j] + f*y[j]).mid()
    return out


# ------------------------------------------------------------------ dressing, x_*, R recursion, Chebyshev
def chan(b):
    s = 2*b + 1; q_b = b*(b + 1); sig = 1 if b == 0 else -1; kb = arb(9)/2 - b*(b - 1); e = 7 - 4*b
    return s, q_b, sig, kb, e


def dress(YA, YD, c, b, mid=False):
    """Y^B_j = [t^j](D_b G_A + sigma_b G_D), j <= len - 3 (exact contraction, FULLPOS4 §8)."""
    s, q_b, sig, kb, e = chan(b)
    d = 1 - c
    out = []
    for j in range(len(YA) - 2):
        v = ((j + 2)*(j + 1)*YA[j + 2] + ((d - c)*j + e*c + 1)*(j + 1)*YA[j + 1] + (-c*d*(j*(j - 1)) + e*c*d*j)*YA[j])/(2*c*d) + kb*YA[j]
        out.append(_m(v + sig*YD[j], mid))
    return out


def xstar(c, s, a=998):
    d = 1 - c; a = arb(a)
    return 1/(a*a*(1 + c + (4*c + s*s*d*d/(a*a)).sqrt()))


def contract(c, xs, s, vecs, M, absnorm=False, mid=False, progress=None):
    """For each vector v in vecs: A_m(v) = sum_j R_{m,j}(c, xs) v_j, m = 0..M.  Returns list of (A, N) with
    N_m = sum_j |R_{m,j}| |v_j| (if absnorm) else None.  Ball arithmetic unless mid."""
    d = 1 - c; d2 = d*d; ic = 1/c
    nv = len(vecs)
    As = [[v[0]] for v in vecs]
    Ns = [[abs(v[0])] for v in vecs] if absnorm else None
    prev = [arb(1)]
    for m in range(1, M + 1):
        top = len(prev) - 1                                   # = 3(m-1) (prev has entries 0..3(m-1))
        Bv = [_m(prev[j + 1] - prev[j]*ic, mid) for j in range(top)] + [_m(-prev[top]*ic, mid)]
        cur = [arb(0)]*(3*m + 1)
        for j in range(0, 3*m):
            acc = d2*(j*(j + s))*Bv[j] if j <= top else arb(0)
            if 1 <= j <= top + 1:
                acc += d*(j*(2*j + s - 1))*Bv[j - 1]
            if 2 <= j <= top + 2:
                acc += (j*(j - 1))*Bv[j - 2]
            cur[j + 1] = _m(xs*acc, mid)
        for i, v in enumerate(vecs):
            tot = arb(0)
            for j in range(1, 3*m + 1):
                tot += cur[j]*v[j]
            As[i].append(_m(tot, mid))
            if absnorm:
                n = arb(0)
                for j in range(1, 3*m + 1):
                    n += abs(cur[j])*abs(v[j])
                Ns[i].append(n.mid())
        prev = cur
        if progress and m % 100 == 0:
            progress(m)
    return [(As[i], Ns[i] if absnorm else None) for i in range(nv)]


def poly_y(Am, Bm, c, xs, b, mid=False):
    """p_m: P(y) = sum A_m y^m + d^2 xs y sum (q_b A_m + B_m) y^m  (degree M+1)."""
    s, q_b, sig, kb, e = chan(b)
    d = 1 - c; dx = d*d*xs; M = len(Am) - 1
    pm = [Am[0]] + [_m(Am[m] + dx*(q_b*Am[m - 1] + Bm[m - 1]), mid) for m in range(1, M + 1)] + [_m(dx*(q_b*Am[M] + Bm[M]), mid)]
    return pm


def cheb(pm, mid=False):
    """Chebyshev coefficients on y in [0,1] (shifted T_l(2y-1)): exact binomial weights, ball arithmetic."""
    D = len(pm) - 1
    # 4^{-m} C(2m, m-l): build per m the row incrementally
    Ph = [arb(0)]*(D + 1)
    for m in range(D + 1):
        w = arb.bin_uiui(2*m, m)/arb(4)**m          # l = 0
        Ph[0] += w*pm[m]
        # C(2m, m-l) = C(2m, m-l+1) (m-l+1)/(m+l)
        for l in range(1, m + 1):
            w = w*(m - l + 1)/(m + l)
            Ph[l] += 2*w*pm[m]
    return [_m(v, mid) for v in Ph]
