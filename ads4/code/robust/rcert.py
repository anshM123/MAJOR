"""Certificates (FLINT/Arb, exact rationals) for the robust functional (ROBUST_FUNCTIONAL.md section 4).

Generalised heavy action (normal form, u = 1/m^2, JJ = J(J+1)):
    FF[g](J,u) = int_0^1 (W0(p) + u W1(p)) D_h(1 - 2 u p^2) dp + (2 aphi JJ - cphi) u^3 .
The plain functional (psi2, psi3, a):  W0 = 2 psi2, W1 = -p^2 psi2 + sigma psi3, aphi = a, cphi = a (c4 - sigma).
Robust margin delta:  FF - delta N >= 0 with N = |B2| + |B3| + |FW|  <=>  for all s2, s3 in {+1,-1}:
    FF[g_{s2,s3}] >= 0,  W0 = 2(1-delta s2) psi2,  W1 = -(1-delta s2) p^2 psi2 + (1-delta s3) sigma psi3,
    aphi = (1-delta sF) a, cphi = (1-delta sF) a (c4 - sigma),  sF = sign(2JJ - c4 + sigma)  (= +1 except +-, J = 4: -1).
Regions (as in FLAT_GRAVITON_BOUND.md section 3, generalised): P (J < JP, all u in [0,1]); C (J >= JP, u >= 1/4);
M (u <= 1/4, bhat >= BA); S (J >= JP, u <= 1/4, bhat <= BA, incl. bhat <= 1).
Usage: python3 rcert.py NAME DELTA REGION [args]   (REGION in P, C, M, S++, S+-, low)
"""
import sys, os, json, time, importlib.util
from fractions import Fraction as Fr
from math import factorial, comb
import flint

HERE = os.path.dirname(os.path.abspath(__file__))
FB = os.path.join(HERE, '..', 'flatbound')
def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod
fbc = _load('fb_common', os.path.join(FB, 'common.py'))          # cg_abs_bound, cg_coeffs, to_fmpz_poly (unchanged, imported)
A = flint.arb
SIGMA = {'++': 1, '+-': -1}
C4 = {'++': 1, '+-': 41}

# ------------------------------------------------------------------ exact polynomials (Fraction lists)
def pmul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else Fr(0)) + (b[i] if i < len(b) else Fr(0)) for i in range(n)]
def pscale(a, c): return [c * x for x in a]
def mom(a, j): return sum(Fr(c) / (i + j + 1) for i, c in enumerate(a))
def pval(a, x):
    r = 0
    for c in reversed(a): r = r * x + c
    return r
def pder(a):
    return [i * a[i] for i in range(1, len(a))] or [Fr(0)]

CERTIFIED = os.path.join(HERE, 'functionals_certified.json')   # tracked; never written by the LP stage

def fingerprint(psi2, psi3, a):
    import hashlib
    canon = 'psi2=' + ','.join(str(c) for c in psi2) + ';psi3=' + ','.join(str(c) for c in psi3) + ';a=' + str(a)
    return hashlib.sha256(canon.encode()).hexdigest()

def load_functional(name):
    """Exact coefficients: from functionals_certified.json (tracked) if NAME is there, else from cache/functionals.json.
    Prints the source and a sha256 fingerprint of the exact coefficients."""
    src = CERTIFIED
    D = json.load(open(CERTIFIED))
    if name not in D:
        src = os.path.join(HERE, 'cache', 'functionals.json'); D = json.load(open(src))
    d = D[name]
    F = ([Fr(s) for s in d['psi2_exact']], [Fr(s) for s in d['psi3_exact']], Fr(d['a_exact']))
    # consistency of psi2/psi3 with x, y when present
    if 'x' in d:
        omp = [Fr(1), Fr(-3), Fr(3), Fr(-1)]
        assert pmul(omp, [Fr(0)] + [Fr(s) for s in d['x']]) == F[0] and pmul(omp, [Fr(s) for s in d['y']]) == F[1]
    print('FUNCTIONAL %s from %s ; sha256(exact coefficients) = %s' % (name, os.path.relpath(src, HERE), fingerprint(*F)), flush=True)
    return F

class Gen:
    def __init__(self, h, W0, W1, aphi, cphi, label):
        self.h, self.W0, self.W1, self.aphi, self.cphi, self.label = h, W0, W1, aphi, cphi, label

def gens(psi2, psi3, a, delta, h, sF):
    s = SIGMA[h]; out = []
    pats = [(1, 1), (1, -1), (-1, 1), (-1, -1)] if delta != 0 else [(0, 0)]
    for s2, s3 in pats:
        W0 = pscale(psi2, 2 * (1 - delta * s2))
        W1 = padd(pscale([Fr(0), Fr(0)] + psi2, -(1 - delta * s2)), pscale(psi3, s * (1 - delta * s3)))
        aphi = (1 - delta * sF) * a
        out.append(Gen(h, W0, W1, aphi, aphi * (C4[h] - s), '%s s2=%+d s3=%+d sF=%+d' % (h, s2, s3, sF)))
    return out

def e_h(h, i): return i * (i - 1) if h == '++' else (i + 3) * (i + 4)

# ------------------------------------------------------------------ region P
def poly_u(g, J):
    JJ = J * (J + 1); kmax = J if g.h == '++' else J - 4
    co = {}; prod = Fr(1)
    for k in range(0, kmax + 1):
        if k > 0: prod *= Fr(-(JJ - e_h(g.h, k)), k * k)
        co[k] = co.get(k, 0) + prod * mom(g.W0, 2 * k)
        co[k + 1] = co.get(k + 1, 0) + prod * mom(g.W1, 2 * k)
    co[3] = co.get(3, 0) + 2 * g.aphi * JJ - g.cphi
    return [Fr(co.get(k, 0)) for k in range(max(co) + 1)]

def sign_var(coeffs):
    s = [c for c in coeffs if c != 0]
    return sum(1 for a, b in zip(s, s[1:]) if (a > 0) != (b > 0))
def interval_poly(P, k, s):
    d = P.degree(); c = P.coeffs()
    x = flint.fmpz_poly([k, 1]); Q = flint.fmpz_poly([c[d]])
    for j in range(d - 1, -1, -1):
        Q = Q * x + flint.fmpz_poly([c[j] * (2 ** (s * (d - j)))])
    return Q
def descartes01(Q):
    d = Q.degree(); c = list(Q.coeffs()) + [0] * (d + 1 - len(Q.coeffs()))
    R = flint.fmpz_poly(list(reversed(c)))(flint.fmpz_poly([1, 1]))
    return sign_var(R.coeffs())
def positive_on_01(P, maxdepth=60):
    """(copied from flatbound/cert_P.py) exact Descartes/bisection: P > 0 on [0,1]."""
    if P(0) <= 0: return False, 0
    stack = [(0, 0)]; n = 0
    while stack:
        k, s = stack.pop(); n += 1
        Q = interval_poly(P, k, s)
        if Q(0) == 0 or Q(1) == 0: return False, n
        if descartes01(Q) == 0: continue
        if s >= maxdepth: return False, n
        stack += [(2 * k, s + 1), (2 * k + 1, s + 1)]
    return True, n

def region_P(F, delta, JP):
    psi2, psi3, a = F
    ok_all = True; cnt = 0; t0 = time.time()
    for h in ['++', '+-']:
        Js = range(0, JP, 2) if h == '++' else range(4, JP)
        for J in Js:
            sF = 1 if (2 * J * (J + 1) - C4[h] + SIGMA[h]) >= 0 else -1
            for g in gens(psi2, psi3, a, delta, h, sF):
                P, L = fbc.to_fmpz_poly(poly_u(g, J))
                ok, n = positive_on_01(P)
                cnt += 1
                if not ok or J < 8 or J % 50 == 0:
                    print('P %s J=%d deg=%d intervals=%d %s' % (g.label, J, P.degree(), n, 'PASS' if ok else 'FAIL'), flush=True)
                ok_all &= ok
    print('VERDICT rcert_P (delta=%s, JP=%d, %d polynomials, %.0fs): %s' % (delta, JP, cnt, time.time() - t0, 'PASS' if ok_all else 'FAIL'))
    return ok_all

# ------------------------------------------------------------------ arb helpers
def fq(x):
    x = Fr(x); return A(flint.fmpq(x.numerator, x.denominator))
def iv(lo, hi): return fq(lo).union(fq(hi))
def ev(poly, x):
    r = A(0)
    for c in reversed(poly): r = r * x + fq(c)
    return r
def integral_upper(f, a, b, ncell):
    tot = A(0)
    for i in range(ncell):
        lo = Fr(a) + (Fr(b) - Fr(a)) * i / ncell; hi = Fr(a) + (Fr(b) - Fr(a)) * (i + 1) / ncell
        tot += A(f(iv(lo, hi)).upper()) * fq(hi - lo)
    return tot.upper()
def sup_on(f, a, b, ncell):
    best = None
    for i in range(ncell):
        lo = Fr(a) + (Fr(b) - Fr(a)) * i / ncell; hi = Fr(a) + (Fr(b) - Fr(a)) * (i + 1) / ncell
        v = f(iv(lo, hi)).upper()
        best = v if best is None else best.max(v)
    return best
def div_omp3(poly):
    """poly / (1-p)^3 (exact; asserts divisibility)."""
    q = list(poly); out = []
    # synthetic division by (1-p) three times: poly = (1-p) r  <=>  r_k = sum_{i<=k} poly_i
    for _ in range(3):
        r = []; acc = Fr(0)
        for c in q[:-1]:
            acc += c; r.append(acc)
        assert acc + q[-1] == 0, 'not divisible by (1-p)'
        q = r
    return q

# ------------------------------------------------------------------ region C
def region_C(F, delta, JP):
    flint.ctx.prec = 200
    psi2, psi3, a = F
    ok_all = True
    for h in ['++', '+-']:
        for g in gens(psi2, psi3, a, delta, h, 1):
            G0 = div_omp3(g.W0); G1 = div_omp3(g.W1)
            gh = lambda x: A(ev(G0, x).abs_upper()) + A(ev(G1, x).abs_upper())     # |W0|+|W1| <= (1-p)^3 gh
            Phi = 2 * fq(g.aphi) * JP * (JP + 1) - fq(g.cphi)
            if h == '++':
                tot = integral_upper(lambda x: gh(x) * (1 - x) ** 3, 0, 1, 4000)
                ok = bool(Phi > 64 * A(tot)) and bool(4 * fq(g.aphi) * JP > 0)
                print('C %s J>=%d: 2aJJ-c = %s > 64 int(|W0|+|W1|) <= %s : %s' % (g.label, JP, Phi.str(6), (64 * A(tot)).str(6), ok))
            else:
                c_a = integral_upper(lambda x: gh(x) / ((1 - x) * (1 + x) ** 4), 0, Fr(1, 2), 2000)
                gam = sup_on(lambda x: gh(x) / (1 + x) ** 4, Fr(1, 2), 1, 4000)
                g1 = sup_on(gh, Fr(1, 2), 1, 4000)
                Bn = A(comb(JP + 4, 8)); assert bool(Bn > A(4) ** 4 / A(3) ** 4)
                U = A(c_a) + A(gam) * Bn.log() / 4 + A(g1) / (4 * (A(3) / 2) ** 4)
                ok1 = bool(Phi > 64 * U)
                mono = bool(2 * fq(g.aphi) * (2 * JP + 1) > 16 * A(gam) * 8 / A(JP - 3))
                ok = ok1 and mono
                print('C %s J>=%d: U<=%s, 2aJJ-c=%s > 64U=%s: %s ; monotone %s' % (g.label, JP, U.str(5), Phi.str(6), (64 * U).str(6), ok1, mono))
            ok_all &= ok
    print('VERDICT rcert_C (delta=%s): %s' % (delta, 'PASS' if ok_all else 'FAIL'))
    return ok_all

# ------------------------------------------------------------------ region M (Lemma M of FLAT_GRAVITON_BOUND.md, unchanged)
ORD = 5
class Jet:
    __slots__ = ['c']
    def __init__(self, c): self.c = list(c) + [A(0)] * (ORD - len(c))
    def __add__(self, o):
        o = o if isinstance(o, Jet) else Jet([o]); return Jet([x + y for x, y in zip(self.c, o.c)])
    __radd__ = __add__
    def __sub__(self, o):
        o = o if isinstance(o, Jet) else Jet([o]); return Jet([x - y for x, y in zip(self.c, o.c)])
    def __rsub__(self, o): return Jet([o]) - self
    def __mul__(self, o):
        if not isinstance(o, Jet): return Jet([x * o for x in self.c])
        r = [A(0)] * ORD
        for i in range(ORD):
            for j in range(ORD - i): r[i + j] += self.c[i] * o.c[j]
        return Jet(r)
    __rmul__ = __mul__
    def inv(self):
        a0 = self.c[0]; r = [A(0)] * ORD; r[0] = 1 / a0
        for n in range(1, ORD):
            s = A(0)
            for k in range(1, n + 1): s += self.c[k] * r[n - k]
            r[n] = -s / a0
        return Jet(r)
    def der(self, k):
        f = 1
        for i in range(2, k + 1): f *= i
        return self.c[k] * f
def poly_jet(coeffs, pj):
    r = Jet([A(0)])
    for c in reversed(coeffs): r = r * pj + Jet([fq(c)])
    return r
def weight_jet(g, p0, u):
    pj = Jet([p0, A(1)])
    w = poly_jet(g.W0, pj) + poly_jet(g.W1, pj) * u
    if g.h == '+-':
        q = (1 - u * pj * pj).inv(); q2 = q * q; q4 = q2 * q2; w = w * (q4 * q4)
    return w
def CR_bound(g, ulo, uhi, NP):
    """Lemma M remainder constant C_R(w,u), sup over u in [ulo,uhi] (interval arithmetic, p-cells 1/NP)."""
    u = fq(ulo).union(fq(uhi))
    n_w4 = A(0); n_pw3 = A(0); n_w2 = A(0); s2 = A(0); s1 = A(0)
    width = fq(Fr(1, NP))
    for i in range(NP):
        p0 = iv(Fr(i, NP), Fr(i + 1, NP))
        w = weight_jet(g, p0, u)
        wd = [w.der(k) for k in range(5)]
        n_w4 += A(wd[4].abs_upper()) * width
        n_pw3 += A((p0 * wd[3]).abs_upper()) * width
        n_w2 += A(wd[2].abs_upper()) * width
        s2 = s2.max(A(wd[2].abs_upper())); s1 = s1.max(A(wd[1].abs_upper()))
    w3at1 = A(weight_jet(g, A(1), u).der(3).abs_upper())
    uu = fq(uhi)
    CR = 2 / A.pi() * (A(2) / 3 * (w3at1 + n_w4) + 6 * uu * n_pw3 + 4 * uu * n_w2) + uu * uu.sqrt() * (3 * A.pi() / 4 * s2 + s1)
    return A(CR.upper())

def region_M(F, delta, BA, u1, NP):
    flint.ctx.prec = 128
    psi2, psi3, a = F
    b = fq(BA); ok_all = True
    neg = lambda x: (-x) if x < 0 else Fr(0)
    for h in ['++', '+-']:
        for g in gens(psi2, psi3, a, delta, h, 1):
            assert g.W0[0] == 0
            alpha = g.W1[0]
            q0 = 2 * (g.W0[2] if len(g.W0) > 2 else 0)
            q1 = 2 * (g.W1[2] if len(g.W1) > 2 else 0) + (16 * g.W0[0] if h == '+-' else 0)
            q2 = 16 * g.W1[0] if h == '+-' else Fr(0)
            gam0, gam1, gam2 = -q0 / 2, -q1 / 2, -q2 / 2
            assert gam0 > 0, 'leading 1/b^3 coefficient not positive'
            CR1 = CR_bound(g, 0, u1, NP); CRm = CR_bound(g, u1, Fr(1, 4), NP)
            # u in [0, u1]
            dl = 8 * fq(u1).sqrt() / b if h == '+-' else A(0)
            bm = b * (1 - dl)
            pi1 = fq(neg(alpha)) + fq(neg(gam1)) / bm ** 2
            aph = fq(g.aphi)
            eta = (fq(max(g.aphi / 2 + g.cphi, 0)) * fq(u1) + fq(neg(gam2)) / bm ** 3) / (aph / 2 * b ** 2)
            val = fq(gam0) * b / (1 + dl) ** 3 - CR1 / (1 - dl) ** 4 - pi1 ** 2 / (2 * aph * (1 - eta) * (1 - dl) ** 2)
            c1 = bool(val > 0) and bool(eta < 1)
            # u in [u1, 1/4]: b_- >= b - 4 (+-)
            bm2 = b - 4 if h == '+-' else b
            pi1m = fq(neg(alpha)) + fq(neg(gam1)) / bm2 ** 2
            etam = (fq(max(g.aphi / 2 + g.cphi, 0)) / 4 + fq(neg(gam2)) / bm2 ** 3) / (aph / 2 * b ** 2)
            A2 = aph / 2 * b ** 2 * (1 - etam); A1 = pi1m / bm2; A0 = CRm / bm2 ** 4
            uu = fq(u1)
            c2 = bool(2 * A2 * uu > A1) and bool(A2 * uu ** 2 - A1 * uu - A0 > 0) and bool(etam < 1)
            ok = c1 and c2
            print('M %s: gamma0=%.5g alpha=%.5g gamma1=%.5g gamma2=%.4g aphi=%.4g; C_R(u<=u1)<=%s C_R(u1..1/4)<=%s; b^4 FF >= %s (u<=u1) %s ; u in [u1,1/4] %s'
                  % (g.label, float(gam0), float(alpha), float(gam1), float(gam2), float(g.aphi), CR1.str(5), CRm.str(5), val.str(5), c1, c2), flush=True)
            ok_all &= ok
    print('VERDICT rcert_M (delta=%s, BA=%s, u1=%s): %s' % (delta, BA, u1, 'PASS' if ok_all else 'FAIL'))
    return ok_all

# ------------------------------------------------------------------ region S (adapted from flatbound/cert_S.py)
PREC_S = 768; NMAX = 420; KS = 12; IS = 9
def build_ak(g):
    lam = [None] + [Fr(1, 4) + e_h(g.h, i) for i in range(1, NMAX + 1)]
    E = [Fr(1)] + [Fr(0)] * KS
    coef = [dict() for _ in range(KS)]
    t0 = [mom(g.W0, 2 * n) for n in range(NMAX + 1)]
    t1 = [mom(g.W1, 2 * n) for n in range(NMAX + 1)]
    fact2 = Fr(1)
    for n in range(0, NMAX + 1):
        if n > 0:
            fact2 *= n * n
            for k in range(min(n, KS), 0, -1): E[k] = E[k] + lam[n] * E[k - 1]
        sgn = -1 if n % 2 else 1
        for k in range(KS):
            if k <= n:
                coef[k][n - k] = coef[k].get(n - k, 0) + Fr(sgn * (-1) ** k) * E[k] * t0[n] / fact2
            if 1 <= k <= n + 1:
                coef[k][n - k + 1] = coef[k].get(n - k + 1, 0) + Fr(sgn * (-1) ** (k - 1)) * E[k - 1] * t1[n] / fact2
    polys = []
    for k in range(KS):
        deg = 2 * max(coef[k]); cl = [A(0)] * (deg + 1)
        for j, c in coef[k].items(): cl[2 * j] = fq(c / 4 ** j)
        polys.append(flint.arb_poly(cl))
    return polys

def tmax_ok(g, Tmax):
    i0 = integral_upper(lambda x: A(ev(g.W0, x).abs_upper()), 0, 1, 2000)
    i1 = integral_upper(lambda x: A(ev(g.W1, x).abs_upper()), 0, 1, 2000)
    return bool(A(i0) <= Tmax) and bool(A(i1) <= Tmax)

def tail_bound(g, Xbar, Tmax):
    lamsum = lambda n: sum(Fr(1, 4) + e_h(g.h, i) for i in range(1, n + 1))
    def term(n):
        Ln = fq(lamsum(n))
        return Tmax * (fq(Xbar) + 1) ** (n + 1) * 2 * (Ln + 1) ** KS / A(factorial(n)) ** 2
    n0 = NMAX + 1
    ratio = (fq(Xbar) + 1) * ((fq(lamsum(n0 + 1)) + 1) / (fq(lamsum(n0)) + 1)) ** KS / A(n0 + 1) ** 2
    assert bool(ratio < A(1) / 2), 'tail ratio'
    return 2 * term(n0)

def s_fun(x):
    if bool(x == 0): return A(1)
    sx = x.sqrt(); return sx.asin() / sx

def cauchy_M(g, bc, r, rho, NPC=160):
    tot = A(0)
    for i in range(NPC):
        plo, phi = Fr(i, NPC), Fr(i + 1, NPC)
        p = iv(plo, phi); ph = fq(phi)
        Wabs = A(ev(g.W0, p).abs_upper()) + rho * A(ev(g.W1, p).abs_upper())
        x = rho * ph * ph; sx = s_fun(x)
        Y = r * ph * sx + (bc + r) * ph * (sx - 1)
        th = 2 * (ph * rho.sqrt()).asin(); fac = th / th.sin()
        if g.h == '+-':
            Wabs = Wabs / (1 - x) ** 8; Y = Y + 4 * th
        tot += A((Wabs * Y.cosh() * fac).upper()) / NPC
    if g.h == '+-':
        Zmin = (bc - r) / (2 * rho.sqrt())
        Zlo = Fr(int((Zmin.lower() * 1000).floor().unique_fmpz()), 1000)
        assert Zlo > 5
        tot = tot * fq(fbc.cg_abs_bound(Zlo))
    return A(tot.upper())

def eval2(gm, t, u):
    r = A(0)
    for i in reversed(range(len(gm))):
        ri = A(0)
        for k in reversed(range(len(gm[i]))): ri = ri * u + gm[i][k]
        r = r * t + ri
    return r
def deriv_t(gm): return [[gm[i][k] * i for k in range(len(gm[i]))] for i in range(1, len(gm))] or [[A(0)]]
def deriv_u(gm): return [[gm[i][k] * k for k in range(1, len(gm[i]))] or [A(0)] for i in range(len(gm))]
def certify_box(gm, gt, gu, tail_fn, t0, t1, u0, u1, depth, stats):
    tm = (t0 + t1) / 2; um = (u0 + u1) / 2
    T = iv(t0, t1); U = iv(u0, u1)
    val = eval2(gm, fq(tm), fq(um)); dT = eval2(gt, T, U); dU = eval2(gu, T, U)
    low = val - A(dT.abs_upper()) * fq((t1 - t0) / 2) - A(dU.abs_upper()) * fq((u1 - u0) / 2)
    if bool(low - tail_fn(u1) > 0):
        stats['boxes'] += 1; return True
    if depth > 26: return False
    lt = dT.abs_upper() * fq((t1 - t0) / 2); lu = dU.abs_upper() * fq((u1 - u0) / 2)
    if bool(A(lt) > A(lu)):
        return certify_box(gm, gt, gu, tail_fn, t0, tm, u0, u1, depth + 1, stats) and certify_box(gm, gt, gu, tail_fn, tm, t1, u0, u1, depth + 1, stats)
    return certify_box(gm, gt, gu, tail_fn, t0, t1, u0, um, depth + 1, stats) and certify_box(gm, gt, gu, tail_fn, t0, t1, um, u1, depth + 1, stats)

def cell(g, ak, bc, hh, JP, Tmax, kappa=3):
    r = min(Fr(2), bc / 2); assert hh <= r / 10
    u0 = ((bc + hh) / (2 * JP + 1)) ** 2
    rho = fq(kappa * u0); rr = fq(r)
    M = cauchy_M(g, fq(bc), rr, rho)
    tb = tail_bound(g, ((bc + r) ** 2) / 4, Tmax)
    gm = [[None] * KS for _ in range(IS)]
    for k in range(KS):
        pk = ak[k]; fi = A(1)
        for i in range(IS):
            err = A((tb / rr ** i).upper())
            gm[i][k] = pk(fq(bc)) / fi + (-err).union(err)
            pk = pk.derivative(); fi = fi * (i + 1)
    # Phi part: (2 aphi JJ - cphi) u^3 = 2 aphi X u^2 - (aphi/2 + cphi) u^3,  X = (bc + t)^2/4
    a2 = 2 * g.aphi
    gm[0][2] += fq(a2 * bc * bc / 4); gm[1][2] += fq(a2 * 2 * bc / 4); gm[2][2] += fq(a2 / 4)
    gm[0][3] -= fq(g.aphi / 2 + g.cphi)
    q = fq(hh) / rr
    def tail_fn(ub):
        x = fq(ub) / rho
        return M * (q ** IS + x ** KS) / ((1 - q) * (1 - x))
    stats = {'boxes': 0}
    ok = certify_box(gm, deriv_t(gm), deriv_u(gm), tail_fn, -hh, hh, Fr(0), u0, 0, stats)
    return ok, stats['boxes'], M, u0

def small_b(g, BS, JP, Tmax):
    X = fq(Fr(BS) ** 2 / 4); u = fq((Fr(BS) / (2 * JP + 1)) ** 2)
    low = fq(mom(g.W0, 0)) - fq(abs(mom(g.W1, 0))) * u
    prod = A(1); tail = A(0)
    for n in range(1, 400):
        prod = prod * (X + fq(Fr(1, 4) + e_h(g.h, n)) * u) / (n * n)
        tail += prod * (fq(abs(mom(g.W0, 2 * n))) + fq(abs(mom(g.W1, 2 * n))) * u)
    for n in range(400, 2000):
        assert bool((X + fq(Fr(1, 4) + e_h(g.h, n + 1)) * u) / (n + 1) ** 2 < A(1) / 2)
    tail += 2 * prod * Tmax
    low = low - tail - fq(max(g.aphi / 2 + g.cphi, 0)) * u ** 3
    return bool(low > 0), low

def region_S(F, delta, h, BS, BA, JP, Tmax=Fr(100)):
    flint.ctx.prec = PREC_S
    psi2, psi3, a = F
    ok_all = True; t0 = time.time(); ncell = 0
    for g in gens(psi2, psi3, a, delta, h, 1):
        flint.ctx.prec = 200
        assert tmax_ok(g, fq(Tmax)), 'Tmax'
        flint.ctx.prec = PREC_S
        ok, low = small_b(g, BS, JP, fq(Tmax))
        print('S %s small-b (b<=%s): lower bound %s %s' % (g.label, BS, low.str(6), 'PASS' if ok else 'FAIL'), flush=True)
        ok_all &= ok
        ak = build_ak(g)
        b = Fr(BS); fails = 0
        while b < BA:
            r = min(Fr(2), (b + Fr(1, 10)) / 2); hh = min(r / 10, Fr(1, 10)); bc = b + hh
            while hh > min(Fr(2), bc / 2) / 10:
                hh = hh / 2; bc = b + hh
            if bc + hh > BA:
                hh = (BA - b) / 2; bc = b + hh
            ok, nbox, M, u0 = cell(g, ak, bc, hh, JP, fq(Tmax)); ncell += 1
            if not ok:
                fails += 1
                print('S %s cell b in [%s, %s] umax=%.3e M<=%s FAIL' % (g.label, float(b), float(b + 2 * hh), float(u0), M.str(4)), flush=True)
            ok_all &= ok
            b = b + 2 * hh
        print('S %s: cells done, fails=%d (%.0fs)' % (g.label, fails, time.time() - t0), flush=True)
    print('VERDICT rcert_S %s (delta=%s, b in [%s,%s], %d cells, %.0fs): %s' % (h, delta, BS, BA, ncell, time.time() - t0, 'PASS' if ok_all else 'FAIL'))
    return ok_all

# ------------------------------------------------------------------ low-energy side (exact)
def low_side(F):
    """k, A, B exact; matter coefficients <= 0 on (0,1] (exact Sturm via sympy)."""
    import sympy as sp
    psi2, psi3, a = F
    k = -2 * mom(psi2, 6) + 2 * mom(psi3, 4) + 4 * a
    x0 = psi2[1]
    c = list(psi2); c[1] -= x0
    assert c[0] == 0
    Lfin = sum(c[i] / (i - 1) for i in range(2, len(c)))
    Aex = 16 * x0 / k; Bex = 16 * Lfin / k
    m = sp.Symbol('m', positive=True)
    S = lambda q: sp.Rational(q.numerator, q.denominator)
    C0 = S(mom(psi2, 2) - mom(psi3, 0)) - 2 * m ** 2 * S(mom(psi2, 0))
    C2m4 = (sum(S(cc) * sp.integrate(sp.Symbol('p') ** i * (sp.Symbol('p') ** 2 - 2 * m ** 2) * (m ** 4 - 6 * m ** 2 * sp.Symbol('p') ** 2 + 6 * sp.Symbol('p') ** 4), (sp.Symbol('p'), 0, 1)) for i, cc in enumerate(psi2))
            - sum(S(cc) * sp.integrate(sp.Symbol('p') ** i * (m ** 4 - 6 * m ** 2 * sp.Symbol('p') ** 2 + 6 * sp.Symbol('p') ** 4), (sp.Symbol('p'), 0, 1)) for i, cc in enumerate(psi3))
            - 12 * S(a))
    ok = True
    for nm, e in [('C0', C0), ('m^4 C2', C2m4)]:
        P = sp.Poly(sp.expand(e), m)
        nr = sp.count_roots(P, 0, 1)
        sgn = sp.sign(e.subs(m, sp.Rational(1, 2)))
        neg = nr == 0 and sgn < 0 and P.eval(0) <= 0 and P.eval(1) < 0
        print('low: %s = %s ; roots in [0,1]: %d ; negative on (0,1]: %s' % (nm, sp.expand(e), nr, neg))
        ok &= bool(neg)
    print('low: k = %s = %.10f ; A = 16 x0/k = %.6f ; B = 16 Lfin/k = %.6f ; B_AdS = B - A log(3/2) = %.6f' % (
        k, float(k), float(Aex), float(Bex), float(Bex) - float(Aex) * 0.4054651081081644))
    print('VERDICT rcert_low: %s' % ('PASS' if ok else 'FAIL'))
    return ok

if __name__ == '__main__':
    name, delta, region = sys.argv[1], Fr(sys.argv[2]), sys.argv[3]
    F = load_functional(name)
    if region == 'P': region_P(F, delta, int(sys.argv[4]) if len(sys.argv) > 4 else 400)
    elif region == 'C': region_C(F, delta, int(sys.argv[4]) if len(sys.argv) > 4 else 400)
    elif region == 'M': region_M(F, delta, Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(120), Fr(sys.argv[5]) if len(sys.argv) > 5 else Fr(1, 100),
                                 int(sys.argv[6]) if len(sys.argv) > 6 else 4000)
    elif region.startswith('S'): region_S(F, delta, region[1:], Fr(1), Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(120), int(sys.argv[5]) if len(sys.argv) > 5 else 400)
    elif region == 'low': low_side(F)
