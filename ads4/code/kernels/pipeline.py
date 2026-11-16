"""
Energy-flux matrix elements  <O(p),m| E(theta) |T(p),h>  from position-space three-point structures, d = 3.

Chain (all exact, polynomial arithmetic over Q[i] with python-flint):
  structure S(X1,Z1;X3,Z3;X4,Z4)  --(detector limit)-->  X1 = (0,-alpha,z), Z1 = (0,-1,0), z = (1,n):
       E(n) = lim r int dt n^i T^0_i = 2 int dalpha S(...)                      [= 2 L[T](infty,z), KKSZ (3.18)]
  --(alpha integral, Wightman i eps)-->
       int dalpha X14^{-a} X13^{-b} = -2 pi i e^{i pi a} Gamma(a+b-1)/(Gamma(a)Gamma(b)) (2 z.x + i eps)^{1-a-b},
       X14 = -(alpha + 2 z.x + i eps), X13 = -(alpha - i eps),  x = x4 - x3
  --(Fourier transform)-->
       int d^3x e^{-ip.x} (2 z.x + i eps)^{-A} (x^2)^{-B} = e^{-i pi A/2} pi^{5/2} 2^{4-2B}/(Gamma(B)Gamma(A+B-1/2))
                                                         * (-p^2)^{A+B-3/2} (-4 z.p)^{-A}     (p in forward cone)
       and  x^mu -> i d/dp_mu  for polynomial prefactors.
  Rest frame p = (1,0,0);  bra polarisation zeta' = conj zeta(a',b') = zeta(b',a') with a' = 1, b' = t, so that the
  charge-m component of <O| is the coefficient of t^{J-m}.

Normalisation: M_m(theta) := int d^3x e^{-ip.x} <0| O_m(x)^dagger E(n) T_h(0) |0>,
               <O_m(p')|E(n)|T_h(p)> = (2 pi)^3 delta^3(p-p') M_m(theta),   M_m(theta) = e^{i(h-m) theta} M_m(0).
"""
from flint import fmpq, fmpq_mpoly_ctx, Ordering
import math

CTX = fmpq_mpoly_ctx.get(('t', 'd', 's', 'I'), Ordering.lex)
t, dd, s, I = CTX.gens()
C = lambda n, m=1: CTX.constant(fmpq(n, m))
ONE, ZERO = C(1), C(0)

def redI(p):
    """reduce I^2 = -1"""
    out = {}
    for (et, ed, es, eI), c in p.to_dict().items():
        sign = -1 if (eI // 2) % 2 else 1
        key = (et, ed, es, eI % 2)
        out[key] = out.get(key, fmpq(0)) + sign*c
    return CTX.from_dict({k: v for k, v in out.items() if v != 0})

def dot(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def det3(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))
def covec(v):          # covector of a vector
    return [-v[0], v[1], v[2]]
def eps_covec(a, b):   # w_mu with eps(a,b,x) = w_mu x^mu
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def pair(w1, w2):      # eta^{mu nu} w1_mu w2_nu
    return -w1[0]*w2[0] + w1[1]*w2[1] + w1[2]*w2[2]
def cov_on(w, v):      # w_mu v^mu
    return w[0]*v[0] + w[1]*v[1] + w[2]*v[2]

E_P = [ZERO, ONE, I]; E_M = [ZERO, ONE, -I]; E_0 = [ONE, ZERO, ZERO]
def zeta_vec(a, b):
    return [a*b, (a*a + b*b)*C(1, 2), I*(a*a - b*b)*C(1, 2)]     # a^2 e_+/2 + b^2 e_-/2 + a b e_0

class Kin:
    """rest-frame kinematics. source: 'plus' (z3 = e_+), 'minus' (z3 = e_-), 'generic' (z3 = zeta(1,d));
       n = (c, sn) detector direction (rational point of the circle)."""
    def __init__(self, source='plus', n=(1, 0)):
        self.p = [ONE, ZERO, ZERO]
        self.z = [ONE, C(*_fr(n[0])), C(*_fr(n[1]))]
        self.z3 = {'plus': E_P, 'minus': E_M, 'generic': zeta_vec(ONE, dd)}[source]
        self.zeta = zeta_vec(t, ONE)                      # zeta(b',a') with a' = 1, b' = t
        self.zz3 = redI(dot(self.z, self.z3)); self.zze = redI(dot(self.z, self.zeta))
        self.z3ze = redI(dot(self.z3, self.zeta)); self.E0 = redI(det3(self.z3, self.zeta, self.z))
        # x-linear forms: 1: z3.x  2: zeta.x  pseudo: 'e1' eps(z3,z,x), 'e2' eps(zeta,z,x), 'e3' eps(z3,zeta,x)
        self.forms = {1: covec(self.z3), 2: covec(self.zeta), 'e1': eps_covec(self.z3, self.z),
                      'e2': eps_covec(self.zeta, self.z), 'e3': eps_covec(self.z3, self.zeta)}
        self._tab = {}

    # ---- scri values of the building blocks, as {(i,j,k,l,n1,n2,eps): coeff}:
    #      P3^i P4^j (x^2)^k (z.x)^l (z3.x)^n1 (zeta.x)^n2 [eps form], P3 = alpha, P4 = alpha + 2 z.x
    def blocks(self):
        return dict(
            V1={(0, 0, -1, 1, 0, 0, 0): C(-1)},
            V3={(1, -1, 0, 0, 1, 0, 0): ONE, (0, -1, 1, 0, 0, 0, 0): self.zz3},
            V4={(-1, 1, 0, 0, 0, 1, 0): ONE, (-1, 0, 1, 0, 0, 0, 0): -self.zze},
            H13={(0, 0, 0, 0, 0, 0, 0): self.zz3},
            H14={(0, 0, 0, 0, 0, 0, 0): self.zze},
            H34={(0, 0, 1, 0, 0, 0, 0): self.z3ze, (0, 0, 0, 0, 1, 1, 0): C(-2)},
            e13={(1, 0, 0, 0, 0, 0, 1): C(-1)},                                   # eps13 * X13
            e14={(0, 1, 0, 0, 0, 0, 2): C(-1)},                                   # eps14 * X14
            e34={(0, 0, 1, 0, 0, 1, 1): C(2), (1, 0, 1, 0, 0, 0, 3): ONE, (0, 0, 2, 0, 0, 0, 0): self.E0})  # eps34*X34

    # ---- derivative table: C_{j,i} = [tau^n] (-X)^j M^i,  X = 2 sum tau_k L_k(p) + sum tau_k tau_l <L_k,L_l>,
    #      M = sum tau_k L_k(z);   u(p+W) = 1 - X, v(p+W) = 4 (1 - M)
    def table(self, n1, n2, eps):
        key = (n1, n2, eps)
        if key in self._tab:
            return self._tab[key]
        names = [1, 2] + (['e%d' % eps] if eps else [])
        n = (n1, n2) + ((1,) if eps else ())
        K = len(names)
        W = [self.forms[nm] for nm in names]
        Lp = [redI(cov_on(w, self.p)) for w in W]; Lz = [redI(cov_on(w, self.z)) for w in W]
        unit = lambda k: tuple(1 if q == k else 0 for q in range(K))
        Xp, Mp = {}, {}
        for k in range(K):
            if n[k] >= 1:
                Xp[unit(k)] = -2*Lp[k]; Mp[unit(k)] = Lz[k]
            for l in range(k, K):
                e = tuple(a + b for a, b in zip(unit(k), unit(l)))
                if all(e[q] <= n[q] for q in range(K)):
                    g = redI(pair(W[k], W[l]))
                    Xp[e] = Xp.get(e, ZERO) - (g if k == l else 2*g)
        Xp = {e: c for e, c in Xp.items() if c != 0}; Mp = {e: c for e, c in Mp.items() if c != 0}
        def mul(A, B):
            out = {}
            for ea, ca in A.items():
                for eb, cb in B.items():
                    e = tuple(a + b for a, b in zip(ea, eb))
                    if all(e[q] <= n[q] for q in range(K)):
                        out[e] = out.get(e, ZERO) + ca*cb
            return {e: redI(c) for e, c in out.items()}
        ntot = sum(n)
        zero_e = tuple(0 for _ in range(K))
        Xpow = [{zero_e: ONE}]; Mpow = [{zero_e: ONE}]
        for _ in range(ntot):
            Xpow.append(mul(Xpow[-1], Xp)); Mpow.append(mul(Mpow[-1], Mp))
        tab = {}
        for j in range(ntot + 1):
            for i in range(ntot + 1):
                acc = ZERO
                for e, c in Xpow[j].items():
                    e2 = tuple(a - b for a, b in zip(n, e))
                    if e2 in Mpow[i]:
                        acc = acc + c*Mpow[i][e2]
                acc = redI(acc)
                if acc != 0:
                    tab[(j, i)] = acc
        self._tab[key] = tab
        return tab

def _fr(q):
    from fractions import Fraction
    f = Fraction(q)
    return f.numerator, f.denominator

def emul(A, B):
    out = {}
    for ka, ca in A.items():
        for kb, cb in B.items():
            if ka[6] and kb[6]:
                raise ValueError('two epsilon forms')
            k = tuple(x + y for x, y in zip(ka[:6], kb[:6])) + (ka[6] or kb[6],)
            out[k] = out.get(k, ZERO) + ca*cb
    return {k: redI(c) for k, c in out.items() if c != 0}
def epow(A, n):
    R = {(0, 0, 0, 0, 0, 0, 0): ONE}
    for _ in range(n):
        R = emul(R, A)
    return R
def eadd(A, B, cB=ONE):
    out = dict(A)
    for k, c in B.items():
        out[k] = out.get(k, ZERO) + cB*c
    return {k: c for k, c in out.items() if c != 0}

def scri_mono(mono, kin):
    """detector-limit expression of a monomial (m1,m3,m4,n13,n14,n34,eps) of structures.py"""
    B = kin.blocks()
    m1, m3, m4, n13, n14, n34, eps = mono
    R = {(0, 0, 0, 0, 0, 0, 0): ONE}
    for nm, pw in (('V1', m1), ('V3', m3), ('H13', n13), ('H14', n14), ('H34', n34), ('V4', m4)):
        if pw:
            R = emul(R, epow(B[nm], pw))
    if eps:
        R = emul(R, B['e' + eps])
    return R

def pow2(n):
    return C(2**n) if n >= 0 else C(1, 2**(-n))

def poch(x, n):
    assert n >= 0
    r = ONE
    for q in range(n):
        r = r*(x + q)
    return r

def amplitude(expr, kin, parity, sig=None, K=None):
    """expr: scri expression (dict).  Returns (poly, K) with
         M(theta) = G * poly(t, d, s, I),
         G = -4 pi i e^{i pi e14} pi^{5/2} 2^{4-2 e34} / (Gamma(e14+Ka) Gamma(e13+Kb) Gamma(e34+KB) Gamma(e34-1/2+Kg)),
       e13 = 5 - s (+1/2 if odd), e14 = e34 = s (+1/2 if odd),  s = (Delta+J)/2,  K = (Ka,Kb,KB,Kg).
       sig: None (symbolic s) or a rational number (then poly has no s)."""
    S = s if sig is None else C(*_fr(sig))
    h = C(1, 2) if parity == 'odd' else ZERO
    e13, e14, e34 = 5 - S + h, S + h, S + h
    esum = 6 if parity == 'odd' else 5
    if K is None:
        K = Kbounds([expr], parity)
    Ka, Kb, KB, Kg = K
    total = ZERO
    for (i, j, k, l, n1, n2, eps), cf in expr.items():
        nfac = esum - 2 - i - j
        if nfac < 0:
            raise ValueError('divergent alpha integral')
        A = esum - 1 - i - j - l
        a, b, Bx = e14 - j, e13 - i, e34 - k
        g = Bx + A - C(3, 2)
        pref = cf*(-1 if i % 2 else 1)*math.factorial(nfac)*poch(a, Ka + j)*poch(b, Kb + i)*poch(Bx, KB + k) \
            * poch(Bx + A - C(1, 2), Kg - (A - k))
        pref = pref*pow2(2*k - l)
        pref = pref*(-I)**(A % 4) if A >= 0 else pref*(I)**((-A) % 4)
        ntot = n1 + n2 + (1 if eps else 0)
        tab = kin.table(n1, n2, eps)
        D = ZERO
        binom = ONE
        binoms = [ONE]
        for jj in range(1, ntot + 1):
            binom = binom*(g - (jj - 1))*C(1, jj)
            binoms.append(binom)
        for (jj, ii), cji in tab.items():
            # 4^{-A} (A)_ii / ii!
            num = 1
            for q in range(ii):
                num *= (A + q)
            D = D + binoms[jj]*cji*C(num, math.factorial(ii))
        D = D*pow2(-2*A)
        D = D*(I**(ntot % 4))*math.factorial(n1)*math.factorial(n2)
        total = total + redI(pref*D)
    return redI(total), K

def Kbounds(exprs, parity):
    esum = 6 if parity == 'odd' else 5
    Ka = max(-j for e in exprs for (i, j, k, l, n1, n2, eps) in e)
    Kb = max(-i for e in exprs for (i, j, k, l, n1, n2, eps) in e)
    KB = max(-k for e in exprs for (i, j, k, l, n1, n2, eps) in e)
    Kg = max((esum - 1 - i - j - l) - k for e in exprs for (i, j, k, l, n1, n2, eps) in e)
    return (Ka, Kb, KB, Kg)

def tcoeffs(poly, J):
    """{m: polynomial in (d, s, I)} from the coefficient of t^{J-m}"""
    out = {}
    for (et, ed, es, eI), c in poly.to_dict().items():
        m = J - et
        out[m] = out.get(m, ZERO) + CTX.from_dict({(0, ed, es, eI): c})
    return {m: v for m, v in out.items() if v != 0}

# ---------------------------------------------------------------- two-point norms
def norm_poly(J, sig=None):
    """N_{O,m} = <O_m(p)|O_m(p)>/(2 pi)^3 delta = pi^{5/2} 2^{4-2 Delta-2J} / (Gamma(Delta+J) Gamma(Delta+J-1/2)) * n_m(Delta),
       returns {m: n_m} with Delta = 2 s - J (polynomials in s, or numbers if sig given).
       <O(x,z')O(0,z)> = (z'.z x^2 - 2 z'.x z.x)^J / (x^2)^{Delta+J}."""
    S = s if sig is None else C(*_fr(sig))
    Delta = 2*S - J
    zb = zeta_vec(ONE, ONE)            # bra: zeta(b',a') with a' = b' = 1
    zk = zeta_vec(ONE, t)              # ket: zeta(a,b) = zeta(1,t)
    p = [ONE, ZERO, ZERO]
    zz = redI(dot(zb, zk)); pb = redI(dot(p, zb)); pk = redI(dot(p, zk))
    total = ZERO
    for k in range(J + 1):
        # (zb.x)^k (zk.x)^k (x^2)^{-Delta-k}: FT = c_B (i zb.d)^k (i zk.d)^k u^g, g = Delta + k - 3/2, u = -p^2
        g = Delta + k - C(3, 2)
        # generating function: u(p + a zb + b zk) = 1 - 2 a pb - 2 b pk - 2 a b zz
        # coefficient of a^k b^k in (1 - X)^g: sum_j binom(g,j) (-1)^j [a^k b^k] X^j, X = 2a pb + 2b pk + 2ab zz
        acc = ZERO
        binom = ONE
        for j in range(0, 2*k + 1):
            if j > 0:
                binom = binom*(g - (j - 1))*C(1, j)
            r = 2*k - j                         # number of (ab zz) factors; j = (k-r) + (k-r) + r
            if r < 0 or r > k:
                continue
            mult = math.factorial(j)//(math.factorial(r)*math.factorial(k - r)**2)
            acc = acc + binom*((-1)**j)*mult*(2**j)*pb**(k - r)*pk**(k - r)*zz**r
        deriv = acc*math.factorial(k)**2*((-1)**k)          # i^{2k}
        coef = math.comb(J, k)*((-2)**k)*zz**(J - k)
        # 2^{4-2B}/(Gamma(B)Gamma(B-1/2)), B = Delta + k, relative to B = Delta + J
        gam = poch(Delta + k, J - k)*poch(Delta + k - C(1, 2), J - k)*C(4**(J - k))
        total = total + coef*deriv*gam
    total = redI(total)
    out = {}
    for (et, ed, es, eI), c in total.to_dict().items():
        m = J - et
        out[m] = out.get(m, ZERO) + CTX.from_dict({(0, 0, es, eI): c})
    return out
