"""Referee re-derivation for FLAT_UNSTRIPPED.md (independent of ads4/code/flatunstripped/).
R1 low sides of B'_k from (2.7a) via Laurent coefficients (own code), sign convention fixed against paper (2.31a)
R2 kernel identities B_k <-> B'_j (k=2,3 forward; B'_2, B'_3 inverse), as rational functions of s on f(s,u) alone
R3 heavy densities from the cut discontinuities with explicit Wigner functions at exact rational points
R4 F_g3 rewritten in B' basis: contacts, k(0); F_U low-side contact coefficients
R5 F_U(++, J=2, m=213/200 and 6/5) by mpmath from the definitional densities (numerical derivative for forward term)
"""
import sympy as sp, mpmath as mp
s, q, e, p = sp.symbols('s q e p')
G, g3, g4, g5, g6, g6p = sp.symbols('G g3 g4 g5 g6 g6p')
ok = True
def rep(name, c):
    global ok; ok &= bool(c); print(('PASS ' if c else 'FAIL ') + name)

def f_low(S, U):
    T = -S - U
    return 8*sp.pi*G/(S*T*U) + 2*sp.pi*G*S*U/T*g3 + g4 + g5*T + g6*T**2 - g6p*S*U

def kern_s(k):  # stripped, acts on f
    return (2*s - q)*(s*(s - q))**sp.Rational(2 - k, 2) if k % 2 == 0 else (s*(s - q))**sp.Rational(3 - k, 2)
def kern_u(k):  # unstripped, acts on t^4 f
    return (2*s - q)*(s*(s - q))**sp.Rational(-(k + 2), 2) if k % 2 == 0 else (s*(s - q))**sp.Rational(-(k + 1), 2)

def integrand(k, unstripped):
    A = lambda S: ((q - S)**4 if unstripped else 1)*f_low(S, -q)   # t = -S-u = q - S
    # second term M(q-s,u): t(q-s) = s
    B = lambda S: (S**4 if unstripped else 1)*f_low(q - S, -q)
    br = A(s) + B(s) if k % 2 == 0 else -A(s) + B(s)
    K = kern_u(k) if unstripped else kern_s(k)
    return K*br

def half_res(expr, poles):
    tot = 0
    for a in poles:
        ser = sp.series(sp.together(expr.subs(s, a + e)), e, 0, 1).removeO()
        tot += sp.expand(ser).coeff(e, -1)
    return sp.expand(sp.simplify(tot/2))

# R1
mB2 = half_res(integrand(2, False), [0, q])
rep('sign convention: 1/2 sum Res (stripped k=2) = %s  (paper -B2 = 16piG/p^2 + 2piG|g3|^2 p^6)' % mB2,
    sp.simplify(mB2 - (16*sp.pi*G/q + 2*sp.pi*G*g3*q**3)) == 0)
LU = {k: half_res(integrand(k, True), [0, q]) for k in range(2, 10)}
LS = {k: half_res(integrand(k, False), [0, q]) for k in range(2, 8)}
claim = {2: 16*sp.pi*G/q + 2*sp.pi*G*g3*q**3 + 4*g4*q**2 + 5*g5*q**3 + (6*g6 - g6p)*q**4,
         3: 2*sp.pi*G*g3*q**2 + 2*g4*q + 3*g5*q**2 + (4*g6 - g6p)*q**3,
         4: 2*g4 + (4*sp.pi*G*g3 + 5*g5)*q + (9*g6 - 3*g6p)*q**2,
         5: g5 + (3*g6 - g6p)*q, 6: 2*g6, 7: sp.Integer(0)}
for k in range(2, 8):
    rep("-B'_%d|low = %s" % (k, LU[k]), sp.simplify(LU[k] - claim[k]) == 0)

# R2 kernel identities on f(s) alone: reflect second term
def eff(k, un):
    K = kern_u(k) if un else kern_s(k)
    t4 = (q - s)**4 if un else 1
    Kref = K.subs(s, q - s)*(s**4 if un else 1)   # weight of f(q-s) term, then s->q-s, ds->-ds, orientation kept
    # \oint K(s) X(q-s) ds = \oint K(q-s') X(s') ds'  (s=q-s', circle->circle same orientation, ds=-ds' and limits swap)
    sec = Kref.subs(s, q - s)  # weight multiplying f(s') coming from second term
    # careful: \oint_C g(s) ds with s = q - s': C maps to C' traversed same orientation; ds = -ds'; so = -\oint g(q-s') ds' ... check numerically below
    return K*t4, sec
# rather than trust orientation algebra, verify identities directly on the full integrands' residues AND on rational kernels
def full_int(k, un):
    K = kern_u(k) if un else kern_s(k)
    return K
Fs, Fu = sp.Function('F')(s), None
# Generic test: replace f by a generic rational test function with poles off {0,q}; identity must hold on ALL residues
a1, a2 = sp.Rational(7, 3), sp.Rational(-11, 5)
def ftest(S, U):
    T = -S - U
    return 1/((S - a1)*(T - a2)) + S**2*U/(T - a1)**3 + 1/((S - a2)**2*(T - a1))   # not s<->t symmetric, arbitrary
def total(k, un, fn, qv):
    A = lambda S: ((qv - S)**4 if un else 1)*fn(S, -qv)
    B = lambda S: (S**4 if un else 1)*fn(qv - S, -qv)
    br = A(s) + B(s) if k % 2 == 0 else -A(s) + B(s)
    K = (kern_u(k) if un else kern_s(k)).subs(q, qv)
    ex = sp.together(K*br)
    # sum over all finite poles = - residue at infinity; compute residue at infinity directly
    return sp.simplify(-sp.series(ex.subs(s, 1/e)/e**2, e, 0, 1).removeO().coeff(e, -1))
qv = sp.Rational(3, 7)
comb = {2: {2: 1, 4: 2*q**2, 6: q**4/2, 3: -4*q, 5: -3*q**3, 7: -q**5/2},
        3: {3: 1, 5: 2*q**2, 7: q**4/2, 4: -q, 6: -q**3/2}}
inv = {2: {2: 1, 3: 4*q, 4: 2*q**2, 5: 3*q**3, 6: q**4/2, 7: q**5/2},
       3: {3: 1, 4: q, 5: 2*q**2, 6: q**3/2, 7: q**4/2}}
# Kernel identity: integrand equality as functions of s, after reflection-symmetrization: compare Keff via symmetrised integrand
def symm_int(k, un, fn):
    A = lambda S: ((q - S)**4 if un else 1)*fn(S, -q)
    B = lambda S: (S**4 if un else 1)*fn(q - S, -q)
    br = A(s) + B(s) if k % 2 == 0 else -A(s) + B(s)
    return (kern_u(k) if un else kern_s(k))*br
fgen = sp.Function('f')
# exact identity of integrands modulo reflection: I(s) + I(q-s) must agree (contour integral of g over circle = that of g(q-s))
for k in (2, 3):
    lhs = symm_int(k, False, ftest); rhs = sum(c*symm_int(j, True, ftest) for j, c in comb[k].items())
    d1 = sp.simplify(sp.together((lhs - rhs) - (lhs - rhs).subs(s, q - s)))
    rep('kernel identity B_%d = sum c_j B\'_j (symmetrised integrands, generic f)' % k, d1 == 0)
    lhs = symm_int(k, True, ftest); rhs = sum(c*symm_int(j, False, ftest) for j, c in inv[k].items())
    d2 = sp.simplify(sp.together((lhs - rhs) - (lhs - rhs).subs(s, q - s)))
    rep("inverse identity B'_%d = sum e_i B_i (symmetrised integrands, generic f)" % k, d2 == 0)
# low-side identities with the EFT
for k in (2, 3):
    rep('low-side identity k=%d (EFT)' % k, sp.simplify(LS[k] - sum(c*LU[j] for j, c in comb[k].items())) == 0)
    rep("low-side inverse k=%d (EFT)" % k, sp.simplify(LU[k] - sum(c*LS[j] for j, c in inv[k].items())) == 0)

# R3 heavy densities with explicit Wigner functions at exact rational points
def dt44(J, y):  # 2F1(4-J, J+5; 1; y)
    return sum(sp.rf(4 - J, n)*sp.rf(J + 5, n)/sp.factorial(n)**2*y**n for n in range(0, J - 3))
def PJ(J, y): return sp.legendre(J, 1 - 2*y)
def heavy(k, un, w, qq, J, cpp, cpm):
    y = qq/w
    Ims = cpm*dt44(J, y)/w**4 if J >= 4 else 0     # Im f on s-cut
    Imt = cpp*PJ(J, y)/w**4                         # Im f on t-cut (t = w)
    K = (kern_u(k) if un else kern_s(k)).subs(q, qq)
    # density  ~ Keff(w) Im_s - Keff(q-w) Im_t, Keff(s) = +-2 K(s) * (t^4 if un), t = q - s
    sg = 1 if k % 2 == 0 else -1
    keff = lambda S: sg*K.subs(s, S)*((qq - S)**4 if un else 1)
    return sp.Rational(keff(w)*Ims - keff(qq - w)*Imt)
w0, q0 = sp.Rational(37, 10), sp.Rational(3, 5)
good = True
for J in (4, 6, 9):
    for k in (2, 3):
        lhs = heavy(k, False, w0, q0, J, 2, 3)
        rhs = sum(c.subs(q, q0)*heavy(j, True, w0, q0, J, 2, 3) if hasattr(c, 'subs') else c*heavy(j, True, w0, q0, J, 2, 3) for j, c in comb[k].items())
        good &= sp.simplify(lhs - rhs) == 0
rep('heavy-density identity at exact rational (w,q), J=4,6,9, k=2,3', good)
# U1-type: unstripped densities = common prefactor x (P +- d44 full)
J = 6; y = q0/w0
d44 = (1 - y)**4*dt44(J, y)
h2 = heavy(2, True, w0, q0, J, 1, 1); Z = 1/(w0*(w0 - q0))
rep("B'_2 density = (2w-q)Z^2 (P + d44full)", sp.simplify(h2 - (2*w0 - q0)*Z**2*(PJ(J, y) + d44)) == 0)
h3 = heavy(3, True, w0, q0, J, 1, 1)
rep("B'_3 density = Z^2 (P - d44full)", sp.simplify(h3 - Z**2*(PJ(J, y) - d44)) == 0)

# R4 functionals
psi2 = (1 - p)**3*p*(65 + 155*p + 47*p**2); psi3 = 2*(1 - p)**3*(81 - 1927*p + 4601*p**2)
phi = {2: psi2, 3: psi3 - 4*p**2*psi2, 4: 2*p**4*psi2 - p**2*psi3, 5: -3*p**6*psi2 + 2*p**4*psi3,
       6: (p**8*psi2 - p**6*psi3)/2, 7: (-p**10*psi2 + p**8*psi3)/2}
eps = sp.Symbol('eps', positive=True)
def lowfun(wts, fwd4, fwd5):
    tot = sp.integrate(sp.expand(sum(wts[j]*LU[j].subs(q, p**2) for j in wts)), (p, eps, 1))
    tot += -fwd4*sp.diff(LU[4], q).subs(q, 0) + fwd5*LU[5].subs(q, 0)
    return sp.expand(tot)
Lg3 = lowfun(phi, 5, 25)
rep('F_g3 in B\' basis (phi_j, -5 dB\'4 + 25 B\'5): contacts %s' % [sp.limit(Lg3.coeff(c), eps, 0) for c in (g4, g5, g6, g6p)],
    all(sp.simplify(Lg3.coeff(c)) == 0 for c in (g4, g5, g6, g6p)))
k0 = sp.limit(-Lg3.coeff(g3)/(sp.pi*G), eps, 0)
rep('k(0) = %s' % k0, k0 == sp.Rational(330151, 12012))
LUf = lowfun({2: psi2, 3: psi3}, 5, 5)
cc = {str(c): sp.limit(LUf.coeff(c), eps, 0) for c in (g4, g5, g6, g6p)}
rep('F_U contact coeffs %s' % cc, cc['g4'] == sp.Rational(60526, 3465) and cc['g5'] == sp.Rational(-2746189, 360360)
    and cc['g6'] == sp.Rational(729983, 90090) and cc['g6p'] == sp.Rational(-358681, 180180))
print('   F_U k(0) =', sp.limit(-LUf.coeff(g3)/(sp.pi*G), eps, 0), ' F_g3 k(0) =', k0)

# R5 F_U(++, J=2) from definitional densities, mpmath, dw-density x w^3
mp.mp.dps = 30
def hpp(k, w, qq):  # ++ unstripped density per |c++|^2 (dropping the factor 2), J=2
    y = qq/w; P = (3*(1 - 2*y)**2 - 1)/2
    Z = 1/(w*(w - qq))
    return ((2*w - qq)*Z**((k + 2)/mp.mpf(2)) if k % 2 == 0 else Z**((k + 1)/mp.mpf(2)))*P
def FU(m):
    w = mp.mpf(m)**2
    ps2 = lambda x: (1 - x)**3*x*(65 + 155*x + 47*x**2)
    ps3 = lambda x: 2*(1 - x)**3*(81 - 1927*x + 4601*x**2)
    I = mp.quad(lambda x: ps2(x)*hpp(2, w, x**2) + ps3(x)*hpp(3, w, x**2), [0, 0.5, 1])
    fw = -5*mp.diff(lambda qq: hpp(4, w, qq), 0) + 5*hpp(5, w, 0)
    return w**3*(I + fw)
for m, ref in [(mp.mpf(213)/200, -3.33384834423), (mp.mpf(6)/5, -1.12774359442)]:
    v = FU(m); rep('m^4 F_U(++,J=2,m=%s) = %s (doc %s)' % (m, mp.nstr(v, 14), ref), abs(v - ref) < 1e-10)
print('VERDICT', 'PASS' if ok else 'FAIL')
