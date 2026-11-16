"""Fixed-u dispersive sum rules for the UNSTRIPPED helicity amplitude  Mu(s,u) := t^4 f(s,u)  (FLAT_UNSTRIPPED.md).

All statements here are exact (sympy, rational arithmetic); each prints PASS/FAIL, final VERDICT.
  U1  heavy densities of B'_k from the cut discontinuities of Mu (full Wigner d_00, d_44, common prefactor);
  U2  kernel identities  B_k^(1) (stripped) = sum_j c_kj(p^2) B'_j  (k = 2..5)  and the inverse  B'_j = sum_i e_ji(p^2) B_i^(1);
  U3  low-energy sides  -B'_k|low (k = 2..10) from (2.7a) + matter (B.1); contacts g4, g5, g6, g6' appear already at k = 2;
      the identities of U2 hold on the low sides as well (independent residue computation);
  U4  F_g3 (completed), M01, M02 rewritten in the unstripped basis; their low sides recomputed FROM the unstripped low sides:
      contacts cancel, k, A, B as in the stripped computation;
  U5  the naive uplift's flat limit F_U = int[psi2 B'_2 + psi3 B'_3] - 5 d B'_4|0 + 5 B'_5|0: low side keeps contacts (g4 ...);
  U6  F' (dtilde44 -> d44 in F_g3) is not any combination of fixed-u sum rules (Puiseux obstruction).
"""
import json, os, sys
import sympy as sp

s, q, p, w = sp.symbols('s q p w', positive=True)          # q = p^2, w = m^2 (heavy), units M = 1
P, d = sp.symbols('P d')                                     # P = d^J_00 = P_J(x), d = d^J_44(x) (full), x = 1-2q/w
G, g3, g4, g5, g6, g6p, ml, a0, a2, eps = sp.symbols('G g3 g4 g5 g6 g6p m_l a0 a2 epsilon', positive=True)
OK = True
def check(name, cond):
    global OK
    OK &= bool(cond)
    print('%-70s %s' % (name, 'PASS' if cond else 'FAIL'))

t = q - s                                                    # fixed u = -q:  t = -s-u = q - s
def Ks(k):   # stripped B^(1)_k kernel (acts on f(s) +/- f(t))
    return (2*s - q)/(s*(s - q))**sp.Rational(k - 2, 2) if k % 2 == 0 else 1/(s*(s - q))**sp.Rational(k - 3, 2)
def Ku(j):   # unstripped B'_j kernel (acts on Mu(s,u) +/- Mu(q-s,u) = t^4 f(s,u) +/- s^4 f(t,u))
    return (2*s - q)/(s*(s - q))**sp.Rational(j + 2, 2) if j % 2 == 0 else 1/(s*(s - q))**sp.Rational(j + 1, 2)
# effective kernel on f(s,u) alone (reflection s -> q-s of the second term; see .md section 1)
def eff_s(k): return 2*Ks(k) if k % 2 == 0 else -2*Ks(k)
def eff_u(j): return 2*Ku(j)*t**4 if j % 2 == 0 else -2*Ku(j)*t**4

# ---------------------------------------------------------------- U1 heavy densities
# s-cut s = w: Im f(w,-q) = (16pi/w^4) dtilde44 = (16pi/w^4) d (w/(w-q))^4 ; Im Mu = t^4 Im f, t = q - w  -> 16pi d
# t-cut (second term at s = w): Im f(q-w,-q) [t-channel, t = w] = (16pi/w^4) P ; Mu(q-s) = s^4 f(t,u) -> 16pi P
ImMu_s = (q - w)**4 * (w/(w - q))**4 * d / w**4
ImMu_t = w**4 * P / w**4
check('U1 Im Mu on s-cut = 16pi d44 (full), on t-cut = 16pi d00 (prefactors cancel)',
      sp.simplify(ImMu_s - d) == 0 and sp.simplify(ImMu_t - P) == 0)
def heavy_u(j):   # per-state raw density of B'_j (common positive factor 16(2J+1)/pi dropped)
    Z = 1/(w*(w - q))
    return sp.simplify((2*w - q)*Z**sp.Rational(j + 2, 2)*(P + d)) if j % 2 == 0 else sp.simplify(Z**sp.Rational(j + 1, 2)*(P - d))
def heavy_s(k):   # stripped B^(1)_k, (2.33)/(2.35a), with dtilde44 = d w^4/(w-q)^4
    Z = 1/(w*(w - q)); a4 = 1/w**4; b4 = 1/(w - q)**4
    return (2*w - q)*Z**sp.Rational(k - 2, 2)*(a4*P + b4*d) if k % 2 == 0 else Z**sp.Rational(k - 3, 2)*(a4*P - b4*d)

# ---------------------------------------------------------------- U2 identities
COMB = {}
for k in [2, 3, 4, 5]:
    if k % 2 == 0:
        COMB[k] = {k: 1, k + 2: 2*q**2, k + 4: q**4/2, k + 1: -4*q, k + 3: -3*q**3, k + 5: -q**5/2}
    else:
        COMB[k] = {k: 1, k + 2: 2*q**2, k + 4: q**4/2, k + 1: -q, k + 3: -q**3/2}
    ker = sp.simplify(eff_s(k) - sum(c*eff_u(j) for j, c in COMB[k].items()))
    hv = sp.simplify(heavy_s(k) - sum(c*heavy_u(j) for j, c in COMB[k].items()))
    check('U2 B^(1)_%d = %s  (kernels, heavy densities)' % (k, ' + '.join('(%s)B\'_%d' % (c, j) for j, c in COMB[k].items())),
          ker == 0 and hv == 0)
# inverse: B'_j = sum_{i=j}^{j+5} e_i(q) B^(1)_i
for j in [2, 3]:
    es = sp.symbols('e0:6')
    expr = sp.together(eff_u(j) - sum(es[i]*eff_s(j + i) for i in range(6)))
    num = sp.Poly(sp.numer(expr), s)
    sol = sp.solve(num.coeffs(), es, dict=True)
    okinv = len(sol) == 1 and sp.simplify(sp.together(eff_u(j) - sum(sol[0].get(es[i], 0)*eff_s(j + i) for i in range(6)))) == 0
    check("U2 inverse: B'_%d = %s" % (j, ' + '.join('(%s)B_%d' % (sp.factor(sol[0][es[i]]), j + i) for i in range(6) if sol[0].get(es[i], 0) != 0)), okinv)
# convergence: kernel x Mu ~ s^{-k-1} s^{J0} ; contour at infinity vanishes iff J0 < k
for j in [2, 3]:
    deg = sp.degree(sp.numer(sp.together(Ku(j))), s) - sp.degree(sp.denom(sp.together(Ku(j))), s)
    check("U2 Regge: deg_s Ku(%d) = %d = -%d-1 -> B'_%d converges iff J0 < %d (J0 < 2 assumed, (2.23)/(2.25))" % (j, deg, j, j, j), deg == -j - 1)

# ---------------------------------------------------------------- U3 low sides
def flow(S, U):
    T = -S - U
    return 8*sp.pi*G/(S*T*U) + 2*sp.pi*G*S*U/T*g3 + g4 + g5*T + g6*T**2 - g6p*S*U
def fmat(S, U):
    T = -S - U
    return a0/(ml**2 - T) + a2*(1 - 6*S*U/ml**4)/(ml**2 - T)
def low(kernel_fn, k, f, poles, unstripped):
    if unstripped:
        A1, A2 = (q - s)**4*f(s, -q), s**4*f(q - s, -q)
    else:
        A1, A2 = f(s, -q), f(q - s, -q)
    br = (A1 + A2) if k % 2 == 0 else (-A1 + A2)
    e = kernel_fn(k)*br
    return sp.expand(sp.simplify(sum(sp.residue(e, s, z) for z in poles)/2))
LOWU = {}; LOWUM = {}
for k in range(2, 11):
    LOWU[k] = low(Ku, k, flow, [0, q], True)
    LOWUM[k] = sp.factor(low(Ku, k, fmat, [0, q, ml**2, q - ml**2], True))
for k in range(2, 8):
    print("    -B'_%d|low (EFT) = %s" % (k, sp.collect(LOWU[k], [G, g4, g5, g6, g6p])))
for k in range(2, 4):
    print("    -B'_%d|low (matter) = %s" % (k, LOWUM[k]))
check("U3 contacts already in -B'_2|low (coefficient of g4 = %s)" % LOWU[2].coeff(g4), LOWU[2].coeff(g4) != 0)
LOWS = {k: low(Ks, k, flow, [0, q], False) for k in [2, 3, 4, 5]}
LOWSM = {k: sp.factor(low(Ks, k, fmat, [0, q, ml**2, q - ml**2], False)) for k in [2, 3, 4, 5]}
for k in [2, 3, 4, 5]:
    rhs = sum(c*LOWU[j] for j, c in COMB[k].items()); rhsm = sum(c*LOWUM[j] for j, c in COMB[k].items())
    check('U3 identity on low sides, k=%d (EFT and matter)' % k,
          sp.simplify(LOWS[k] - rhs) == 0 and sp.simplify(LOWSM[k] - rhsm) == 0)

# ---------------------------------------------------------------- U4 F_g3, M01, M02 in the unstripped basis
x = sp.Symbol('x')
def unstripped_weights(psi2, psi3):
    """phi_j(p) such that int psi2 B_2 + psi3 B_3 = int sum_j phi_j B'_j  (pointwise in p)."""
    phi = {}
    for k, psi in [(2, psi2), (3, psi3)]:
        for j, c in COMB[k].items():
            phi[j] = sp.expand(phi.get(j, 0) + psi*sp.sympify(c).subs(q, p**2))
    return phi
def low_from_unstripped(phi, a, eps_):
    integrand = sp.expand(sum(phi[j]*LOWU[j].subs(q, p**2) for j in phi))
    integrand_m = sum(phi[j]*LOWUM[j].subs(q, p**2) for j in phi)
    tot = sp.integrate(integrand, (p, eps_, 1))
    # forward: -a dB^(1)_4|0 + a B^(1)_5|0 = -a dB'_4|0 + 5a B'_5|0   (from U2 at p -> 0)
    fw = -a*sp.diff(LOWU[4], q).subs(q, 0) + 5*a*LOWU[5].subs(q, 0)
    fwm = -a*sp.diff(LOWUM[4], q).subs(q, 0) + 5*a*LOWUM[5].subs(q, 0)
    return sp.expand(tot + fw), integrand_m, fwm
fwd_check = sp.simplify(sp.diff(sum(c*LOWU[j] for j, c in COMB[4].items()), q).subs(q, 0) - (sp.diff(LOWU[4], q).subs(q, 0) - 4*LOWU[5].subs(q, 0)))
check("U4 forward map: d_q B^(1)_4|0 = d_q B'_4|0 - 4 B'_5|0 ; B^(1)_5|0 = B'_5|0", fwd_check == 0 and sp.simplify(LOWS[5].subs(q, 0) - LOWU[5].subs(q, 0)) == 0)
here = os.path.dirname(os.path.abspath(__file__))
funcs = {'F_g3': (sp.expand((1 - p)**3*p*(65 + 155*p + 47*p**2)), sp.expand(2*(1 - p)**3*(81 - 1927*p + 4601*p**2)), sp.Integer(5))}
jpath = os.path.join(here, '..', 'robust', 'functionals_certified.json')
if os.path.exists(jpath):
    J = json.load(open(jpath))
    for name in ['M01', 'M02']:
        ps2 = sum(sp.Rational(c)*p**i for i, c in enumerate(J[name]['psi2_exact']))
        ps3 = sum(sp.Rational(c)*p**i for i, c in enumerate(J[name]['psi3_exact']))
        funcs[name] = (sp.expand(ps2), sp.expand(ps3), sp.Rational(J[name]['a_exact']))
for name, (ps2, ps3, a) in funcs.items():
    phi = unstripped_weights(ps2, ps3)
    lowtot, im, fwm = low_from_unstripped(phi, a, eps)
    contacts = [lowtot.coeff(c) for c in (g4, g5, g6, g6p)]
    kk = sp.nsimplify(-lowtot.coeff(g3)/(sp.pi*G))
    k0 = sp.simplify(kk.subs(eps, 0))
    grav = sp.expand(lowtot.coeff(G).subs(g3, 0)/(16*sp.pi))       # = L(eps)
    x0 = -sp.limit(eps*sp.diff(grav, eps), eps, 0)
    Lfin = sp.limit(grav - x0*sp.log(1/eps), eps, 0)
    A = 16*x0/k0; B = 16*Lfin/k0
    # matter: pointwise integrand vs stripped
    ims = ps2*LOWSM[2].subs(q, p**2) + ps3*LOWSM[3].subs(q, p**2)
    okm = sp.simplify(sp.together(im - ims)) == 0
    print('    %s unstripped weights: %s' % (name, ', '.join("B'_%d: deg %d" % (j, sp.degree(phi[j], p)) for j in sorted(phi))))
    check('U4 %s from unstripped sum rules: contacts g4,g5,g6,g6\' -> %s ; k(0) = %s ; A = %.6f ; B = %.6f ; matter = stripped'
          % (name, contacts, k0, float(A), float(B), ), all(c == 0 for c in contacts) and okm)

# ---------------------------------------------------------------- U5 naive uplift flat limit F_U
ps2, ps3, a = funcs['F_g3']
lowU = sp.expand(sp.integrate(sp.expand(ps2*LOWU[2].subs(q, p**2) + ps3*LOWU[3].subs(q, p**2)), (p, eps, 1))
                 - 5*sp.diff(LOWU[4], q).subs(q, 0) + 5*LOWU[5].subs(q, 0))
cU = {str(c): sp.nsimplify(sp.limit(lowU.coeff(c), eps, 0)) for c in (g4, g5, g6, g6p)}
print('    F_U low side contact coefficients (eps->0):', cU)
check('U5 F_U (naive transcription of Omega_g3 with C_k -> B\'_k) keeps contact couplings', any(v != 0 for v in cU.values()))
# heavy action of F_U in normal form (per |c|^2, u = 1/m^2):
#   ++ : int [psi2 (2-uq)/(1-uq)^2 + u psi3/(1-uq)^2] P_J + (10JJ-20) u^3 ;  +- : int [psi2 (2-uq) - u psi3]/(1-uq)^2 d44 + (10JJ-190) u^3
JJ = sp.Symbol('JJ'); uu = sp.Symbol('u')
dq_P = -JJ/w; dq_d = -(JJ - 16)/w            # d/dq of P_J(1-2q/w), d44(1-2q/w) at q = 0
fw_pp = sp.simplify(w**3*(-5*(sp.diff((2*w - q)/(w*(w - q))**3, q).subs(q, 0) + 2*w/(w**2)**3*dq_P) + 5/(w**2)**3))
fw_pm = sp.simplify(w**3*(-5*(sp.diff((2*w - q)/(w*(w - q))**3, q).subs(q, 0) + 2*w/(w**2)**3*dq_d) - 5/(w**2)**3))
print('    F_U forward terms (normal form): ++ %s ; +- %s' % (sp.expand(fw_pp*w**3), sp.expand(fw_pm*w**3)), ' (times u^3)')

# ---------------------------------------------------------------- U6 F' is not a sum-rule combination
r = sp.Symbol('r', positive=True)                # Z = r^2
wZ = (q + sp.sqrt(q**2 + 4/r**2))/2              # w as function of Z = 1/(w(w-q)), w > q
ser = sp.series(1/wZ**4, r, 0, 7).removeO()
odd = sp.expand(ser).coeff(r, 5)
check("U6 1/w^4 = %s + ... as series in sqrt(Z): odd power present -> F' needs a^4 = 1/w^4 in span{Z^n}: impossible" % sp.expand(ser),
      sp.simplify(odd) != 0)

print('VERDICT sumrules: %s' % ('PASS' if OK else 'FAIL'))
