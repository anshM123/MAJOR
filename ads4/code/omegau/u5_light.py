"""u5: light side of Omega^U (M02 primary, M01, F_g3), and the conditional AdS bound vs the heavy threshold n.

(a) [D, exact sympy] flat low sides -B'_k|low, k = 2..7, of the unstripped sum rules (same formula as flatunstripped/sumrules.py U3,
    re-implemented): Einstein only in B'_2 (16 pi G/q); contacts g4, g5, g6, g6' in B'_2..B'_6; |g3|^2 (= |ghat3|^2) in B'_2..B'_4.
(b) [D] flat Omega^U: contacts cancel, k(0), A, B, B_AdS (reproduces FLAT_UNSTRIPPED.md / ROBUST_FUNCTIONAL.md).
(c) [D given (F), (P)] AdS: the light side of each moment is a polynomial in nu^2 whose top coefficient is the flat one; the lower ones
    a_{k,c,j} (AdS corrections, not computed for <TTTT>) leave residues  sum_k lev_{k,c} a_{k,c} / Lam^2 + O(Lam^-4)  in -Omega^U_light.
    lev_{k,c} = (flat top coefficient) x int phi_k p^{2d-2}  (+ the forward pieces).  Tabulated, with the induced shift of the bound on
    |ghat3|^2 Lam^8 for illustrative sizes (|a| = 36 = 4 Df^2 - 1/4 at Df = 3, the scalar c3 STU value; g_f = G_f log Lam).
    Same for the |ghat3|^2 coefficient k_AdS = k(0) + O(Lam^-2).
(d) [D given counting rule + RB] Einstein: no term in C_k, k >= 3, iff J0(G^E) < k; leverage of a test residue 16 pi G (nu^2+9/4)^{-k}
    in each k-term, relative to the leading 16 pi G Lam^2 L_AdS (as in adslight/t4 for Omega).
(e) conditional bound |ghat3|^2 Lam^8 <= n^8 16 L_AdS(Lam/n)/k  and  sqrt(t4^2 + alpha_T^2) <= 240 sqrt(.)/Lam^4  for thresholds n.
Output: cache/u5.log
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp, json
from fractions import Fraction as Fr

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

s, q, p = sp.symbols('s q p', positive=True)
G, g3, g4, g5, g6, g6p, eps = sp.symbols('G g3 g4 g5 g6 g6p epsilon', positive=True)
def Ku(j):
    return (2*s - q)/(s*(s - q))**sp.Rational(j + 2, 2) if j % 2 == 0 else 1/(s*(s - q))**sp.Rational(j + 1, 2)
def flow(S, U):
    T = -S - U
    return 8*sp.pi*G/(S*T*U) + 2*sp.pi*G*S*U/T*g3 + g4 + g5*T + g6*T**2 - g6p*S*U
def low(k):
    A1, A2 = (q - s)**4*flow(s, -q), s**4*flow(q - s, -q)
    br = (A1 + A2) if k % 2 == 0 else (-A1 + A2)
    e = Ku(k)*br
    return sp.expand(sp.simplify(sum(sp.residue(e, s, z) for z in (0, q))/2))
LOW = {k: low(k) for k in range(2, 9)}
pr("(a) flat -B'_k|low (EFT, through g6, g6'):")
for k in range(2, 9):
    pr(f"    k={k}: {sp.collect(LOW[k], [G, g4, g5, g6, g6p])}")
ok = all(sp.simplify(LOW[k].coeff(G).subs(g3, 0)) == 0 for k in range(3, 9)) and sp.simplify(LOW[2].coeff(G).subs(g3, 0) - 16*sp.pi/q) == 0
pr(f"    Einstein term only in B'_2 (16 pi G/q), none in B'_3..B'_8: {ok}")

# ---------------------------------------------------------------- functionals
here = os.path.dirname(os.path.abspath(__file__))
Jd = json.load(open(os.path.join(here, '..', 'robust', 'functionals_certified.json')))
FUN = {'F_g3': (sp.expand((1 - p)**3*p*(65 + 155*p + 47*p**2)), sp.expand(2*(1 - p)**3*(81 - 1927*p + 4601*p**2)), sp.Integer(5))}
for name in ('M02', 'M01'):
    FUN[name] = (sum(sp.Rational(c)*p**i for i, c in enumerate(Jd[name]['psi2_exact'])),
                 sum(sp.Rational(c)*p**i for i, c in enumerate(Jd[name]['psi3_exact'])), sp.Rational(Jd[name]['a_exact']))
def phis(ps2, ps3):
    return {2: ps2, 3: sp.expand(ps3 - 4*p**2*ps2), 4: sp.expand(2*p**4*ps2 - p**2*ps3), 5: sp.expand(-3*p**6*ps2 + 2*p**4*ps3),
            6: sp.expand((p**8*ps2 - p**6*ps3)/2), 7: sp.expand((-p**10*ps2 + p**8*ps3)/2)}
mu = lambda f, j: sp.integrate(sp.expand(f*p**j), (p, 0, 1))

CONST = {}
for name, (ps2, ps3, a) in FUN.items():
    PH = phis(ps2, ps3)
    integrand = sp.expand(sum(PH[k]*LOW[k].subs(q, p**2) for k in PH))
    tot = sp.integrate(integrand, (p, eps, 1))
    fw = -a*sp.diff(LOW[4], q).subs(q, 0) + 5*a*LOW[5].subs(q, 0)
    lowtot = sp.expand(tot + fw)
    contacts = {str(c): sp.limit(lowtot.coeff(c), eps, 0) for c in (g4, g5, g6, g6p)}
    k0 = sp.nsimplify(sp.limit(-lowtot.coeff(g3)/(sp.pi*G), eps, 0))
    x0 = ps2.coeff(p, 1)
    Lfin = sp.integrate(sp.expand((ps2 - x0*p)/p**2), (p, 0, 1))
    A = 16*x0/k0; B = 16*Lfin/k0; BA = B - A*sp.log(sp.Rational(3, 2))
    CONST[name] = dict(k0=k0, A=float(A), B=float(B), BA=float(BA), ps2=ps2, PH=PH, a=a)
    pr(f"(b) {name}: contacts in flat Omega^U low side {contacts}; k(0) = {k0} = {float(k0):.6f}; A = {float(A):.6f}; B = {float(B):.6f}; B_AdS = {float(BA):.6f}")

# ---------------------------------------------------------------- (c) AdS lower-coefficient leverage
pr("(c) AdS residues from lower-degree coefficients of the moment light sides (leading 1/Lam^2 unknowns), normalised functionals (k(0) = 1):")
pr("    lev_{k,c} = (top flat coeff of -B'_k in coupling c) x int phi_k p^{2d-2} / k(0);  forward: -a d/dq of the 1/Lam^2 term is 0 for d <= 1")
LEV = {}
for name, C in CONST.items():
    PH, a, k0 = C['PH'], C['a'], C['k0']
    rows = []
    for c in (g4, g5, g6, g6p, g3):
        levs = {}
        for k in range(2, 8):
            poly = sp.Poly(sp.expand(LOW[k].coeff(c).subs(G, 1)), q) if LOW[k].coeff(c) != 0 else None
            if poly is None:
                continue
            d = poly.degree(); top = poly.coeff_monomial(q**d)
            if d >= 1:
                levs[k] = float(top*mu(PH[k], 2*d - 2)/k0)
        # forward pieces: -a d_{nu^2}(-C4) at 0 picks the coefficient of nu^2 of the AdS polynomial: for top degree d4 >= 2 the 1/Lam^2
        # unknown multiplies nu^{2(d4-1)} -> contributes at q-degree d4-1 >= 1 -> derivative at 0 only if d4 - 1 == 1
        p4 = sp.Poly(sp.expand(LOW[4].coeff(c).subs(G, 1)), q) if LOW[4].coeff(c) != 0 else None
        if p4 is not None and p4.degree() == 2:
            levs['fw4'] = float(-a*p4.coeff_monomial(q**2)/k0)
        p5 = sp.Poly(sp.expand(LOW[5].coeff(c).subs(G, 1)), q) if LOW[5].coeff(c) != 0 else None
        if p5 is not None and p5.degree() >= 1:
            levs['fw5'] = float(5*a*p5.coeff_monomial(q**p5.degree())/k0) if p5.degree() == 1 else 0.0
        LEV[(name, str(c))] = levs
        rows.append(f"      {str(c):4s}: " + ", ".join(f"{kk}: {v:+.3g}" for kk, v in levs.items()) + f"   (sum |lev| = {sum(abs(v) for v in levs.values()):.3g})")
    pr(f"  {name}:"); [pr(r) for r in rows]
pr("    size of the induced shift of the bound on |ghat3|^2 Lam^8 (normalised k = 1): delta ~ sum_k lev a/(pi Lam^2) x (g_f/G_f);")
pr("    illustrative |a| = 36 (scalar: -C_2[c3 STU] = nu^2 - (2Df+1/2)(2Df-1/2), Df = 3), g_f/G_f = log Lam (flat two-sided bound size):")
for name in ('M02', 'F_g3'):
    for Dg in (20, 100, 1000):
        Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
        tot = sum(sum(abs(v) for v in LEV[(name, c)].values()) for c in ('g4', 'g5', 'g6', 'g6p'))
        pr(f"    {name} Dg={Dg}: contacts sum|lev| = {tot:.3g} -> delta(|ghat3|^2 Lam^8) <~ {tot*36*np.log(Lam)/(np.pi*Lam**2):.3g}"
           f"  vs 16 L_AdS/k ~ {CONST[name]['A']*np.log(Lam) + CONST[name]['BA']:.3g};  k_AdS - k(0): sum|lev(g3)| x 36/Lam^2 = "
           f"{sum(abs(v) for v in LEV[(name, 'g3')].values())*36/Lam**2:.3g} (relative)")

# ---------------------------------------------------------------- (d) Einstein test residues
pr("(d) Einstein residues: counting rule => C_k|_Einstein = 0 iff J0(G^E) < k (k = 3..7 all need only J0 <= 2 = RB_2) [I].")
pr("    leverage of a test residue 16 pi G (nu^2+9/4)^{-k} in the k-term of Omega^U, relative to 16 pi G Lam^2 L_AdS:")
from scipy.integrate import quad
for name in ('M02', 'F_g3'):
    PH = {k: sp.lambdify(p, CONST[name]['PH'][k]) for k in range(2, 8)}
    ps2 = sp.lambdify(p, CONST[name]['ps2']); a = float(CONST[name]['a'])
    for Dg in (20, 100, 1000):
        Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
        L = quad(lambda t: ps2(t)/(t*t + 2.25/Lam**2), 0, 1, limit=400, points=[1.5/Lam])[0]
        lev = []
        for k in range(3, 8):
            ek = 1.0 if k % 2 == 0 else 0.5
            v = quad(lambda t: PH[k](t)*ek*Lam**(2*k)*(Lam**2*t*t + 2.25)**(-k), 0, 1, limit=400, points=[1.5/Lam, 10/Lam])[0]
            lev.append(abs(v)/(Lam**2*L))
        f4 = a*Lam**10*4*(2.25)**(-5)/(Lam**2*L)       # -a Lam^10 d_{nu^2}(nu^2+9/4)^{-4}|_0
        f5 = 2.5*a*Lam**10*(2.25)**(-5)/(Lam**2*L)
        pr(f"    {name} Dg={Dg}: k=3..7: " + " ".join(f"{v:.2g}" for v in lev) + f";  forward C4: {f4:.2g}, C5: {f5:.2g}")

# ---------------------------------------------------------------- (e) bound vs threshold
pr("(e) conditional bound: |ghat3|^2 Lam^8 <= n^8 16 L_AdS(Lam/n)/k(0)  [needs heavy positivity for m >= n Lam with smearing Lam/n, and all inputs]")
def LAdS(ps2f, Lam):
    return quad(lambda t: ps2f(t)/(t*t + 2.25/Lam**2), 0, 1, limit=400, points=[1.5/Lam])[0]
for name in ('M02', 'F_g3', 'M01'):
    ps2f = sp.lambdify(p, CONST[name]['ps2']); k0 = float(CONST[name]['k0'])
    for n in (1.0, 1.2, 1.3, 1.5, 2.0, 2.357, 3.0):
        vals = []
        for Dg in (20, 50, 100, 300, 1000):
            Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
            X = n**8*16*LAdS(ps2f, Lam/n)/k0
            vals.append(f"{Dg}: {X:.3g} (sqrt(t4^2+aT^2) <= {240*np.sqrt(max(X, 0))/Lam**4:.2g})")
        pr(f"    {name} n={n}: " + "; ".join(vals))
open(os.path.join(HERE, 'cache', 'u5.log'), 'w').write("\n".join(out) + "\n")
