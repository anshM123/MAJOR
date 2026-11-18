"""nf5: light side of Omega^{NF,SC} (NF02) vs Omega^U (M02).
(a) [D, exact] flat: contacts cancel in sum_k phi_k B'_k (no forward terms); k(0), A, B, B_AdS of NF02.
(b) [D given (F), (P)] AdS lower-coefficient leverages lev_{k,g} (same definition as u5_light.py (c); no forward entries now), and the induced
    shift of the bound at n = 1, 1.3, 1.5 (|a| = d x 36, g_f = c_g G_f log(Lam/n), c_g = 1, 30), relative to 16 L_AdS/k.
(c) leverage of a test residue 16 pi G (nu^2+9/4)^{-k} / N_k(nu) in each k-term (the superconvergent replacements are normalised by N_k; for the
    Einstein term the residue is exactly zero given RB_2 -- this measures what a failure of (P) for the replacements would cost), relative to
    the leading 16 pi G Lam^2 L_AdS; Omega^U(M02) forward pieces for comparison.
(d) residual one-pair flip of the even replacements Psi_4/N^ff_4, Psi_6/N^ff_6 with NF02 weights (c2 model of superconv/, leading Regge) / N_s.
Output: cache/nf5.log
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'superconv'))
import numpy as np, sympy as sp
from fractions import Fraction as Fr
from scipy.integrate import quad
from nf_common import NF, KS, EPS, SIGMA, state
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
    e = Ku(k)*((A1 + A2) if k % 2 == 0 else (-A1 + A2))
    return sp.expand(sp.simplify(sum(sp.residue(e, s, z) for z in (0, q))/2))
LOW = {k: low(k) for k in range(2, 8)}

F = NF()
PH = {k: sum(sp.Rational(c.numerator, c.denominator)*p**i for i, c in enumerate(F.PHx[k])) for k in KS}   # normalised (k(0) = 1)
integrand = sp.expand(sum(PH[k]*LOW[k].subs(q, p**2) for k in KS))
lowtot = sp.expand(sp.integrate(integrand, (p, eps, 1)))
contacts = {str(c): sp.limit(lowtot.coeff(c), eps, 0) for c in (g4, g5, g6, g6p)}
k0 = sp.nsimplify(sp.limit(-lowtot.coeff(g3)/(sp.pi*G), eps, 0))
x0 = PH[2].coeff(p, 1); Lfin = sp.integrate(sp.expand((PH[2] - x0*p)/p**2), (p, 0, 1))
A = 16*x0/k0; B = 16*Lfin/k0; BA = B - A*sp.log(sp.Rational(3, 2))
pr(f"(a) NF02 (normalised): flat contacts in Omega^NF low side {contacts}; k(0) = {k0} (=1 by normalisation; unnormalised {F.k0}); "
   f"A = {float(A):.4f}; B = {float(B):.4f}; B_AdS = {float(BA):.4f}")
mu = lambda f, j: sp.integrate(sp.expand(f*p**j), (p, 0, 1))
LEV, DEG = {}, {}
for c in (g4, g5, g6, g6p, g3):
    levs = {}
    for k in range(2, 8):
        if LOW[k].coeff(c) == 0:
            continue
        poly = sp.Poly(sp.expand(LOW[k].coeff(c).subs(G, 1)), q); d = poly.degree()
        if d >= 1:
            levs[k] = float(poly.coeff_monomial(q**d)*mu(PH[k], 2*d - 2)/(sp.pi if c == g3 else 1)); DEG[(str(c), k)] = d
    LEV[str(c)] = levs
    pr(f"(b) lev {str(c):4s}: " + ", ".join(f"k={kk}: {v:+.4g}" for kk, v in levs.items()) + f"   sum|lev| = {sum(abs(v) for v in levs.values()):.4g}")
pr("    (Omega^U, M02, from u5_light.py: sum|lev| = 0.64 (g4), 2.21 (g5), 7.03 (g6, incl. fw 3.94), 2.29 (g6'), 1.63 (|g3|^2))")
def LAdS(Ls):
    f = sp.lambdify(p, PH[2]); return quad(lambda t: f(t)/(t*t + 2.25/Ls**2), 0, 1, limit=400, points=[1.5/Ls])[0]
for Dg in (50, 100, 300, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    for n in (1.0, 1.5, 2.1, 2.9):
        Ls = Lam/n; X = 16*LAdS(Ls)
        tot = sum(abs(v)*DEG.get((c, k), 1) for c in ('g4', 'g5', 'g6', 'g6p') for k, v in LEV[c].items())
        dk = sum(abs(v)*DEG.get(('g3', k), 1) for k, v in LEV['g3'].items())*36/Ls**2
        pr(f"    Dg={Dg} n={n}: contact shift / (16 L_AdS/k) = {tot*36*np.log(Ls)/(np.pi*Ls**2)/X:.2e} (c_g=1), {30*tot*36*np.log(Ls)/(np.pi*Ls**2)/X:.2e} (c_g=30);"
           f"  16L_AdS/k = {X:.4g};  k_AdS - k(0) relative ~ {dk:.2g}")

# ---------------------------------------------------------------- (c) normalised test residues
def Nk(k, nu):
    a = k + 5
    if k % 2:
        return abs(-(32/(a - 1))*(2*nu + 1j)*(4*nu**2 + 2*k**2 - 9)/(4*nu**2 + 8*(a - 3)**2 + 1))
    c = {4: (216, 263, 81, 393), 6: (536, 2023, 169, 649)}[k]
    return 4*(16*nu**4 + c[0]*nu**2 + c[1])/((4*nu**2 + c[2])*(4*nu**2 + c[3]))
pr("(c) |test residue 16piG(nu^2+9/4)^{-k}/N_k(nu)| in the k-term, relative to 16piG Lam^2 L_AdS (NF02; no forward terms):")
phf = {k: (lambda kk: (lambda t: float(F.phi(kk, t))))(k) for k in KS}
for Dg in (20, 100, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); L = LAdS(Lam)
    levs = []
    for k in range(3, 8):
        v = quad(lambda t: abs(phf[k](t))*EPS[k]*Lam**(2*k)*(Lam**2*t*t + 2.25)**(-k)/Nk(k, Lam*t), 0, 1, limit=400, points=[1.5/Lam, 10/Lam])[0]
        levs.append(v/(Lam**2*L))
    # the same for Omega^U (M02) with the same N_k normalisation, plus its forward pieces
    from common import load as loadU
    from rcommon import peval
    FM, PM = loadU('M02')
    psi2m = lambda t: float(peval(FM.psi2, t)); LM = quad(lambda t: psi2m(t)/(t*t + 2.25/Lam**2), 0, 1, limit=400, points=[1.5/Lam])[0]
    levm = [quad(lambda t: abs(float(peval(PM[k], t)))*EPS[k]*Lam**(2*k)*(Lam**2*t*t + 2.25)**(-k)/Nk(k, Lam*t), 0, 1, limit=400,
                 points=[1.5/Lam, 10/Lam])[0]/(Lam**2*LM) for k in range(3, 8)]
    hh = 1e-4; g4 = lambda nu2: (nu2 + 2.25)**(-4)/Nk(4, np.sqrt(nu2))
    f4 = FM.a*Lam**10*abs((g4(hh) - g4(0.0))/hh)/(Lam**2*LM); f5 = 2.5*FM.a*Lam**10*(2.25)**(-5)/Nk(5, 0.0)/(Lam**2*LM)
    pr(f"    Dg={Dg}: NF02 k=3..7: " + " ".join(f"{v:.2g}" for v in levs) + "  (no forward)  |  Omega^U(M02) k=3..7: " + " ".join(f"{v:.2g}" for v in levm)
       + f", forward C4: {f4:.2g}, C5: {f5:.2g}")
    pr(f"      L_AdS(NF02)/L_AdS(M02) = {L/LM:.3g}  (the per-own-bound normalisation: a ratio R of leverages = (L_NF/L_M02) x (shape factor R L_M02/L_NF);"
       f" shape factors k=3..7: " + " ".join(f"{vm/v/(L/LM):.2g}" for v, vm in zip(levs, levm)) + ")")

# ---------------------------------------------------------------- (d) residual flip of the even replacements, NF02 weights
from c1_flip_contamination import K
def rho_sc(nu, k):
    c = {4: (210, 216, 263), 6: (1890, 536, 2023)}[k]
    return c[0]*(2*nu - 5j)*(2*nu - 1j)/((2*nu + 3j)*(2*nu + 7j)*(16*nu**4 + c[1]*nu**2 + c[2]))
from numpy.polynomial.legendre import leggauss
xg, wg = leggauss(48); pn = (xg + 1)/2; pw = wg/2
from nf_common import flat_parts_nf
pr("(d) residual one-pair flip of Psi_4/N^ff_4, Psi_6/N^ff_6 (NF02 weights; superconv/c2 model, leading Regge; bound / N_s):")
pf, wf = leggauss(600); pf = (pf + 1)/2; wf = wf/2
worst = 0.0
for Dg in (100, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); worst = 0.0
    for h in ('++', '+-'):
        for mfl in (1.05, 1.2, 1.4, 2.0):
            u = 1/mfl**2
            Js, B2, B3, parts = flat_parts_nf(F, h, u, 12, pf, wf)
            for J in ((0, 2, 4, 6) if h == '++' else (4, 6, 8)):
                st = state(Dg, Lam, mfl, J)
                if st is None:
                    continue
                i = int(np.searchsorted(Js, J)); Ns = abs(B2[i]) + abs(B3[i])
                fl = sum(pw[j]*u**(k - 2)*abs(float(F.phi(k, pn[j])))*2*EPS[k]*abs(rho_sc(Lam*pn[j], k))*abs(K(Lam*pn[j], st['beta'], 4, 0))/2
                         for j in range(len(pn)) for k in (4, 6))
                worst = max(worst, fl/Ns)
    pr(f"    Dg={Dg}: max smeared flip / N_s over m/Lam in {{1.05,1.2,1.4,2}}, small J: {worst:.2g}")
open(os.path.join(HERE, 'cache', 'nf5.log'), 'w').write("\n".join(out) + "\n")
