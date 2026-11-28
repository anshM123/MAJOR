"""v4: the 3-state configuration A=(b0,J=2), B=(b4,J=4), C=(b4,J=6), all at twist exactly g, in closed form.
Derived by hand (README), checked here symbolically:  |c|^2 = rho_C r1(C) [75/(T_C-30) - 25/(T_A-2)],
rho_A = rho_C r1(C)/r1(A), rho_B = rho_A - rho_C, rho's fixed by the endpoint budget 3.
Also: the paper's k=2-only bound and the light/contact sensitivity factors."""
import sympy as sp

g = sp.symbols('g', positive=True)


def T(D): return (D - 1) * (D - 2)


def Wt(b, D): return (T(D) - b * (b - 1)) / sp.Integer(3) + sp.Rational(3, 2)


TA, TB, TC = T(g + 2), T(g + 4), T(g + 6)
r1A = 16 * sp.Integer(6) / TA                      # b=0, J=2: s=0
r2A = r1A * sp.Rational(25, 4) * (6 - 2) / (TA - 2)  # s=1
r1C = 16 * sp.Integer(42 - 20) / (TC - 20)          # b=4, J=6: s=4
r2C = r1C * sp.Rational(25, 4) * (42 - 30) / (TC - 30)
rC = 1
rA = r1C / r1A
rB = rA - rC
s = 3 / (rA * Wt(0, g + 2) + rB * Wt(4, g + 4) + rC * Wt(4, g + 6))
c2 = sp.simplify(s * (rC * r2C - rA * r2A))
print("|c|^2(g) =", sp.factor(c2))
ser = sp.series(c2 * g ** 6, g, sp.oo, 3)
print("|c|^2 g^6 =", ser)
lead = sp.limit(c2 * g ** 6, g, sp.oo)
print("lim |c|^2 g^6 =", lead, " -> tau g^3 -> 4 sqrt(lead) =", sp.N(4 * sp.sqrt(lead), 8))
print("rho_B > 0 for all g>0:", sp.simplify(rB), ">0 ?", [sp.N(rB.subs(g, x)) for x in (1, 5, 10, 100)])
for gv in [10, 20, 40, 80, 100, 1000]:
    K = sp.Rational(729, 4800) * (gv - 2) ** 2 * (1 - sp.Rational(49, (2 * gv + 9) ** 2))
    print(f"g={gv}: 3-state tau = {sp.N(4*sp.sqrt(c2.subs(g, gv)), 6)}  (tau g^3 = {sp.N(4*sp.sqrt(c2.subs(g, gv))*gv**3, 6)});"
          f"  paper k=2-only tau_* = sqrt(48/K) = {sp.N(sp.sqrt(48/K), 6)}")
