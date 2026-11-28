#!/usr/bin/env python3
"""Exact verifier for the angular Racah dual's two-ladder feasibility proof.

Theorem domain: g >= 4, n = floor(g/2), Delta >= J+g real;
b=0: J even >=0; b=4: J any integer >=4.

Run: python verify_angular_racah_feasibility.py
Requires Python 3 and SymPy. Every mathematical assertion uses exact symbolic
or integer arithmetic. No spin/twist scan, floating point, or asymptotic
approximation is used. The generated ledger is a complete finite certificate
for an inequality on a continuous unbounded domain.

Source conventions pinned to anshM123/MAJOR commit
01b01e6811825e446632d23ee8b33aab15b52b59, ads4/code/angular_opt/.
"""
from pathlib import Path
import json
import sympy as S

OUT = Path(__file__).resolve().parent
CHECKS = []


def check(label, expression):
    assert S.cancel(expression) == 0, label
    CHECKS.append(label)
    print("PASS", label, flush=True)


def positive_coefficients(label, polynomial):
    assert all(c > 0 for c in polynomial.coeffs()), label
    CHECKS.append(label)
    print("PASS", label, "terms", len(polynomial.terms()), flush=True)


# 1. Exact all-degree central-factorial identity, for arbitrary s and k.
x, s, k = S.symbols("x s k")
bx = (x+1)*(x+s)/((2*x+s)*(2*x+s+1))
dx = x*(x+s-1)/((2*x+s-1)*(2*x+s))
plus_ratio = (x+1)*(x+s+k)/((x-k+1)*(x+s))
minus_ratio = (x-k)*(x+s-1)/(x*(x+s+k-1))
check("arbitrary-degree central-factorial identity",
      bx*(plus_ratio-1)+dx*(minus_ratio-1)
      - k*k/((x-k+1)*(x+s+k-1)))

# These ratios are for P_k(x(x+s)) = x_down_k (x+s)_k.
# Multiplying by P_k removes apparent poles and proves a polynomial identity.
X, Y, q = S.symbols("X Y q")
r = (k-1)*(k+s-1)
check("ladder telescoping factor", (Y-X)/(Y-r) - (1-(X-r)/(Y-r)))
check("coefficient telescoping eigenvalue", ((k*k-q)-k*k)+q)


def finite_numerator(j, helicity_s=9):
    """Polynomial numerator of the finite sum, in q=n^2 and Y=N(N+s)."""
    den = S.prod(Y-i*(i+helicity_s) for i in range(j))
    total, ck, ak = den, S.Integer(1), S.Integer(1)
    for ell in range(1, j+1):
        ck *= ((ell-1)**2-q)/S.Integer(ell*ell)
        ak *= j*(j+helicity_s)-(ell-1)*(ell+helicity_s-1)
        total += ck*ak*S.prod(Y-i*(i+helicity_s) for i in range(ell, j))
    return S.expand(total), S.expand(den)


def finite_y(j, helicity_s=9):
    num, den = finite_numerator(j, helicity_s)
    return S.cancel(num/den)


# 2. Helicity-zero energy and physical weight comparisons.
c, u, v, rr = S.symbols("c u v rr")
energy0 = lambda aa, bb: (aa*aa+bb*bb-2*c*aa*bb)/(1-c)
w = (1+rr)*c*v-rr*u
check("b0 transfer contraction identity",
      energy0(u, v)-energy0(v, w)-(1-rr*rr)*(u-c*v)**2/(1-c))
check("b0 energy parameter derivative",
      S.diff(energy0(u, v), c)-(u-v)**2/(1-c)**2)
T, C, g, J, n = S.symbols("T C g J n")
c0, cJ = 1-2*q/T, 1-2*q/(T-C)
check("b0 final energy ratio",
      (1+c0)/(1+cJ)-1-q*C/(T*(T-C-q)))
B = g*g+5*g-S.Rational(3, 2)
Tmin = (J+g-1)*(J+g-2)
check("b0 denominator lower bound identity",
      Tmin-J*(J+1)-((2*g-4)*J+(g-1)*(g-2)))
check("b0 first weight slack",
      Tmin-B-(J-4)*(J+2*g+1)-S.Rational(15, 2))
check("b0 second weight slack",
      2*(Tmin+S.Rational(9, 2)-B)-J*(J+1)
      -(J-4)*(J+4*g-3)-4)
check("b0 spin-two polynomial",
      finite_y(2, 1)-(1-6*q/Y+6*q*(q-1)/(Y*(Y-2))))
check("b0 spin-two weight margin",
      (g*(g+1)+S.Rational(9, 2))/3-B/6
      -(g*g-3*g+S.Rational(21, 2))/6)

# 3. General Sonin energy algebra. q and q1 here are recurrence coefficients.
yy, vv, p, hh, hh1, qq, qq1 = S.symbols("yy vv p hh hh1 qq qq1")
yp = yy+vv/p
vp = vv-qq1*yp
E = yy*yy+hh*yy*vv+hh/qq*vv*vv
Ep = yp*yp+hh1*yp*vp+hh1/qq1*vp*vp
dec = p*p*(hh1/qq1-hh/qq)+p*(hh-hh1)/2
check("b4 exact Sonin energy decrement",
      (Ep-E-dec*vv*vv/(p*p)).subs(hh1, 2/p-hh))
check("b4 energy square completion",
      E-(1-hh*qq/4)*yy*yy-hh/qq*(vv+qq*yy/2)**2)
z = S.symbols("z")
Hbound = 2/(x+1)+2/(x+9)-1/(x+2)-1/(x+10)
positive_coercivity = 2*z**4+81*z**3+1075*z*z+4806*z+3528
check("b4 uniform coercivity polynomial",
      (S.Rational(4, 3)-(2*x+9)*Hbound/4).subs(x, z+5)
      - positive_coercivity/(6*(z+6)*(z+7)*(z+14)*(z+15)))
positive_coefficients("b4 coercivity numerator positive", S.Poly(positive_coercivity, z))
check("b4 B_j >= 1/8 certificate",
      ((x+1)*(x+9)/((2*x+9)*(2*x+10))-S.Rational(1, 8)).subs(x, z+1)
      -(2*z*z+25*z+14)/(8*(z+6)*(2*z+11)))
check("b4 spectral denominator margin",
      (2*n-2)*(2*n+7)-4*n*n-(10*n-14))
check("b4 h5 elementary upper bound", Hbound.subs(x, 5)-S.Rational(4, 15))

# 4. Independently construct y_j through the spin recurrence.
ys = [S.Integer(1), 1-10*q/Y]
for j in range(1, 6):
    pj = S.Rational((j+1)*(j+9), 2*j+10)
    pm = S.Rational(j*(j+8), 2*j+8)
    qj = (2*j+9)*q/(Y-j*(j+9))
    ys.append(S.cancel(ys[-1]+(pm*(ys[-1]-ys[-2])-qj*ys[-1])/pj))
for j in range(7):
    check(f"finite sum agrees with independent recurrence at x={j}", ys[j]-finite_y(j))

# 5. Exact low-spin b4 certificates x=0,1,2,3,4.
check("b4 x=0", ys[0]-1)
check("b4 x=1", 1-ys[1]-10*q/Y)
check("b4 x=2", 1-ys[2]-22*q*(Y-3*q-7)/(Y*(Y-10)))
B3 = 36*(Y-10)*(Y-22)-234*(q-1)*(Y-22)+364*(q-1)*(q-4)
check("b4 x=3", 1-ys[3]-q*B3/(Y*(Y-10)*(Y-22)))
Y3 = (2*n+1)*(2*n+10)
check("b4 x=3 initial numerator",
      B3.subs({q: n*n, Y: Y3})
      -4*(n**4+297*n**3+4405*n*n-1089*n-338))
check("b4 x=3 initial derivative",
      S.diff(B3, Y).subs({q: n*n, Y: Y3})-18*(3*n*n+88*n-11))
positive_coefficients("b4 x=3 initial numerator positive",
    S.Poly((n**4+297*n**3+4405*n*n-1089*n-338).subs(n, z+1), z))
check("b4 x=3 second derivative", S.diff(B3, Y, 2)-72)
Q4 = (2*Y**3-21*Y**2*q-115*Y**2+70*Y*q**2+868*Y*q+1806*Y
      -70*q**3-1540*q**2-7462*q-6768)
D4 = Y*(Y-10)*(Y-22)*(Y-36)
check("b4 x=4 exact formula", 1-ys[4]-26*q*Q4/D4)
m, t = S.symbols("m t")
R4 = (m**4*(38*t*t+462*t+1578)
      +m**3*(152*t**3+3332*t*t+23036*t+52080)
      +m*m*(99*t**4+3486*t**3+40027*t*t+182448*t+281098)
      +m*(24*t**5+1176*t**4+20368*t**3+150612*t*t+480668*t+546588)
      +2*t**6+126*t**5+2987*t**4+32844*t**3+169559*t*t+402738*t+356184)
NN = 2*(m+2)+2+t
check("b4 x=4 sum-of-squares identity",
      Q4.subs({q:(m+2)**2, Y:NN*(NN+9)})-2*m**4*(m-2*t-7)**2-R4)
positive_coefficients("b4 x=4 remaining coefficients positive", S.Poly(R4, m, t))

# 6. Continuous-domain initial-energy certificate.
D5 = Y*(Y-10)*(Y-22)*(Y-36)*(Y-52)
A = S.cancel(D5*ys[5])
P = S.cancel(D5*(ys[4]-ys[5])/(18*q))
assert S.denom(A) == S.denom(P) == 1
A, P = S.expand(A), S.expand(P)
B5 = S.expand(4225*q*(Y-70)*P*P+1235*q*A*P)
h = S.symbols("h")
E5 = (ys[5]**2+h*S.Rational(4225, 324)*(Y-70)/(19*q)*(ys[5]-ys[4])**2
      -h*S.Rational(65, 18)*ys[5]*(ys[5]-ys[4]))
check("b4 E5 cleared rational identity", E5-(19*A*A+h*B5)/(19*D5*D5))
u, t = S.symbols("u t")
certificates = []
for H in (0, 4):
    F = 190*(2*Y+25)**2*D5**2-(8*n*n+36*n+25)**2*(285*A*A+H*B5)
    gap = (S.Rational(2, 3)*(Y+S.Rational(25, 2))**2*19*D5**2
           -(4*n*n+18*n+S.Rational(25, 2))**2*(19*A*A+S.Rational(H, 15)*B5))
    check(f"b4 certificate normalization H={H}", F-60*gap)
    raw = S.Poly(S.expand(F), q, Y, n)
    ns = S.Poly(u+1, u, t)
    Ys = S.Poly(4*(u+1)**2+30*(u+1)+36+t, u, t)
    out = S.Poly(0, u, t)
    for (iq, iy, inn), coefficient in raw.terms():
        out += coefficient*ns**(2*iq+inn)*Ys**iy
    support = {(a, b) for b in range(13) for a in range(25-2*b)}
    assert set(out.monoms()) == support
    assert len(support) == 169
    positive_coefficients(f"b4 continuous-domain certificate H={H}", out)
    assert min(out.coeffs()) == 760
    rows = []
    for b in range(13):
        row = {a: cc for (a, bb), cc in out.terms() if bb == b}
        assert min(row.values()) == row[24-2*b]
        rows.append({"t_power": b, "minimum_coefficient": str(min(row.values())),
                     "minimum_at_u_power": 24-2*b})
    certificates.append({"H": H, "h_denominator": 15,
        "substitutions": {"n": "u+1", "Y": "4*n^2+30*n+36+t", "q": "n^2"},
        "domain": "u>=0, t>=0", "term_count": 169,
        "minimum_coefficient": "760", "row_minima": rows,
        "terms": [{"u": mon[0], "t": mon[1], "coefficient": str(cc)}
                  for mon, cc in out.terms()]})

# 7. Gap conversion, objective, and the polynomial degree restriction.
check("minimum Y at tail entry", (2*n+3)*(2*n+12)-(4*n*n+30*n+36))
check("maximum gap weight denominator",
      B.subs(g, 2*n+2)-(4*n*n+18*n+S.Rational(25, 2)))
D2sq = (S.rf(4, 2)/S.factorial(2))**2
check("dual objective coefficient", (n*n*(n*n-1)/4)/D2sq-n*n*(n*n-1)/400)

payload = {"theorem": "Two-ladder feasibility for real g>=4, integer n=floor(g/2), and real Delta>=J+g.",
    "source_repository": "https://github.com/anshM123/MAJOR",
    "source_commit": "01b01e6811825e446632d23ee8b33aab15b52b59",
    "proof_method": "Exact recurrence energies plus integer polynomial certificates.",
    "energy_lemma": "E5 < (2/3)*[(Y+25/2)/(4*n^2+18*n+25/2)]^2 for n>=1, Y>=4*n^2+30*n+36, 0<=h<=4/15.",
    "D5": str(D5), "A": str(A), "P": str(P),
    "B5_definition": "4225*q*(Y-70)*P^2 + 1235*q*A*P",
    "energy_numerator": "19*A^2+h*B5", "energy_denominator": "19*D5^2",
    "certificate_definition": "190*(2*Y+25)^2*D5^2-(8*n^2+36*n+25)^2*(285*A^2+H*B5)",
    "certificates": certificates,
    "verification": {"exact_checks_passed": len(CHECKS), "checks": CHECKS,
                     "sympy_version": S.__version__}}
(OUT / "coefficient_ledger.json").write_text(json.dumps(payload, indent=2)+"\n")
print(f"VERIFIED: {len(CHECKS)} exact checks; 338 positive integer coefficients.")
print("The coefficient ledger is coefficient_ledger.json.")
