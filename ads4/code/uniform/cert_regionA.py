"""
cert_regionA: Theorem-1-type positivity of the B2 Regge term of M02 (and M01) at large curved impact parameter.

Claim [certified, given (R2) of UNIFORM_POSITIVITY.md and correct mpmath.iv enclosures of + - * / sin cos]:
  for psi2 of M02 (resp. M01) and every M > 0, beta > 0 with M beta >= TAU0 (= 23.11 for M02, 23.49 for M01):
      F_M(beta) := int_0^1 psi2(p) P_{-1/2+iMp}(cosh beta) dp > 0.
Proof = Corollary A of UNIFORM_POSITIVITY.md with f -> psi2:
  (i)  phi = psi2/p > 0 and strictly decreasing on [0,1)            [exact Sturm counts, below]
       => Lemma 2: kernel K = S(Mt) + R_M(t), R_M > 0 (Lemma 1 applied to phi(p) * 2p/(e^{2 pi M p} - 1)).
  (ii) S(tau) = int psi2 sin(p tau) dp > 0 for tau >= TAU0           [interval cells on [TAU0, T*], analytic tail beyond T*]
  (iii) (2.1): F_M(beta) = (sqrt2/pi) int_beta^oo K(t)(cosh t - cosh beta)^{-1/2} dt, and M t >= M beta >= TAU0.
S closed form (integration by parts, deg psi2 = 8):  S(tau) = Im sum_j (a_j - e^{i tau} b_j) i^{j+1}/tau^{j+1},  a_j = psi2^{(j)}(0), b_j = psi2^{(j)}(1).
Cell lemma: |S'| <= L := int p psi2 (psi2 >= 0); if S(c) >= s > 0 then S > 0 on (c - s/L, c + s/L).
Output: cache/cert_regionA.log
"""
import os, sys, json
from fractions import Fraction as Fr
import sympy as sp
from mpmath import iv, mp
HERE = os.path.dirname(os.path.abspath(__file__))
iv.prec = 120
lines = []
def pr(s):
    print(s, flush=True); lines.append(s)

D = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_certified.json')))
p = sp.Symbol('p')
ok_all = True
for name, TAU0 in (('M02', Fr(2311, 100)), ('M01', Fr(2349, 100))):
    c = [sp.Rational(str(Fr(s))) for s in D[name]['psi2_exact']]
    f = sum(ci*p**i for i, ci in enumerate(c))
    x = [sp.Rational(str(Fr(s))) for s in D[name]['x']]
    g = sum(ci*p**i for i, ci in enumerate(x))
    assert sp.expand(f - (1 - p)**3*p*g) == 0
    n_g = sp.Poly(g, p).count_roots(0, 1)
    num = sp.expand(-3*g + (1 - p)*sp.diff(g, p))
    n_num = sp.Poly(num, p).count_roots(0, 1)
    pr(f"{name}: psi2 = (1-p)^3 p g(p); roots of g on [0,1]: {n_g} (g(0) = {g.subs(p, 0)} > 0); phi' = (1-p)^2 N(p), roots of N on [0,1]: {n_num} (N(0) = {num.subs(p, 0)} < 0)")
    okA = (n_g == 0 and g.subs(p, 0) > 0 and n_num == 0 and num.subs(p, 0) < 0)
    deg = sp.degree(f, p)
    a = [sp.diff(f, p, j).subs(p, 0) for j in range(deg + 1)]
    b = [sp.diff(f, p, j).subs(p, 1) for j in range(deg + 1)]
    L = sp.integrate(p*f, (p, 0, 1))
    pr(f"   a_j = {a}\n   b_j = {b}\n   L = int p psi2 = {L} = {float(L):.6f}")
    cyc = [(0, 1), (-1, 0), (0, -1), (1, 0)]       # i^{j+1} = cr + i ci for j = 0,1,2,3 (mod 4)
    q = lambda v: iv.mpf(int(sp.Rational(v).p))/iv.mpf(int(sp.Rational(v).q))
    aq = [q(v) for v in a]
    bq = [q(v) for v in b]
    def S_iv(t):
        t = iv.mpf(t)
        co, si = iv.cos(t), iv.sin(t)
        s = iv.mpf(0)
        for j in range(deg + 1):
            cr, ci = cyc[j % 4]
            s += ((aq[j] - bq[j]*co)*ci - bq[j]*si*cr)/t**(j + 1)
        return s
    # analytic tail: tau^{deg+1} S >= poly(tau) := non-oscillatory part - sum 2|b_j| tau^{deg-j}
    tau = sp.Symbol('tau', positive=True)
    nonosc = sum(a[j]*cyc[j % 4][1]*tau**(deg - j) for j in range(deg + 1))
    osc = sum(2*abs(b[j])*tau**(deg - j) for j in range(deg + 1))
    cpoly = sp.Poly(sp.expand(nonosc - sum(abs(a[j]*cyc[j % 4][1])*tau**(deg - j) for j in range(deg + 1) if j != 2) - osc
                              + 0*tau), tau)
    # (the j = 2 non-oscillatory term -a_2/tau^3 is the positive leading term; all others are bounded in absolute value)
    lead = -a[2]
    cpoly = sp.Poly(lead*tau**(deg - 2) - sum(abs(a[j]*cyc[j % 4][1])*tau**(deg - j) for j in range(deg + 1) if j != 2) - osc, tau)
    Tstar = 1
    while not (cpoly.eval(Tstar) > 0 and sp.Poly(cpoly, tau).count_roots(Tstar, None) == 0):
        Tstar += 1
    pr(f"   tail: tau^{deg+1} S(tau) >= c(tau) = {cpoly.as_expr()};  c > 0 on [T*, oo) with T* = {Tstar} (exact root count)")
    # cells on [TAU0, T*]
    Lf = q(L)
    cpt = mp.mpf(TAU0.numerator)/TAU0.denominator
    left = cpt; ncell = 0; okcells = True
    # first cell must contain TAU0: start at c = TAU0 + small, require c - s/L <= TAU0
    c0 = cpt
    while True:
        s = S_iv(c0)
        if not (s.a > 0):
            okcells = False; pr(f"   FAIL at tau = {c0}: S in {s}"); break
        r = mp.mpf((s.a/Lf.b).a)       # lower endpoint of the interval quotient (conservative)
        if ncell == 0 and not (c0 - r <= cpt):
            okcells = False; pr("   FAIL: first cell does not contain TAU0"); break
        if c0 + r >= Tstar:
            ncell += 1; break
        c0 = c0 + r*mp.mpf(0.98)
        ncell += 1
    pr(f"   cells on [{float(TAU0)}, {Tstar}]: {ncell}, all S(c) > 0 with radius s/L: {okcells}")
    okn = okA and okcells
    ok_all = ok_all and okn
    pr(f"VERDICT {name}: region (A) [M beta >= {float(TAU0)}]: {'PASS' if okn else 'FAIL'}")
pr(f"VERDICT ALL: {'PASS' if ok_all else 'FAIL'}")
open(os.path.join(HERE, 'cache', 'cert_regionA.log'), 'w').write("\n".join(lines) + "\n")
