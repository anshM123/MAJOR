"""ANGULAR_MATCH item 2: the uniform constant K(g).

(a) W = |a_2|^2/|a_0|^2 on the b = 4 ladder from the ANGULAR_ROUTE §8 closed form and the §4 norms (symbolic J, Delta).
    a_k := a_{m=4+k}; a_2 is ladder index k = 2, i.e. O-helicity m = 6, harmonic h - m = -4 for source h = +2.
(b) The chain  R(J,Delta) := [(T-12)/3]/W  >=  K(g) = (729/4800)(g-2)^2 (1 - 49/(2g+9)^2)  for J >= 6, Delta - J >= g,
    proved symbolically step by step, then checked in exact rational arithmetic on a J-scan for several g.
(c) True minima (with and without the +3/2 of the endpoint weight) and the resulting bounds 4 sqrt(3/K).
Output: cache/am2_constant.txt
"""
import sys, os, time
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp

out = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'am2_constant.txt'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

t0 = time.time()
J, D, k = sp.symbols('J Delta k', positive=True)
# ---- (a) W from closed forms -------------------------------------------------------------------------------
F = lambda kk: sp.binomial(J-4, kk)*sp.rf(4, kk)/sp.rf(D-5-kk, kk)          # F_{4+k}/F_4  (ANGULAR_ROUTE §8)
N = lambda m: sp.binomial(2*J, J+m)*sp.gamma(D+m-1)*sp.gamma(D-m-1)        # N_m up to m-independent factor (§4)
def ratio2(kk):
    r = F(kk)**2*N(4)/N(4+kk)
    return sp.factor(sp.simplify(sp.combsimp(sp.gammasimp(sp.expand_func(r)))))
W1, W2 = ratio2(1), ratio2(2)
C = J*(J+1); T = (D-1)*(D-2)
W_claim = 100*(C-20)*(C-30)/((T-20)*(T-30))
P('(a) |a_1/a_0|^2 =', W1)
P('    |a_2/a_0|^2 =', W2)
P('    W - 100(C-20)(C-30)/((T-20)(T-30)) =', sp.simplify(W2 - W_claim))
P('    W - §6 form 100(J-5)(J-4)(J+5)(J+6)/((D-7)(D-6)(D+3)(D+4)) =',
  sp.simplify(W2 - 100*(J-5)*(J-4)*(J+5)*(J+6)/((D-7)*(D-6)*(D+3)*(D+4))))
# ladder-ratio form of §8: (a_{k+1}/a_k)^2 = ((k+4)/(k+1))^2 X_k/Y_k, s = 4+k
lad = 1
for kk in range(2):
    s = 4+kk; lad *= sp.Rational(kk+4, kk+1)**2*(C-s*(s+1))/(T-s*(s+1))
P('    same from the §8 normalised ladder ratio:', sp.simplify(lad - W_claim) == 0)

# ---- (b) the chain, symbolically ---------------------------------------------------------------------------
x, y, a, v, A = sp.symbols('x y a v A', positive=True)
# R = (T-12)(T-20)(T-30)/(300 (C-20)(C-30)); with x = J+1/2, y = Delta-3/2: C = x^2-1/4, T = y^2-1/4
Rxy = (y**2-sp.Rational(49, 4))*(y**2-sp.Rational(81, 4))*(y**2-sp.Rational(121, 4))/(300*(x**2-sp.Rational(81, 4))*(x**2-sp.Rational(121, 4)))
Rdirect = ((T-12)*(T-20)*(T-30)/(300*(C-20)*(C-30))).subs({J: x-sp.Rational(1, 2), D: y+sp.Rational(3, 2)})
P('(b) R = [(T-12)/3]/W in (x,y) variables, identity check:', sp.simplify(sp.together(Rxy - Rdirect)) == 0)
# step 1: (y^2-a)/(x^2-a) - y^2/x^2 = a(y^2-x^2)/(x^2(x^2-a)) >= 0 for y > x > sqrt(a)
P('    step 1: (y^2-a)/(x^2-a) - y^2/x^2 =', sp.factor((y**2-a)/(x**2-a) - y**2/x**2), ' (>= 0 for y>x>sqrt a; a = 81/4, 121/4 need x > 11/2, i.e. J >= 5)')
# step 2: 1 - 49/(4y^2) >= 1 - 49/(2g+9)^2 since 2y = 2(J+g)-3 >= 2g+9 for J >= 6
P('    step 2: 2y = 2J + 2g - 3 >= 2g + 9  <=>  J >= 6')
# step 3: min_{v>0} (v+A)^6/v^4 at v = 2A, value 729 A^2/16, with v = x, A = g-2 (y = x + g - 2)
fv = (v+A)**6/v**4
crit = sp.solve(sp.diff(sp.log(fv), v), v)
P('    step 3: critical point of (v+A)^6/v^4:', crit, ' value:', sp.simplify(fv.subs(v, crit[0])), '; second log-derivative there:',
  sp.simplify(sp.diff(sp.log(fv), v, 2).subs(v, crit[0])))
P('    => R >= (1/300)(1-49/(2g+9)^2)(729/16)(g-2)^2 = K(g)  for all real J >= 6, Delta-J >= g > 2  [all steps exact]')
# monotonicity in Delta: d/dT of (T-12)(T-20)(T-30) > 0 for T > 30
Ts = sp.symbols('T_s')
dp = sp.diff((Ts-12)*(Ts-20)*(Ts-30), Ts)
P('    monotone in Delta: d/dT (T-12)(T-20)(T-30) =', sp.expand(dp), '; real roots:', [sp.N(r, 6) for r in sp.solve(dp, Ts)], '(all < 30)')

def K(g):
    g = Fr(g); return Fr(729, 4800)*(g-2)**2*(1 - Fr(49)/(2*g+9)**2)
def Rex(Jv, g, plus=False):
    Dv = Fr(Jv) + Fr(g); Tv = (Dv-1)*(Dv-2); Cv = Fr(Jv*(Jv+1))
    Wv = 100*(Cv-20)*(Cv-30)/((Tv-20)*(Tv-30))
    num = (Tv-12)/3 + (Fr(3, 2) if plus else 0)
    return num/Wv
P('\n(b\') exact-rational scan, worst case Delta = J + g, J = 6 .. Jmax (beyond Jmax the step-3 bound is increasing in x)')
P('   g      Jmax     K(g)               min_J R (argmin)          min_J R+ (argmin)        R/K-1        all R>=K')
res = {}
for g in [10, 30, 100, 300, 1000, 3000]:
    Jmax = 12*g; Kg = K(g)
    best = (None, None); bestp = (None, None); ok = True
    for Jv in range(6, Jmax+1):
        r = Rex(Jv, g)
        if r < Kg: ok = False
        if best[0] is None or r < best[0]: best = (r, Jv)
        if Jv % 1 == 0 and (abs(Jv - 2*(g-2)) < 0.2*g + 20 or Jv < 40):
            rp = Rex(Jv, g, True)
            if bestp[0] is None or rp < bestp[0]: bestp = (rp, Jv)
    # beyond Jmax: lower bound (y^6/x^4)(1-49/(4y^2))/300 is increasing for x > 2(g-2) and exceeds K there
    res[g] = (Kg, best, bestp)
    P(f'  {g:5d} {Jmax:7d}  {float(Kg):16.6f}  {float(best[0]):16.6f} ({best[1]:5d})  {float(bestp[0]):16.6f} ({bestp[1]:5d})  {float(best[0]/Kg-1):.3e}  {ok}')
P('\n(c) bounds sqrt(t4^2+alpha_T^2) <= 4 sqrt(3/K) at D = L_g = 0')
mp.mp.dps = 30
for g in [100, 300, 1000, 3000]:
    Kg, best, bestp = res[g]
    P(f'  g={g:5d}:  4sqrt(3/K(g)) = {mp.nstr(4*mp.sqrt(3/mp.mpf(Kg.numerator)*Kg.denominator), 10)}'
      f'   with exact min R: {mp.nstr(4*mp.sqrt(3/mp.mpf(best[0].numerator)*best[0].denominator), 10)}'
      f'   with exact min incl. +3/2: {mp.nstr(4*mp.sqrt(3/mp.mpf(bestp[0].numerator)*bestp[0].denominator), 10)}'
      f'   asymptotic 160/(9g) = {mp.nstr(mp.mpf(160)/(9*g), 10)}')
Kg = res[1000][0]
P('  K(1000) exact =', Kg, '=', mp.nstr(mp.mpf(Kg.numerator)/Kg.denominator, 15))
P('  §9 constant 151250 gives', mp.nstr(4*mp.sqrt(mp.mpf(3)/151250), 10))
P('\ndone in %.1f s' % (time.time()-t0))
