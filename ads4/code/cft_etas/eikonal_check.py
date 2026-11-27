"""Leading-eikonal (high-momentum) graviton propagation in a type-D background for 16piG L = R + 6 + F(I, J).
Self-dual bivector block: background W = q diag(-2,1,1) (type D; the brane has q = -mu/(2 r^3)-type fixed sign outside the horizon),
fluctuation of a null plane wave at leading eikonal order: dW = v v^T with v.v = 0 (type N).  Second variation:
  d2F = F_I d2I + F_J d2J + F_II dI^2 + 2 F_IJ dI dJ + F_JJ dJ^2,
  dI = 2 tr(W dW),  dJ = 3 tr(W^2 dW),  and the k^4 parts of d2I, d2J are 2 tr(dW^2) and 3 tr(W dW^2), which vanish (dW^2 = 0).
So the k^4 (momentum-enhanced) term is  K(q) * [polarisation factor]^2 * k^4  with  K = 4 F_II tr(W dW)^2/... ; computed below."""
import sympy as sp
q, a_, b_, c_, g2, g3, M = sp.symbols('q a b c g2 g3 M', real=True)
v1, v2, v3 = sp.symbols('v1 v2 v3')
W = q*sp.diag(-2, 1, 1)
v = sp.Matrix([v1, v2, v3])
null = {v3: sp.sqrt(-v1**2 - v2**2)}
dW = v*v.T
print('tr(dW^2) on null v =', sp.simplify((dW*dW).trace().subs(null)), '; tr(W dW^2) =', sp.simplify((W*dW*dW).trace().subs(null)))
dI = sp.simplify((2*(W*dW).trace()).subs(null))
dJ = sp.simplify((3*(W*W*dW).trace()).subs(null))
print('dI =', sp.factor(dI), '  dJ =', sp.factor(dJ))
FII, FIJ, FJJ = sp.symbols('F_II F_IJ F_JJ')
K4 = sp.factor(sp.expand(FII*dI**2 + 2*FIJ*dI*dJ + FJJ*dJ**2))
print('k^4 coefficient (self-dual block) =', K4)
Kq = sp.factor(sp.simplify(K4/(36*q**2*v1**4)))
print('  => K(q) =', Kq, '   [times 36 q^2 v1^4; real waves: v1^4 -> 2 Re(v1^2)^2-type square >= 0]')
# operators: F = a I^2 + b IJ + c J^2  -> F_II = 2a, F_IJ = b, F_JJ = 2c
Kop = sp.expand(Kq.subs({FII: 2*a_, FIJ: b_, FJJ: 2*c_}))
print('F = a I^2 + b IJ + c J^2 :  K(q) =', Kop)
print('pure beta5 IJ (a=c=0):   K(q) =', Kop.subs({a_: 0, c_: 0}), '  -> one sign for all polarisations at fixed q-sign')
ex = sp.factor(Kop.subs({a_: g2**2/(2*M**2), b_: g2*g3/M**2, c_: g3**2/(2*M**2)}))
print('heavy scalar exchange (g2 X I + g3 X J):  K(q) =', ex, '  (perfect square: causal)')
# causality K(q) >= 0 for all q in the exterior range <=> b bounded by a, c:
qs = sp.symbols('qs', positive=True)
Kpos = Kop.subs(q, qs)   # take the sign convention where exterior q>0 after orientation; minimise over qs
qmin = sp.solve(sp.diff(Kpos, qs), qs)
print('min over q>0 at q* =', qmin, ' ; K(q*) =', [sp.simplify(Kpos.subs(qs, s)) for s in qmin])
print('=> K(q) >= 0 for all q>0  <=>  b <= (2/3)*... ; at the optimum the condition reads b^2 <= 4ac (Remark 4) when q* lies in range')
# Brane orientation: q(r) = -1/(2 r^3) (Remark 7: I = 48 q^2, horizon q = -1/2), so the exterior is q in [-1/2, 0); s = |q|.
s = sp.symbols('s', positive=True)
Kb = sp.expand(Kop.subs(q, -s))
print('on the brane exterior (s = |q| in (0, 1/2]):  K =', Kb)
print('  pure IJ: K = 3 b s  => causal iff b >= -O(slack): a LOWER bound only; b > 0 (which lowers eta/s by -1944 b) is unconstrained')
print('  with a, c: interior minimum at s* = -b/(3c) gives b^2 <= 4ac; endpoint s = 1/2 gives b >= -(4/3) a - (3/4) c')
print('  endpoint check:', sp.solve(sp.Eq(Kb.subs(s, sp.Rational(1, 2)), 0), b_))
# gradient operator Y = grad I . grad J : leading eikonal part of its second variation is 2 (k.k) dI dJ = 0 for null k
