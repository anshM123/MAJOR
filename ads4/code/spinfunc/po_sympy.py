"""
po_sympy (PARITY_ODD_CLOSE.md sect. 1.3): symbolic (sympy) checks of the algebra that turns the single-S^0-point structure of the 3x3
[E1, E2, O] action into Gamma_OO == Gamma_++, and its numerical instantiation on the oo_p3 K2 rows.
 L1  helicity parametrisation: same-helicity structures H(++), H(--) (parity swaps them), E_s = H++ + H--, O = i c (H++ - H--)
     (12 side; the 34 side carries the conjugate phase).  For ANY bilinear functional W (complex entries) and its parity image,
     the S^0-summed OO and E_s E_s entries are  OO = 2|c|^2 (D - X),  E_sE_s = 2 (D + X),  D = W++,++ + W--,--,  X = W++,-- + W--,++.
     => OO / E_sE_s is a nu-independent constant  <=>  X = 0  (no ++ <-> -- flip), given D, X independent functions.
 L2  order 0: if M0_EE = l1 u1 u1'^T + l2 u2 u2'^T (rank-1 helicity decomposition) and the O row of M0 is prop. to u1'^T (odd structure
     couples only to channel 1), the row null vector w (M0_{O,.} = w^T M0_{E,.}) satisfies w.u2 = 0, i.e. w is the dual vector of channel 1.
 L3  if w = alpha t1*, v = beta t1* (t1* = row 1 of T0^-1) and M_{OO} = w^T M_EE v at every nu and order, then
     M_OO = alpha beta (T0^-1 M_EE T0^-T)_{11}, and with Gamma = m^2[G/(a env) - D] the fitted a absorbs alpha beta (state/reference):
     Gamma_OO == Gamma_++ and a_O/a_++ = kappa(Delta)/kappa(Delta_ref).
Usage: python3 po_sympy.py  -> cache/po_sympy.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import numpy as np

out = []


def pr(s=''):
    print(s, flush=True); out.append(s)


# ---------------- L1
Wpp, Wpm, Wmp, Wmm, c = sp.symbols('W_pp_pp W_pp_mm W_mm_pp W_mm_mm c')
cc = sp.conjugate(c)
# single n0 bilinear form on (H++, H--) x (H++, H--); parity image swaps ++ <-> --
W1 = sp.Matrix([[Wpp, Wpm], [Wmp, Wmm]])
W2 = sp.Matrix([[Wmm, Wmp], [Wpm, Wpp]])
Es = sp.Matrix([1, 1]); O12 = sp.I*c*sp.Matrix([1, -1]); O34 = -sp.I*cc*sp.Matrix([1, -1])
EE = sp.expand((Es.T*(W1 + W2)*Es)[0]); OO = sp.expand((O12.T*(W1 + W2)*O34)[0])
D, X = Wpp + Wmm, Wpm + Wmp
pr('L1: S0-summed  E_sE_s = ' + str(sp.factor(EE)) + ' ;  OO = ' + str(sp.factor(OO)))
assert sp.simplify(EE - 2*(D + X)) == 0 and sp.simplify(OO - 2*c*cc*(D - X)) == 0
pr('    E_sE_s = 2(D + X), OO = 2|c|^2 (D - X)  [checked]; OO - |c|^2 E_sE_s = ' + str(sp.factor(sp.expand(OO - c*cc*EE))))
# single-n0 rank deficiency (W-- rows/cols = 0 at order 0) => X = 0 at order 0
pr('    order 0 (single n0 sees only ++: W--,-- = W++,-- = W--,++ = 0 at n0 = +1): X = ' + str(X.subs({Wpm: 0, Wmp: 0})))

# ---------------- L2
l1, l2 = sp.symbols('lambda1 lambda2', nonzero=True)
u1 = sp.Matrix(sp.symbols('u11 u12')); u2 = sp.Matrix(sp.symbols('u21 u22'))
u1p = sp.Matrix(sp.symbols('v11 v12')); u2p = sp.Matrix(sp.symbols('v21 v22'))
g = sp.symbols('gamma', nonzero=True)
M0EE = l1*u1*u1p.T + l2*u2*u2p.T
rowO = l1*g*u1p.T
w = sp.Matrix(sp.symbols('w1 w2'))
sol = sp.solve(list(w.T*M0EE - rowO), list(w), dict=True)[0]
wv = sp.Matrix([sol[w[0]], sol[w[1]]])
pr('L2: w solving w^T M0_EE = O-row:  w.u2 = ' + str(sp.simplify((wv.T*u2)[0])) + ' ;  w.u1 = ' + str(sp.simplify((wv.T*u1)[0])))
assert sp.simplify((wv.T*u2)[0]) == 0

# ---------------- L3
T = sp.Matrix(2, 2, sp.symbols('T11 T12 T21 T22')); Ti = T.inv()
A = sp.Matrix(2, 2, sp.symbols('A11 A12 A21 A22'))
al, be = sp.symbols('alpha beta')
wv = al*Ti[0, :].T; vv = be*Ti[0, :].T
lhs = (wv.T*A*vv)[0]; rhs = al*be*(Ti*A*Ti.T)[0, 0]
pr('L3: w^T A v - alpha beta (T^-1 A T^-T)_11 = ' + str(sp.simplify(lhs - rhs)))
assert sp.simplify(lhs - rhs) == 0
ks, kr, Gpp, m, env, Dc, a = sp.symbols('kappa_s kappa_r G_pp m env D a', positive=True)
GOO = ks/kr*Gpp
aO = a*ks/kr
GamO = m**2*(GOO/(aO*env) - Dc); Gampp = m**2*(Gpp/(a*env) - Dc)
pr('    Gamma_OO - Gamma_++ with a_O = a kappa_s/kappa_r: ' + str(sp.simplify(GamO - Gampp)))

# ---------------- numerical instantiation on the rows (po_identity output)
pr('Instantiation: see cache/po_identity.txt (kappa nu-constancy and kappa_s/kappa_r = a_O/a_++ on every oo_p3 row).')
open(os.path.join(HERE, 'cache', 'po_sympy.txt'), 'w').write("\n".join(out) + "\n")
