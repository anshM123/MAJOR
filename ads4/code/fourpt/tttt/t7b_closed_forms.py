"""
Closed forms [ID] of the exact Fock MFT data of t7_fock_l2 (cache/fock_l1.pkl, cache/fock_l2.pkl), in CHL's normalisation:
   T^{E,h} = C_l(Delta,J) X_h,  C_l = (-1)^l 2^{3l-1/2} 2^{J/2} ((2-Delta+J)/2)_{2l} ((Delta+J+1-4l)/2)_{2l}   (chl_xspace; all pi's cancel),
   claim:  p^MFT[T^{E,h}] = A_h P_FK(1,1; Delta, J) / K_l,   K_1 = 16, K_2 = 2304 = 48^2,   for n >= 1
           (n = 0: the two helicity structures coincide, one primary; the diagonal entries are then 2x the formula),
A_opp = (-J)_{2l}/(J+1)_{2l}, A_same = (Delta-1)_{2l}/(2-Delta)_{2l};  P_FK(1,1) = MFT coefficient of two Delta = 1 scalars (Fitzpatrick-Kaplan).
At double twist C_l has a simple zero for n < 2l ((2-Delta+J)/2 = -n) that cancels the simple pole of X_h (order eps^-1, t7); the limit is
taken with Delta = 6 + N - 2 eps (l = 2) exactly as in t7.
Output: out_t7b_closed_forms.txt
"""
import os, sys, pickle
import sympy as sp
from chl_xspace import EPS, _lead
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
def PFK(D1, D2, Delta, J):
    n = sp.Rational(Delta - D1 - D2 - J, 2); h = sp.Rational(3, 2); rf = sp.rf
    return (rf(D1-h+1, n)*rf(D2-h+1, n)*rf(D1, J+n)*rf(D2, J+n)/(sp.factorial(n)*sp.factorial(J)*rf(J+h, n)*rf(D1+D2+n-2*h+1, n)
            *rf(D1+D2+2*n+J-1, J)*rf(D1+D2+n+J-h, n)))
K = {1: 16, 2: 2304}
pr("# closed form p[T^{E,h}] = A_h P_FK(1,1)/K_l (fourpt/tttt/t7b_closed_forms.py); entries: (l, N, J, h): p_exact / formula")
for l in (1, 2):
    res = pickle.load(open(os.path.join(HERE, 'cache', 'fock_l%d.pkl' % l), 'rb'))
    allok = True; cnt = 0
    for key, r in sorted(res.items()):
        if r is None: continue
        _, N, Jv = key; Dv = 2*(l + 1) + N; n = sp.Rational(2*l + N - Jv, 2)
        for h in r['names']:
            # T = C_l X_h,  X_h = eps^k S_h (leading),  C_l ~ c eps^{kc}:  p[T] = p[S]/|c|^2 (orders must cancel)
            D = Dv - 2*EPS
            kc, cc = _lead(sp.expand_func(sp.rf((2 - D + Jv)/2, 2*l))*sp.rf(sp.Rational(Dv + Jv + 1 - 4*l, 2), 2*l))
            assert kc + r['k'][h] == 0, (key, h)
            c2 = 2**(6*l - 1)*2**Jv*cc**2
            pT = r['p'][h]/c2
            A = sp.rf(-Jv, 2*l)/sp.rf(Jv + 1, 2*l) if h == 'opp' else sp.rf(Dv - 1, 2*l)/sp.rf(2 - Dv, 2*l)
            ratio = sp.nsimplify(pT/(A*PFK(1, 1, Dv, Jv)/K[l]))
            expect = 2 if n == 0 else 1
            ok = ratio == expect; allok &= ok; cnt += 1
            pr(f"  l={l} N={N} J={Jv} n={n} {h}: p[T^E] = {pT};  ratio to formula = {ratio} {'(n = 0: sum of both helicities)' if n == 0 else ''} {'OK' if ok else 'FAIL'}")
    pr(f"  l={l}: {cnt} entries, all consistent with the closed form: {allok}")
open(os.path.join(HERE, 'out_t7b_closed_forms.txt'), 'w').write("\n".join(out) + "\n")
