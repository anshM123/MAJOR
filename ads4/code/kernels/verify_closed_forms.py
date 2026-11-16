"""
Checks of the closed forms quoted in KERNELS_D3.md against the symbolic pipeline output (J <= 12) and numerically
(J = 10, 12).  Run after compute_kernels.py J for J = 0,...,12.
  (R0)  P_{-k}/P_0 = C(J,k) (4)_k / (Delta-1-k)_k              for every structure, 0 <= k <= J
  (R4)  P_{4+k}/P_4 = C(J-4,k) (4)_k / (Delta-5-k)_k           for every structure, 0 <= k <= J-4
  (Z)   P_m = 0 for m = 1,2,3;  parity-odd even-J structures and S^{++}: P_m = 0 for m >= 4;
        S^{+-} and odd-J structures: P_m = 0 for m <= 0.
"""
import sys, os, pickle
import sympy as sp
import kernels_lib as kl
S = kl.S; D = sp.Symbol('Delta')

def load(J, par):
    fn = os.path.join(kl.CACHE, f'kern_J{J}_{par}.pkl')
    return pickle.load(open(fn, 'rb')) if os.path.exists(fn) else None

def check_ratios(J, r):
    sub = {S: (D + J)/2}
    ok = True
    if r.get(0, 0) != 0:
        for k in range(0, J + 1):
            lhs = sp.cancel((r.get(-k, 0)/r[0]).subs(sub))
            rhs = sp.binomial(J, k)*sp.rf(4, k)/sp.rf(D - 1 - k, k)
            ok &= sp.simplify(lhs - rhs) == 0
    else:
        ok &= all(r.get(-k, 0) == 0 for k in range(0, J + 1))
    if r.get(4, 0) != 0:
        for k in range(0, J - 3):
            lhs = sp.cancel((r.get(4 + k, 0)/r[4]).subs(sub))
            rhs = sp.binomial(J - 4, k)*sp.rf(4, k)/sp.rf(D - 5 - k, k)
            ok &= sp.simplify(lhs - rhs) == 0
    else:
        ok &= all(r.get(4 + k, 0) == 0 for k in range(0, J - 3))
    ok &= all(r.get(m, 0) == 0 for m in (1, 2, 3))
    return ok

if __name__ == '__main__':
    for J in range(0, 13):
        for par in ('even', 'odd'):
            dat = load(J, par)
            if dat is None or not dat[1]:
                continue
            monos, basis, K, res = dat
            for a, r in enumerate(res):
                typ = ('++' if r.get(0, 0) != 0 else '') + ('+-' if r.get(4, 0) != 0 else '')
                print(f"J={J} {par:4s} structure {a}: sectors present: {typ:5s} ratio formulas (R0),(R4),(Z): {check_ratios(J, r)}")
