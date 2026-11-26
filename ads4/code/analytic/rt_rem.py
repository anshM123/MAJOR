"""
rt_rem.py: Lemma Loc (TAIL_CERT §4) remainder constants for Route T' at moderate L, both kernels, several truncation orders J0.
Computes, with TAIL_CERT's certified K-functions (tc_remainder.logK_real / logK_small, imported read-only):
  k = 2 (dressed, full realised f):  lI2(J0) >= log int_[0.4, 6000] |f| K_3(nu) dnu,  small-nu part lK0 + L + log(0.4 sup|f|)
  k = 3 (psi3):                       lI3(J0) >= log int_0^1 K_4(Lam p) dp (p-cells of 1/1000, nu >= 0.4),  small part lK0_4 + L
and then  log E_rho2 <= 3 log d + log(2/Lam) + log2 + log_loc(m, L) + logaddexp(lI2, small2),             m = J0 - 4
          log E_rho3 <= 3 log d + log(d^2/(3.9 a^2)) + 2 log Lam + log(sum|psi3 coeffs|) + log2 + log_loc(m3, L) + logaddexp(lI3, small3),  m3 = J0 - 5
(log_loc(m, L) = log 2.53 + (m-1)(log 1.0006 - L) + Q_m + log(1 + 40*1.5^m), valid for 2 <= m <= a - 3 and L >= ~6 (a/N <= 0.5003 w)).
Output cache/rt_rem.json {J0: {lI2, lK0_2, lI3, lK0_3}}.
"""
import os, sys, json, math
sys.path.insert(0, '/home/user/MAJOR/ads4/code/tailcert')
sys.path.insert(0, '/home/user/MAJOR/ads4/code/contour')
sys.dont_write_bytecode = True
import tc_remainder as RM
from an_common import CACHE
LAM = math.sqrt(998*999)

def compute(J0):
    m = J0 - 4; m3 = J0 - 5
    lK0, cL = RM.logK_small(3, J0, m, 7.0)
    edges = [0.4, 1.0] + [float(k) for k in range(2, 1301)]
    ls = []
    for n1, n2 in zip(edges[:-1], edges[1:]):
        ls.append(math.log(21.33815*(n2 - n1)) + RM.logK_real(n1, n2, 3, J0, m))
    import cb_common as KC, cb_lines as CL
    W = KC.load()
    for T in range(1300, 6000, 10):
        intf, _ = CL.tail(W, 0.0, T)
        ls.append(RM.lup(intf) + RM.logK_real(T, T + 10, 3, J0, m))
    mx = max(ls); lI2 = mx + math.log(sum(math.exp(x - mx) for x in ls))
    lK0_3, _ = RM.logK_small(4, J0, m3, 7.0)
    ls3 = []
    for k in range(1, 1000):
        n1, n2 = LAM*k/1000, LAM*(k + 1)/1000
        if n1 < 0.4: continue
        ls3.append(RM.logK_real(n1, n2, 4, J0, m3) + math.log(1/1000))
    mx = max(ls3); lI3 = mx + math.log(sum(math.exp(x - mx) for x in ls3))
    return dict(lI2=lI2, lK0_2=lK0, lI3=lI3, lK0_3=lK0_3)

if __name__ == '__main__':
    fn = os.path.join(CACHE, 'rt_rem.json')
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    for J0 in [int(x) for x in sys.argv[1:]]:
        if str(J0) in out: continue
        out[str(J0)] = compute(J0)
        print(J0, out[str(J0)], flush=True)
        json.dump(out, open(fn, 'w'))
