"""
p2spin/s6_nu6cost.py (P2_SPINNING.md §7.1): Lemma-I action of the nu^6 structure that CHMRS's O(1/m^4) series contains
((eta^2-1) nu^6 P, coefficient -1/72 in the seed, validated in s5_order4), relative to Theorem W's K4 representative 4a = eta^2 nu^4 P,
at every sampled H_LJ(1000) state (lsw/cache/costs_1000.json; lsw_common imported read-only).  -> cache/s6_nu6cost.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lsw'))
import numpy as np
import lsw_common as W
R = json.load(open(os.path.join(HERE, '..', 'lsw', 'cache', 'costs_1000.json')))
out = []
for r in R:
    b = r['beta']; P, _ = W.P_and_dP(b); P = P.real
    ch = np.cosh(b); n2 = W.NU**2
    s4 = W.act(ch**2*n2**2*P); s6 = W.act((ch**2 - 1)*n2**3*P)
    out.append('J=%3d m/L=%.3f beta=%.4f : |act(nu^6 (eta^2-1) P)| / |act(eta^2 nu^4 P)|  R1 %.3g  R2 %.3g ;  (1/72)*ratio(R2) = %.3g'
               % (r['J'], r['mfl'], b, abs(s6[0]/s4[0]), abs(s6[1]/s4[1]), abs(s6[1]/s4[1])/72))
print('\n'.join(out)); open(os.path.join(HERE, 'cache', 's6_nu6cost.txt'), 'w').write('\n'.join(out) + '\n')
