"""t0check/chk_consistency.py (INDEP_CHECK_T0 §2): consistency of the objects used by T0's certificate.
 1. NF02 exact coefficients: prl/data/nf02_psi.json (normalised by k0) vs robust/functionals_noforward.json (Part F's input).
 2. psi used by the float hybrid evaluator (ls_core.NFw(Dg)) vs the exact psi, at Dg = 300, 1000, 5000 (Dg-independence).
 3. Part F's object: rcert.poly_u (delta = 0 and sign patterns) / k0  vs  flat_parts_nf (B2, B3) used in H (ls_core.hyb_fast).
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'robust')); sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
from fractions import Fraction as Fr
import numpy as np
import t0_common as T
import nf_cert, rcert as R
import ls_core as LS

d = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'nf02_psi.json')))
k0 = Fr(d['k0'])
F = nf_cert.load('NF02')
print('1. k0 =', k0, '; nf02_psi psi2 == functionals_noforward psi2:', [Fr(c) for c in d['psi2_exact']] == F[0],
      '; psi3:', [Fr(c) for c in d['psi3_exact']] == F[1], '; t0 PSI == unnorm/k0:', T.PSI[2] == [c/k0 for c in F[0]] and T.PSI[3] == [c/k0 for c in F[1]])
ps = np.linspace(0, 1, 101)
for Dg in (300, 1000, 5000):
    NF = LS.NFw(Dg)
    e2 = max(abs(NF.psi2(p) - float(T.pval(T.PSI[2], Fr(p).limit_denominator(10**9)))) for p in ps)
    e3 = max(abs(NF.psi3(p) - float(T.pval(T.PSI[3], Fr(p).limit_denominator(10**9)))) for p in ps)
    print(f'2. Dg={Dg}: max|psi2_float - psi2_exact| = {e2:.2e}, psi3: {e3:.2e}')
Q = LS.Q
worst = 0
for h in ('++', '+-'):
    Js = [0, 2, 4, 8, 20, 64] if h == '++' else [4, 5, 6, 9, 21, 63]
    for J in Js:
        for u in (Fr(1), Fr(1, 2), Fr(1, 7), Fr(1, 50)):
            pf, wf = Q.gl01(max(400, 2*J + 100))
            Js_, B2, B3, _ = Q.flat_parts_nf(LS.NFw(1000), h, float(u), J, pf, wf); i = int(np.searchsorted(Js_, J))
            g = R.gens(F[0], F[1], F[2], 0, h, 1)[0]
            P = R.poly_u(g, J)
            ex = float(R.pval(P, u)/k0)
            fl = B2[i] + B3[i]
            worst = max(worst, abs(ex - fl)/max(1e-12, abs(B2[i]) + abs(B3[i])))
            # sign patterns: FF_g(c) == (1 - c s2) B2 + (1 - c s3) B3
            for g in R.gens(F[0], F[1], F[2], Fr(199, 1000), h, 1):
                s2 = int(g.label.split('s2=')[1][:2]); s3 = int(g.label.split('s3=')[1][:2])
                ex = float(R.pval(R.poly_u(g, J), u)/k0)
                fl = (1 - 0.199*s2)*B2[i] + (1 - 0.199*s3)*B3[i]
                worst = max(worst, abs(ex - fl)/max(1e-12, abs(B2[i]) + abs(B3[i])))
print(f'3. max relative |poly_u/k0 - flat_parts_nf| over (h, J, u, pattern): {worst:.2e}')
