"""
p2spin/s5b_dphi.py (P2_SPINNING.md §7.1): Delta_phi-(in)dependence of the top-nu coefficients of CHMRS's O(1/m2^2) Regge-moment
series (ancillary regge_moments_order8.m, path argv[1]); k = 2, alpha = 2 Delta_phi + k - 1, d = 3; Legendre-reduced to (P, P').
-> cache/s5b_dphi.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import s5_order4 as S5, s3_mech as M
out = []
eta = sp.Symbol('eta')
for Dphi in (1, 2, 3, 5, 8):
    t2, nu, et = S5.chmrs_order2(sys.argv[1], Dphi=Dphi, alpha=2*Dphi + 1)
    sol = {0: (1, 0), 1: (0, 1)}
    for k in range(0, 8):
        c0 = -(2*(k + 1)*et*sol[k + 1][0] + (k*(k + 1) + M.nu2 + sp.Rational(1, 4))*sol[k][0])/(et**2 - 1)
        c1 = -(2*(k + 1)*et*sol[k + 1][1] + (k*(k + 1) + M.nu2 + sp.Rational(1, 4))*sol[k][1])/(et**2 - 1)
        sol[k + 2] = (c0, c1)
    t2 = sp.expand(t2.subs(nu, sp.sqrt(M.nu2)))
    r0 = sp.expand(sp.cancel(sum(t2.coeff(sp.Symbol(f'G{k}'))*sol[k][0] for k in range(9))))
    r1 = sp.expand(sp.cancel(sum(t2.coeff(sp.Symbol(f'G{k}'))*sol[k][1] for k in range(9))))
    out.append(f'Dphi={Dphi}: nu^6 P: {sp.factor(r0.coeff(M.nu2, 3))};  nu^4 P: {sp.factor(r0.coeff(M.nu2, 2))};  nu^4 P\': {sp.factor(r1.coeff(M.nu2, 2))}')
    print(out[-1], flush=True)
open(os.path.join(HERE, 'cache', 's5b_dphi.txt'), 'w').write('\n'.join(out) + '\n')
