"""p2con/pc_p4b.py (P2_CONTRACT.md §2.2): K4 of the BEYOND-HYBRID O(1/m^4) remainder (scalar seed, Delta_phi = 3):
    R4(nu) = A4^CHMRS(nu, beta)  -  nu^4 g2(bhat p)       (g2 = the flat-scaling part of A4, pc_flat; flat kinematics as the hybrid)
costed exactly (de-shift form; A4 and nu^4 g2 entire) per sample state, in units of the K4 representative 4a = eta^2 nu^4 P, and
compared with the raw cost of the flat nu^6 term alone (P2_SPINNING §7.1).  Also the cost of the +3/16 J0 y^2 flat discrepancy (pc_flat).
-> cache/p4b.txt, cache/p4b.json.   Usage: python3 pc_p4b.py <regge_moments_order8.m>"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'lsw')); sys.path.insert(0, os.path.join(HERE, '..', 'p2spin'))
import numpy as np, sympy as sp
from scipy.special import j0, j1
import lsw_common as W
import s5_order4 as S5, s3_mech as M
import pc_flat as PF
out = []
pr = lambda s: (print(s, flush=True), out.append(s))
t2, nu_s, et = S5.chmrs_order2(sys.argv[1], Dphi=3, alpha=7)
t2 = sp.expand(t2.subs(nu_s, sp.sqrt(M.nu2)))
r0, r1 = PF.reduce(t2, et)
g2 = PF.flat_part(r0, r1, et, 2)
f0 = sp.lambdify((et, M.nu2), r0, 'numpy'); f1 = sp.lambdify((et, M.nu2), r1, 'numpy')
g2f = sp.lambdify((PF.bsym, PF.J0s, PF.J1s), g2, 'numpy')
R = json.load(open(os.path.join(HERE, '..', 'lsw', 'cache', 'costs_1000.json')))
nu = W.NU; rows = []
pr(__doc__.split('->')[0])
pr('  J   m/L    beta  | K4(R4: beyond-hybrid remainder) | K4(raw full A4) | K4(raw nu^6 term) | K4(3/16 J0 y^2) ')
for r in R:
    b = r['beta']; J = r['J']; mfl = r['mfl']
    P, dPb = W.P_and_dP(b); P = P.real; dPe = dPb.real/np.sinh(b)
    eta = np.cosh(b)
    A4 = f0(eta, nu**2)*P + f1(eta, nu**2)*dPe
    bb = (2*J + 1)/mfl*nu/W.L
    flat = nu**4*g2f(bb, j0(bb), j1(bb))
    c4a = abs(W.act(eta**2*nu**4*P)[0])
    row = dict(J=J, mfl=mfl, beta=b, K_rem=abs(W.act(A4 - flat)[0])/c4a, K_full=abs(W.act(A4)[0])/c4a,
               K_nu6=abs(W.act(-(eta**2 - 1)*nu**6*P/72)[0])/c4a, K_316=abs(W.act((3/16)*nu**4*j0(bb))[0])/c4a)
    rows.append(row)
    pr(f"  {J:3d} {mfl:.3f} {b:.4f} | {row['K_rem']:10.3g} | {row['K_full']:10.3g} | {row['K_nu6']:10.3g} | {row['K_316']:8.3g}")
pr('max over states: K4(remainder) = %.3g; for J >= 50: %.3g' % (max(x['K_rem'] for x in rows), max(x['K_rem'] for x in rows if x['J'] >= 50)))
json.dump(rows, open(os.path.join(HERE, 'cache', 'p4b.json'), 'w'))
open(os.path.join(HERE, 'cache', 'p4b.txt'), 'w').write('\n'.join(out) + '\n')
