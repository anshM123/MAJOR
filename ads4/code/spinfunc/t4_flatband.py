"""
t4_flatband (RESULT_D_CONDITIONAL.md §4): the flat D = 4 IR band of NF02 with a sharp cutoff p > eps = m_IR/M (CHLPS 2201.06602 (3.6)
analogue), scaling limit F_eps(b) = 2 int_eps^1 psi2(p) J0(b p) dp (both helicities; psi3 is u-suppressed), against the AdS band of
LARGEJ_OFFDIAG with m_IR := 3/(2 R_AdS), i.e. eps = 1.5/Lam.  Prints the negative b-hat interval(s) and b_c/eps^(-2/3).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
import r3_redesign as R3

v = R3.nf_vector()
x, w = np.polynomial.legendre.leggauss(64)
for Dg in (100, 300, 1000):
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0)); eps = 1.5/Lm
    ed = np.unique(np.concatenate([np.geomspace(eps, 1, 40), np.linspace(eps, 1, 400)]))
    p = np.concatenate([(b - a)*(x + 1)/2 + a for a, b in zip(ed[:-1], ed[1:])])
    pw = np.concatenate([(b - a)*w/2 for a, b in zip(ed[:-1], ed[1:])])
    psi2 = v[:R3.D + 1] @ R3.basis_fn(p)
    bh = np.geomspace(5, 20*Lm, 4000)
    F = np.array([2*np.sum(pw*psi2*j0(b*p)) for b in bh])
    neg = F < 0
    edges = np.where(np.diff(neg.astype(int)) != 0)[0]
    ivs = []
    st = None
    for i in range(len(bh)):
        if neg[i] and st is None: st = bh[i]
        if (not neg[i]) and st is not None: ivs.append((st, bh[i - 1])); st = None
    if st is not None: ivs.append((st, bh[-1]))
    print(f"Delta_g={Dg} Lam={Lm:.1f} eps={eps:.4g}: flat-cutoff NF02 negative b-hat intervals " +
          ", ".join(f"[{a:.1f}, {b:.1f}]" for a, b in ivs[:4]) + (f" (+{len(ivs) - 4} more)" if len(ivs) > 4 else "") +
          f"; b_c = {ivs[0][0]:.1f} = {ivs[0][0]/Lm**(2/3):.2f} Lam^(2/3) = {ivs[0][0]*eps**(2/3):.2f} eps^(-2/3); first end {ivs[0][1]/Lm:.2f} Lam")
    for b in (20, 30, 50, 100, 200):
        i = np.argmin(abs(bh - b)); print(f"     b-hat {b}: b^3 F = {bh[i]**3*F[i]:+.1f}")
