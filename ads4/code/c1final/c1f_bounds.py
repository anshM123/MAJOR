#!/usr/bin/env python3
"""C1_FINAL step 4: the headline bounds with c1 = 0, from the exact formulas of C1_COMPUTE §4 (imports ../c1/c1_quadratic.py and,
through it, ../dgscan/bounds.py read-only; writes only to ./cache).  Also prints what a hypothetical |c1| would have done (for the
record of the removed error terms).
Output: cache/c1f_bounds.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'c1'))
import mpmath as mp
import c1_quadratic as CQ          # read-only: D_bounds, Q_bounds, DK

out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

P("# Headline bounds at c1 = 0 (C1_FINAL) -- exact forms R0 (x + sqrt(1+x^2)) and 2 L_AdS + coef c1^2 evaluated at c1 = 0")
for Dg in (300, 1000):
    d0 = CQ.D_bounds(Dg, 0); q0 = CQ.Q_bounds(Dg, 0)
    dK = CQ.D_bounds(Dg, 0, K=1 + CQ.DK[Dg])
    d1 = CQ.D_bounds(Dg, 1); q1 = CQ.Q_bounds(Dg, 1)
    P(f"Delta_g = {Dg}: Lam = {mp.nstr(d0['Lam'], 9)}")
    P(f"   Result D: sqrt(t4^2+aT^2) <= R0 = {mp.nstr(d0['R0'], 6)}  (K = 1);  with K = k(0)(1+dk): {mp.nstr(dK['R0'], 6)};"
      f"  |g3|^2 Lam^8 <= {mp.nstr(d0['g3'], 7)};  t4 interval symmetric [-R0, R0]")
    P(f"   Result Q: upper end g4~ Lam^6/8piG <= 2 L_AdS = {mp.nstr(q0['L2'], 7)}  (c1^2 term = 0)")
    P(f"   removed caveat (per unit |c1|, for the record): D x = {mp.nstr(d1['x'], 6)} (ratio {mp.nstr(d1['R']/d0['R0'], 7)}),"
      f" Q coefficient {mp.nstr(q1['eps_free'], 6)}")
open(os.path.join(HERE, 'cache', 'c1f_bounds.txt'), 'w').write("\n".join(out) + "\n")
