"""Regge-machinery check of the rescaled functionals (ROBUST_FUNCTIONAL.md section 5.3): smearing scale Lam_s = Lam/n, physical states m >= Lam.
'last negative m/Lam' <= 1 for an envelope means: positive on every physical grid point with that envelope.  [NE]"""
import sys
from check import load
from ads import regge_scan
CASES = [('M02', 2.2), ('M02', 4.0), ('M02', 6.8), ('M01', 3.2), ('T1', 1.75), ('T3', 2.25), ('T10', 4.0)]
for name, n in CASES:
    F = load(name)
    for Dg in [100, 1000]:
        print('== %s n=%.2f Dg=%d' % (name, n, Dg), flush=True)
        regge_scan(F, Dg, n)
