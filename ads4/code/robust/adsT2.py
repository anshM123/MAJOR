"""adsT for the envelope-targeted functionals only (ROBUST_FUNCTIONAL.md section 5.3)."""
from check import load
from ads import regge_scan, hybrid_scan
for name, n in [('T1', 1.75), ('T3', 2.25), ('T10', 4.0)]:
    F = load(name)
    for Dg in [100, 1000]:
        print('== %s n=%.2f Dg=%d' % (name, n, Dg), flush=True)
        regge_scan(F, Dg, n)
