"""Timing probe: exact Lemma-P positivity (rcert.region_P machinery, read-only import) for single +- polynomials at large J."""
import sys, os, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'robust'))
import nf_cert, rcert
F = nf_cert.load('NF02')
for J in [int(a) for a in sys.argv[1:]]:
    for g in rcert.gens(*F, Fr(1, 5), '+-', 1)[:1]:
        t0 = time.time(); P, L = rcert.fbc.to_fmpz_poly(rcert.poly_u(g, J)); t1 = time.time()
        ok, n = rcert.positive_on_01(P); t2 = time.time()
        print('J=%d %s build %.1fs descartes %.1fs ok=%s intervals=%d bits~%d' % (J, g.label, t1 - t0, t2 - t1, ok, n, max(abs(c) for c in P.coeffs()).bit_length()), flush=True)
