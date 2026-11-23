"""NF02 flat certificate, region P extended (CERT_CLOSE.md §1): exact Lemma-P polynomials + Descartes (rcert.poly_u /
rcert.positive_on_01, read-only import from ads4/code/robust), for the +- channel, J0 <= J < J1, all u in [0,1], margin DELTA.
Usage: python3 nf_Pext.py DELTA J0 J1   (e.g. 1/5 400 800)"""
import sys, os, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'robust'))
import nf_cert, rcert           # nf_cert.load prints the sha256 fingerprint of the exact NF02 coefficients
delta, J0, J1 = Fr(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
F = nf_cert.load('NF02')
ok_all = True; cnt = 0; t0 = time.time()
for J in range(J0, J1):
    sF = 1 if (2 * J * (J + 1) - 41 - 1) >= 0 else -1
    for g in rcert.gens(*F, delta, '+-', sF):
        P, L = rcert.fbc.to_fmpz_poly(rcert.poly_u(g, J))
        ok, n = rcert.positive_on_01(P); cnt += 1; ok_all &= ok
        if not ok or J % 50 == 0:
            print('Pext %s J=%d deg=%d intervals=%d %s (%.0fs)' % (g.label, J, P.degree(), n, 'PASS' if ok else 'FAIL', time.time() - t0), flush=True)
print('VERDICT nf_Pext (delta=%s, +-, %d <= J < %d, %d polynomials, %.0fs): %s' % (delta, J0, J1, cnt, time.time() - t0, 'PASS' if ok_all else 'FAIL'))
