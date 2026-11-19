"""
r7_checks (LARGEJ_OFFDIAG.md §1.3, item (c)): bookkeeping checks of the 2x2 normalisation used by q3 / l9 for the off-diagonal.
 (i)  the reference diagonal ratio M~ref_++/M~ref_+- (order 0, C2, Delta_ref = 2e4) entering sq = sqrt|Re(.)| must be REAL POSITIVE
      (otherwise the per-helicity normalisation is not a real congruence and |.| would hide an indefinite leading action);
 (ii) the per-row helicity constants a_++, a_+- (common to C2 and C3H) must be positive;
 (iii) antisymmetric part of the state off-diagonal |M01 - M10| / |M01 + M10| (only the symmetric part enters lambda^T M lambda);
 (iv) basis conditioning: the real helicity basis T is defined only up to scale/sign per vector; check the LO action is diagonal in it
      (|M0_01| / sqrt|M0_00 M0_11| at the state, order 0).
Usage: python3 r7_checks.py -> cache/r7_checks.txt
"""
import os, sys, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import r2_ladder                                 # noqa
import q3_analysis as Q
A = Q.A
L8 = A.L


def main():
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    pr("# r7: normalisation checks (C2, K2 records)")
    for J in (8, 12, 16, 18, 20, 24, 28, 32, 36):
        RC = A.Recs('C2', J, 2)
        Rb = dict(L8.records(J)); Rb.update(A._load('C2', J, '_k2')); Rb.update(A._load('C2', J, ''))
        Bref = L8.basis(J, A.DREF, Rb)
        if Bref is None:
            pr(f"  J={J}: no reference basis"); continue
        Ds = sorted(set(k[0] for k in RC.R if k[0] < A.DREF and k[0] >= 300))
        g = {}
        try:
            g = Q.even_grid(J)
        except Exception:
            pass
        for D in Ds:
            m = A.mphys(D, J)
            B = L8.basis(J, D, Rb)
            if B is None:
                continue
            nus = sorted(k[1] for k in RC.R if k[0] == D)
            rat, asym, lo = [], [], []
            for nu in nus:
                st, rf = RC.get(D, nu), RC.get(A.DREF, nu)
                if st is None or rf is None:
                    continue
                Mr = A.helmat(rf, Bref, 2)[0]
                rat.append(Mr[0, 0]/Mr[1, 1])
                Ms = A.helmat(st, B, 2)
                S = sum(Ms[:3])
                asym.append(abs(S[0, 1] - S[1, 0])/max(abs(S[0, 1] + S[1, 0]), 1e-300))
                lo.append(abs(Ms[0][0, 1])/np.sqrt(abs(Ms[0][0, 0]*Ms[0][1, 1])))
            rat = np.array(rat)
            ah = g[D]['a'] if D in g else {}
            pr(f"  J={J:2d} Delta={D:6.0f} (m={m:.0f}): ref ++/+- ratio Re in [{rat.real.min():.3g}, {rat.real.max():.3g}], "
               f"max |Im/Re| {np.max(np.abs(rat.imag/rat.real)):.1e}; a_h = " + ", ".join(f"{h} {v:.3e}" for h, v in ah.items()) +
               f"; max asym {max(asym):.1e}; max LO off/diag {max(lo):.1e}")
    open(os.path.join(HERE, 'cache', 'r7_checks.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
