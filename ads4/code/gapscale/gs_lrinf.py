"""
gapscale/gs_lrinf.py (GAP_SCALING.md §3.1): re-run the CERTIFIED Lemma LR-inf remainder ratio r(beta) of lrtight/lrt3_rcert.py
at other gaps, by patching Lam (and Delta_g) in lrtight/lrt_common BEFORE importing lrt3_rcert.  Nothing in lrtight/ is written:
CACHE is redirected to gapscale/cache/lr_<Dg>/, into which the Lam-independent inputs (lrt1_constants.json: int|psi|, A1min, A2p;
lrt2_m1.json: certified M1 lower bounds -- M1(beta) has no Lam in it) are copied read-only.
Usage: python3 gs_lrinf.py Dg   -> cache/lr_<Dg>/lrt3_rcert.{json,txt}
Grade: [certified] in the same sense as lrt3 (Arb, 128 bits) for the patched Lam = sqrt((Dg-2)(Dg-1)).
"""
import os, sys, shutil
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
LRT = os.path.join(HERE, '..', 'lrtight')
sys.path.insert(0, LRT)
Dg = int(sys.argv[1])
out = os.path.join(HERE, 'cache', 'lr_%d' % Dg)
os.makedirs(out, exist_ok=True)
for f in ('lrt1_constants.json', 'lrt2_m1.json'):
    shutil.copy(os.path.join(LRT, 'cache', f), os.path.join(out, f))
import lrt_common as LC
LC.TAU = Dg
LC.LAM2 = Fr((Dg - 2)*(Dg - 1))
LC.LAM = LC.A(LC.LAM2).sqrt()
LC.CACHE = out
import lrt3_rcert as L3
assert L3.LAM is LC.LAM and L3.CACHE == out
L3.pr('# gapscale: lrt3_rcert re-run at Dg = %d, Lam = %s' % (Dg, LC.LAM.str(10)))
L3.main()
