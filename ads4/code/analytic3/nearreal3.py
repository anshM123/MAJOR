#!/usr/bin/env python3
"""nearreal3.py La Lb [k0]: bottom and left contour pieces (psi2 dressed, psi3) with the exact state constants of norm_extra.patch_state
(valid for every L; needed below L = 3.05).  Output cache/nearreal3_La_Lb_k{k0}.json in contour3's format (as ratios x F inside)."""
import sys, os, math, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kernel3                      # sets the analytic/analytic2 paths
import norm_extra as NE
import rt_cert as RC
import rt_pp as PPm

if __name__ == '__main__':
    La, Lb = float(sys.argv[1]), float(sys.argv[2]); k0 = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    cst = NE.patch_state(La)
    d = math.exp(-La); t0 = time.time()
    A, B, den, Fi = PPm.rec_step(1500 + 4 + 8 + 1, 9, d)
    assert (A + B)/den <= 2.0 and Fi/den <= 1.0, (A, B, den, Fi)       # REC hypotheses of the near-real geometric tails
    F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
    res = dict(F=F, k0=k0, state_const=cst, rec=[(A + B)/den, Fi/den])
    for (name, p, Q, Qp, dress) in (('psi2', 3, RC.Q2, RC.Q2p, True), ('psi3', 4, RC.Q3, RC.Q3p, False)):
        fac = d**3*(2/RC.LAMf)*2 if name == 'psi2' else d**3*(d*d/RC.TAU_COEF)*RC.LAMup*2
        Y = (p - 0.5) - 0.25 + k0
        res[name] = dict(Y=Y, bottom=fac*RC.piece_bottom(La, p, Q, Qp, dress), left=fac*RC.piece_left(La, p, Q, dress, Y))
        print(name, {k: (v/F if k != 'Y' else v) for k, v in res[name].items()}, f'[{time.time()-t0:.0f}s]', flush=True)
    json.dump(res, open(HERE / 'cache' / f'nearreal3_{La}_{Lb}_k{k0}.json', 'w'))
