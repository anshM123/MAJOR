"""lx_phase3.py: the general lemma (tau up to 1150, phase III active): rt_pp.cross vs the independent re-derivation
lx_common.bound_indep, which repairs the phase-III nits (|lam+| <= y + u + sqrt(y^2-1+u^2+2uy) in the final factor,
|r| >= sqrt(y^2-1-u^2) in the adiabatic integral, first-window-step factor, step bound dy = 2 tau/(a+1)^2).
cross >= indep at every point means the nits are absorbed by the slack already in rt_pp.cross."""
import math, sys
sys.path.insert(0, '/home/user/MAJOR/ads4/code/analytic')
import rt_pp
from lx_common import bound_indep, LAMf
worst = 1e9
for tau in (LAMf, 1000.0, 1031.0, 1080.0, 1150.0):
    for sg in (0.0, 0.25, 2.5, 5.5):
        for c in (0.5, 0.7, 1.0, 1.5, 2.0):
            try:
                x = math.log(rt_pp.cross(sg, tau - 1 if tau > 1000 else tau, tau, c)); y = bound_indep(sg, tau - 1 if tau > 1000 else tau, tau, c)
            except AssertionError:
                continue
            worst = min(worst, x - y)
            print(f'tau={tau:.2f} sigma={sg} c={c}: log cross = {x:.3f}, log indep = {y:.3f}, cross - indep = {x - y:+.3f}')
print('min (cross - indep) =', round(worst, 4))
