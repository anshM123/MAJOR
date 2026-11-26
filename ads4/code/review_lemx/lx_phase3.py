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

# dense scan of the boxes PP0 actually uses (c in {0.7, 1, 1.5}; box [t-1, t]); per box the certificate uses min over c,
# so compare min_c cross with min_c indep as well as each c separately
w_each = 1e9; w_min = 1e9; at = None
for t in range(949, 1151):
    for sg in (0.0, 0.25, 0.5, 1.0, 1.25, 2.5, 3.5, 4.5, 5.5, 6.5):
        xs = [math.log(rt_pp.cross(sg, max(900.0, t - 1.0), float(t), c)) for c in (0.7, 1.0, 1.5)]
        ys = [bound_indep(sg, max(900.0, t - 1.0), float(t), c) for c in (0.7, 1.0, 1.5)]
        e = min(x - y for x, y in zip(xs, ys))
        if e < w_each: w_each, at = e, (t, sg)
        w_min = min(w_min, min(xs) - min(ys))
print(f'dense scan tau_hi in [949, 1150], box [tau_hi - 1, tau_hi], c in (0.7, 1, 1.5): min_c-wise (cross - indep) = {w_each:.4f} at {at};'
      f' min (min_c cross - min_c indep) = {w_min:.4f}')
