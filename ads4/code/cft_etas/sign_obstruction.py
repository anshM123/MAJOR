"""Test of one-signedness of the relaxation-function deviation phi(t) via log-convexity of its Laplace transform.
If phi >= 0 (phi != 0) then s -> log F(s) is convex on s >= 0 (Hoelder); if phi <= 0 then log(-F) is convex.
Horizon-layer deformations a=b=1+Amp*chi_w (monotone coupling, positive depth measure; horizon_bump.py), units 2piT=1.5.
A negative second difference of log F with F > 0 everywhere proves phi is not one-signed."""
import numpy as np
import horizon_bump as hb
for w in (0.8, 0.4, 0.1):
    ss = np.arange(0.0, 4.01, 0.5); ss[0] = 1e-4
    F = np.array([hb.E(s, 0.1, w) - hb.E(s, 0.0, w) for s in ss])
    lf = np.log(F)
    d2 = lf[2:] - 2*lf[1:-1] + lf[:-2]
    print('w=%.2f  F(s), s=0..4 step 0.5:' % w, ' '.join('%.4e' % x for x in F))
    print('        second differences of log F:', ' '.join('%+.4f' % x for x in d2))
    print('        all F>0: %s ; log-concave somewhere: %s  => phi NOT one-signed' % (bool(np.all(F > 0)), bool(np.any(d2 < -1e-6))))
