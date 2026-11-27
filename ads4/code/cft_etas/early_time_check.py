"""Holographic check of the 'analytic data drop out' claim (Lemma R3): a deformation supported at depth u >= 1-w changes
NO large-frequency (OPE) coefficient: F(s) decays exponentially in s, at the rate set by the tortoise depth of the layer
(F(s) ~ exp(-2 r*(1-w) * 1.5 s), r*(u) = Int_0^u du'/f0(u') ), i.e. by bulk causality the early-time commutator is untouched for
t < 2 r*(1-w).  Units r+=1 (2piT = 1.5)."""
import numpy as np
from scipy.integrate import quad
import horizon_bump as hb
for w in (0.6, 0.4):
    rs = quad(lambda u: 1/(1 - u**3), 0, 1 - w)[0]
    print('w=%.2f  tortoise depth r*(1-w) = %.4f  -> predicted decay rate of F(s) in s: 2*1.5*r* = %.4f' % (w, rs, 3*rs))
    prev = None
    for s in (2, 4, 6, 8, 10, 12):
        F = hb.E(s, 0.1, w) - hb.E(s, 0.0, w)
        rate = (np.log(abs(prev)) - np.log(abs(F)))/2 if prev else float('nan')
        print('   s=%2d  F(s)=%+.4e   local decay rate -dlog|F|/ds = %.4f' % (s, F, rate))
        prev = F
