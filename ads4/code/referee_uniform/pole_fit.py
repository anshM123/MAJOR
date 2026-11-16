"""Referee: is the pole model |rho| <= 20 max(1,(1-x^2)^-4)/m^2 an envelope of the controlled skeleton values?"""
import sys, os, pickle, numpy as np
U = '/home/user/MAJOR/ads4/code/uniform'; sys.path.insert(0, U); sys.path.insert(0, U + '/../robust')
from rcommon import gl01
pn, pw = gl01(48)
for Dg in (20, 50, 100, 1000):
    C = pickle.load(open(f'{U}/cache/skel_t3_{Dg}.pkl', 'rb'))
    worst = {}
    for (mfl, J), sk in C.items():
        if sk is None or J > int(sys.argv[1]): continue
        x = pn/mfl; m = sk['m']
        for k in (2, 3):
            rho, err, ctrl = sk[k]
            env = 20*np.maximum(1, (1 - np.minimum(x*x, .999))**-4)/m**2
            r = np.where(ctrl, (np.abs(rho)+err)/env, 0)
            r0 = np.where(ctrl, np.abs(rho)/env, 0)
            i = int(np.argmax(r))
            band = 'near' if mfl < 1.5 else 'far'
            kk = (k, band)
            if r[i] > worst.get(kk, (0,))[0]:
                worst[kk] = (round(r[i],3), round(r0[i],3), mfl, J, round(x[i],3))
    print(Dg, {k: v for k, v in sorted(worst.items())})
