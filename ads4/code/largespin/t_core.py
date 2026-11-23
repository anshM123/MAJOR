import sys, time; sys.dont_write_bytecode=True
import numpy as np, ls_core as C
import mpmath as mp
# check OmLR against direct mpmath at a few t (k=2) using P from hyp2f1 on Gauss nodes
L=C.lam(1000); F=C.NFw()
for t in [0.01,0.1,0.4,1.0]:
    p,w=C.gl01(600)
    z=-np.sinh(t/2)**2
    P=np.array([float(mp.hyp2f1(0.5+1j*L*pp,0.5-1j*L*pp,1,z).real) for pp in p])
    d=(F.psi2(p)*w*P).sum()
    t0=time.time(); o=C.OmLR(2,[t])[0]; t1=time.time()
    print(t, d, o, o/d-1, 'time',t1-t0)
