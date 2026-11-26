# How does the normalised prefix action vary over the physical twist range h in [0,hphys(c)]?
import sys,math,time
from common import *
for c in [float(x) for x in sys.argv[1].split(',')]:
  for b in (0,4):
    K=need_k(c);hp=hphys(c,b);t0=time.time()
    vals=[S0(arb(c),arb(hp*t),b,K) for t in (0,.25,.5,.75,1.0,1.25)]
    print('c %.3f b %d K %d hphys %.3e hK^2 %.3f |'%(c,b,K,hp,hp*K*K),' '.join('%.4e'%float(v.mid()) for v in vals),' (%.1fs)'%(time.time()-t0),flush=True)
