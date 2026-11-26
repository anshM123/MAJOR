# deviation of the prefix from its linear-in-h interpolation, t=h/hphys
import sys
from common import *
c=float(sys.argv[1]);b=int(sys.argv[2]);K=need_k(c);hp=hphys(c,b)
v0=S0(arb(c),arb(0),b,K);v1=S0(arb(c),arb(hp),b,K)
for t in [0.1,0.3,0.5,0.7,0.9,1.0,1.02,1.04,1.06,1.08,1.1,1.12,1.15,1.2]:
    v=S0(arb(c),arb(hp*t),b,K);lin=v0+(v1-v0)*t
    print('t %.2f S %.6e dev %.3e'%(t,float(v.mid()),float((v-lin).mid())),flush=True)
