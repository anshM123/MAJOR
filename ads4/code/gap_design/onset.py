# deviation of the prefix from its linear-in-y interpolant beyond the physical cap y=h/hcap(c) in [1,1.2]
import sys
from common2 import *
b=int(sys.argv[2])
for c in [float(x) for x in sys.argv[1].split(',')]:
    K=needK(c,target_log2=math.log2(1e-7));hp=hcap_tw(c)
    v0=Sh(c,0,b,K);v1=Sh(c,hp,b,K);row=[]
    for t in (1.02,1.04,1.06,1.08,1.10,1.12):
        v=Sh(c,hp*t,b,K);row.append(float((v-v0-(v1-v0)*t).mid())/float(v1.mid()))
    print('c %.2f b %d K %d rel.dev at y=1.02..1.12:'%(c,b,K),' '.join('%.1e'%x for x in row),flush=True)
