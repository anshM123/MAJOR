# Route A (G1 cap, author's triangle driver) cost per unit c, cells sized so that the h-rectangle
# overshoots the physical twist cap at c_hi by at most a factor OV; optional split of h into nh pieces.
import sys,math
from cost_model import best
from common2 import needK,hcap_tw
from mpmath import mpf
OV=float(sys.argv[2]);b=int(sys.argv[3]) if len(sys.argv)>3 else 0
for c in [float(x) for x in sys.argv[1].split(',')]:
    # width with hcap(cl)/hcap(cu)=OV
    lo,hi=0.0,0.2
    for _ in range(50):
        m=(lo+hi)/2
        if hcap_tw(c)/hcap_tw(c+m)<=OV:lo=m
        else:hi=m
    w=lo;cu=c+w;K=needK(cu,target_log2=math.log2(1e-7));H=hcap_tw(c)
    out=[]
    for nh in (1,2,3):
        tot=0;ok=True;det=[]
        for i in range(nh):
            r=best('A',[mpf(c),mpf(cu),mpf(H*i/nh),mpf(H*(i+1)/nh)],b,K)
            if r is None:ok=False;break
            tot+=r['est_s'];det.append((r['ac'],r['ah'],r['Q'],r['log10v']))
        if ok:out.append((nh,round(tot,1),round(tot/w/3600,2),det))
    print('c %.3f width %.4f K %d hK^2 %.3f :'%(c,w,K,H*K*K),out,flush=True)
