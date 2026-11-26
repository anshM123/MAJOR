# compare h-piece partitions of [0,H] for one route-A cell (cost model)
import sys,math
from cost_model import best
from common2 import needK,hcap_tw
from mpmath import mpf
c=float(sys.argv[1]);w=float(sys.argv[2]);b=int(sys.argv[3]) if len(sys.argv)>3 else 0
cu=c+w;K=needK(cu,target_log2=math.log2(1e-7));H=hcap_tw(c)
for parts in [[0,.5,1],[0,1/3,2/3,1],[0,.25,.5,.75,1],[0,.6,.85,1],[0,.55,.8,.93,1],[0,.5,.75,.9,1]]:
    tot=0;det=[]
    for a,bb in zip(parts[:-1],parts[1:]):
        r=best('A',[mpf(c),mpf(cu),mpf(H*a),mpf(H*bb)],b,K)
        if r is None:tot=float('inf');break
        tot+=r['est_s'];det.append((r['ac'],r['ah'],r['Q'],r['log10v']))
    print(parts,round(tot,1),det,flush=True)
