# Route B (termination-preserving x-box, n>=K) majorant/cost on c-cells; x in [0, min(1/(K(K+s)), xphys(cu))]
import sys,math
from cost_model import best
from common2 import needK,xphys
from mpmath import mpf
b=int(sys.argv[3]) if len(sys.argv)>3 else 0;s=2*b+1
for c in [float(x) for x in sys.argv[1].split(',')]:
  for w in [float(x) for x in sys.argv[2].split(',')]:
    cl,cu=c,c+w;K=needK(cu);X=min(1/(K*(K+s)),xphys(cu))
    r=best('B',[mpf(cl),mpf(cu),mpf(0),mpf(X)],b,K)
    print('c %.3f w %.4f K %d X/xK %.3f :'%(c,w,K,X*K*(K+s)),r and (r['ac'],r['ah'],r['Q'],r['log10v'],r['terms'],r['est_s'],'h/unit-c %.2f'%(r['est_s']/w/3600)),flush=True)
