import sys,math
from cost_model import best
from common2 import needK,hcap_tw
from mpmath import mpf
b=0
for c in [float(x) for x in sys.argv[1].split(',')]:
  for w in [float(x) for x in sys.argv[2].split(',')]:
    cl,cu=c,c+w;K=needK(cu)
    rA=best('A',[mpf(cl),mpf(cu),mpf(0),mpf(hcap_tw(cl))],b,K)
    ov=hcap_tw(cl)/hcap_tw(cu)
    rC=best('C',[mpf(math.sqrt(cl)),mpf(math.sqrt(cu)),mpf(0),mpf(1)],b,K)
    print('c %.3f w %.4f K %d | A(overshoot y=%.3f): %s | C: %s'%(c,w,K,ov,
      rA and (rA['ac'],rA['ah'],rA['Q'],rA['log10v'],rA['terms'],rA['est_s']),
      rC and (rC['ac'],rC['ah'],rC['Q'],rC['log10v'],rC['terms'],rC['est_s'])),flush=True)
