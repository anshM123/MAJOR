# predicted L1/centre of route-A rectangles as a function of the twist overshoot omega (box sampling)
import sys,math,subprocess
from common2 import hcap_tw,needK
b=int(sys.argv[1])
for c in [float(x) for x in sys.argv[2].split(',')]:
  for om in [float(x) for x in sys.argv[3].split(',')]:
    lo,hi=0.0,0.2
    for _ in range(50):
        m=(lo+hi)/2
        if hcap_tw(c)/hcap_tw(c+m)<=om:lo=m
        else:hi=m
    cu=round(c+lo,5);K=needK(cu,target_log2=math.log2(1e-7))
    args=['python3','box_cheb_h.py',str(c),str(cu),str(b),str(K),'16','16']
    if c>=0.55:
        for fl,fu in ((0,.5),(.5,1)):print(subprocess.run(args+[str(fl),str(fu)],capture_output=True,text=True).stdout.strip(),flush=True)
    else:print(subprocess.run(args,capture_output=True,text=True).stdout.strip(),flush=True)
