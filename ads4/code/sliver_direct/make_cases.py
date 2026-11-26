from fractions import Fraction as F
from base_bounds import hcap,hi
from check_lowfreq import exp_u
from check_norm import plan,S
from major_rect import major
from contour_direct import ds
from math import ceil,prod
import json,sys
L,H=map(F,sys.argv[1:3]);dl=F(int((1/exp_u(H))*10**7),10**7);du=F(ceil(exp_u(-L)*10**7),10**7)
xc=F(ceil(hcap(du,F(998))/du**2*10**15/2),10**15)
cs=[]
for l,u in [('159','400.1'),('400','700.1'),('700','850.1'),('850','950.1'),('950','1040.1')]:
 l,u=F(l),F(u)
 c=dict(name=f'real_3_{float(l)}_{float(u)}',d0=ds((dl+du)/2),dd=ds((du-dl)/2),h0=ds(xc),dh=ds(xc),nu0=ds((l+u)/2),dnu=ds((u-l)/2),nuImag='0',fullTwist=True,coupled=True,K=60,digits=100,epsilon='.08',channels=[[3,0],[3,4]],rhos=['2','4','8','16','32','64','128','256'])
 while F(major(c,3,4)[1],S)>F('.003'):c['K']+=10
 cs.append(c)
cs.append(dict(cs[-1],name='real_4_950_1040.1',channels=[[4,0],[4,4]]))
for c in cs:
 p,b=c['channels'][-1];axes,t,*_=plan(c,p,b);print(c['name'],c['K'],prod(x['N']+1 for x in axes),[x['m'] for x in axes],float(F(t,S)),flush=True)
json.dump(cs,open(f'boxes_{float(L)}_{float(H)}.json','w'),indent=2)
