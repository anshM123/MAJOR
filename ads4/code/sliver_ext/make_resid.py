from fractions import Fraction as F
from base_bounds import hcap
from check_lowfreq import exp_u
from check_norm import plan,S
from contour_direct import ds
import json,math,sys
Llo,Lhi=map(F,sys.argv[1:3]);dl=F(int((1/exp_u(Lhi))*10**6),10**6);du=F(math.ceil(exp_u(-Llo)*10**6),10**6);hc=F(math.ceil(hcap(du,F(998))*10**15/2),10**15)
assert dl < 1/exp_u(Lhi) and du>=exp_u(-Llo)
cs=[]
for l,u in [('159','400.1'),('400','700.1'),('700','850.1'),('850','950.1'),('950','1040.1')]:
 lo,up=F(l),F(u);c=dict(name=f'real_3_{l}_{u}',d0=ds((dl+du)/2),dd=ds((du-dl)/2),h0=ds(hc),dh=ds(hc),nu0=ds((lo+up)/2),dnu=ds((up-lo)/2),nuImag='0',fullTwist=True,K=80,digits=90,epsilon='.05',channels=[[3,0],[3,4]])
 from check_norm import major
 while True:
  M,t,g=major(c,3,4,[F(1)]*3)
  if F(t,S)<F('.001'):break
  c['K']+=10
  if c['K']>1500:raise ValueError('tail')
 axes,t,*_=plan(c,3,4);n=math.prod(a['N']+1 for a in axes);print(c['name'],c['K'],n,[a['m']for a in axes],float(F(t,S)),flush=True);cs.append(c)
json.dump(cs,open(f'resid_{float(Llo)}_{float(Lhi)}.json','w'),indent=2)
c=dict(cs[-1]);c.update(name='real_4_950_1040.1',channels=[[4,0],[4,4]])
json.dump([c],open(f'p4_{float(Llo)}_{float(Lhi)}.json','w'),indent=2)
