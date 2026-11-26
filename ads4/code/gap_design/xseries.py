# Taylor coefficients in y = x/x_phys(c) (x=h/c) of the UNNORMALISED prefix numerator S_0, at fixed c.
# Non-rigorous diagnostic (arb balls, midpoint table). Prints log10|s_m| for selected m.
import sys,time,math
from common import *
from flint import arb_series
c=float(sys.argv[1]);b=int(sys.argv[2]);M=int(sys.argv[3]);K=int(sys.argv[4]) if len(sys.argv)>4 else need_k(c)
ctx.cap=M+1
A,Bs=parts(b,K);s=2*b+1;q=b*(b+1)
cc=arb(c);hp=hphys(c,b);xp=arb(hp)/cc   # x_phys
t0=time.time()
y=arb_series([0,1])
X=xp*y           # x as series in y
h=cc*X
S=None
for k in range(K,-1,-1):
    Fk=A[k]+h*(q*A[k]+Bs[k])
    if S is None:S=Fk;continue
    p=k*(k+s)
    # r_k = c + (c-1) sum_{m>=1} (c p x_phys)^m y^m
    g=cc*p*xp
    co=[cc]+[(cc-1)*g**m for m in range(1,M+1)]
    S=Fk+arb_series(co)*S
co=S.coeffs()
print('c',c,'b',b,'K',K,'M',M,'time %.1f'%(time.time()-t0))
tot=arb(0)
for m,v in enumerate(co):
    if m<6 or m%25==0 or m==len(co)-1:
        r=v.rad();mid=v.mid()
        print(m,'%.3e'%float(mid),'rad %.1e'%float(r))
import json
json.dump([str(v.mid()) for v in co],open(f'xser_{c}_{b}_{M}.json','w'))
