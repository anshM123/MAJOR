import numpy as np, math, json, sys, os
from scipy.special import logsumexp
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from decimal import Decimal as D, localcontext, ROUND_CEILING
sys.set_int_max_str_digits(0)
COEFF=Path(os.environ.get('SLIVER_COEFF', str(Path(__file__).resolve().parent.parent/'sliver_merge'/'coeff_compact60.txt'))); S=10**60
f=open(COEFF);N,DD=map(int,next(f).split());assert DD==60
AI=[];GI=[];maxrad=0
for l in f:
 a,r,g,t=map(int,l.split());AI.append(abs(a)+r);GI.append(abs(g)+t);maxrad=max(maxrad,r,t)
assert len(AI)==N+1
TI=[((k+1)**2*AI[k+1]+(k*k+8*k+65)*AI[k]+2*GI[k]+1)//2 for k in range(N)]
A=np.array([math.log(a+1)-60*math.log(10) for a in AI[:-1]]);T=np.array([math.log(t+1)-60*math.log(10) for t in TI]);k=np.arange(N);p=k*(k+9)
# For b=0 or 4, |k2|<=k^2+8k+65 includes q.
Rhos=[F(1)+F(j,20) for j in range(1,41)]+[F(j,4) for j in range(13,81)]+[F(j) for j in range(21,81)]

def cap(c):
 S=10**22;sl=F(isqrt(c.numerator*S*S//c.denominator),S);v=(1-sl)**2/998**2
 return F(math.ceil(v*10**18),10**18)

def shortK(l,u,H):
 z=np.logaddexp(A,T+math.log(float(H)))+k*math.log(float(u));tail=np.logaddexp.accumulate(z[::-1])[::-1]
 idx=np.where(tail<math.log(1e-12))[0]
 assert len(idx)>0
 return int(idx[0])

def float_bounds(l,u,Hl,Hu,K,axis):
 out=[];ps=p[:K];# p_j(s=9) not automatically worst for ratio under complex, include both p_j s=1 and 9 endpoints
 for rho in Rhos:
  R=(rho+1/rho)/2
  if axis=='c':
   c0=(l+u)/2;dc=(u-l)/2;H0=(Hl+Hu)/2;dH=(Hu-Hl)/2
   pairs=[(c0+z*dc,F(0)) for z in [-R,R]]+[(c0+z*dc,H0+z*dH) for z in [-R,R]]
  else:
   pairs=[(c,H*(1+z)/2) for c,H in [(l,Hl),(u,Hu)] for z in [-R,R]]
  HH=max(abs(h) for c,h in pairs);norm=max(abs(1-c) for c,h in pairs)**3
  if max(h for c,h in pairs)*(K-1)*(K+8)>=1:continue
  f=np.zeros(K)
  for pp in [ps,k[:K]*(k[:K]+1)]:
   for c,h in pairs:
    f=np.maximum(f,np.abs(float(c)-pp*float(h))/(1-pp*float(h)))
  M=np.r_[0.,np.cumsum(np.log(np.maximum(f,1e-300)))]+np.logaddexp(A[:K+1],T[:K+1]+math.log(float(HH)))+math.log(float(norm))+math.log(2/float(rho-1))
  out.append((rho,M))
 assert out
 rho=np.array([math.log(float(r)) for r,v in out]);M=np.array([v for r,v in out]);
 def err(m):return logsumexp(np.min(M-m*rho[:,None],axis=0))
 lo,hi=1,500
 while err(hi)>math.log(1e-10):hi*=2
 while hi-lo>1:
  m=(hi+lo)//2
  if err(m)>math.log(1e-10):lo=m
  else:hi=m
 # choose optimum radius for each k, only compute selected exact.
 m=hi+2;choice=np.argmin(M-m*rho[:,None],axis=0)
 # reduce to unique radii
 used=sorted(set(choice));rs=[out[i][0] for i in used];ii={v:i for i,v in enumerate(used)}
 return m,rs,[ii[i] for i in choice],err(m)

def exact_check(l,u,Hl,Hu,K,axis,m,rs,choices):
 with localcontext() as ctx:
  ctx.prec=45;ctx.rounding=ROUND_CEILING
  def dec(f):f=F(f);return D(f.numerator)/D(f.denominator)
  total=D(0);bounds=[None]*(K+1)
  for idx,rho in enumerate(rs):
   R=(rho+1/rho)/2
   if axis=='c':
    c0=(l+u)/2;dc=(u-l)/2;H0=(Hl+Hu)/2;dH=(Hu-Hl)/2
    pairs=[(c0+z*dc,F(0)) for z in [-R,R]]+[(c0+z*dc,H0+z*dH) for z in [-R,R]]
   else:pairs=[(c,H*(1+z)/2) for c,H in [(l,Hl),(u,Hu)] for z in [-R,R]]
   H=max(abs(h) for c,h in pairs);norm=max(abs(1-c) for c,h in pairs)**3
   assert max(h for c,h in pairs)*(K-1)*(K+8)<1
   coef=dec(F(2)*rho**(-m)/(rho-1)*norm);prod=D(1);hd=dec(H)
   for j in range(K+1):
    if choices[j]==idx:bounds[j]=prod*(D(AI[j])/D(S)+hd*D(TI[j])/D(S))*coef
    if j<K:
     mx=max(abs(c-h*j*(j+s))/(1-h*j*(j+s)) for c,h in pairs for s in (1,9));prod*=dec(mx)
  assert all(x is not None for x in bounds)
  total=sum(bounds,D(0));assert total<D('1e-10'),(axis,total)
  return str(total)

def check_tail(l,u,H,K):
 with localcontext() as ctx:
  ctx.prec=50;ctx.rounding=ROUND_CEILING
  def dec(x):return D(x.numerator)/D(x.denominator)
  c=dec(u);h=dec(H);v=D(0)
  for j in range(N-1,K,-1):v=v*c+D(AI[j])/D(S)+h*D(TI[j])/D(S)
  power=D(1);base=c;ex=K+1
  while ex:
   if ex&1:power*=base
   base*=base;ex//=2
  v*=power;assert v<D('1e-11'),v
  # beyond N-1 retained ordinary tail at N-1 <= 2^-140 over c<=.9.
  return str(v)

def run(l,u):
 Hl,Hu=cap(l),cap(u);K=shortK(l,u,Hl)
 if Hl*(K-1)*(K+8)>=1: return None
 tails=check_tail(l,u,Hl,K)
 res=dict(l=str(l),u=str(u),Hl=str(Hl),Hu=str(Hu),K=K,physical_tail=tails,maxrad=maxrad)
 for ax in ['c','y']:
  m,rs,choice,err=float_bounds(l,u,Hl,Hu,K,ax);print(ax,m,len(rs),err,flush=True)
  e=exact_check(l,u,Hl,Hu,K,ax,m,rs,choice);res[ax]=dict(deg=m,rhos=[str(r) for r in rs],choices=choice,error=e)
 return res
if __name__=='__main__':
 l,u=map(F,sys.argv[1:3]);print('start',l,u,flush=True);r=run(l,u);assert r;print({k:v if k not in ['c','y'] else {k1:v1 for k1,v1 in v.items() if k1!='choices'} for k,v in r.items()},flush=True);json.dump(r,open(sys.argv[3],'w'))
