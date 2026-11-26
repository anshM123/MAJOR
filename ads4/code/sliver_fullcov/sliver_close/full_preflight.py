from fractions import Fraction as F
from decimal import Decimal as D,localcontext,ROUND_CEILING
from math import isqrt,ceil,log
import json,sys,time
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from coupled_plan import cap,shortK,check_tail,float_bounds,exact_check,AI,TI,S,N,maxrad

def sqrtlo(x):return F(isqrt(x.numerator*10**80//x.denominator),10**40)
def sqrtup(x):return sqrtlo(x)+F(1,10**40)
# Lower bound for log(product_{ell>=K} (p_ell-c p_j)/(p_ell-p_j)).
def lower_exp(j,K,u,s,R=4):
 p=j*(j+s);A=F(p)+F(s*s,4);B=u*p+F(s*s,4);v=F(K)+F(s,2)
 return sum((A**r-B**r)/(r*(2*r-1)*v**(2*r-1)) for r in range(1,R+1))
# actual max absolute numerator coefficients; common across both channels
maxa=[0]*(N+1);maxt=[0]*(N+1)
for j in range(N-1,-1,-1):maxa[j]=max(AI[j],maxa[j+1]);maxt[j]=max(TI[j],maxt[j+1])
LOG2_UP=F('0.69314718055994530941723212145817656807550013436026')+F(1,10**49)
# Verify upper using exact atanh expansion with geometric remainder.
z=F(1,3);J=80;assert 2*sum(z**(2*j+1)/(2*j+1) for j in range(J))+2*z**(2*J+1)/((2*J+1)*(1-z*z))<LOG2_UP

def residue_groups(l,u,H,K,chunk=16):
 groups=[]
 for a in range(1,K,chunk):
  b=min(K-1,a+chunk-1);g=sqrtlo(1-H*b*(b+9));assert g>0
  r=(1-g)/(1+g);assert 0<=r<1
  # Use k=K; lower log integral is decreasing in k and increasing in j.
  L=min(lower_exp(a,K,u,s) for s in (1,9));e=int(L/LOG2_UP)
  M=F(maxa[a+1],S)+F(maxt[a+1],S*a*(a+1))
  # |rho_j| <= (2/j) c^-8 sum_{k>j} (...) e^-L; count at most K-a.
  B=F(2*(K-a)*(b-a+1),a)*l**-8*M*F(1,2**e)*2/(g*(1-r))
  groups.append((r,B))
 return groups

def residue_degree(l,u,H,K,degree=None):
 groups=residue_groups(l,u,H,K)
 # Exact directed decimal arithmetic at positive arguments; for power, explicit repeated squaring.
 with localcontext() as ctx:
  ctx.prec=45;ctx.rounding=ROUND_CEILING
  def dec(x):return D(x.numerator)/D(x.denominator)
  gd=[(dec(r),dec(B)) for r,B in groups]
  def powup(x,n):
   v=D(1)
   while n:
    if n&1:v*=x
    x*=x;n//=2
   return v
  def err(m):return sum((B*powup(r,m+1) for r,B in gd),D(0))
  if degree is not None:return degree,str(err(degree))
  lo,hi=0,1000
  while err(hi)>D('1e-10'):hi*=2
  while hi-lo>1:
   mid=(hi+lo)//2
   if err(mid)>D('1e-10'):lo=mid
   else:hi=mid
  m=hi+2;er=err(m);assert er<D('1e-10')
 return m,str(er)

def make_cell(l,u):
 Hl,Hu=cap(l),cap(u);K=shortK(l,u,Hl);assert Hl*(K-1)*(K+8)<1
 for s in (1,9):
  for k in (0,K-1):
   p=k*(k+s)
   for c,h in [(l,Hl),(u,Hu)]:assert abs(c-p*h)<=1-p*h and 1-p*h>0
 tails=check_tail(l,u,Hl,K);m,rs,choice,_=float_bounds(l,u,Hl,Hu,K,'c');ec=exact_check(l,u,Hl,Hu,K,'c',m,rs,choice)
 my,ey=residue_degree(l,u,Hl,K)
 res=dict(l=str(l),u=str(u),Hl=str(Hl),Hu=str(Hu),K=K,physical_tail=tails,maxrad=maxrad,
          c=dict(deg=m,rhos=[str(r) for r in rs],choices=choice,error=ec),y=dict(deg=my,error=ey,method='residue_r4_group16'))
 return res

def choose(l,lim=F('0.8')):
 lim=min(lim,F(99,100)*(1+l)/2)
 # largest cell endpoints on 10^-4 mesh up to .9, width <= .1
 aa=int(l*10000);bb=min(9000,aa+1000);H=cap(l)
 def ok(b):
  K=shortK(l,F(b,10000),H)
  return H*(K-1)*(K+8)<=lim
 if not ok(bb):
  while bb-aa>1:
   m=(aa+bb)//2
   if ok(m):aa=m
   else:bb=m
  bb=aa
 assert F(bb,10000)>l
 return F(bb,10000)
if __name__=='__main__':
 out=Path(__file__).resolve().parent;rows=[];l=F('.35');st=time.monotonic()
 while l<F('.9'):
  u=choose(l);z=make_cell(l,u);fn=f'cell_{int(l*10000):04d}_{int(u*10000):04d}.json';json.dump(z,open(out/fn,'w'));rows.append(dict(cell=fn,l=z['l'],u=z['u'],K=z['K'],mc=z['c']['deg'],my=z['y']['deg'],ec=z['c']['error'],ey=z['y']['error'],tail=z['physical_tail'],nodes=2*(ceil(z['c']['deg']*1.4)+1)*(ceil(z['y']['deg']*1.4)+1)))
  print(rows[-1],flush=True);l=u
 assert l==F('.9');print('preflight only',len(rows),'cells',sum(z['nodes'] for z in rows),'nominal nodes',time.monotonic()-st,flush=True);json.dump({'status':'PREFLIGHT_ONLY','rows':rows},open(out/'full_preflight.json','w'),indent=2)
