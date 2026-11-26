#!/usr/bin/env python3
"""A rescaled downward-closed Taylor triangle, with rigorous geometric Cauchy tail."""
import sys,subprocess,json,time,math
from math import lcm
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
def up(n,d):return -((-n)//d)
R=F(1999,2000)
def need_k(c):
 t=c/R;p,q=t.numerator,t.denominator
 def ok(k):return (p**(k+1)<<2700)<=q**(k+1)
 lo,hi=-1,1
 while not ok(hi):hi*=2
 while hi-lo>1:
  m=(lo+hi)//2
  if ok(m):hi=m
  else:lo=m
 return hi
l,u=map(F,sys.argv[1:3]);b=int(sys.argv[3]);D=100;S=10**D;K=need_k(u);s=2*b+1;q=b*(b+1)
assert 0<F('.05')<=l<u<=F('.56') and b in (0,4)
with open('direct_coefficients.txt') as f:
 N,DD=map(int,f.readline().split());assert N>=K+1 and DD>=D
 div=10**(DD-D);A=[];Ar=[];G=[];Gr=[]
 for k in range(K+2):
  a,r,g,t=map(int,f.readline().split());d=div<<k;A.append(a//d);Ar.append(up(r,d)+1);G.append(g//d);Gr.append(up(t,d)+1)
h=F(up(u.numerator*10**12,u.denominator*(1001-b)*(1002+b)),10**12)
c0=(l+u)/2;dc=(u-l)/2;h0=dh=h/2;den=lcm(c0.denominator,dc.denominator,h0.denominator,dh.denominator);C0,DC,H0,DH=[int(x*den)for x in (c0,dc,h0,dh)]
assert 1-h*(K-1)*(K-1+s)>0
F0=[];F1=[]
for k in range(K+1):
 k1=(k+1)**2;k2=-k*k+(8-4*b)*k+9-2*b*(b-1)
 t=(2*k1*A[k+1]+k2*A[k])//2+q*A[k]+(G[k]if b==0 else -G[k]);e=up(2*k1*Ar[k+1]+abs(k2)*Ar[k],2)+1+q*Ar[k]+Gr[k]
 F0.append(abs(A[k]+H0*t//den)+Ar[k]+up(H0*e,den)+1);F1.append(up(DH*(abs(t)+e),den))
best=None
for ac in range(1,9):
 for ah in range(1,8):
  rc,rh=2**ac,2**ah
  if den-(H0+rh*DH)*(K-1)*(K-1+s)<=0:continue
  v=0
  for k in range(K,-1,-1):
   p=k*(k+s)
   if k<K:v=up(2*v*(abs(C0-H0*p)+rc*DC+rh*DH*p),den-(H0+rh*DH)*p)
   v+=F0[k]+rh*F1[k]
  # Number of pairs (i,j) with ac*i+ah*j=m <= m+1. Tail sum exactly
  # sum_(m>=Q) (m+1)/2^m = (Q+2)*2^(1-Q).
  target=S//10**10
  Q=max(ac,ah)+1
  while 2*v*(Q+2)>target*(1<<Q):Q+=1
  eps=F(2*v*(Q+2),1<<Q);assert eps<=target
  mc,mh=(Q-1)//ac,(Q-1)//ah
  terms=sum((Q-1-ac*i)//ah+1 for i in range(mc+1))
  if best is None or terms<best[0]:best=(terms,ac,ah,Q,mc,mh,eps,v)
assert best
terms,ac,ah,Q,mc,mh,eps,v=best
name=f'triangle_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}'
print('PREFLIGHT',str(l),str(u),'b',b,'K',K,'cap',float(h),'ac,ah',ac,ah,'Q',Q,'mc,mh',mc,mh,'terms',terms,'analytic tail',float(eps/S),'log10M',math.log10(v)-D,flush=True)
with open(name+'.in','w')as f:
 f.write(f'{K} {ac} {ah} {Q} {b} {D} {den} {C0} {DC} {H0} {DH}\n');f.write(''.join(f'{a} {ar} {g} {gr}\n'for a,ar,g,gr in zip(A,Ar,G,Gr)))
st=time.monotonic()
with open(name+'.in')as f,open(name+'.out','w')as out,open(name+'.log','w')as err:subprocess.run(['./direct_sign_triangle'],stdin=f,stdout=out,stderr=err,check=True)
z={a:int(b)for a,b in map(str.split,open(name+'.out'))};lb=F(z['LOWER'])/S-eps/S
physical=lb*(1-u)**3/(1+q*h)**4-F('1e-18') if lb>=0 else lb-F('1e-18')
res=dict(c_lo=str(l),c_hi=str(u),b=b,K=K,ac=ac,ah=ah,Q=Q,mc=mc,mh=mh,terms=terms,hcap=str(h),analytic_tail=str(eps/S),poly_lower=str(lb),physical_lower=str(physical),center=float(F(z['CENTER'],S)),L1=float(F(z['L1'],S)),roundoff=float(F(z['ERROR'],S)),lower=float(physical),status='PASS_CONDITIONAL' if physical>0 else 'FAIL_INCONCLUSIVE',seconds=time.monotonic()-st)
json.dump(res,open(name+'.json','w'),indent=2);print({k:v for k,v in res.items()if k not in ('analytic_tail','poly_lower','physical_lower')},flush=True)
