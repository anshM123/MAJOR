#!/usr/bin/env python3
"""Exact input transfer, pole gap, analytic tails and final sign ledgers."""
import json,sys,math
from fractions import Fraction as F
from base_bounds import sqrtlo
from tight_spectral_tail import tail
sys.set_int_max_str_digits(0)
def up(n,d):return -((-n)//d)
last=F('.268');out=[]
for lo,hi,src in [('.268','.28','coeff_long_17200_prefix.txt'),('.28','.29','coeff_long_17200_prefix.txt'),('.29','.30','coeff_long_17200_prefix.txt')]:
 l,u=F(lo),F(hi);assert l==last
 results=[]
 for b in (0,4):
  name=f'next_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}'
  r=json.load(open(name+'.json'));assert (F(r['c_lo']),F(r['c_hi']),r['b'])==(l,u,b)
  if r['status']!='PASS_CONDITIONAL':break
  with open(name+'.in') as f:
   K,ac,ah,Q,bb,D,den,C0,DC,H0,DH=map(int,f.readline().split());a=[list(map(int,x.split()))for x in f]
  assert (K,ac,ah,Q,bb)==(r['K'],r['ac'],r['ah'],r['Q'],b) and len(a)==K+2
  assert F(C0-DC,den)==l and F(C0+DC,den)==u and H0==DH
  h=F(H0+DH,den);assert h==F(r['hcap']);s=2*b+1;q=b*(b+1)
  assert h>=min(u/((1001-b)*(1002+b)),(1-sqrtlo(l))**2/998**2)
  assert 1-h*(K-1)*(K-1+s)>0
  assert tail(u,K)<F('1e-7');allowance=F('1e-7')+F('1e-18');assert F(r['physical_allowance'])==allowance and u<=F('.56')
  with open(src)as f:
   N,DD=map(int,f.readline().split());assert N>=K+1 and DD>=D
   for i,row in enumerate(a):
    orig=list(map(int,f.readline().split()));div=10**(DD-D)<<i
    for m,e,m0,e0 in (row[:2]+orig[:2],row[2:]+orig[2:]):assert e>=0 and abs(m*div-m0)+e0<=e*div
  F0=[];F1=[];rc,rh=2**ac,2**ah
  for i in range(K+1):
   A,Ar,G,Gr=a[i];An,Arn,_,_=a[i+1];k1=(i+1)**2;k2=-i*i+(8-4*b)*i+9-2*b*(b-1)
   t=(2*k1*An+k2*A)//2+q*A+(G if b==0 else -G)
   e=up(2*k1*Arn+abs(k2)*Ar,2)+1+q*Ar+Gr
   F0.append(abs(A+H0*t//den)+Ar+up(H0*e,den)+1);F1.append(up(DH*(abs(t)+e),den))
  v=0
  for i in range(K,-1,-1):
   p=i*(i+s)
   if i<K:
    gap=den-(H0+rh*DH)*p;assert gap>0
    v=up(2*v*(abs(C0-H0*p)+rc*DC+rh*DH*p),gap)
   v+=F0[i]+rh*F1[i]
  S=10**D;eps=F(2*v*(Q+2),(1<<Q)*S)
  assert eps==F(r['analytic_tail'])<=F('1e-10')
  z={a:int(v)for a,v in map(str.split,open(name+'.out'))};assert z['LOWER']==z['CENTER']-z['L1']-z['ERROR'] and z['L1']>=0 and z['ERROR']>=0
  lb=F(z['LOWER'],S)-eps;H=lb*(1-u)**3/(1+q*h)**4-allowance
  assert lb==F(r['poly_lower'])>0 and H==F(r['physical_lower'])>0
  results.append(dict(c_lo=str(l),c_hi=str(u),b=b,K=K,lower_exact=str(H),lower=float(H)))
 if len(results)!=2:break
 out.extend(results);last=u
 for r in results:print('PASS',r['c_lo'],r['c_hi'],r['b'],r['lower'])
assert last>=F('.28')
json.dump(dict(status='PARTIAL_SIGN_COVERAGE',c_lo='.268',c_hi=str(last),open_L_lower=-math.log1p(-float(last)),open_L_upper=2.3,rows=out),open('new_sign_ledger.json','w'),indent=2)
print('PASS executed sign ledger; upper endpoint',last)
