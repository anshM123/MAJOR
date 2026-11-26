#!/usr/bin/env python3
"""Verify source-vector enclosure transfer, exact Cauchy tails, spin caps,
and final rational inequalities in the new extension ledgers.
The separate C++ runs are needed to reproduce the signed Taylor/Chebyshev sums.
"""
import json, sys, math
from fractions import Fraction as F
from pathlib import Path
from prefix_improvements import tail,TIGHT_ALLOW,conformal_tail
sys.set_int_max_str_digits(0)
def up(n,d):return -((-n)//d)
def run(include_tight=True):
 chains={0:[('triangle','.258','.262'),('triangle','.262','.264')],4:[('triangle','.258','.262'),('triangle','.262','.263')]}
 with open('direct_coefficients.txt') as f:
  N,DD=map(int,f.readline().split());source=[list(map(int,x.split()))for x in f];assert len(source)==N+1
 out=[];tops={}
 for b in (0,4):
  last=F('.258')
  for kind,ls,us in chains[b]:
   l,u=F(ls),F(us);assert l==last;last=u;tops[b]=u
   name=f'{kind}_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}'
   r=json.load(open(name+'.json'));z=dict((a,int(v))for a,v in map(str.split,open(name+'.out')))
   assert (F(r['c_lo']),F(r['c_hi']),r['b'])==(l,u,b)
   K=r['K'];s=2*b+1;q=b*(b+1);assert K+1<=N
   h=F(r['hcap']);assert h>=u/((1001-b)*(1002+b)) and 1-h*(K-1)*(K-1+s)>0
   c=u/F(1999,2000)
   if kind!='tight':
    assert (c.numerator**(K+1)<<2700)<=c.denominator**(K+1)
    allowance=F('1e-18')
   else:
    assert tail(u,K)<TIGHT_ALLOW;allowance=TIGHT_ALLOW+F('1e-18')
    assert F(r['physical_allowance'])==allowance
   with open(name+'.in')as f:
    head=list(map(int,f.readline().split()));a=[list(map(int,x.split())) for x in f]
   if kind=='direct':
    k,mc,mh,bb,D,den,C0,DC,H0,DH=head;sc=1;rc,rh=r['rc'],r['rh'];assert (mc,mh)==(r['mc'],r['mh'])
   else:
    k,ac,ah,Q,bb,D,den,C0,DC,H0,DH=head;sc=2;rc,rh=2**ac,2**ah
    assert (ac,ah,Q)==(r['ac'],r['ah'],r['Q']);assert sum((Q-1-ac*i)//ah+1 for i in range((Q-1)//ac+1))==r['terms']
   assert (k,bb)==(K,b) and len(a)==K+2 and 80<=D<=DD
   assert F(C0-DC,den)==l and F(C0+DC,den)==u and H0==DH and F(H0+DH,den)==h
   for i,(row,orig) in enumerate(zip(a,source)):
    div=10**(DD-D)*sc**i
    for m,e,m0,e0 in (row[:2]+orig[:2],row[2:]+orig[2:]):
     assert e>=0 and abs(m*div-m0)+e0<=e*div
   F0=[];F1=[]
   for i in range(K+1):
    A,Ar,G,Gr=a[i];An,Arn,_,_=a[i+1];k1=(i+1)**2;k2=-i*i+(8-4*b)*i+9-2*b*(b-1)
    t=(sc*k1*An+k2*A)//2+q*A+(G if b==0 else -G)
    e=up(sc*k1*Arn+abs(k2)*Ar,2)+1+q*Ar+Gr
    F0.append(abs(A+H0*t//den)+Ar+up(H0*e,den)+1);F1.append(up(DH*(abs(t)+e),den))
   v=0
   for i in range(K,-1,-1):
    p=i*(i+s)
    if i<K:
     gap=den-(H0+rh*DH)*p;assert gap>0
     v=up(sc*v*(abs(C0-H0*p)+rc*DC+rh*DH*p),gap)
    v+=F0[i]+rh*F1[i]
   S=10**D
   if kind=='direct':eps=F(v*rc*rh,(rc-1)*(rh-1))*(F(1,rc**(mc+1))+F(1,rh**(mh+1)))/S
   else:eps=F(2*v*(Q+2),(1<<Q)*S)
   assert eps==F(r['analytic_tail'])<=F('1e-10')
   assert z['LOWER']==z['CENTER']-z['L1']-z['ERROR'] and z['L1']>=0 and z['ERROR']>=0
   lb=F(z['LOWER'],S)-eps;assert lb==F(r['poly_lower'])>0
   H=lb*(1-u)**3/(1+q*h)**4-allowance
   assert H==F(r['physical_lower']) and H>0 and r['status']=='PASS_CONDITIONAL'
   print(kind,l,u,b,'K',K,'lower',float(H),flush=True)
   out.append(dict(kind=kind,c_lo=str(l),c_hi=str(u),b=b,K=K,lower_exact=str(H),lower_approx=float(H)))
 assert conformal_tail(5333)<F(1,2**110)
 cm=min(tops.values());json.dump(dict(status='PASS_CONDITIONAL',c_lo='.258',c_hi_b0=str(tops[0]),c_hi_b4=str(tops[4]),c_hi_both=str(cm),L_hi_both=-math.log1p(-float(cm)),rows=out),open('ext_ledger.json','w'),indent=2)
 print('PASS: exact source transfer, spin caps, pole gaps, Cauchy majorants/tails and final inequalities;',len(out),'channel cells')
if __name__=='__main__':run()
