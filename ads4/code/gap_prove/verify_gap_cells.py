#!/usr/bin/env python3
"""Verify G1-capped direct sign cells (modelled on the author's verify_extension.py).
For every piece of every cell listed in cells.json (c-chain per b, contiguous from START):
 * exact transfer of every input coefficient ball from the K=3299 table (scale 10^D, A_k,G_k * 2^-k);
 * prefix condition (c_hi/R)^(K+1) 2^2700 <= 1, R=1999/2000; spin+twist cap (G1) recomputed independently,
   hcap <= 998^-2, real pole gap 1-hcap(K-1)(K-1+s)>0;
 * h-pieces tile [0,hcap] exactly; header encodes exactly [c_lo,c_hi] x [h_lo,h_hi];
 * Cauchy majorant v (complex pole gaps den-(H0+rh DH)p>0) and triangle tail eps=2v(Q+2)/2^Q recomputed;
 * LOWER == CENTER-L1-ERROR, L1,ERROR>=0, poly lower = LOWER/S-eps > 0;
 * physical lower = poly_lower (1-c_hi)^3/(1+q h_hi)^4 - 1e-18 > 0, equal to the json value.
Writes ledger.json (in this script's directory). The C++ sign runs themselves are not re-executed here.
Usage: verify_gap_cells.py RUNDIR TABLE
"""
import json,sys,math,os
from math import isqrt
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
def up(n,d):return -((-n)//d)
HERE=os.path.dirname(os.path.abspath(__file__))
START=F('.28')
def caps(l,u,b):
 c1=F(up(u.numerator*10**15,u.denominator*(1001-b)*(1002+b)),10**15)
 slo=F(isqrt((l.numerator*10**40)//l.denominator),10**20);assert slo*slo<=l and 0<slo<1
 v=(1-slo)**2/998**2;c2=F(up(v.numerator*10**15,v.denominator),10**15)
 assert c1>=u/((1001-b)*(1002+b))
 # independent exact check that c2 >= (1-sqrt(l))^2/998^2: exhibit r<=sqrt(c2*998^2) with (1-r)^2<=l
 w=c2*998**2;r=F(isqrt((w.numerator*10**60)//w.denominator),10**30);assert r*r<=w and r<1 and (1-r)**2<=l
 return min(c1,c2)
def run(rundir,table):
 cells=json.load(open(os.path.join(HERE,'cells.json')))
 with open(table) as f:
  N,DD=map(int,f.readline().split());source=[list(map(int,x.split()))for x in f];assert len(source)==N+1
 out=[];tops={}
 for b in (0,4):
  last=START
  for cell in [c for c in cells if c['b']==b]:
   l,u,m,D=F(cell['c_lo']),F(cell['c_hi']),cell['m'],cell['D'];assert l==last and l<u<=F('.56');last=u;tops[b]=u
   s=2*b+1;q=b*(b+1);S=10**D;hcap=caps(l,u,b);assert hcap<=F(1,998**2)
   lows=[];hprev=F(0);secs=0
   for pi in range(m):
    name=f'g1tri_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_m{m}_p{pi}_D{D}'
    r=json.load(open(os.path.join(rundir,name+'.json')));z=dict((a,int(v))for a,v in map(str.split,open(os.path.join(rundir,name+'.out'))))
    assert (F(r['c_lo']),F(r['c_hi']),r['b'],r['m'],r['piece'],r['D'])==(l,u,b,m,pi,D) and F(r['hcap'])==hcap
    K=r['K'];assert K+1<=N
    c=u/F(1999,2000);assert (c.numerator**(K+1)<<2700)<=c.denominator**(K+1)
    assert 1-hcap*(K-1)*(K-1+s)>0
    ha,hb=F(r['h_lo']),F(r['h_hi']);assert ha==hprev and hb==hcap*(pi+1)/m;hprev=hb
    with open(os.path.join(rundir,name+'.in'))as f:
     head=list(map(int,f.readline().split()));a=[list(map(int,x.split())) for x in f]
    k,ac,ah,Q,bb,DD2,den,C0,DC,H0,DH=head;rc,rh=2**ac,2**ah
    assert (k,bb,DD2)==(K,b,D) and len(a)==K+2 and 80<=D<=DD and (ac,ah,Q)==(r['ac'],r['ah'],r['Q'])
    assert sum((Q-1-ac*i)//ah+1 for i in range((Q-1)//ac+1))==r['terms']
    assert F(C0-DC,den)==l and F(C0+DC,den)==u and F(H0-DH,den)==ha and F(H0+DH,den)==hb and H0>=DH>0 and DC>0
    for i,(row,orig) in enumerate(zip(a,source)):
     div=10**(DD-D)*2**i
     for mm,e,m0,e0 in (row[:2]+orig[:2],row[2:]+orig[2:]):
      assert e>=0 and abs(mm*div-m0)+e0<=e*div
    F0=[];F1=[]
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
    eps=F(2*v*(Q+2),(1<<Q)*S)
    assert eps==F(r['analytic_tail'])<=F('1e-10')
    assert z['LOWER']==z['CENTER']-z['L1']-z['ERROR'] and z['L1']>=0 and z['ERROR']>=0
    lb=F(z['LOWER'],S)-eps;assert lb==F(r['poly_lower'])>0
    H=lb*(1-u)**3/(1+q*hb)**4-F('1e-18')
    assert H==F(r['physical_lower']) and H>0 and r['status']=='PASS_CONDITIONAL'
    lows.append(H);secs+=r['seconds']
   assert hprev==hcap
   H=min(lows)
   print('cell',l,u,'b',b,'K',K,'m',m,'D',D,'lower',float(H),'sec',round(secs),flush=True)
   out.append(dict(c_lo=str(l),c_hi=str(u),b=b,K=K,m=m,D=D,hcap=str(hcap),lower_exact=str(H),lower_approx=float(H),cpu_seconds=round(secs,1)))
 cm=min(tops.values()) if len(tops)==2 else START
 json.dump(dict(status='PASS_CONDITIONAL',note='G1-capped direct sign cells; conditional on the retained inputs of the direct method (seed allowance 1e-18 valid for c<=0.56)',c_lo=str(START),c_hi_b0=str(tops.get(0,START)),c_hi_b4=str(tops.get(4,START)),c_hi_both=str(cm),L_hi_both=-math.log1p(-float(cm)),rows=out),open(os.path.join(HERE,'ledger.json'),'w'),indent=2)
 print('PASS: transfer, prefix K, G1 caps, h-tiling, pole gaps, Cauchy majorants/tails, final inequalities;',len(out),'channel cells; top',{b:str(t) for b,t in tops.items()})
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
