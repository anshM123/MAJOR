#!/usr/bin/env python3
"""G1-capped direct sign cell (prover copy of the author's direct_sign_triangle.py logic).

Differences from code/sliver_direct_cont/direct_sign_triangle.py (author's file is NOT modified):
 * h-cap (Lemma G1, GAP_LEMMAS.md):  hcap = min( up12(c_hi/((1001-b)(1002+b))), up12((1-s_lo)^2/998^2) ),
   where s_lo = isqrt(floor(c_lo*10^40))/10^20 <= sqrt(c_lo), so (1-s_lo)^2 >= (1-sqrt(c_lo))^2 (rigorous upward rounding),
   both caps rounded upward to denominator HDEN=10^15.
 * optional split of [0,hcap] into m equal h-pieces [hcap*i/m, hcap*(i+1)/m] (smaller effective h-radius:
   the complex pole condition becomes 1-(h0+2^ah*dh)p>0 with a smaller dh); each piece is an ordinary
   run of the author's unchanged C++ binary direct_sign_triangle (H0>=DH>0 holds for every piece).
 * the fixed-point scale D is an input (default 100); the C++ requires D>=80. Raise D when 2c>1 amplifies rounding.
The cell passes iff every h-piece passes; the reported lower bound is the minimum over pieces.
Usage: gap_triangle.py c_lo c_hi b [m] [D] [i_first i_last]
"""
import sys,subprocess,json,time,math,os
from math import lcm,isqrt
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
def up(n,d):return -((-n)//d)
R=F(1999,2000)
HDEN=10**15
TABLE=os.environ.get('GAP_TABLE','direct_coefficients.txt')
BIN=os.environ.get('GAP_BIN','./direct_sign_triangle')
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
def g1_cap(l,u,b):
 cap1=F(up(u.numerator*HDEN,u.denominator*(1001-b)*(1002+b)),HDEN)
 slo=F(isqrt((l.numerator*10**40)//l.denominator),10**20)
 assert slo*slo<=l and 0<slo<1
 v=(1-slo)**2/998**2
 cap2=F(up(v.numerator*HDEN,v.denominator),HDEN)
 return min(cap1,cap2),cap1,cap2
def main():
 l,u=map(F,sys.argv[1:3]);b=int(sys.argv[3]);m=int(sys.argv[4]) if len(sys.argv)>4 else 1
 D=int(sys.argv[5]) if len(sys.argv)>5 else 100
 i0,i1=(int(sys.argv[6]),int(sys.argv[7])) if len(sys.argv)>7 else (0,m-1)
 S=10**D;K=need_k(u);s=2*b+1;q=b*(b+1)
 assert F('.05')<=l<u<=F('.56') and b in (0,4) and m>=1 and D>=80
 hcap,cap1,cap2=g1_cap(l,u,b)
 assert 1-hcap*(K-1)*(K-1+s)>0
 with open(TABLE) as f:
  N,DD=map(int,f.readline().split());assert N>=K+1 and DD>=D
  div=10**(DD-D);A=[];Ar=[];G=[];Gr=[]
  for k in range(K+2):
   a,r,g,t=map(int,f.readline().split());d=div<<k;A.append(a//d);Ar.append(up(r,d)+1);G.append(g//d);Gr.append(up(t,d)+1)
 rows=[]
 for pi in range(i0,i1+1):
  ha,hb=hcap*pi/m,hcap*(pi+1)/m
  c0=(l+u)/2;dc=(u-l)/2;h0=(ha+hb)/2;dh=(hb-ha)/2;den=lcm(c0.denominator,dc.denominator,h0.denominator,dh.denominator);C0,DC,H0,DH=[int(x*den)for x in (c0,dc,h0,dh)]
  assert H0>=DH>0
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
    target=S//10**10
    Q=max(ac,ah)+1
    while 2*v*(Q+2)>target*(1<<Q):Q+=1
    eps=F(2*v*(Q+2),1<<Q);assert eps<=target
    mc,mh=(Q-1)//ac,(Q-1)//ah
    terms=sum((Q-1-ac*i)//ah+1 for i in range(mc+1))
    if best is None or terms<best[0]:best=(terms,ac,ah,Q,mc,mh,eps,v)
  assert best,'no admissible radius (complex pole): increase m'
  terms,ac,ah,Q,mc,mh,eps,v=best
  name=f'g1tri_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_m{m}_p{pi}_D{D}'
  print('PREFLIGHT',name,str(l),str(u),'b',b,'K',K,'hcap',float(hcap),'cap1',float(cap1),'cap2',float(cap2),'piece',float(ha),float(hb),'ac,ah',ac,ah,'Q',Q,'mc,mh',mc,mh,'terms',terms,'analytic tail',float(eps/S),'log10M',math.log10(v)-D,flush=True)
  with open(name+'.in','w')as f:
   f.write(f'{K} {ac} {ah} {Q} {b} {D} {den} {C0} {DC} {H0} {DH}\n');f.write(''.join(f'{a} {ar} {g} {gr}\n'for a,ar,g,gr in zip(A,Ar,G,Gr)))
  st=time.monotonic()
  with open(name+'.in')as f,open(name+'.out','w')as out,open(name+'.log','w')as err:subprocess.run([BIN],stdin=f,stdout=out,stderr=err,check=True)
  z={a:int(bb)for a,bb in map(str.split,open(name+'.out'))};lb=F(z['LOWER'])/S-eps/S
  physical=lb*(1-u)**3/(1+q*hb)**4-F('1e-18') if lb>=0 else lb-F('1e-18')
  res=dict(name=name,c_lo=str(l),c_hi=str(u),b=b,K=K,D=D,m=m,piece=pi,h_lo=str(ha),h_hi=str(hb),hcap=str(hcap),cap1=str(cap1),cap2=str(cap2),ac=ac,ah=ah,Q=Q,mc=mc,mh=mh,terms=terms,analytic_tail=str(eps/S),poly_lower=str(lb),physical_lower=str(physical),center=float(F(z['CENTER'],S)),L1=float(F(z['L1'],S)),roundoff=float(F(z['ERROR'],S)),lower=float(physical),status='PASS_CONDITIONAL' if physical>0 else 'FAIL_INCONCLUSIVE',seconds=time.monotonic()-st)
  json.dump(res,open(name+'.json','w'),indent=2);print({k:v for k,v in res.items()if k not in ('analytic_tail','poly_lower','physical_lower','h_lo','h_hi','hcap','cap1','cap2')},flush=True)
  rows.append(res)
  if res['status']!='PASS_CONDITIONAL':break  # cell already failed; skip remaining pieces
 print('CELL',str(l),str(u),'b',b,'m',m,'min lower',min(r['lower'] for r in rows),'ALLPASS',all(r['status']=='PASS_CONDITIONAL' for r in rows),flush=True)
if __name__=='__main__':main()
