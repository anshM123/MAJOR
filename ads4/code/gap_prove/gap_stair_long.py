#!/usr/bin/env python3
"""LONG-TABLE G1-capped STAIRCASE direct sign cell for 0.56<=c<=0.9 (Lemma G2) (variant of gap_triangle.py; prover copy of the author's direct_sign_triangle.py logic).

Differences from code/sliver_direct_cont/direct_sign_triangle.py (author's file is NOT modified):
 * h-cap (Lemma G1, GAP_LEMMAS.md):  hcap = min( up12(c_hi/((1001-b)(1002+b))), up12((1-s_lo)^2/998^2) ),
   where s_lo = isqrt(floor(c_lo*10^40))/10^20 <= sqrt(c_lo), so (1-s_lo)^2 >= (1-sqrt(c_lo))^2 (rigorous upward rounding),
   both caps rounded upward to denominator HDEN=10^15.
 * optional split of [0,hcap] into m equal h-pieces [hcap*i/m, hcap*(i+1)/m] (smaller effective h-radius:
   the complex pole condition becomes 1-(h0+2^ah*dh)p>0 with a smaller dh); each piece is an ordinary
   run of the author's unchanged C++ binary direct_sign_triangle (H0>=DH>0 holds for every piece).
 * the fixed-point scale D is an input (default 100); the C++ requires D>=80. Raise D when 2c>1 amplifies rounding.
The cell passes iff every h-piece passes; the reported lower bound is the minimum over pieces.
Pieces: h-boundaries 0=h_0<...<h_mb=hu<...<h_{mb+mt}=hcap, equal splits of [0,hu] (mb pieces) and [hu,hcap] (mt pieces),
hu=min(hcap, (1-sqrt(c_hi))^2/998^2 rounded to 1e-15) (hu is only a partition point; no rigour needed).
Staircase: piece i (h in [h_i,h_{i+1}]) uses the c-range [c_lo, c_i], c_i=min(c_hi, up6((1-998 r_i)^2)),
r_i=isqrt(floor(h_i 10^40))/10^20<=sqrt(h_i). By Lemma G1 every physical state with h>=h_i has sqrt(c)<=1-998 sqrt(h_i), so
c<=(1-998 sqrt(h_i))^2<=(1-998 r_i)^2<=c_i: the pieces cover all physical states of the cell. K=k_tail(c_i) per piece: least K with tail(c_i,K)<1e-7 (envelope 2^2337, tight_spectral_tail.py, copied below);
allowance 1e-7+1e-14 (Lemma G2). Table: coeff_long_17200.txt (env GAP_TABLE).
Usage: gap_stair.py c_lo c_hi b mb mt D
"""
import sys,subprocess,json,time,math,os
from math import lcm,isqrt
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
def up(n,d):return -((-n)//d)
R=F(1999,2000)
HDEN=10**15
TABLE=os.environ.get('GAP_TABLE','coeff_long_17200.txt')
TAILE=2337;RT=F(199,200);HT=F(1,998**2);ALLOW=F('1e-7')+F('1e-14')
def tail(c,K):
 t=F(c)/RT;assert 0<t<1
 return 2**TAILE*t**(K+1)*((1+20*HT)/(1-t)+2*HT*((K+5)**2/(1-t)+2*(K+5)*t/(1-t)**2+t*(1+t)/(1-t)**3))
def k_tail(c):
 lo,hi=1,2
 while tail(c,hi)>=F('1e-7'):hi*=2
 while hi-lo>1:
  m=(lo+hi)//2
  if tail(c,m)<F('1e-7'):hi=m
  else:lo=m
 return hi
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
def stair_c(ha,u):
 if ha==0:return u
 r=F(isqrt((ha.numerator*10**40)//ha.denominator),10**20);assert r*r<=ha and 998*r<1
 v=(1-998*r)**2;return min(u,F(up(v.numerator*10**6,v.denominator),10**6))
def hgrid(l,u,b,mb,mt):
 hcap=g1_cap(l,u,b)[0]
 hu=min(hcap,F(round(float((1-math.sqrt(float(u)))**2/998**2)*10**15),10**15)) if mt>0 else hcap
 hb_list=[hu*i/mb for i in range(mb+1)]+[hu+(hcap-hu)*j/mt for j in range(1,mt+1)]
 assert hb_list[0]==0 and hb_list[-1]==hcap and all(x<y for x,y in zip(hb_list,hb_list[1:]))
 return hb_list
def main():
 l,u=map(F,sys.argv[1:3]);b=int(sys.argv[3]);mb=int(sys.argv[4]);mt=int(sys.argv[5]);m=mb+mt
 D=int(sys.argv[6])
 i0,i1=0,m-1
 S=10**D;K=k_tail(u);s=2*b+1;q=b*(b+1)
 assert F('.56')<=l<u<=F('.9') and b in (0,4) and mb>=1 and mt>=0 and D>=80
 hcap,cap1,cap2=g1_cap(l,u,b)
 assert 1-hcap*(K-1)*(K-1+s)>0
 with open(TABLE) as f:
  N,DD=map(int,f.readline().split());assert N>=K+1 and DD>=D
  div=10**(DD-D);A=[];Ar=[];G=[];Gr=[]
  for k in range(K+2):
   a,r,g,t=map(int,f.readline().split());d=div<<k;A.append(a//d);Ar.append(up(r,d)+1);G.append(g//d);Gr.append(up(t,d)+1)
 rows=[]
 hb_list=hgrid(l,u,b,mb,mt)
 for pi in range(i0,i1+1):
  ha,hb=hb_list[pi],hb_list[pi+1]
  ui=stair_c(ha,u);K=k_tail(ui);assert tail(ui,K)<F('1e-7') and l<ui<=u and 1-hcap*(K-1)*(K-1+s)>0
  c0=(l+ui)/2;dc=(ui-l)/2;h0=(ha+hb)/2;dh=(hb-ha)/2;den=lcm(c0.denominator,dc.denominator,h0.denominator,dh.denominator);C0,DC,H0,DH=[int(x*den)for x in (c0,dc,h0,dh)]
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
  name=f'g2st_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_m{mb}x{mt}_p{pi}_D{D}'
  print('PREFLIGHT',name,str(l),str(u),'piece c_hi',float(ui),'b',b,'K',K,'hcap',float(hcap),'cap1',float(cap1),'cap2',float(cap2),'piece',float(ha),float(hb),'ac,ah',ac,ah,'Q',Q,'mc,mh',mc,mh,'terms',terms,'analytic tail',float(eps/S),'log10M',math.log10(v)-D,flush=True)
  with open(name+'.in','w')as f:
   f.write(f'{K} {ac} {ah} {Q} {b} {D} {den} {C0} {DC} {H0} {DH}\n');f.write(''.join(f'{a} {ar} {g} {gr}\n'for a,ar,g,gr in zip(A[:K+2],Ar[:K+2],G[:K+2],Gr[:K+2])))
  st=time.monotonic()
  with open(name+'.in')as f,open(name+'.out','w')as out,open(name+'.log','w')as err:subprocess.run([BIN],stdin=f,stdout=out,stderr=err,check=True)
  z={a:int(bb)for a,bb in map(str.split,open(name+'.out'))};lb=F(z['LOWER'])/S-eps/S
  physical=lb*(1-ui)**3/(1+q*hb)**4-ALLOW if lb>=0 else lb-ALLOW
  res=dict(name=name,c_lo=str(l),c_hi=str(u),c_hi_piece=str(ui),b=b,K=K,D=D,mb=mb,mt=mt,piece=pi,h_lo=str(ha),h_hi=str(hb),hcap=str(hcap),cap1=str(cap1),cap2=str(cap2),ac=ac,ah=ah,Q=Q,mc=mc,mh=mh,terms=terms,analytic_tail=str(eps/S),poly_lower=str(lb),physical_lower=str(physical),center=float(F(z['CENTER'],S)),L1=float(F(z['L1'],S)),roundoff=float(F(z['ERROR'],S)),lower=float(physical),allowance=str(ALLOW),status='PASS_CONDITIONAL' if physical>0 else 'FAIL_INCONCLUSIVE',seconds=time.monotonic()-st)
  json.dump(res,open(name+'.json','w'),indent=2);print({k:v for k,v in res.items()if k not in ('analytic_tail','poly_lower','physical_lower','h_lo','h_hi','hcap','cap1','cap2')},flush=True)
  rows.append(res)
  if res['status']!='PASS_CONDITIONAL':break  # cell already failed; skip remaining pieces
 print('CELL',str(l),str(u),'b',b,'mb,mt',mb,mt,'min lower',min(r['lower'] for r in rows),'ALLPASS',all(r['status']=='PASS_CONDITIONAL' for r in rows),flush=True)
if __name__=='__main__':main()
