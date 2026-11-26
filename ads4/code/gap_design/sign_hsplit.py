#!/usr/bin/env python3
"""GAP-DESIGN variant of the author's sign_next.py (sliver_full_band_reduction): identical checker, but the
h-range is a sub-interval [hlo_frac*H, hhi_frac*H] of [0,H], H = min(spin cap(c_hi), twist cap(c_lo)) (G1),
so that cells with 1.5 H K^2 >= 1 can be covered by several h-pieces. Usage:
  python sign_hsplit.py c_lo c_hi b K coeff_file hlo_frac hhi_frac   (fractions as exact decimals, e.g. 0 .5)
"""
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
l,u=map(F,sys.argv[1:3]);b=int(sys.argv[3]);D=100;S=10**D;K=1251;s=2*b+1;q=b*(b+1)
from tight_spectral_tail import tail
from base_bounds import sqrtlo
TIGHT_ALLOW=F("1e-7")
assert 0<F(".05")<=l<u<=F(".9") and b in (0,4)
K=int(sys.argv[4]) if len(sys.argv)>4 else 1251
assert tail(u,K)<TIGHT_ALLOW
allowance=TIGHT_ALLOW+(F('1e-18') if u<=F('.56') else F('1e-14'))
h=min(F(up(u.numerator*10**12,u.denominator*(1001-b)*(1002+b)),10**12),F(math.ceil((1-sqrtlo(l))**2*10**15/998**2),10**15))
hl_f,hu_f=(F(sys.argv[6]),F(sys.argv[7])) if len(sys.argv)>7 else (F(0),F(1));assert 0<=hl_f<hu_f<=1
hlo,hhi=h*hl_f,h*hu_f
c0=(l+u)/2;dc=(u-l)/2;h0=(hlo+hhi)/2;dh=(hhi-hlo)/2;den=lcm(c0.denominator,dc.denominator,h0.denominator,dh.denominator);C0,DC,H0,DH=[int(x*den)for x in (c0,dc,h0,dh)]
assert 1-hhi*(K-1)*(K-1+s)>0
# GAP-DESIGN: precision guard. The C++ l1 rounding-error recurrence multiplies E by
# g_k=(|a|+dc+2r)/(d-r) per step (2^-k storage doubles the ratio); its product can exceed 10^100 for c>0.5,
# which destroys the D=100 working precision. Choose D so that 10^-D * prod g_k stays tiny.
lg=0.0
for k in range(K-1,-1,-1):
 p=k*(k+s);a2=2*abs(C0-p*H0);gk=(a2+2*DC+2*p*DH)/(den-p*H0-p*DH);lg=max(0.0,lg+math.log10(gk)) if gk>0 else lg
D=max(D,100+int(lg)+30);S=10**D
print('PRECISION growth log10',round(lg,1),'D',D,flush=True)
with open(sys.argv[5] if len(sys.argv)>5 else 'old_coefficients.txt') as f:
 N,DD=map(int,f.readline().split());assert N>=K+1 and DD>=D
 div=10**(DD-D);A=[];Ar=[];G=[];Gr=[]
 for k in range(K+2):
  a,r,g,t=map(int,f.readline().split());d=div<<k;A.append(a//d);Ar.append(up(r,d)+1);G.append(g//d);Gr.append(up(t,d)+1)
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
name=f'hs_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_{hl_f.numerator}_{hl_f.denominator}_{hu_f.numerator}_{hu_f.denominator}'
print('PREFLIGHT',str(l),str(u),'b',b,'K',K,'cap',float(h),'ac,ah',ac,ah,'Q',Q,'mc,mh',mc,mh,'terms',terms,'analytic tail',float(eps/S),'log10M',math.log10(v)-D,flush=True)
with open(name+'.in','w')as f:
 f.write(f'{K} {ac} {ah} {Q} {b} {D} {den} {C0} {DC} {H0} {DH}\n');f.write(''.join(f'{a} {ar} {g} {gr}\n'for a,ar,g,gr in zip(A,Ar,G,Gr)))
st=time.monotonic()
with open(name+'.in')as f,open(name+'.out','w')as out,open(name+'.log','w')as err:subprocess.run(['./direct_sign_triangle'],stdin=f,stdout=out,stderr=err,check=True)
z={a:int(b)for a,b in map(str.split,open(name+'.out'))};lb=F(z['LOWER'])/S-eps/S
physical=lb*(1-u)**3/(1+q*h)**4-allowance if lb>=0 else lb-allowance
res=dict(c_lo=str(l),c_hi=str(u),b=b,K=K,h_lo=str(hlo),h_hi=str(hhi),ac=ac,ah=ah,Q=Q,mc=mc,mh=mh,terms=terms,hcap=str(h),analytic_tail=str(eps/S),poly_lower=str(lb),physical_lower=str(physical),center=float(F(z['CENTER'],S)),L1=float(F(z['L1'],S)),roundoff=float(F(z['ERROR'],S)),lower=float(physical),physical_allowance=str(allowance),status='PASS_CONDITIONAL' if physical>0 else 'FAIL_INCONCLUSIVE',seconds=time.monotonic()-st)
json.dump(res,open(name+'.json','w'),indent=2);print({k:v for k,v in res.items()if k not in ('analytic_tail','poly_lower','physical_lower')},flush=True)
