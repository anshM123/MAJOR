#!/usr/bin/env python3
"""EXECUTOR-2 copy of code/gap_prove/verify_gap_cells.py (all checks identical), adapted to:
 * the executor-2 range: cells.json / ledger.json in code/gap_exec2/, chain per b contiguous from START=0.67 and required to reach END=0.777;
 * kind 'longu' (gap_exec2/gap_stair_long_u.py, unscaled checker direct_sign_triangle_u): transfer check WITHOUT the 2^-k factor,
   B-coefficient with k1*A_(k+1) and Cauchy-majorant recurrence without the factor 2 (exactly the unscaled C++ recurrence);
 * extra checks: SHA-256 of the long table against code/sliver_fullband/LARGE_FILES_SHA256.txt; the designer's overshoot
   omega=(1-sqrt(c_lo))^2/(1-sqrt(c_hi))^2 <= 1.08 (design rule only; rigour comes from G1 + staircase/h-tiling).
Original docstring follows.
Verify G1-capped direct sign cells (modelled on the author's verify_extension.py).
For every piece of every cell listed in cells.json (c-chain per b, contiguous from START):
 * exact transfer of every input coefficient ball from the K=3299 table (scale 10^D, A_k,G_k * 2^-k);
 * prefix condition (c_hi/R)^(K+1) 2^2700 <= 1, R=1999/2000; spin+twist cap (G1) recomputed independently,
   hcap <= 998^-2, real pole gap 1-hcap(K-1)(K-1+s)>0;
 * h-pieces tile [0,hcap] exactly; header encodes exactly [c_lo,c_hi(piece)] x [h_lo,h_hi];
 * staircase cells (gap_stair.py): piece c_hi >= (1-998 sqrt(h_lo))^2 checked exactly (G1), prefix K per piece;
 * Cauchy majorant v (complex pole gaps den-(H0+rh DH)p>0) and triangle tail eps=2v(Q+2)/2^Q recomputed;
 * LOWER == CENTER-L1-ERROR, L1,ERROR>=0, poly lower = LOWER/S-eps > 0;
 * physical lower = poly_lower (1-c_hi)^3/(1+q h_hi)^4 - 1e-18 > 0, equal to the json value.
Writes ledger.json (in this script's directory). The C++ sign runs themselves are not re-executed here.
Long cells (kind 'long', gap_stair_long.py, c>=0.56, Lemma G2): table coeff_long_17200.txt (LONGTABLE), prefix condition
tail(c_hi_piece,K)<1e-7 with the exact envelope-2^2337 tail (re-implemented here), allowance 1e-7+1e-14.
Usage: verify_gap_cells.py RUNDIR TABLE [LONGTABLE]
"""
import json,sys,math,os
from math import isqrt
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
def up(n,d):return -((-n)//d)
HERE=os.path.dirname(os.path.abspath(__file__))
START=F('.67');END=F('.777')
import hashlib
def sha_ok(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for blk in iter(lambda:f.read(1<<24),b''):h.update(blk)
 ref=[l.split()[0] for l in open(os.path.join(HERE,'..','sliver_fullband','LARGE_FILES_SHA256.txt')) if l.split()[1]=='coeff_long_17200.txt'][0]
 assert h.hexdigest()==ref,'long table SHA-256 mismatch';print('long table sha256 OK',ref,flush=True)
def omega(l,u):return ((1-math.sqrt(float(l)))/(1-math.sqrt(float(u))))**2
def caps(l,u,b):
 c1=F(up(u.numerator*10**15,u.denominator*(1001-b)*(1002+b)),10**15)
 slo=F(isqrt((l.numerator*10**40)//l.denominator),10**20);assert slo*slo<=l and 0<slo<1
 v=(1-slo)**2/998**2;c2=F(up(v.numerator*10**15,v.denominator),10**15)
 assert c1>=u/((1001-b)*(1002+b))
 # independent exact check that c2 >= (1-sqrt(l))^2/998^2: exhibit r<=sqrt(c2*998^2) with (1-r)^2<=l
 w=c2*998**2;r=F(isqrt((w.numerator*10**60)//w.denominator),10**30);assert r*r<=w and r<1 and (1-r)**2<=l
 return min(c1,c2)
TAILE=2337;RT=F(199,200);HT=F(1,998**2);ALLOW=F('1e-7')+F('1e-14')
def tail(c,K):
 t=F(c)/RT;assert 0<t<1
 return 2**TAILE*t**(K+1)*((1+20*HT)/(1-t)+2*HT*((K+5)**2/(1-t)+2*(K+5)*t/(1-t)**2+t*(1+t)/(1-t)**3))
def run(rundir,table,longtable=None):
 cells=json.load(open(os.path.join(HERE,'cells.json')))
 sha_ok(longtable)
 tabs={}
 def load(path,kmax):
  if path not in tabs or len(tabs[path][2])<kmax+2:
   with open(path) as f:
    N,DD=map(int,f.readline().split());rows=[]
    for i in range(min(N+1,kmax+2)):rows.append(list(map(int,f.readline().split())))
   tabs[path]=(N,DD,rows)
  return tabs[path]
 kmax_long=max([c['K'] for c in cells if c.get('kind') in ('long','longu')],default=0)+2
 out=[];tops={}
 for b in (0,4):
  last=START
  for cell in [c for c in cells if c['b']==b]:
   l,u,m,D=F(cell['c_lo']),F(cell['c_hi']),cell['m'],cell['D'];assert l==last and l<u<=F('.9');last=u;tops[b]=u
   kind=cell.get('kind','tri');stair=kind in ('stair','long','longu');lng=kind in ('long','longu');uns=kind=='longu'
   assert lng and omega(l,u)<=1.08,('omega',float(l),float(u),omega(l,u))
   N,DD,source=load(longtable,kmax_long) if lng else load(table,3400)
   s=2*b+1;q=b*(b+1);S=10**D;hcap=caps(l,u,b);assert hcap<=F(1,998**2)
   lows=[];hprev=F(0);secs=0;Kmax=0
   for pi in range(m):
    if uns:name=f'g2su_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_m{cell["mb"]}x{cell["mt"]}_p{pi}_D{D}'
    elif lng:name=f'g2st_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_m{cell["mb"]}x{cell["mt"]}_p{pi}_D{D}'
    elif stair:name=f'g1st_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_m{cell["mb"]}x{cell["mt"]}_p{pi}_D{D}'
    else:name=f'g1tri_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_m{m}_p{pi}_D{D}'
    r=json.load(open(os.path.join(rundir,name+'.json')));z=dict((a,int(v))for a,v in map(str.split,open(os.path.join(rundir,name+'.out'))))
    assert (F(r['c_lo']),F(r['c_hi']),r['b'],r.get('m',r.get('mb',0)+r.get('mt',0)),r['piece'],r['D'])==(l,u,b,m,pi,D) and F(r['hcap'])==hcap
    if stair:assert (r['mb'],r['mt'])==(cell['mb'],cell['mt'])
    K=r['K'];assert K+1<=N;Kmax=max(Kmax,K)
    ha,hb=F(r['h_lo']),F(r['h_hi']);assert ha==hprev and hb>ha;hprev=hb
    ui=F(r['c_hi_piece']) if stair else u
    if stair:
     assert l<ui<=u
     if ha>0 and ui<u:
      # staircase: physical states with h>=ha have c<=(1-998 sqrt(ha))^2 (G1); exhibit r<=sqrt(ha) with (1-998r)^2<=ui
      rr=F(isqrt((ha.numerator*10**60)//ha.denominator),10**30);assert rr*rr<=ha and 0<=998*rr<=1 and (1-998*rr)**2<=ui
    else:assert hb==hcap*(pi+1)/m
    if lng:assert tail(ui,K)<F('1e-7') and F(r['allowance'])==ALLOW and F('.56')<=l
    else:c=ui/F(1999,2000);assert (c.numerator**(K+1)<<2700)<=c.denominator**(K+1) and u<=F('.56')
    assert 1-hcap*(K-1)*(K-1+s)>0
    with open(os.path.join(rundir,name+'.in'))as f:
     head=list(map(int,f.readline().split()));a=[list(map(int,x.split())) for x in f]
    k,ac,ah,Q,bb,DD2,den,C0,DC,H0,DH=head;rc,rh=2**ac,2**ah
    assert (k,bb,DD2)==(K,b,D) and len(a)==K+2 and 80<=D<=DD and (ac,ah,Q)==(r['ac'],r['ah'],r['Q'])
    assert sum((Q-1-ac*i)//ah+1 for i in range((Q-1)//ac+1))==r['terms']
    assert F(C0-DC,den)==l and F(C0+DC,den)==ui and F(H0-DH,den)==ha and F(H0+DH,den)==hb and H0>=DH>0 and DC>0
    for i,(row,orig) in enumerate(zip(a,source)):
     div=10**(DD-D)*(1 if uns else 2**i)
     for mm,e,m0,e0 in (row[:2]+orig[:2],row[2:]+orig[2:]):
      assert e>=0 and abs(mm*div-m0)+e0<=e*div
    F0=[];F1=[]
    for i in range(K+1):
     A,Ar,G,Gr=a[i];An,Arn,_,_=a[i+1];k1=(i+1)**2;k2=-i*i+(8-4*b)*i+9-2*b*(b-1)
     t=((1 if uns else 2)*k1*An+k2*A)//2+q*A+(G if b==0 else -G)
     e=up((1 if uns else 2)*k1*Arn+abs(k2)*Ar,2)+1+q*Ar+Gr
     F0.append(abs(A+H0*t//den)+Ar+up(H0*e,den)+1);F1.append(up(DH*(abs(t)+e),den))
    v=0
    for i in range(K,-1,-1):
     p=i*(i+s)
     if i<K:
      gap=den-(H0+rh*DH)*p;assert gap>0
      v=up((1 if uns else 2)*v*(abs(C0-H0*p)+rc*DC+rh*DH*p),gap)
     v+=F0[i]+rh*F1[i]
    eps=F(2*v*(Q+2),(1<<Q)*S)
    assert eps==F(r['analytic_tail'])<=F('1e-10')
    assert z['LOWER']==z['CENTER']-z['L1']-z['ERROR'] and z['L1']>=0 and z['ERROR']>=0
    lb=F(z['LOWER'],S)-eps;assert lb==F(r['poly_lower'])>0
    H=lb*(1-ui)**3/(1+q*hb)**4-(ALLOW if lng else F('1e-18'))
    assert H==F(r['physical_lower']) and H>0 and r['status']=='PASS_CONDITIONAL'
    lows.append(H);secs+=r['seconds']
   assert hprev==hcap
   if uns:assert r.get('checker')=='unscaled'
   H=min(lows)
   print('cell',l,u,'b',b,'K',Kmax,'m',m,'D',D,'lower',float(H),'sec',round(secs),flush=True)
   out.append(dict(kind=cell.get('kind','tri'),c_lo=str(l),c_hi=str(u),b=b,K=Kmax,m=m,mb=cell['mb'],mt=cell['mt'],D=D,omega=round(omega(l,u),4),hcap=str(hcap),lower_exact=str(H),lower_approx=float(H),cpu_seconds=round(secs,1)))
 assert tops.get(0)==END and tops.get(4)==END,('coverage incomplete',tops)
 cm=min(tops.values()) if len(tops)==2 else START
 json.dump(dict(status='PASS_CONDITIONAL',note='Executor 2: G1-capped long-table staircase direct sign cells on [0.67,0.777]; conditional on the retained inputs of the direct method: allowance 1e-18 (K=3299 table, c<=0.56) and, for long cells, Lemma G2 allowance 1e-7+1e-14 (coeff_long_17200, c<=0.9)',c_lo=str(START),c_hi_b0=str(tops.get(0,START)),c_hi_b4=str(tops.get(4,START)),c_hi_both=str(cm),L_hi_both=-math.log1p(-float(cm)),rows=out),open(os.path.join(HERE,'ledger.json'),'w'),indent=2)
 print('PASS: long-table sha256, transfer, prefix K, G1 caps, omega<=1.08, h-tiling, pole gaps, Cauchy majorants/tails, final inequalities;',len(out),'channel cells; top',{b:str(t) for b,t in tops.items()})
if __name__=='__main__':run(sys.argv[1],sys.argv[2],sys.argv[3] if len(sys.argv)>3 else None)
