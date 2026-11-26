#!/usr/bin/env python3
"""Adaptive chain of LONG-TABLE staircase cells (gap_stair_long.py, c>=0.56, Lemma G2) for one b (sequential; one process at a time).
Usage: run_chain2.py b c_start c_stop w0 wmax RUNDIR mt
On PASS the width grows by 1.25x (capped at wmax, endpoints on the 1/2000 grid); on FAIL it halves (min 1/2000) and wmax drops to 0.8x the failed width.
Choose m (h-pieces) minimal with complex-pole admissibility for ah=1, and D from the a-priori rounding growth.
Appends passing cells to RUNDIR/chain_b.json.
"""
import sys,os,json,math,subprocess,time
from fractions import Fraction as F
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from gap_triangle import need_k,g1_cap
from gap_stair import stair_c,hgrid
from gap_stair_long import k_tail as need_k
b=int(sys.argv[1]);c=F(sys.argv[2]);stop=F(sys.argv[3]);w=F(sys.argv[4]);wmax=F(sys.argv[5]);rd=sys.argv[6]
G=F(1,2000)
def snap(x):return F(math.floor(x/G))*G
drv=os.path.join(os.path.dirname(os.path.abspath(__file__)),'gap_stair_long.py');mtop=int(sys.argv[7])
cf=os.path.join(rd,f'chain_{b}.json');done=json.load(open(cf)) if os.path.exists(cf) else []  # shared with run_chain.py (never run both for one b); entries carry kind
s=2*b+1
while c<stop:
 u=min(stop,c+max(G,snap(w)));K=need_k(u);hcap=g1_cap(c,u,b)[0];p=(K-1)*(K-1+s)
 mt=mtop;mb=2  # mb=1 (single bottom piece) is slower: nearer the complex pole, larger majorant
 while True:
  hl=hgrid(c,u,b,mb,mt);okp=True
  for ha,hb in zip(hl,hl[1:]):
   Ki=need_k(stair_c(ha,u));pi_=(Ki-1)*(Ki-1+s)
   if 1-(hb+(hb-ha)/2)*pi_<=0:okp=False
  if okp:break
  mb+=1
 m=mb+mt
 # error injected at Horner step k is amplified by prod_{j<k} ratio_j, ratio_j=2(|c0-p h0|+dc+p dh)/(1-p(h0+dh)) (C++ l1 recurrence)
 mx=0.0
 for ha,hb in zip(hl,hl[1:]):
  h0=(ha+hb)/2;dh=(hb-ha)/2;lg=0.0;ui=stair_c(ha,u);Ki=need_k(ui)
  for k in range(Ki):
   pk=k*(k+s);lg+=math.log10(float(2*(abs((c+ui)/2-pk*h0)+(ui-c)/2+pk*dh)/(1-pk*(h0+dh))));mx=max(mx,lg)
 D=100+(math.ceil(mx)+10 if mx>0 else 0)
 st=time.time()
 subprocess.run(['python3',drv,str(c),str(u),str(b),str(mb),str(mt),str(D)],cwd=rd,stdout=open(os.path.join(rd,(f'lgcell_{c}_{u}_{b}_m{mb}x{mt}_D{D}.log').replace('/','o')),'w'),stderr=subprocess.STDOUT)
 ok=True;low=None
 for pi in range(m):
  name=f'g2st_{c.numerator}_{c.denominator}_{u.numerator}_{u.denominator}_{b}_m{mb}x{mt}_p{pi}_D{D}.json'
  try:
   r=json.load(open(os.path.join(rd,name)));ok&=r['status']=='PASS_CONDITIONAL';low=r['lower'] if low is None else min(low,r['lower'])
  except Exception:ok=False
 print(time.strftime('%H:%M:%S'),'b',b,'cell',float(c),float(u),'K',K,'mb,mt',mb,mt,'D',D,'PASS' if ok else 'FAIL','lower',low,'wall',round(time.time()-st),flush=True)
 if ok:
  done.append(dict(kind='long',c_lo=str(c),c_hi=str(u),b=b,m=m,mb=mb,mt=mt,D=D,K=K,lower=low,wall=round(time.time()-st,1)));json.dump(done,open(cf,'w'),indent=1)
  c=u;w=min(wmax,w*F(5,4))
 else:
  if w<=G:print('STOP: fails at minimal width',flush=True);break
  wmax=min(wmax,snap(w*F(4,5)));w=max(G,w/2)  # never retry the failed width
