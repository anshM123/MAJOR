#!/usr/bin/env python3
"""Executor-2 runner: G1-capped long-table staircase cells (Lemma G2 table coeff_long_17200.txt) on c in [0.67, 0.777], one b,
sequential (one C++ process at a time per invocation).  Mirrors code/gap_prove/run_chain3.py (same mb rule: minimal mb with
complex-pole admissibility for ah=1 on every piece; same K rule k_tail; same caps/staircase from gap_prove), except:
 * the cell list is fixed (CELLS json: list of [c_lo, c_hi] strings), taken from code/gap_design/plan_highc_056_0777.json;
 * a FAIL_INCONCLUSIVE cell is split in two halves in c (on the 1/4000 grid) and both halves are rerun;
 * mb (bottom h-pieces): among the pole-admissible mb in [mb_min, mb_min+3], the one minimising the predicted cost
   sum_pieces K*terms (driver preflight with PREP_ONLY=1; only the choice of pieces depends on it, not the rigour);
 * checker: mode 'u' = unscaled (gap_exec2/gap_stair_long_u.py + direct_sign_triangle_u), D = 100 + growth of the
   C++ l1 recurrence factor (|a|+dc+2r)/(d-r) (+10); mode 's' = the prover's scaled driver gap_prove/gap_stair_long.py
   with the prover's D rule.
Usage: run_exec2.py b CELLS RUNDIR mode mt
Appends passing cells to RUNDIR/chain_b.json (kind 'longu' or 'long').
"""
import sys,os,json,math,subprocess,time
from fractions import Fraction as F
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'..','gap_prove'))
from gap_triangle import g1_cap
from gap_stair import stair_c,hgrid
from gap_stair_long import k_tail as need_k
b=int(sys.argv[1]);cells=[tuple(map(F,x)) for x in json.load(open(sys.argv[2]))];rd=sys.argv[3];mode=sys.argv[4];mtop=int(sys.argv[5])
drv=os.path.join(HERE,'gap_stair_long_u.py') if mode=='u' else os.path.join(HERE,'..','gap_prove','gap_stair_long.py')
pref='g2su' if mode=='u' else 'g2st'
cf=os.path.join(rd,f'chain_{b}.json');done=json.load(open(cf)) if os.path.exists(cf) else []
s=2*b+1
def plan(c,u,mb=None):
 mt=mtop
 if mb is not None:return mb,mt,Dof(c,u,mb,mt)
 mb=2
 while True:
  hl=hgrid(c,u,b,mb,mt);okp=True
  for ha,hb in zip(hl,hl[1:]):
   Ki=need_k(stair_c(ha,u));pi_=(Ki-1)*(Ki-1+s)
   if 1-(hb+(hb-ha)/2)*pi_<=0:okp=False
  if okp:break
  mb+=1
 return mb,mt,Dof(c,u,mb,mt)
def Dof(c,u,mb,mt):
 hl=hgrid(c,u,b,mb,mt);mx=0.0
 for ha,hb in zip(hl,hl[1:]):
  h0=(ha+hb)/2;dh=(hb-ha)/2;lg=0.0;ui=stair_c(ha,u);Ki=need_k(ui)
  for k in range(Ki):
   pk=k*(k+s)
   if mode=='u':g=(abs((c+ui)/2-pk*h0)+(ui-c)/2+2*pk*dh)/(1-pk*h0-pk*dh)
   else:g=2*(abs((c+ui)/2-pk*h0)+(ui-c)/2+pk*dh)/(1-pk*(h0+dh))
   lg+=math.log10(float(g));mx=max(mx,lg)
 return 100+(math.ceil(mx)+10 if mx>0 else 0)
queue=list(cells)
G=F(1,4000)
while queue:
 c,u=queue.pop(0)
 if any(F(r['c_lo'])==c and F(r['c_hi'])==u for r in done):continue
 mb0,mt,D=plan(c,u)
 if mode=='u':
  best=None;pd=os.path.join(rd,f'prep_{b}');os.makedirs(pd,exist_ok=True)
  for mbx in range(mb0,mb0+4):
   _,_,Dx=plan(c,u,mbx)
   o=subprocess.run(['python3',drv,str(c),str(u),str(b),str(mbx),str(mt),str(Dx)],cwd=pd,env=dict(os.environ,PREP_ONLY='1'),capture_output=True,text=True).stdout
   for fn in os.listdir(pd):os.remove(os.path.join(pd,fn))
   pf=[ln.split() for ln in o.splitlines() if ln.startswith('PREFLIGHT')]
   cost=sum(int(t[t.index('K')+1])*int(t[t.index('terms')+1]) for t in pf) if len(pf)==mbx+mt else None
   print('  prep mb',mbx,'D',Dx,'cost K*terms',cost,'terms',[int(t[t.index('terms')+1]) for t in pf],flush=True)
   if cost is not None and (best is None or cost<best[0]):best=(cost,mbx,Dx)
  _,mb,D=best
 else:mb=mb0
 m=mb+mt;K=need_k(u)
 st=time.time()
 subprocess.run(['python3',drv,str(c),str(u),str(b),str(mb),str(mt),str(D)],cwd=rd,stdout=open(os.path.join(rd,(f'x2cell_{mode}_{c}_{u}_{b}_m{mb}x{mt}_D{D}.log').replace('/','o')),'w'),stderr=subprocess.STDOUT)
 ok=True;low=None;cpu=0.0
 for pi in range(m):
  name=f'{pref}_{c.numerator}_{c.denominator}_{u.numerator}_{u.denominator}_{b}_m{mb}x{mt}_p{pi}_D{D}.json'
  try:
   r=json.load(open(os.path.join(rd,name)));ok&=r['status']=='PASS_CONDITIONAL';low=r['lower'] if low is None else min(low,r['lower']);cpu+=r['seconds']
  except Exception:ok=False
 print(time.strftime('%H:%M:%S'),'b',b,'mode',mode,'cell',str(c),str(u),float(c),float(u),'K',K,'mb,mt',mb,mt,'D',D,'PASS' if ok else 'FAIL','lower',low,'wall',round(time.time()-st),flush=True)
 if ok:
  done.append(dict(kind='longu' if mode=='u' else 'long',c_lo=str(c),c_hi=str(u),b=b,m=m,mb=mb,mt=mt,D=D,K=K,lower=low,wall=round(time.time()-st,1),cpu=round(cpu,1)));json.dump(done,open(cf,'w'),indent=1)
 else:
  mid=F(math.floor((c+u)/2/G))*G
  if not (c<mid<u):print('STOP: cannot split further',str(c),str(u),flush=True);break
  print('SPLIT',str(c),str(mid),str(u),flush=True)
  queue[:0]=[(c,mid),(mid,u)]
