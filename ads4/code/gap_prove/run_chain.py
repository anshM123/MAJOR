#!/usr/bin/env python3
"""Adaptive chain of G1-capped cells for one b (sequential; one process at a time).
Usage: run_chain.py b c_start c_stop w0 wmax RUNDIR
On PASS the width grows by 1.5x (capped at wmax, endpoints on the 1/2000 grid); on FAIL it halves (min 1/2000).
Choose m (h-pieces) minimal with complex-pole admissibility for ah=1, and D from the a-priori rounding growth.
Appends passing cells to RUNDIR/chain_b.json.
"""
import sys,os,json,math,subprocess,time
from fractions import Fraction as F
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from gap_triangle import need_k,g1_cap
b=int(sys.argv[1]);c=F(sys.argv[2]);stop=F(sys.argv[3]);w=F(sys.argv[4]);wmax=F(sys.argv[5]);rd=sys.argv[6]
G=F(1,2000)
def snap(x):return F(math.floor(x/G))*G
drv=os.path.join(os.path.dirname(os.path.abspath(__file__)),'gap_triangle.py')
cf=os.path.join(rd,f'chain_{b}.json');done=json.load(open(cf)) if os.path.exists(cf) else []
s=2*b+1
while c<stop:
 u=min(stop,c+max(G,snap(w)));K=need_k(u);hcap=g1_cap(c,u,b)[0];p=(K-1)*(K-1+s)
 m=1
 while 1-hcap*(1+F(1,2*m))*p<=0:m+=1
 # a-priori rounding growth of the Horner l1 error (ratio of the C++ recurrence, real point bound)
 lg=0.0
 for k in range(K):
  pk=k*(k+s);r=float(2*(u+pk*hcap)/(1-pk*hcap))
  if r>1:lg+=math.log10(r)
 D=100+max(0,math.ceil(lg)+10)
 st=time.time()
 subprocess.run(['python3',drv,str(c),str(u),str(b),str(m),str(D)],cwd=rd,stdout=open(os.path.join(rd,f'cell_{c}_{u}_{b}_m{m}_D{D}.log'.replace('/','o')),'w'),stderr=subprocess.STDOUT)
 ok=True;low=None
 for pi in range(m):
  name=f'g1tri_{c.numerator}_{c.denominator}_{u.numerator}_{u.denominator}_{b}_m{m}_p{pi}_D{D}.json'
  try:
   r=json.load(open(os.path.join(rd,name)));ok&=r['status']=='PASS_CONDITIONAL';low=r['lower'] if low is None else min(low,r['lower'])
  except Exception:ok=False
 print(time.strftime('%H:%M:%S'),'b',b,'cell',float(c),float(u),'K',K,'m',m,'D',D,'PASS' if ok else 'FAIL','lower',low,'wall',round(time.time()-st),flush=True)
 if ok:
  done.append(dict(c_lo=str(c),c_hi=str(u),b=b,m=m,D=D,K=K,lower=low,wall=round(time.time()-st,1)));json.dump(done,open(cf,'w'),indent=1)
  c=u;w=min(wmax,w*F(3,2))
 else:
  if w<=G:print('STOP: fails at minimal width',flush=True);break
  w=max(G,w/2)
