#!/usr/bin/env python3
"""Append an already-computed long/stair cell (all pieces PASS) to RUNDIR/chain_b.json (used when a chain wrapper was restarted).
Usage: adopt_cell.py RUNDIR prefix c_lo c_hi b mb mt D K kind   -> prints next c_lo if adopted, else c_lo"""
import sys,os,json
from fractions import Fraction as F
rd,pre,l,u,b,mb,mt,D,K,kind=sys.argv[1:11];l,u=F(l),F(u);b,mb,mt,D,K=map(int,(b,mb,mt,D,K))
low=None;ok=True;wall=0
for pi in range(mb+mt):
 p=os.path.join(rd,f'{pre}_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}_m{mb}x{mt}_p{pi}_D{D}.json')
 if not os.path.exists(p):ok=False;break
 r=json.load(open(p));ok&=r['status']=='PASS_CONDITIONAL';low=r['lower'] if low is None else min(low,r['lower']);wall+=r['seconds']
cf=os.path.join(rd,f'chain_{b}.json');done=json.load(open(cf))
if ok and not any(F(d['c_lo'])==l and F(d['c_hi'])==u for d in done):
 done.append(dict(kind=kind,c_lo=str(l),c_hi=str(u),b=b,m=mb+mt,mb=mb,mt=mt,D=D,K=K,lower=low,wall=round(wall,1)));json.dump(done,open(cf,'w'),indent=1)
print(str(u) if ok else str(l))
