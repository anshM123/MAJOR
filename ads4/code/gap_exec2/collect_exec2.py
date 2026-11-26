#!/usr/bin/env python3
"""Collect passing cells from RUNDIR/chain_{0,4}.json into code/gap_exec2/cells.json (input of verify_exec2_cells.py).
Same format as code/gap_prove/collect_cells.py; cells are sorted by c_lo, and cells superseded by a split are not present
(the runner only records PASS cells)."""
import json,sys,os
from fractions import Fraction as F
rd=sys.argv[1];HERE=os.path.dirname(os.path.abspath(__file__));out=[]
for b in (0,4):
 rows=sorted(json.load(open(os.path.join(rd,f'chain_{b}.json'))),key=lambda r:F(r['c_lo']))
 out+=[dict(kind=r['kind'],c_lo=r['c_lo'],c_hi=r['c_hi'],b=b,m=r['m'],mb=r['mb'],mt=r['mt'],D=r['D'],K=r['K'],lower=r['lower'],wall=r['wall'],cpu=r.get('cpu')) for r in rows]
json.dump(out,open(os.path.join(HERE,'cells.json'),'w'),indent=1);print(len(out),'cells')
