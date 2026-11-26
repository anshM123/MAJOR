"""Fail-closed exact union of the previously rechecked and newly rechecked ledgers."""
from fractions import Fraction as F
from pathlib import Path
from math import factorial,log1p
import json
from verify_merged import union
ROOT=Path(__file__).resolve().parent
old=json.load(open(ROOT/'verified_previous.json'))['rows']
new=json.load(open(ROOT/'verified_new.json'))['new_rows']
rows=old+new
U=union(rows);print('union',U)
# No float comparison participates in coverage or positivity.
# This test covers both direct-action L and L_phys=L+log(1+qh).
t=F(3,10)-F(20,998**2)
assert t>0 and sum(t**j/factorial(j) for j in range(4))>1/(1-F('.258'))
t=F(23,10);m=30;partial=sum(t**j/factorial(j) for j in range(m+1));tail=t**(m+1)/factorial(m+1)/(1-t/(m+2));assert partial+tail<10
expect={'0':[['129/500','9/10']],'4':[['129/500','9/10']]}
if U!=expect:
    print('INCOMPLETE: no full-coverage certificate.');raise SystemExit(1)
assert all(F(r['lower_exact'])>0 for r in rows)
mins={str(b):min(F(r['lower_exact']) for r in rows if r['b']==b) for b in (0,4)}
# The common simple rounded floor is selected downwards using exact arithmetic.
floor=min(mins.values());floor_dec=F(floor*10**8//1,10**8);assert 0<floor_dec<=floor
new_count=sum(json.load(open(p))['count'] for p in list(ROOT.glob('*_ref_ledger.json'))+list(ROOT.glob('*_group_ledger.json')))
result=dict(status='FULL_SLIVER_PASS_CONDITIONAL',scope='Delta_g=1000, a>=998, J>=1001, b=0,4; retained angular identity and spectral inputs',c_union=U,L_target=['0.3','2.3'],L_direct_endpoints_approx=[-log1p(-.258),-log1p(-.9)],floor_exact_by_channel={b:str(v) for b,v in mins.items()},common_floor=str(floor_dec),new_scalar_nodes=new_count,rows=rows)
json.dump(result,open(ROOT/'FULL_COVERAGE_LEDGER.json','w'),indent=2)
print('FULL SLIVER PASS_CONDITIONAL, common floor',float(floor_dec),'new nodes',new_count,'cells',len(rows)//2)
