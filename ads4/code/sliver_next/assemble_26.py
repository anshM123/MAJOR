from fractions import Fraction as F
import json
from check_norm import parameters
from check_lowfreq import exp_u
from base_bounds import hi
r=json.load(open('results_2.6_2.9.json'));assert len(r)==12 and all(x['status']=='ENCLOSED' for x in r)
d,dd,*_=parameters(r[0]['config']);dl,du=d-dd,d+dd
assert dl<=1/exp_u(F('2.9')) and du>=exp_u(-F('2.6'))
assert all(parameters(x['config'])[:2]==[d,dd] for x in r)
# floor from 16 retained contiguous LRINF beta cells [3.9,4.3]. Verify actual beta span via cosh^2(beta/2)=exp(L).
# If beta=3.9, cosh^2(beta/2)<exp(2.6); beta=4.3 gives cosh^2(beta/2)>exp(2.9).
assert ((exp_u(F('1.95'))+exp_u(F('-1.95')))/2)**2 < 1/exp_u(F('-2.6'))
assert ((1/exp_u(F('-2.15'))+1/exp_u(F('2.15')))/2)**2 > exp_u(F('2.9'))
rows=[];floor=F('2.5e-6')
for r in json.load(open('p3_2.6_2.9_budget.json')):
 b=r['b'];assert F(r['dlo'])==dl and F(r['dhi'])==du
 E=F(r['total_exact']);ct=F(0)
 for p,fn in [(3,f'ct_2.6_2.9_3_{b}.json'),(4,f'ct_2.6_2.9_4_{b}_high.json')]:
  c=json.load(open(fn))[0];assert c['p']==p and c['b']==b and F(c['dlo'])==dl and F(c['dhi'])==du
  ct+=F(c['total'])
 E+=ct
 row=dict(L=['2.6','2.9'],b=b,floor=str(floor),error_exact=str(E),contour_approx=float(ct),ratio_upper=float(hi(E/floor)),lower_approx=float(floor-E),status='PASS_CONDITIONAL' if E<floor else 'FAIL')
 rows.append(row);print(json.dumps(row))
json.dump(rows,open('ledger_26.json','w'),indent=2)
