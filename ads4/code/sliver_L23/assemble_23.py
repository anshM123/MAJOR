from fractions import Fraction as F
from pathlib import Path
import json
from check_norm import parameters,cap_check
from check_lowfreq import exp_u
from base_bounds import hi
floor=F('3e-6');dl=F('.0742735');du=F('.1002589')
assert dl<=1/exp_u(F('2.6')) and du>=exp_u(-F('2.3'))
# Retained LRINF contiguous beta cells [3.6,4.0], min floor 3.4606212402066455e-6.
# Certified beta range follows cosh^2(beta/2)=exp(L).
assert ((exp_u(F('1.8'))+exp_u(F('-1.8')))/2)**2 < 1/exp_u(F('-2.3'))
assert ((1/exp_u(F('-2'))+1/exp_u(F('2')))/2)**2 > exp_u(F('2.6'))
res=[]
for b in (0,4):
    path=Path('results_23_p3_budget.json')
    if not path.exists():
        if b:continue
        path=Path('partial23_p3_budget.json')
    r=next(x for x in json.load(path.open()) if x['b']==b)
    names=[Path(f'ct_23_{p}_{b}.json')for p in (3,4)]
    if not all(p.exists() for p in names):continue
    assert F(r['dlo'])==dl and F(r['dhi'])==du
    R=F(r['total_exact']);ct=F(0)
    for p,fn in zip((3,4),names):
        c=json.load(fn.open())[0];assert c['p']==p and c['b']==b and F(c['dlo'])==dl and F(c['dhi'])==du
        ct+=F(c['total'])
    E=R+ct
    row=dict(L=['2.3','2.6'],b=b,floor=str(floor),error_exact=str(E),contour_approx=float(ct),residual_budget_approx=float(R),ratio_upper=float(hi(E/floor)),lower_approx=float(floor-E),status='PASS_CONDITIONAL' if E<floor else 'FAIL')
    print(json.dumps(row));res.append(row)
json.dump(res,open('ledger_23.json','w'),indent=2)
