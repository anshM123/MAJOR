from fractions import Fraction as F
from pathlib import Path
import subprocess,json,sys,time
sys.set_int_max_str_digits(0)
rows=json.load(open('actual_subset.json'));rows=[rows[i] for i in [0,3,5,8,9,10,11]]
K=3253;D=2500;S=10**D
lines=[f'{K} {D} {len(rows)}\n']
for _,c,s,x in rows:lines.append(' '.join(f'{v.numerator} {v.denominator}' for v in [F(float(c)),F(float(s)),F(float(x))])+'\n')
t=time.perf_counter();proc=subprocess.run(['./gaussian_cut_moments_fast'],input=''.join(lines),text=True,stdout=open('gauss_high_vals.txt','w'),check=True)
vals=[tuple(map(int,l.split())) for l in open('gauss_high_vals.txt')];assert [x[0]for x in vals]==list(range(K+1));w=max(u-l for _,l,u in vals);assert w<10**100
print('PASS actual repository subset of',len(rows),'centers',K+1,'even moments, digits',D,'elapsed',time.perf_counter()-t,'maximum integer width digits',len(str(w)),flush=True)
