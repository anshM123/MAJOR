#!/usr/bin/env python3
"""Actual repo subset, exact IEEE dyadics, directed moment transport.
Independent numerical checks are diagnostics only; the interval output is
rigorous from the recurrence and directed MPFR operations.
"""
from fractions import Fraction as F
from pathlib import Path
import subprocess,json,sys,time,mpmath as mp
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
rows=json.load(open(HERE/'actual_subset.json'))
def inp(rows,K,D):
    r=[f'{K} {D} {len(rows)}\n']
    for _,c,s,x in rows:
        a=[F(float(t)) for t in (c,s,x)];r.append(' '.join(f'{v.numerator} {v.denominator}' for v in a)+'\n')
    return ''.join(r)
def run(rows,K,D):
    t=time.perf_counter();r=subprocess.run([str(HERE/'gaussian_cut_moments')],input=inp(rows,K,D),text=True,capture_output=True,check=True);vals=[tuple(map(int,l.split())) for l in r.stdout.splitlines()];assert [j for j,_,_ in vals]==list(range(K+1));width=max(u-l for _,l,u in vals);print('K',K,'digits',D,'centers',len(rows),'elapsed',time.perf_counter()-t,'max radius numerator',str(width)[:30],'digits of width',len(str(width)),flush=True);return vals
mp.mp.dps=75
# Independent shifted-center quadrature; relevant tiny cutoff variation at the last center.
for row in [rows[0],rows[4],rows[7],rows[-1]]:
    vals=run([row],8,55);_,c,s,x=row;c,s,x=[mp.mpf(F(float(t)).numerator)/F(float(t)).denominator for t in (c,s,x)];L=mp.sqrt(997002)+100
    for j,l,u in vals:
        f=lambda v:mp.exp(-v*v/2)*((c+s*v)/1000)**(2*j)*(1-(mp.erfc((L-c-s*v)/10)+mp.erfc((L+c+s*v)/10))/2)
        q=mp.quad(f,[-30,-12,0,12,30])*s*x/(2 if c==0 else 1);q*=mp.mpf(10)**55
        # quadrature not certified; compare to interval within 1e-10 scaled units
        assert l-mp.mpf('1e-10')<=q<=u+mp.mpf('1e-10'),(row[0],j,str(q-l),str(q-u))
print('36 independent high-precision quadrature diagnostics agree')
vals=run(rows,3253,2500);assert max(u-l for _,l,u in vals)<10**100
# actual sample spectrum only, not full repaired weight.
print('PASS 12 actual repository centers through 6506 powers; no full-weight claim')
