#!/usr/bin/env python3
"""Check the grouped cutoff mass against the exact binary64 repaired-weight inventory.
Run with path to ads4/prl/data/k2/k2_continuum_D_1000_rep.json.
"""
from fractions import Fraction as F
from collections import defaultdict
from math import ceil
import json,sys
v=json.load(open(sys.argv[1]));S=10**80;groups=defaultdict(int)
for (c,s),x in zip(v['centres_sigma'],v['x_re']):
    c,s,x=F(float(c)),F(float(s)),F(float(x));assert 0<=c<=1029 and F('.449')<s<F('.901')
    groups[ceil(c)]+=ceil(abs(x)*S)
assert len(v['x_re'])==len(v['centres_sigma'])==3564
# sqrt(2 pi) sigma < 2.27, D=100+2 sigma^2<102,
# nu_c > 1098.49. Complete the square and use erfc(t)<=exp(-t^2).
assert 2*F(22,7)*F('.901')**2<F('2.27')**2
assert 100+2*F('.901')**2<102 and F('998.49')**2<997002
# log(2)<.694, e^-y<=2^-floor(y/.694)
from base_bounds import logu
assert logu(2)<F('.694')
B=0
for c,v in groups.items():
    n=F((109849-100*c)**2,10000*102)/F('.694');m=n.numerator//n.denominator
    B+=ceil(F(227*v,100*(1<<m)))
assert F(B,S)<F('7.1e-24')
print('PASS',len(groups),'groups; exact upper',str(F(B,S)),'decimal',float(F(B,S)))
