#!/usr/bin/env python3
"""Exact pole-free finite reduction on all remaining c; not a sign certificate."""
from fractions import Fraction as F
from base_bounds import sqrtlo
from tight_spectral_tail import tail
from math import ceil
import json
start=F('.268');end=F('.9');last=start;out=[]
old=json.load(open('inherited_prefix_grid.json'))['rows']
for row in old:
 l=max(start,F(row['c_lo']));u=F(row['c_hi'])
 if u<=start:continue
 assert l==last and l<u
 # The tail bound is checked at the selected integer, not inferred by floating point.
 # Choose via a binary search; correctness only uses the final exact inequality.
 a,b=0,min(17199,row['K'])
 assert tail(u,b)<F(1,2**110)
 while b-a>1:
  k=(a+b)//2
  if tail(u,k)<F(1,2**110):b=k
  else:a=k
 K=b;assert tail(u,K)<F(1,2**110) and K+1<=17200
 h=F(ceil((1-sqrtlo(l))**2*10**15/998**2),10**15)
 gap=1-h*(K-1)*(K+8);assert gap>0
 out.append(dict(c_lo=str(l),c_hi=str(u),K=K,hcap=str(h),pole_gap=str(gap),pole_gap_decimal=float(gap),sign='NOT_RUN'))
 last=u
assert last==end
print('PASS finite reduction of whole remaining band:',len(out),'contiguous rectangles, Kmax',max(x['K']for x in out),'min pole gap',min(x['pole_gap_decimal']for x in out))
json.dump(dict(status='PASS_PREFIX_AND_POLE_GAPS_ONLY',c_lo=str(start),c_hi=str(end),rows=out),open('full_grid_results.json','w'),indent=2)
