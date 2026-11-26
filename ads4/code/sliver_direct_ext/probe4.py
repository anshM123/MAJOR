import sys,mpmath as mp
from probe import S0
from fractions import Fraction as F
from prefix_improvements import tail,TIGHT_ALLOW
def kt(c):
 lo,hi=500,4000
 while hi-lo>1:
  m=(lo+hi)//2
  if tail(c,m)<TIGHT_ALLOW: hi=m
  else: lo=m
 return hi
b=int(sys.argv[1])
for cl,cu in [('.262','.266'),('.266','.27'),('.27','.274')]:
 K=kt(F(cu));hc=mp.mpf(cu)/((1001-b)*(1002+b))
 v=[float(S0(mp.mpf(c),hc*t,b,K)) for c in (cl,cu) for t in (.9,1)]
 print(cl,cu,K,' '.join('%.2e'%x for x in v),flush=True)
