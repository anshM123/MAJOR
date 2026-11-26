# float probe, unchanged cap1 rectangle (as in direct_sign_triangle.py), K=need_k(c_hi)
import sys,mpmath as mp
from probe import S0
from fractions import Fraction as F
R=F(1999,2000)
def need_k(c):
 t=c/R;p,q=t.numerator,t.denominator
 def ok(k):return (p**(k+1)<<2700)<=q**(k+1)
 lo,hi=-1,1
 while not ok(hi):hi*=2
 while hi-lo>1:
  m=(lo+hi)//2
  if ok(m):hi=m
  else:lo=m
 return hi
b=int(sys.argv[1])
for cl,cu in [('.258','.262'),('.258','.265'),('.262','.266'),('.266','.27'),('.27','.274'),('.274','.278'),('.278','.282')]:
 K=need_k(F(cu));hc=mp.mpf(cu)/((1001-b)*(1002+b))
 v=[float(S0(mp.mpf(c),hc*t,b,K)) for c in (cl,cu) for t in (.9,1)]
 print(cl,cu,K,' '.join('%.2e'%x for x in v),flush=True)
