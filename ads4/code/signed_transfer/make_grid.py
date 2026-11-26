#!/usr/bin/env python3
"""Generate candidate (NOT certified) full-twist boxes. Standard library only."""
import argparse, json, math
from fractions import Fraction as F
from decimal import Decimal
S=10**24
def ceil(x):return -(-x.numerator//x.denominator)
def ds(n):return format(Decimal(n)/Decimal(S),'f')
def main():
 p=argparse.ArgumentParser();p.add_argument('--d-min',default='.02599');p.add_argument('--d-max',default='.44933');p.add_argument('--ratio',default='1.040816326530612244');p.add_argument('--nu-min',type=int,default=160);p.add_argument('--nu-max',type=int,default=1150);p.add_argument('--nu-width',type=int,default=20);p.add_argument('--m',type=int,default=28);p.add_argument('--K',type=int,default=200);p.add_argument('--R',type=int,default=8);p.add_argument('--h-splits',type=int,default=16);p.add_argument('--output',default='candidate_grid.json');args=p.parse_args()
 lo=int(F(args.d_min)*S);top=ceil(F(args.d_max)*S);ratio=F(args.ratio);assert 0<lo<top<S and ratio>1 and args.nu_width>0 and args.nu_min>=0 and args.nu_min<args.nu_max
 out=[];i=0
 while lo<top:
  hi=min(top,ceil(lo*ratio));assert hi>lo
  if (hi-lo)%2:hi+=1 # preserve exact midpoint/halfwidth and cover the requested range
  root=math.isqrt((S-hi)*S);hcap=ceil(F(hi*hi*S,998**2*(S+root)**2));hcap+=hcap%2
  for v in range(args.nu_min,args.nu_max,args.nu_width):
   u=min(v+args.nu_width,args.nu_max)
   out.append(dict(name=f'd{i:04d}_nu{v:04d}',d0=ds((lo+hi)//2),dd=ds((hi-lo)//2),h0=ds(hcap//2),dh=ds(hcap//2),nu0=str(Decimal(v+u)/2),dnu=str(Decimal(u-v)/2),m=args.m,K=args.K,R=args.R,hSplits=args.h_splits,fullTwist=True,maxRemainder='0.001',maxTransfer='10'))
  lo=hi;i+=1
 with open(args.output,'w') as f:json.dump(out,f,indent=2)
 print(f'{len(out)} UNVERIFIED candidate boxes ({i} d slabs). Use preflight / adaptive subdivision before a full run.')
if __name__=='__main__':main()
