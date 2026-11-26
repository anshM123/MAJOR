#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
import json,sys,struct,math
sys.set_int_max_str_digits(0)
# Verify exact inventory transfer and numerical seed enclosure widths.
w=json.load(open('weights.json'));assert len(w['x_re'])==len(w['centres_sigma'])==3564
h=2166136261
for cs,x in zip(w['centres_sigma'],w['x_re']):
 for a in struct.pack('<ddd',*cs,x):h=((h^a)*16777619)&0xffffffff
assert h==1930466845
ss={k:(int(l),int(u))for k,l,u in map(str.split,open('seed_values.txt'))};keys=('F1','F9','F49','K0','LF9','LF49')
assert all(0<=ss[k][1]-ss[k][0]<10**40 for k in keys)
# Independently verify the Chebyshev transform's norm-one property.
for n in (369,241,371,242,433,325,435,326):
 for k in range(n+1):
  assert sum((1 if j==0 else 2)*math.comb(k,(k-j)//2)for j in range(k,-1,-2))==2**k
# Reassemble exact moment intervals from shards.
K,D=1254,1600
rows=[list(map(lambda l:tuple(map(int,l.split())),open(f'g{K}_{D}_{i}.txt')))for i in range(4)]
for j,line in enumerate(open(f'S_{K}_{D}.txt')):
 k,l,u=map(int,line.split());assert k==j and l==sum(a[j][1]for a in rows) and u==sum(a[j][2]for a in rows) and l<=u
# Check the computed coefficient-vector radii, as distinct from the coherent seed error.
f=open('direct_coefficients.txt');N,D=map(int,f.readline().split());assert (N,D)==(1252,1600);rad=0
for k,line in enumerate(f):a,e,g,t=map(int,line.split());assert e>=0 and t>=0;rad=max(rad,e,t)
assert k==N and 2*rad<10**(D-375)
# Correlated pole seed numerical error budget and cutoff bound.
rr=json.load(open('pole_bound_results.json'));rrfull=json.load(open('pole_bound_full_results.json'));coef={1:F(-96),9:F(5184,5),49:F(-4704,5)}
for rrs,budget in [(rr,F('2.4e-20')),(rrfull,F('1.8e-15'))]:
 for b in (0,4):
  e=sum(abs(coef[r['a']])*F(r['action_bound'])*(F('7.1e-24')/r['a']+F('1e-40'))/998 for r in rrs if r['b']==b)
  assert e<budget
# verify exact sign ledger and contiguous c coverage for both channels
cells=[(F('.05'),F('.065')),(F('.065'),F('.09'))];out=[]
for l,u in cells:
 for b in (0,4):
  name=f'direct_{l.numerator}_{l.denominator}_{u.numerator}_{u.denominator}_{b}';r=json.load(open(name+'.json'));z={a:int(v)for a,v in map(str.split,open(name+'.out'))}
  assert r['status']=='PASS_CONDITIONAL' and (F(r['c_lo']),F(r['c_hi']),r['b'])==(l,u,b)
  assert z['LOWER']==z['CENTER']-z['L1']-z['ERROR']
  h=F(r['hcap']);assert h>=u/((1001-b)*(1002+b));K=r['K'];assert 1-h*(K-1)*(K+2*b)>0
  c=u/F(1999,2000);assert (c.numerator**(K+1)<<2700)<=c.denominator**(K+1)
  lb=F(z['LOWER'],10**1100)-F(r['analytic_tail']);assert lb==F(r['poly_lower'])>0
  H=lb*(1-u)**3/(1+b*(b+1)*h)**4-F('1e-18');assert H==F(r['physical_lower']) and H>0
  out.append(dict(c_lo=str(l),c_hi=str(u),b=b,lower=str(H),approx=float(H)))
print('PASS: full inventory, six seed widths, exact moment sum, coefficient radii, Chebyshev transform, seed budget, four direct sign ledgers')
for r in out:print(r['c_lo'],r['c_hi'],r['b'],r['approx'])
json.dump(dict(status='PASS_CONDITIONAL',c_lo='.05',c_hi='.09',rows=out,prior_high_L_lower='2.3',open_L_lower_approx=-math.log(.91),open_L_upper='2.3'),open('direct_ledger.json','w'),indent=2)
