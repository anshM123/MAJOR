#!/usr/bin/env python3
"""Regression tests for signed integer transport and exact truncated spectral moments."""
from fractions import Fraction as F
from math import comb,log10
import sys,subprocess,json
sys.set_int_max_str_digits(0)
# Deliberately include small spectral degree, so the degree cutoff is exercised.
for kind in (0,4):
 for M in (1,3,9,13):
  V=13;K=29;D=140;S=10**D
  # Synthetic signed measures with rational points and nontrivial truncation.
  atoms=[(F(1,7),F(3,5)),(F(5,11),F(-1,6)),(F(2,3),F(9,10))]
  T=[sum(w*x**j for x,w in atoms)for j in range(M+1)]+[F(0)]
  vals=[(int(v*S),1)for v in T[:-1]]
  inp=f'{kind} {M} {K} {V}\n'+'\n'.join(f'{v} {e}'for v,e in vals)+'\n'
  out=subprocess.run(['./long_transform'],input=inp,text=True,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,check=True).stdout
  rows={int(j):(int(a),int(e))for j,a,e in map(str.split,out.splitlines())}
  I=[[F(0)]*(M+2)for _ in range(K+1)]
  P=[T[0]]
  for m in range(1,K+1):
   for j in range(M+1):I[m][j]=(m*T[j]-sum((m-l)*(V*I[l][j+1]+I[l][j]/4)for l in range(1,m)))/m**2
   P.append(-V*I[m][1]-I[m][0]/4)
  if kind==0:expected={m:I[m][0]for m in range(1,K+1)}
  else:expected={m:sum(P[l]*comb(m-l+3,3)for l in range(m+1))for m in range(K+1)}
  assert rows.keys()==expected.keys()
  for m,v in expected.items():
   a,e=rows[m];assert e>=0 and abs(v*S-a)<=e,(kind,M,m)
print('PASS 236 exact transport checks, including both channels and active spectral truncation')
# Exact independent bounds of the ideal weight, using the Bernstein basis.
r=json.load(open('nf02_exact.json'));k0=F(r['k0']);assert k0>0
for key,mult in [('psi2_exact',2),('psi3_exact',997002)]:
 a=list(map(F,r[key]));n=len(a)-1;be=[sum(a[j]*F(comb(i,j),comb(n,j))for j in range(i+1))for i in range(n+1)]
 bound=mult*sum(map(abs,be))/((n+1)*k0);assert bound<60000000
 print('PASS ideal L1 bound',key,float(bound))
# Check all previously enclosed ideal coefficients overlap the new vector.
path=sys.argv[1]if len(sys.argv)>1 else 'coeff_long_1800.txt'
with open('old_coefficients.txt')as a,open(path)as b:
 n,d=map(int,next(a).split());m,e=map(int,next(b).split());assert m>=n and e>=d;fac=10**(e-d);count=0;diff=0
 for i,(x,y) in enumerate(zip(a,b)):
  A,ar,G,gr=map(int,x.split());B,br,H,hr=map(int,y.split());assert abs(G*fac-H)<=gr*fac+hr
  diff=max(diff,abs(A*fac-B));count+=1
 assert count==n+1 and diff<10**(e-107)
 print('PASS old/new ideal enclosure overlaps',count,'raw max absolute center difference <1e-107')
