from fractions import Fraction as F
from math import comb
import sys,os
sys.set_int_max_str_digits(0)
N=int(sys.argv[1]);D=3000;S=10**D
# coefficient enclosure arithmetic at fixed decimal scale.
def up(n,d):return -((-n)//d)
def mul(x,a):
 a=F(a);x=x if a>=0 else x[::-1];return x[0]*a.numerator//a.denominator,up(x[1]*a.numerator,a.denominator)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def mr(x):m=(x[0]+x[1])//2;return m,max(m-x[0],x[1]-m)
fs=[(0,0)]
with open('raw_long.txt')as f:
 for i in range(1,N+4):
  j,v,e=map(int,f.readline().split());assert j==i and e>=0;fs.append((v-e,v+e))
Z=[sub(add(mul(fs[k+2],12),mul(fs[k+3],k+3)),mul(fs[k+1],k+1))for k in range(N+1)]
A=[mul(sub(Z[k],mul(Z[k-1] if k else (0,0),F(k,k+3))),F(comb(k+3,3)**2,k+3))for k in range(N+1)]
G=[]
with open('ideal_long_corrected.txt')as f:
 for i in range(N+1):
  j,v,e=map(int,f.readline().split());assert j==i and e>=0;G.append((v-e,v+e))
assert F('20.771893')<F(A[0][0],S)<F(A[0][1],S)<F('20.771894')
with open(f'coeff_long_{N}.txt.tmp','w')as f:
 f.write(f'{N} {D}\n')
 for a,g in zip(A,G):f.write(' '.join(map(str,(*mr(a),*mr(g))))+'\n')
os.replace(f'coeff_long_{N}.txt.tmp',f'coeff_long_{N}.txt')
import math
print('PASS vector',N+1,'log10 max abs enclosure width',math.log10(max(u-l for l,u in A+G))-D,flush=True)
