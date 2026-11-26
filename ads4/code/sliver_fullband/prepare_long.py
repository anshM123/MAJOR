#!/usr/bin/env python3
import sys,json,subprocess,time
from fractions import Fraction as F
from math import isqrt,comb
sys.set_int_max_str_digits(0)
K=6000;D=3000;S=10**D

def up(n,d):return -((-n)//d)
def ri(x):q=F(x);return (q.numerator*S//q.denominator,up(q.numerator*S,q.denominator))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[1],-x[0]
def sub(x,y):return add(x,neg(y))
def mul(x,a):
 a=F(a);x=x if a>=0 else x[::-1];return x[0]*a.numerator//a.denominator,up(x[1]*a.numerator,a.denominator)
def divi(x,y):
 assert 0<y[0]<=y[1];v=[F(a*S,b)for a in x for b in y];return min(v).__floor__(),max(v).__ceil__()
def mr(x):m=(x[0]+x[1])//2;return m,max(m-x[0],x[1]-m)
def mrstr(x):return ' '.join(map(str,mr(x)))
# All six low-precision seed approximations are treated as exact numbers;
# their coherent discrepancies are budgeted at the action level.
seed={k:(int(l),int(u))for k,l,u in map(str.split,open('seed_values.txt'))}
assert all(seed[k][1]-seed[k][0]<10**40 for k in ['F1','F9','F49','K0','LF9','LF49'])
seed={k:((l+u)//2)*10**(D-80)for k,(l,u) in seed.items()}
rows=[tuple(map(int,line.split()))for line in open(f'S_{K}_{D}.txt')];assert len(rows)==K+1 and all(i==j and l<=u for i,(j,l,u)in enumerate(rows))
plain=[(l,u)for _,l,u in rows]
g0,e0=[tuple(map(int,l.split()))for l in open('log_constants_3000.txt')];G=[g0];E=[e0];M=[(seed['K0'],seed['K0'])]
for j in range(K+4):
 G.append(mul(G[-1],F(2*j+1,10**6)));E.append(mul(sub(G[-2],E[-1]),F(1,10**6)));M.append(add(mul(M[-1],F(2*j+1,10**6)),mul(E[-1],2)))
y=list(map(F,json.load(open('weights.json'))['y_re']));L=[]
for j in range(K+1):
 v=(0,0)
 for i,a in enumerate(y):v=add(v,mul(add(mul(M[j+i],10**(6*i)),mul(M[j+i+1],10**(6*(i+1)))),a))
 L.append(v)
T=[(0,0)]*(K+1)
for a,coef in [(1,F(-96)),(9,F(5184,5)),(49,F(-4704,5)),(9,F(-9216,5)),(49,F(55296,5))]:
 gaussian=coef in (F(-96),F(5184,5),F(-4704,5));v=seed[('F' if gaussian else 'LF')+str(a)];q=(v,v);mom=plain if gaussian else L
 for j in range(K+1):
  T[j]=add(T[j],mul(q,coef));q=mul(sub(mom[j],mul(q,a)),F(1,4000000))
lam=isqrt(997002*S*S);T=[neg(divi(t,(lam,lam+1)))for t in T]
with open('raw_moments.in','w')as f:f.write(f'0 {K} 17203 1000000\n');f.write('\n'.join(map(mrstr,T))+'\n')

r=json.load(open('nf02_exact.json'));co=list(map(F,r['psi3_exact']));k0=F(r['k0']);J=K
mom=[ri(997002*sum(a/F(2*j+i+1)for i,a in enumerate(co))/k0)for j in range(J+1)]
with open('ideal4.in','w')as f:f.write(f'4 {J} 17201 997002\n');f.write('\n'.join(map(mrstr,mom))+'\n')
print('PASS moment inputs',K,D,flush=True)
