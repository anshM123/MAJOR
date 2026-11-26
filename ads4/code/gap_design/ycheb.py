# Chebyshev coefficients of the normalised prefix in y=h/hphys(c) on [0,Y] at fixed c (sampling at Chebyshev nodes).
import sys,math,time
from common import *
c=float(sys.argv[1]);b=int(sys.argv[2]);n=int(sys.argv[3]);Y=float(sys.argv[4]) if len(sys.argv)>4 else 1.0
K=need_k(c);hp=hphys(c,b)
pi=arb.pi()
nodes=[(arb(1)+((2*j+1)*pi/(2*n)).cos())/2*Y for j in range(n)]
vals=[S0(arb(c),arb(hp)*t,b,K,norm=False) for t in nodes]
co=[]
for m in range(n):
    s=sum((vals[j]*((m*(2*j+1))*pi/(2*n)).cos() for j in range(n)),arb(0))*2/n
    if m==0:s/=2
    co.append(s)
print('c',c,'b',b,'Y',Y,'K',K)
for m in range(n):
    if m<6 or m%25==0 or m==n-1:print(m,'%.3e'%float(co[m].mid()))
