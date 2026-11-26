# 2D Chebyshev coefficients (sampling, non-rigorous) of the unnormalised prefix S on a termination-
# preserving box c in [cl,cu], x in [0,xhi] (default xhi=1/(K(K+s))). Reports centre, L1 of the rest,
# and the L1 tail outside degree (mc,mx) to estimate certificate degrees.
import sys,time
from common2 import *
cl,cu=float(sys.argv[1]),float(sys.argv[2]);b=int(sys.argv[3]);K=int(sys.argv[4])
nc=int(sys.argv[5]);nx=int(sys.argv[6]);s=2*b+1
xhi=float(sys.argv[7]) if len(sys.argv)>7 else 1/(K*(K+s))
pi=arb.pi();t0=time.time()
def nodes(n,lo,hi):return [arb(lo)+(arb(hi)-arb(lo))*(1+((2*j+1)*pi/(2*n)).cos())/2 for j in range(n)]
C=nodes(nc,cl,cu);X=nodes(nx,0,xhi)
V=[[Sx(c,x,b,K) for x in X] for c in C]
def dct(vals,n):
    out=[]
    for m in range(n):
        s=sum((vals[j]*((m*(2*j+1))*pi/(2*n)).cos() for j in range(n)),arb(0))*2/n
        out.append(s/2 if m==0 else s)
    return out
T1=[dct(row,nx) for row in V]
T=[dct([T1[i][j] for i in range(nc)],nc) for j in range(nx)]  # T[j][i]
co=[[float(T[j][i].mid()) for j in range(nx)] for i in range(nc)]
L1=sum(abs(co[i][j]) for i in range(nc) for j in range(nx))-abs(co[0][0])
print('box c[%g,%g] b %d K %d xhi %.4e (y_hi at cu %.3f, nmin(cu) %d) time %.0fs'%(cl,cu,b,K,xhi,xhi/xphys(cu),nmin(cu,b),time.time()-t0))
print(' centre %.4e  L1rest %.4e  ratio %.3f'%(co[0][0],L1,L1/co[0][0]))
print(' |c10| %.3e |c01| %.3e |c11| %.3e |c20| %.3e |c02| %.3e'%(abs(co[1][0]),abs(co[0][1]),abs(co[1][1]),abs(co[2][0]),abs(co[0][2])))
for m in (4,8,16,24,32,48):
    if m<nc:print(' c-degree tail >=%d: %.3e'%(m,sum(abs(co[i][j]) for i in range(m,nc) for j in range(nx))))
for m in (2,4,8,16):
    if m<nx:print(' x-degree tail >=%d: %.3e'%(m,sum(abs(co[i][j]) for i in range(nc) for j in range(m,nx))))
