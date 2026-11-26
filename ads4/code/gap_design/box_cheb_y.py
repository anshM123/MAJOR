# 2D Chebyshev coefficients (sampling, non-rigorous) of the unnormalised prefix on the curvilinear
# physical cell c in [cl,cu], h = y*(1-sqrt c)^2/998^2, y in [0,Y]. Predicts PASS/L1 of a tensor certificate.
import sys,time
from common2 import *
cl,cu=float(sys.argv[1]),float(sys.argv[2]);b=int(sys.argv[3])
nc=int(sys.argv[4]);ny=int(sys.argv[5]);Y=float(sys.argv[6]) if len(sys.argv)>6 else 1.0
K=int(sys.argv[7]) if len(sys.argv)>7 else needK(cu)
pi=arb.pi();t0=time.time()
def nodes(n,lo,hi):return [arb(lo)+(arb(hi)-arb(lo))*(1+((2*j+1)*pi/(2*n)).cos())/2 for j in range(n)]
C=nodes(nc,cl,cu);Yn=nodes(ny,0,Y)
V=[[Sh(c,y*(1-c.sqrt())**2/998**2,b,K) for y in Yn] for c in C]
def dct(vals,n):
    out=[]
    for m in range(n):
        s=sum((vals[j]*((m*(2*j+1))*pi/(2*n)).cos() for j in range(n)),arb(0))*2/n
        out.append(s/2 if m==0 else s)
    return out
T1=[dct(row,ny) for row in V]
T=[dct([T1[i][j] for i in range(nc)],nc) for j in range(ny)]
co=[[float(T[j][i].mid()) for j in range(ny)] for i in range(nc)]
L1=sum(abs(co[i][j]) for i in range(nc) for j in range(ny))-abs(co[0][0])
mn=min(float(v.mid()) for row in V for v in row)
print('ycell c[%g,%g] b %d K %d Y %g: centre %.4e L1rest %.4e ratio %.3f  min sample %.4e  |c10| %.2e |c01| %.2e  ctail>=8 %.1e ytail>=4 %.1e  (%.0fs)'%(cl,cu,b,K,Y,co[0][0],L1,L1/co[0][0],mn,abs(co[1][0]),abs(co[0][1]),
  sum(abs(co[i][j]) for i in range(8,nc) for j in range(ny)),sum(abs(co[i][j]) for i in range(nc) for j in range(4,ny)),time.time()-t0),flush=True)
