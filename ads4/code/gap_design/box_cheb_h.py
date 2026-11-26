# 2D Chebyshev sampling (non-rigorous) of the unnormalised prefix on a route-A rectangle
# c in [cl,cu], h in [hl,hu]*H with H = twist cap at cl (G1); predicts centre/L1 of the triangle certificate.
import sys,time
from common2 import *
cl,cu=float(sys.argv[1]),float(sys.argv[2]);b=int(sys.argv[3]);K=int(sys.argv[4])
nc=int(sys.argv[5]);nh=int(sys.argv[6]);fl,fu=(float(sys.argv[7]),float(sys.argv[8])) if len(sys.argv)>8 else (0,1)
H=min(hcap_tw(cl),cu/((1001-b)*(1002+b)))
pi=arb.pi();t0=time.time()
def nodes(n,lo,hi):return [arb(lo)+(arb(hi)-arb(lo))*(1+((2*j+1)*pi/(2*n)).cos())/2 for j in range(n)]
C=nodes(nc,cl,cu);Hn=nodes(nh,H*fl,H*fu)
V=[[Sh(c,h,b,K) for h in Hn] for c in C]
def dct(vals,n):
    out=[]
    for m in range(n):
        s=sum((vals[j]*((m*(2*j+1))*pi/(2*n)).cos() for j in range(n)),arb(0))*2/n
        out.append(s/2 if m==0 else s)
    return out
T1=[dct(row,nh) for row in V]
T=[dct([T1[i][j] for i in range(nc)],nc) for j in range(nh)]
co=[[float(T[j][i].mid()) for j in range(nh)] for i in range(nc)]
L1=sum(abs(co[i][j]) for i in range(nc) for j in range(nh))-abs(co[0][0])
print('A-rect c[%g,%g] b %d K %d h[%g,%g]H overshoot %.3f: centre %.4e L1rest %.4e ratio %.3f  ctail>=8 %.1e htail>=4 %.1e (%.0fs)'%(cl,cu,b,K,fl,fu,H/hcap_tw(cu),co[0][0],L1,L1/co[0][0],
  sum(abs(co[i][j]) for i in range(8,nc) for j in range(nh)),sum(abs(co[i][j]) for i in range(nc) for j in range(4,nh)),time.time()-t0),flush=True)
