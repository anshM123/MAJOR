# Shared non-rigorous (arb midpoint) tools for the gap-design numerics.
# Loads the realised direct coefficient table (midpoints only) and evaluates the
# prefix S_0(c,h) = sum_k F_k(h) U_k(c,h) exactly as in sliver_direct/direct_sign_cheb.py.
import sys,math,os
from flint import arb,ctx
from fractions import Fraction as Fr
sys.set_int_max_str_digits(0)
ctx.prec=int(os.environ.get('PREC','4800'))
TABLE=os.environ.get('TABLE','/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/gd/direct_coefficients.txt')
_A=None
def load():
    global _A,_G,_N
    if _A is not None:return
    f=open(TABLE);N,D=map(int,f.readline().split());sc=arb(10)**D
    A=[];G=[]
    for l in f:
        a,r,g,t=l.split();A.append(arb(a)/sc);G.append(arb(g)/sc)
    _A,_G,_N=A,G,N
def coeffs():
    load();return _A,_G
def need_k(c):
    # as in the author's drivers: (c/R)^(K+1) <= 2^-2700, R=1999/2000
    t=float(c)/(1999/2000);return int(math.ceil(2700*math.log(2)/(-math.log(t))))-1
def Fk_parts(b,K):
    """returns lists A_k, Bsig_k=(B_k+sigma G_k) for k<=K"""
    A,G=coeffs();s=2*b+1;q=b*(b+1);sg=1 if b==0 else -1;Bs=[]
    for k in range(K+1):
        k1=(k+1)**2;k2=-k*k+(8-4*b)*k+9-2*b*(b-1)
        Bs.append((k1*A[k+1]+k2*A[k])/2+sg*G[k])
    return A[:K+1],Bs
_cache={}
def parts(b,K):
    if (b,K) not in _cache:_cache[(b,K)]=Fk_parts(b,K)
    return _cache[(b,K)]
def S0(c,h,b,K,norm=True):
    A,Bs=parts(b,K);s=2*b+1;q=b*(b+1);c=arb(c);h=arb(h);S=arb(0)
    for k in range(K,-1,-1):
        Fk=(1+q*h)*A[k]+h*Bs[k]
        if k<K:
            p=k*(k+s);S=Fk+(c-p*h)/(1-p*h)*S
        else:S=Fk
    return S*(1-c)**3/(1+q*h)**4 if norm else S
def hphys(c,b):
    c=float(c);return min(c/((1001-b)*(1002+b)),(1-math.sqrt(c))**2/998**2)
