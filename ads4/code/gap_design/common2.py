# Non-rigorous (arb, midpoint) evaluation tools on the author's 17,201-pair table (scale 10^3000).
# Termination-preserving coordinates: x=h/c, U_k = c^k prod_{j<k}(1-x p_j)/(1-c x p_j), p_j=j(j+s).
import sys,math,os
from flint import arb,ctx
sys.set_int_max_str_digits(0)
ctx.prec=int(os.environ.get('PREC','3400'))
TABLE=os.environ.get('TABLE17','/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/fullband/x/sliver_full_band_reduction/coeff_long_17200.txt')
_A=[];_G=[]
def load(K):
    if len(_A)>=K+2:return
    f=open(TABLE);N,D=map(int,f.readline().split());sc=arb(10)**D
    _A.clear();_G.clear()
    for i,l in enumerate(f):
        if i>=K+2:break
        a,r,g,t=l.split();_A.append(arb(a)/sc);_G.append(arb(g)/sc)
_cache={}
def parts(b,K):
    if (b,K) not in _cache:
        load(K);A,G=_A,_G;sg=1 if b==0 else -1;Bs=[]
        for k in range(K+1):
            k1=(k+1)**2;k2=-k*k+(8-4*b)*k+9-2*b*(b-1)
            Bs.append((k1*A[k+1]+k2*A[k])/2+sg*G[k])
        _cache[(b,K)]=(A[:K+1],Bs)
    return _cache[(b,K)]
def Sx(c,x,b,K,norm=False):
    """prefix numerator in (c,x); h=c x"""
    A,Bs=parts(b,K);s=2*b+1;q=b*(b+1);c=arb(c);x=arb(x);h=c*x;S=None
    for k in range(K,-1,-1):
        Fk=A[k]+h*(q*A[k]+Bs[k])
        if S is None:S=Fk;continue
        p=k*(k+s);S=Fk+c*(1-x*p)/(1-c*x*p)*S
    return S*(1-c)**3/(1+q*h)**4 if norm else S
def xphys(c):  # twist cap in x (binding for c>=0.26)
    c=float(c);return (1-math.sqrt(c))**2/(998**2*c)
def nmin(c,b):
    """smallest n allowed by the twist cap at c: c<=n(n+s)/((n+998)(n+998+s))"""
    s=2*b+1;n=1001-b
    lo,hi=n,10**6
    f=lambda n:n*(n+s)/((n+998)*(n+998+s))>=c
    if f(lo):return lo
    while hi-lo>1:
        m=(lo+hi)//2
        if f(m):hi=m
        else:lo=m
    return hi
def Sh(c,h,b,K,norm=False):
    """prefix numerator in (c,h)"""
    c=arb(c);h=arb(h);return Sx(c,h/c,b,K,norm)
def hcap_tw(c):return (1-math.sqrt(float(c)))**2/998**2
def needK(c,E=2337,target_log2=-110):
    # author's tight envelope tail (tight_spectral_tail.py): 2^E t^(K+1)[...] < 2^target ; float version
    t=float(c)/(199/200);H=1/998**2
    def tl(K):return E*math.log(2)+(K+1)*math.log(t)+math.log((1+20*H)/(1-t)+2*H*((K+5)**2/(1-t)+2*(K+5)*t/(1-t)**2+t*(1+t)/(1-t)**3))
    K=10
    while tl(K)>=target_log2*math.log(2):K+=1
    return K
