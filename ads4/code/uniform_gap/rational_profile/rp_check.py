# Independent exact checks of Theorem RP's load-bearing steps, written from the paper's Eqs. fp-U, fp-V and the note's definitions.
from fractions import Fraction as Fr
import random, math, sympy as sp
random.seed(7)
def U(b,J,Delta,k):           # paper Eq. fp-U: (b-J)_k (J+b+1)_k / ((b+2-Delta)_k (Delta+b-1)_k), zero for k>J-b
    if k>J-b: return Fr(0)
    num=Fr(1); den=Fr(1)
    for j in range(k):
        num*= (b-J+j)*(J+b+1+j); den*=(b+2-Delta+j)*(Delta+b-1+j)
    return num/den
def Lb(poly,b,J,Delta):       # L(z^k) := U_k
    return sum(c*U(b,J,Delta,k) for k,c in enumerate(poly) if c)
def mul(p,q):
    r=[Fr(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        if a:
            for j,c in enumerate(q):
                if c: r[i+j]+=a*c
    return r
def binpoly(m):               # (1-z)^m
    return [Fr((-1)**l*math.comb(m,l)) for l in range(m+1)]
def Dop(poly,b):              # D_b = z(1-z)/2 d^2 + [(7/2-2b)z+1/2] d + 9/2-b(b-1)
    n=len(poly); out=[Fr(0)]*(n+1)
    for k,c in enumerate(poly):
        if not c: continue
        if k>=2: out[k-1]+=c*Fr(k*(k-1),2); out[k]-=c*Fr(k*(k-1),2)
        if k>=1: out[k]+=c*(Fr(7,2)-2*b)*k; out[k-1]+=c*Fr(k,2)
        out[k]+=c*(Fr(9,2)-b*(b-1))
    return out
def Vb(poly,b,J,Delta):
    T=(Delta-1)*(Delta-2)
    d=Dop(poly,b); full=[ (poly[k] if k<len(poly) else 0)+ d[k]/T for k in range(len(d))]
    return Lb(full,b,J,Delta)
# 1. positive basis lemma L_b[z^k(1-z)^m0] >= 0 for a>=M, m0<=M+1 (small M to test exhaustively-ish, with fractional twists)
bad=0; tests=0
for M in (3,6,10):
    for _ in range(400):
        b=random.choice((0,4)); J=b+random.randint(0,14); a=Fr(random.randint(M*7,M*7+60),7)
        if a<M: continue
        Delta=J+a+2; k=random.randint(0,J-b+2); m0=random.randint(0,M+1)
        v=Lb(mul([Fr(0)]*k+[Fr(1)],binpoly(m0)),b,J,Delta); tests+=1
        if v<0 or (k<=J-b and v==0): bad+=1
print('lemma tests',tests,'failures',bad)
# negative control: the lemma should FAIL for twists well below M (a<M) sometimes
fails=0
for _ in range(400):
    b=0; J=random.randint(2,14); a=Fr(random.randint(1,20),10); Delta=J+a+2; m0=10; k=0
    if Lb(binpoly(m0),b,J,Delta)<0: fails+=1
print('control (a<<M, m0=10): negative values found',fails,'of 400 (expected > 0)')
# 2. full action for the actual strict profile at M=998 on many heavy states (exact rationals)
M=998; m=M+6; eps=Fr(1,(M+1)**8); rho=[Fr(j,100*m*m) for j in range(1,m+1)]
def Fstrict(N):                # Taylor coefficients to degree N
    R=[Fr(1)]+[Fr(0)]*N
    for r in rho:              # multiply by 1/(1-r z)
        for i in range(1,N+1): R[i]+=r*R[i-1]
    A=mul(binpoly(M)[:N+1],R)[:N+1]
    num=[Fr(0),eps,Fr(0),Fr(1)]
    return mul(num,A)[:N+1]
import time; t0=time.time()
F=Fstrict(70)
mins={}
for b in (0,4):
    for J in list(range(b,b+12))+[20,33,50,64]:
        for a in (Fr(998),Fr(9985,10),Fr(1200),Fr(5000),Fr(10**5)):
            Delta=J+a+2; v=Vb(F,b,J,Delta)
            key=(b,J); mins[key]=min(mins.get(key,v),v)
neg=[k for k,v in mins.items() if v<=0]
print('full-action states tested',len(mins)*5,'nonpositive:',neg,' smallest:',min(float(v) for v in mins.values()), 'time',round(time.time()-t0,1))
# 3. light identity Q = -100 g_{2,3} - 36 g_{3,3}
y=sp.symbols('y'); p=3; g=[sp.Integer(1)]; gm1=sp.Integer(0)
for k in range(0,3):
    nxt=sp.expand(((8*k**2+8*p*k+4*p-1-4*y)*g[k]-4*(k-1+p)**2*(g[k-1] if k>=1 else gm1))/(4*(k+1)**2)); g.append(nxt)
nu=sp.symbols('nu'); Q=nu**6-sp.Rational(245,4)*nu**4+sp.Rational(9619,16)*nu**2-sp.Rational(51575,64)
print('Q + 100 g2 + 36 g3 =',sp.simplify(Q+100*g[2].subs(y,nu**2)+36*g[3].subs(y,nu**2)))
# 4. the constant 280 at M=998 (strict version)
CI=Fr(100,99)*(6*(1+Fr(4,M+5))+2*(1+Fr(2,M+3))/Fr(M+1)**6)
print('460800*C_I/35.99 =',float(460800*CI/Fr(3599,100)),'< 280^2 =',280**2, 460800*CI/Fr(3599,100)<280**2)
for Dg in (1000,10**4,10**5,10**6): print('tau_RP(%g) ='%Dg, 280*math.sqrt(math.log((Dg+2)/3))/(Dg-2)**2)
