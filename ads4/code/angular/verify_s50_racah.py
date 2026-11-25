"""ANGULAR_ROUTE §50 checks (exact rationals): ladder difference factorisation, Jacobi->Racah transform,
Racah orthogonality with weights omega_n and the author's norm formula (s = 1 and 9)."""
from fractions import Fraction as Fr
import random
def poch(a,k):
    r=Fr(1)
    for i in range(k): r*=a+i
    return r
def fact(k): return poch(1,k)
def U(n,N,s,k): return poch(-n,k)*poch(n+s,k)/(poch(-N,k)*poch(N+s,k))
def Lfun(coeffs,n,N,s): return sum(cf*U(n,N,s,k) for k,cf in enumerate(coeffs) if k<=n)
def Qj(j):  # 2F1(-j,j+2;1;z) coefficients
    return [poch(-j,i)*poch(j+2,i)/(fact(i)**2) for i in range(j+1)]
def R(j,n,N,s):
    return sum(poch(-j,i)*poch(j+2,i)*poch(-n,i)*poch(n+s,i)/(poch(1,i)*poch(1-N,i)*poch(N+s+1,i)*fact(i)) for i in range(j+1))
random.seed(5); bad=0
for _ in range(60):
    s=random.choice([1,9]); N=random.randint(4,30); n=random.randint(0,N-1); k=random.randint(0,n)
    c=Fr(n*(n+s),N*(N+s))
    if U(n,N,s,k)-U(n,N,s,k+1)!=(1-c)*poch(-n,k)*poch(n+s,k)/(poch(1-N,k)*poch(N+s+1,k)): bad+=1
    j=random.randint(0,N-1); q=Qj(j); oq=[ (q[i] if i<len(q) else 0)-(q[i-1] if 0<i<=len(q) else 0) for i in range(len(q)+1)]
    if Lfun(oq,n,N,s)!=(1-c)*R(j,n,N,s): bad+=1
print('factorisation + Jacobi->Racah: mismatches',bad)
bad=0
for s in (1,9):
    for N in (6,11,17):
        w=[(1-Fr(n*(n+s),N*(N+s)))*Fr(2*n+s,s) for n in range(N)]
        for j in range(N):
            for k in range(j,N):
                ip=sum(w[n]*R(j,n,N,s)*R(k,n,N,s) for n in range(N))
                if j!=k and ip!=0: bad+=1
                if j==k:
                    h=Fr((N+1)*(N+s-1),2*s*(j+1))
                    for l in range(j): h*=Fr((N+l+2)*(N+s-2-l),(N+s+1+l)*(N-1-l))
                    if ip!=h: bad+=1
print('orthogonality + norms (s=1,9; N=6,11,17): mismatches',bad)
