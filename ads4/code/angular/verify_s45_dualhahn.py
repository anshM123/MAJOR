from fractions import Fraction as Fr
import random
def poch(a,k):
    r=Fr(1)
    for i in range(k): r*=a+i
    return r
def fact(k): return poch(Fr(1),k)
def pochpm(A,x,j):  # (A+i nu)_j (A-i nu)_j with x=nu^2
    r=Fr(1)
    for i in range(j): r*=(A+i)**2+x
    return r
def gk(p,x,k):
    A=Fr(2*p-1,2)
    s=sum(poch(-k,j)*pochpm(A,x,j)/(poch(p,j)**2*fact(j)) for j in range(k+1))
    return poch(p,k)**2/fact(k)**2*s
def U(n,N,b,k):
    db=2*b+1; return poch(-n,k)*poch(n+db,k)/(poch(-N,k)*poch(N+db,k))
def Kdirect(p,x,n,N,b): return sum(U(n,N,b,k)*gk(p,x,k) for k in range(n+1))
def R(m,n,a,b,p):
    return sum(poch(-m,l)*poch(a+1,l)*poch(-2*n-a-2*b,l)/(poch(-2*n-2*b,l)*poch(1-p-n,l)*fact(l)) for l in range(m+1))
def Kclosed(p,x,n,N,b):
    a=N-n
    s=sum(poch(-n,j)*pochpm(Fr(1,2),x,j)/(poch(p,j)*fact(j)**2)*R(n-j,n,a,b,p) for j in range(n+1))
    return U(n,N,b,n)*poch(p,n)/fact(n)*s
random.seed(3); bad=0
for _ in range(40):
    p=random.choice([2,3,4]); b=random.choice([0,4]); n=random.randint(1,14); N=n+random.randint(2,40); x=Fr(random.randint(0,400),random.randint(1,7))
    if Kdirect(p,x,n,N,b)!=Kclosed(p,x,n,N,b): bad+=1
print('closed form vs direct: mismatches', bad, '/40')
# five-point equation
def K(p,x,n,N,b): return Kdirect(p,x,n,N,b)
bad=0
for _ in range(40):
    p=random.choice([2,3,4]); b=random.choice([0,4]); n=random.randint(1,12); N=n+random.randint(3,30); x=Fr(random.randint(0,300),random.randint(1,5))
    db=2*b+1
    Ap=Fr((n+1)*(n+db),(2*n+db)*(2*n+db+1)); Am=Fr(n*(n+db-1),(2*n+db)*(2*n+db-1))
    Bm=Fr((N+p)**2*(N+db-1),(2*N+db)*N*(2*N+db-1)); Bp=Fr((N+db-p)**2*(N+1),(2*N+db)*(N+db)*(2*N+db+1))
    k0=K(p,x,n,N,b)
    dn=Ap*(K(p,x,n+1,N,b)-k0)+Am*(K(p,x,n-1,N,b)-k0)
    dN=Bp*(K(p,x,n,N+1,b)-k0)+Bm*(K(p,x,n,N-1,b)-k0)
    m2=(N-n)*(N+n+db)
    lhs=(x+Fr(2*p-1,2)**2)*k0; rhs=m2*((dN-dn)+Fr(p*p,N*(N+db))*k0)
    if lhs!=rhs: bad+=1
print('five-point equation: mismatches', bad, '/40')
