# Independent referee check: Hall viscosity from lam*Q, Q = C C Ctilde, on AdS4 brane, ingoing EF.
import sympy as sp, itertools, time
t0=time.time()
v,r,x,y,e=sp.symbols('v r x y e')
X=[v,r,x,y]
p=sp.Function('p')(v,r); c=sp.Function('c')(v,r)
f=r**2-1/r
H=sp.Matrix([[p,c],[c,-p]])
S=(sp.eye(2)+e*H+e**2*H*H/2)*r**2
g=sp.zeros(4)
g[0,0]=-f; g[0,1]=g[1,0]=1
for i in range(2):
    for j in range(2): g[2+i,2+j]=S[i,j]
def tr(ex,n=2):
    ex=sp.expand(ex)
    return sum(ex.coeff(e,k)*e**k for k in range(n+1))
g0=g.subs(e,0); g1=g.diff(e).subs(e,0); g2=(g.diff(e,2)/2).subs(e,0)
gi0=g0.inv()
gi=gi0 - e*gi0*g1*gi0 + e**2*(gi0*g1*gi0*g1*gi0 - gi0*g2*gi0)
gi=gi.applyfunc(tr)
dg=[[[sp.diff(g[a,b],X[k]) for k in range(4)] for b in range(4)] for a in range(4)]
Gam=[[[tr(sum(gi[a,d]*(dg[d][b][cc]+dg[d][cc][b]-dg[b][cc][d]) for d in range(4))/2) for cc in range(4)] for b in range(4)] for a in range(4)]
print('Gamma',time.time()-t0)
Rud={}
for a,b,cc,d in itertools.product(range(4),repeat=4):
    if cc>=d: continue
    val=sp.diff(Gam[a][d][b],X[cc])-sp.diff(Gam[a][cc][b],X[d])+sum(Gam[a][cc][k]*Gam[k][d][b]-Gam[a][d][k]*Gam[k][cc][b] for k in range(4))
    val=tr(val)
    Rud[(a,b,cc,d)]=val; Rud[(a,b,d,cc)]=-val
def Rup(a,b,cc,d): return Rud.get((a,b,cc,d),0)
Rl={}
for a,b,cc,d in itertools.product(range(4),repeat=4):
    Rl[(a,b,cc,d)]=tr(sum(g[a,k]*Rup(k,b,cc,d) for k in range(4)))
Ric=sp.Matrix(4,4,lambda b,d: tr(sum(Rup(a,b,a,d) for a in range(4))))
Rs=tr(sum(gi[a,b]*Ric[a,b] for a in range(4) for b in range(4)))
print('Riemann',time.time()-t0)
C={}
for a,b,cc,d in itertools.product(range(4),repeat=4):
    C[(a,b,cc,d)]=tr(Rl[(a,b,cc,d)]-sp.Rational(1,2)*(g[a,cc]*Ric[d,b]-g[a,d]*Ric[cc,b]-g[b,cc]*Ric[d,a]+g[b,d]*Ric[cc,a])+Rs/6*(g[a,cc]*g[d,b]-g[a,d]*g[cc,b]))
# mixed C^{ab}_{cd}
Cm={}
pairs=[(a,b) for a in range(4) for b in range(4) if a<b]
for (a,b) in pairs:
    for (cc,d) in pairs:
        Cm[(a,b,cc,d)]=tr(sum(gi[a,k]*gi[b,l]*C[(k,l,cc,d)] for k in range(4) for l in range(4)))
print('Weyl',time.time()-t0)
sqrtg=r**2  # det(e^H)=1 to O(e^2)
ORIENT=1
# eps^{abcd} = -ORIENT/sqrtg [abcd]  (eps_{vrxy}=+ORIENT sqrtg)
from sympy.combinatorics import Permutation
def lc(*idx):
    if len(set(idx))<4: return 0
    return Permutation(list(idx)).signature()
# Ctilde^{ab}_{cd} = 1/2 eps^{abef} C_{efcd}
Ct={}
for (a,b) in pairs:
    for (cc,d) in pairs:
        s=0
        for (e1,f1) in pairs:
            s+= lc(a,b,e1,f1)*C[(e1,f1,cc,d)]   # 1/2 * 2 for e<f
        Ct[(a,b,cc,d)]=tr(-ORIENT*s/sqrtg)
# Q = C^{ab}_{cd} C^{cd}_{ef} Ct^{ef}_{ab}; sums over pairs get factor 2 each (antisym) -> factor 8
Q=0
for P1 in pairs:
    for P2 in pairs:
        a1=Cm[P1+P2]
        if a1==0: continue
        for P3 in pairs:
            Q+= a1*Cm[P2+P3]*Ct[P3+P1]
Q=tr(8*Q)
print('Q built',time.time()-t0)
L2Q=sp.expand(Q.coeff(e,2)*sqrtg)
L2E=sp.expand(tr(sqrtg*(Rs+6)).coeff(e,2))
print('Q order0,1:',sp.simplify(Q.coeff(e,0)),sp.simplify(Q.coeff(e,1)))
# chi EOM
from sympy.calculus.euler import euler_equations
def EL(L,fn):
    return sp.expand(euler_equations(L,[fn],[v,r])[0].lhs)
EE=EL(L2E,c); EQ=EL(L2Q,c)
print('Einstein chi EOM:',sp.simplify(EE))
# isolate phi terms in Q's chi EOM
print('Q chi EOM:',sp.simplify(EQ))
EOUT=sp.simplify(EE)
