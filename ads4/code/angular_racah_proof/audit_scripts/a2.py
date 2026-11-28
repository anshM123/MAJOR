import sympy as S
c,u,v,r=S.symbols('c u v r')
E=lambda a,b:(a*a+b*b-2*c*a*b)/(1-c)
w=(1+r)*c*v-r*u
print('E two forms', S.simplify(E(u,v)-((1+c)*v**2+(u-c*v)**2/(1-c))))
print('contract', S.simplify(E(v,w)-((1+c)*v**2+r**2*(u-c*v)**2/(1-c))))
print('dE/dc', S.simplify(S.diff(E(u,v),c)-(u-v)**2/(1-c)**2))
print('E(1,c)', S.simplify(E(1,c)-(1+c)))
T,C,nu,g,J=S.symbols('T C nu g J')
c0=1-2*nu/T; cJ=1-2*nu/(T-C)
print('ratio', S.simplify((1+c0)/(1+cJ)-1-nu*C/(T*(T-C-nu))))
# gap
Tm=(J+g-1)*(J+g-2)
print('T-C', S.expand(Tm-J*(J+1)-((2*g-4)*J+(g-1)*(g-2))))
Bg=g**2+5*g-S.Rational(3,2)
print('slack1', S.expand(Tm-Bg-((J-4)*(J+2*g+1)+S.Rational(15,2))))
h=Tm+S.Rational(9,2)-Bg
print('slack2', S.expand(2*h-J*(J+1)-((J-4)*(J+4*g-3)+4)))
# J=2 F2
Cc=6; 
n=S.Symbol('n')
F2=1+(-n)*n*Cc/T+(-n)*(-n+1)*n*(n+1)/4*Cc*(Cc-2)/(T*(T-2))
print('F2', S.simplify(F2-(1-6*nu/T+6*nu*(nu-1)/(T*(T-2)))).subs(nu,n**2).simplify())
W0=lambda TT: TT/3+S.Rational(3,2)
print('J2 weight', S.expand(W0(g*(g+1))-Bg/6-(g**2-3*g+S.Rational(21,2))/6))
# Numerical brute check of the helicity-0 chain bound on a few values (exact rationals)
from fractions import Fraction as Fr
import random
random.seed(1)
def F(xmax,nn,TT):
    vals=[Fr(1)]
    nu_=nn*nn
    cx=lambda x: 1-Fr(2*nu_)/(TT-x*(x+1))
    vals.append(cx(0))
    for x in range(1,xmax):
        vals.append(((2*x+1)*cx(x)*vals[x]-x*vals[x-1])/(x+1))
    return vals
def R0(J,nn,TT):
    tot=Fr(0);c_=Fr(1);CC=J*(J+1)
    for k in range(0,min(nn,J)+1):
        if k>0:
            i=k-1; c_=c_*Fr((i-nn)*(nn+i),k*k)*(CC-i*(i+1))/(TT-i*(i+1))
        tot+=c_
    return tot
bad=0
for trial in range(300):
    g_=Fr(random.randint(400,4000),100); nn=int(g_//2); J=random.choice(range(0,60,2))
    D=J+g_+Fr(random.randint(0,3000),100)
    TT=(D-1)*(D-2)
    vals=F(J+1,nn,TT)
    assert vals[J]==R0(J,nn,TT)
    if J>=4:
        CC=J*(J+1)
        assert vals[J]**2 <= 1+Fr(nn*nn*CC)/(TT*(TT-CC-nn*nn))
print('b0 recurrence==4F3 and chain bound on 300 random exact pts OK')
