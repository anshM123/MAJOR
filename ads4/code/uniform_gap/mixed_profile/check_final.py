#!/usr/bin/env python3
"""Exact regression/constant checker for the all-gap mixed rational-profile theorem.
Requires Python 3 and SymPy. The analytic proof, not the samples, proves the continuum.
"""
from fractions import Fraction as Q
from math import comb
import sympy as s

def ladder(n,a,spin):
    U=[Q(1)]
    for j in range(n):
        U.append(U[-1]*(n-j)*(n+spin+j)/((n+a-j)*(n+a+spin+j)))
    return U

def L(U,k,m):
    assert k>=0 and m>=0
    if k>=len(U):return Q(0)
    return sum((-1)**j*comb(m,j)*U[k+j] for j in range(min(m,len(U)-k-1)+1))

M=998
assert Q(9,2)*998+Q(15,2)<Q(998**2,100)
assert 1-Q(1,4)-Q(1,100)-Q(7,10)==Q(1,25)
assert Q(28,5)-Q(36,M)>Q(11,2)
assert Q(7,10)*M>Q(9,2)
constant=Q(460800)*Q(2,11)*Q(100,99)*(2*(1+Q(3,M+4))+(1+Q(2,M+3))*Q(1,(M+1)**7))*Q(M+1,M)**2
assert constant<420**2
print("PASS exact constants; squared rigidity coefficient <=",float(constant))

z,nu,k,m,b=s.symbols('z nu k m b')
ld=k/z-m/(1-z)
D=s.Rational(1,2)*z*(1-z)*(ld**2+s.diff(ld,z))+((s.Rational(7,2)-2*b)*z+s.Rational(1,2))*ld+s.Rational(9,2)-b*(b-1)
C=k*(s.Rational(7,2)-2*b-m)+s.Rational(9,2)-b*(b-1)-m/2
assert s.cancel(D-(k*(k-1)*(1-z)/(2*z)+k/(2*z)+C+m*(m-9+4*b)*z/(2*(1-z))))==0
# g_{2,4} from the Legendre equation after z=(x-1)/(x+1).
g0=s.Integer(1);g1=4-s.Rational(1,4)-nu**2
g2=s.expand(((2+8+4-s.Rational(1,4)-nu**2)*g1-16*g0)/4)
# Correct recurrence: (k+1)^2 g_{k+1} =
# [2 k^2 + 2 p k + p - 1/4 - nu^2]g_k - (k+p-1)^2 g_{k-1}
assert s.expand(4*g2-(nu**4-s.Rational(35,2)*nu**2+s.Rational(569,16)))==0

count=0;minimum=Q(100)
for ag in (Q(998),Q(1000),Q(2001,2),Q(100000)):
    m=int(ag);amp=Q(7,10)*m*m
    for a in (ag,ag+Q(3,7),ag+1000):
        for b in (0,4):
            for n in (0,1,2,3,5,13,31,95):
                J=n+b;T=(J+a)*(J+a+1);U=ladder(n,a,2*b+1)
                for k in (1,2,3,4,8):
                    base=L(U,k,m)
                    C=k*(Q(7,2)-2*b-m)+Q(9,2)-b*(b-1)-Q(m,2)
                    dressed=(1+C/T)*base+Q(k,2)/T*L(U,k-1,m)+Q(m*(m-9+4*b),2)/T*L(U,k+1,m-1)
                    if k>=2:dressed+=Q(k*(k-1),2)/T*L(U,k-1,m+1)
                    direct=Q(0)
                    for j in range(min(m,n+1-k)+1):
                        ell=k+j;ul=U[ell] if ell<len(U) else Q(0);um=U[ell-1] if 1<=ell<=len(U) else Q(0)
                        gamma=Q(-ell*ell+(8-4*b)*ell+9-2*b*(b-1),2)
                        direct+=(-1)**j*comb(m,j)*((1+gamma/T)*ul+Q(ell*ell,2)*um/T)
                    assert direct==dressed
                    assert dressed-amp/T*base>=Q(1,25)*base
                    if k==1:assert dressed>0
                    if base>0:minimum=min(minimum,(dressed-amp/T*base)/base)
                    count+=1
print("PASS",count,"exact dressing and mixed-action regression tests; sampled minimum ratio",float(minimum))

# Low-spin check: exact formula for the long spin-2 representation.
T=s.symbols('T')
U1=6/T;U2=24/(T*(T-2))
assert s.cancel(T*((1+8/T)*U1+1/(2*T))/U2-(T-2)*(13*T/48+2))==0
assert s.cancel(T*((1+s.Rational(21,2)/T)*U2+2*U1/T)/U2-(s.Rational(3,2)*T+s.Rational(19,2)))==0
print("PASS low-spin reduction and sign constants")

# Why the original positive basis cannot simply be used at quadratic degree.
a=Q(998);n=1;spin=9;M=998**2//4
v=L(ladder(n,a,spin),0,M)
assert v==Q(-741509,503496)
print("PASS quadratic-degree basis counterexample:",v)
