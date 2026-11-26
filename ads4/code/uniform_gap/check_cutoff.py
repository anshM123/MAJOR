#!/usr/bin/env python3
"""Symbolic identities and exact elementary constants for the new cutoff family.
This checks arithmetic in the analytic proof, not swappability of CFT sum rules.
"""
from check_uniform import *
import sympy as s
from math import factorial
x=s.symbols('x');I=s.I
Rm=-2*(2*x-I)*(2*x+I)**2*(2*x+3*I)*(2*x+5*I)*(2*x+11*I)/(x*(2*x-7*I)*(2*x-3*I)*(2*x+7*I));Rp=Rm.subs(x,-x)
P=(x*x+s.Rational(1,4))**2*(x*x+s.Rational(9,4))*(x*x+s.Rational(25,4))*(x*x+s.Rational(81,4))
Bpol=-x*(2*x-7*I)**2*(2*x-3*I)**2*(2*x+3*I)*(2*x+7*I)**2/2048
assert s.cancel(P.subs(x,x+I)/Rm-Bpol)==0
# ratio of the cosh and Gamma factors from Gamma recurrence is (z+3i/2)/(z-3i/2).
assert s.cancel(P.subs(x,x+I)/P.subs(x,x-I)*(x+3*I/2)/(x-3*I/2)-Rm/Rp)==0
Q=x**6-s.Rational(245,4)*x**4+s.Rational(9619,16)*x**2-s.Rational(51575,64)
N10=1024*x**10-195328*x**8+8229504*x**6-90780512*x**4+251603732*x**2-114141375
raw=-(4*x*x+1)*N10/(8*(4*x*x+9)*(4*x*x+49))
assert s.cancel(Rm*Q.subs(x,x+I)+Rp*Q.subs(x,x-I)-raw)==0
Eraw=-32*(16*x**4-440*x**2+273)/((4*x**2+49)*(4*x**2+9))
assert s.cancel(Rm/((x+I)**2+s.Rational(9,4))+Rp/((x-I)**2+s.Rational(9,4))-Eraw)==0
# elementary e bounds by power-series and geometric upper tail
assert sum(F(5)**k/factorial(k) for k in range(10))>128 # e^5 > 2^7
assert sum(F(3)**k/factorial(k) for k in range(13))+F(3)**13/factorial(13)/(1-F(3,14))<21
assert sum(F(2,9)**k/factorial(k) for k in range(8))+F(2,9)**8/factorial(8)/(1-F(2,81))<F(4,3)
# F lower polynomial ratios, gamma modulus lower constant
p_ratio=F(1,5)**2*F(9,13)*F(25,29)*F(81,85)
assert p_ratio>F(1,50)
assert F(1,50)*F(1,4)*F(1,7)>F(1,2048)
# F shifted upper, B upper from P shift and gamma ratio
assert F(5,4)**2*F(13,4)*F(29,4)*F(85,4)<2**10
assert 2*F(22,7)*F(4,3)<9
assert F(7,4)**5<F(17,8)**4 # (7/4)^(5/4) < 17/8
assert F(3,4)<F(15,16)**4
assert sum(F(10,63)**k/factorial(k) for k in range(8))+F(10,63)**8/factorial(8)/(1-F(10,567))<F(6,5)
assert F(22,7)*2*F(6,5)*F(17,8)*F(15,16)<16
assert 16*F(7**4*3**3,2048)<512
# polynomial-times-exponential envelope
assert 24**12<2**60
bc=bern_1(c2);cells=[bc]
for _ in range(4):cells=[p for z in cells for p in split1(z)]
assert max(abs(z) for p in cells for z in p)<32
assert 2**(25+5+60)<=2**90 # |chi R±| bound; kernel 1/2, integral <=2
assert 2**(20+5+60)<=2**90 # |chi| bound
# raw light integrands (absolute rational majorants)
assert 32*440/16<2**10
assert sum(abs(int(c)) for c in s.Poly(N10,x).all_coeffs())<2**29
assert F(2**29,8)<2**30
# finite-spin polynomial recurrence bounds, scaled freq y=(nu+2)/Lambda
# for k<=64, (k+1)^2 ghat_{k+1} = ...; after multiplying by (k!)^2,
# majorant by (1+y^2)^k, use (y²+a)/(1+y²)+b/(1+y²)^2<=1 if a+b<=1.
a=F(2*64**2+2*3*64+3,D);b=F((64+2)**2*64**2,D*D)
assert a+b<F(1,100)
assert F((64+10)**2,D)<F(1,100)
# bounded impact: Cauchy radius z=1600/Lambda^2 for bhat<=40.
assert F(80)/(1-F(40,998))+F(1600,998)/(1-F(40,998)**2)<86
assert (1-F(1600,D))**3>F(1,2)
assert F(16)/(1-F(64,D))<17
# beta<=17/998<.018; exp(beta)-1<=2 beta by geometric series
assert F(17,998)<F(1,2)
assert (1-F(64,D))**3>F(1,2)
r=F(225,16*64)
assert r<F(1,4)
# geometric sum at r=1/4: sum[1+t0(k+10)^2+(k+1)^2/(2R)]r^k <4
r=F(1,4);s0=1/(1-r);s1=r/(1-r)**2;s2=r*(1+r)/(1-r)**3
assert s0+(s2+20*s1+100*s0)/D+(s2+2*s1+s0)/128<4
# common kernel envelope 16(1+y)^2 exp(140(1+y)).
assert F(140,998)<F(1,2)
assert 2+F(1002,998)<4
assert 26<32
# angular cutoff difference: interior <=2^102 e^-9440; exterior <=2^105 e^-440
# e^5>2^7 -> tails below 2^-500 after adding.
assert 105-7*(440//5)<=-511
assert 102-7*(9440//5)<-511
# light E error <= 2^104(3e^-10000+e^-1000); K <=2^121(3^7 e^-10000+3^6*2^12 e^-1000)
assert sum(comb(6,i)*factorial(i) for i in range(7))<2**12
assert 104+2-7*(1000//5)<-501
assert 121+10+12-7*(1000//5)<-501
# fixed-spin normalized lower bounds from exact ledgers
r1=json.load(open(Path(__file__).with_name('results_J32.json')))['rows']
rows=r1+json.load(open(Path(__file__).with_name('physical_33_48.json')))['rows']+json.load(open(Path(__file__).with_name('physical_49_64.json')))['rows']
assert len(rows)==94 and all(z['pass'] for z in rows)
lo=[]
for z in rows:
 J=z['J'];b=z['b'];M=J-b+1
 if 'normalized_bound' in z:lo.append(F(z['normalized_bound']))
 else:lo.append(F(z['min_Bernstein'])/(1+F(J*(J+1),D))**(M+3))
assert min(lo)>F(1,10**8)
print('PASS: Gamma factorization, light polynomial identity, uniform cutoff constants, all 94 spin-channel lower bounds > 1e-8; angular and light cutoff errors < 2^-500.')
print('minimum normalized ideal low-spin bound',float(min(lo)))
rr=json.load(open(Path(__file__).with_name('impact_tensor_X400_J65_results.json')))
assert rr['X']==400 and rr['J_min']==65 and all(z['pass'] and F(z['H_lower'])>F(1519,100000)+F(1,2**500) for z in rr['rows'])
print('PASS: cutoff transfer to all-spin bhat<=40 rectangle')
