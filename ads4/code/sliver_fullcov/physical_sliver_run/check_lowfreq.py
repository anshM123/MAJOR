#!/usr/bin/env python3
"""Verify constants in direct beta-Cauchy low-frequency bound and integrate logged EPSW masses.
All comparisons use integer/Fraction arithmetic; float is not used for certificates.
The EPSW masses and floor are inherited, not recomputed here.
"""
from fractions import Fraction as F
from pathlib import Path
import csv, json
S=10**40

def cd(n,d): return -(-n//d)
def up(x): x=F(x); return cd(x.numerator*S,x.denominator)
def exp_pos(x):
    assert x>=0
    k=0
    while x>S//2: x=cd(x,2); k+=1
    t=s=S
    for j in range(1,71):
        t=cd(t*x,S*j); s+=t
    # Since x/S <= 1/2, the true omitted tail is less than 1/S.
    s+=2
    for _ in range(k): s=cd(s*s,S)
    return s

def exp_lower_pos(x):
    assert x>=0
    k=0
    while x>S//2: x=x//2; k+=1
    t=s=S
    for j in range(1,71): t=t*x//(S*j);s+=t
    for _ in range(k):s=s*s//S
    return s

def exp_u(x):
    x=F(x)
    if x>=0:return F(exp_pos(up(x)),S)
    a=-x;lo=(a.numerator*S)//a.denominator
    return F(cd(S*S,exp_lower_pos(lo)),S)

def log_unit(z):
    assert 1<=z<=2
    v=(z-1)/(z+1);s=2*sum((v**(2*j+1)/ (2*j+1) for j in range(80)),F(0))
    err=2*v**161/(161*(1-v*v))
    return s,s+err

def log_box(x):
    x=F(x); assert x>0
    k=0
    while x<1:x*=2;k-=1
    while x>2:x/=2;k+=1
    lo,hi=log_unit(x);l2,u2=log_unit(F(2));
    return lo+k*(l2 if k>=0 else u2),hi+k*(u2 if k>=0 else l2)

def verify_constants():
    a=F(998);tmax=F('.26');tmin=F('.012');pi=F(22,7)
    # K0 <=1.6 via Robbins/Stirling on the physical range n>=a.
    assert (pi/2)**2*(1+2/a)*exp_u(F(1,3)/a)<F('1.6')**2
    # Dressing Cauchy radius t/200.
    assert 1+(80000+1000*tmax+8*tmax*tmax)/a**2<F('1.081')
    assert (F('.0065')*F('1.55'))**2*(F('.25')-tmax/200)>F('.005')**2
    assert (F('.0101')*F('.995'))**2*(F('.25')-tmax/200)>F('.005')**2
    assert F('1.555') <= 1+F('.75')*(1-tmax)
    # Central angular measure. Jacobian and half-Gaussian ratio.
    assert F('1.004')**2*(1-tmax*tmax/(16*(1-tmax)))>1
    assert F('2.16')**2*F('.999')**2*F(7,8)>(2*F('1.004'))**2
    # log(1+X^2) >= 7 X^2/8 on 0<=X<=1/2, by log(1+u)>=u-u^2/2.
    assert 2*F(7,8)*999*F('1.55')**2>4200
    # Normalization lower bound: integrate Gaussian to X=1; omitted relative tail <= exp(-999/2).
    assert exp_u(-F(999,2))<F('.001')
    # beta low and high Chernoff rates.
    x=F('.752');assert F('.75')*(1+9*tmax/a)<x
    assert a*(-log_box(x)[1]-(1-x))>36
    v=F('.662');m=tmax;assert 1-F('.25')/(1-tmax)>v
    rate=log_box(m/v)[0]+(1/m-1)*log_box((1-m)/(1-v))[0]
    assert a*rate>1000
    # Angular tail X>=1/2, derivative negative (k in [1,1.555], nu<=200).
    assert 200*2*F('1.555')<999
    pref=pi*32/tmin*exp_u(F(999,1996));assert pref<14000
    damp=F(999,2)*log_box(F(5,4))[0]
    assert pref*exp_u(F(100,1)/F('1.555')-damp)<F('1e-14')
    assert pref*exp_u(100-damp)<1
    # Combine 3 regions. q=1-r after accounting for enlargement; d<=2t.
    C=F('1.6')*F('1.081')*pi/2
    principal=C*(2/F('1.55'))**3*(F('2.16')+F('1e-14'))
    lower=C*(2/F('.995'))**3*exp_u(-36+F('2.02'))*(F('2.16')*exp_u(F(200**2,1)/(2*F(7,8)*999))+1)
    outer=C*(2/F('.995'))**3*exp_u(-1000+pi*100)
    assert principal+lower+outer <14
    return {'principal_constant_upper':str(float(principal)), 'excluded_beta_constant_upper':str(float(lower)), 'outer_constant_upper':str(float(outer)), 'combined_lt_14':True}

def K(v):return 14*exp_u(F('0.0065')*v+v*v/4200)

def main():
    print(json.dumps(verify_constants()))
    rows=list(csv.DictReader((Path(__file__).parent/'residual_masses.csv').open()))
    assert len(rows)==199
    # [0,.2] + [.2,1] EPSW masses: 7.8938e-7+1.1063e-6 < 1.9e-6.
    assert F('7.8938e-7')+F('1.1063e-6')<F('1.9e-6') and K(F(1))<15
    # Inherited global cutoff L1 bound (TAIL_AUDIT's weighted bound is stronger).
    cutoff=F('1.34e-14');floor=F('1.87e-6');total=15*F('1.9e-6');out=[]
    for i,row in enumerate(rows,1):
        lo=F(row['lo']);hi=F(row['frequency_ceiling']);mass=F(row['logged_mass']);assert lo==i and hi>=i+1 and mass>=0
        total+=F('1.002')*mass*K(hi)
        if i in (39,79,119,139,159,179,199):
            # Last interval only needs nu<=i+1, so avoid extending theorem beyond 200.
            E=F(2,998)*(total+cutoff*K(min(hi,F(200))))
            # For i=199, K(200.1) is a monotone upper envelope for K(200), not a claim beyond 200.
            out.append({'frequency_end':i+1,'error_upper_scaled_1e40':str(up(E)),'error_upper':str(float(F(up(E),S))), 'ratio_upper':str(float(F(up(E/floor),S)))})
            if i==159:assert E<F('1.46e-7') and E<F('.08')*floor
    print(json.dumps(out,indent=2));print('PASS: conditional on the retained EPSW residual masses, cutoff integral, and F_cert>=1.87e-6.')
if __name__=='__main__':main()
