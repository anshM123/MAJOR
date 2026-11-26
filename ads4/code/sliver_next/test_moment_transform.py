#!/usr/bin/env python3
"""Exact ideal-weight moment inputs, global transport, and pointwise ideal-action checks.
Realized-weight moments are NOT generated here; a radius for the repaired moments
must be obtained independently before this can certify the full low-L sliver.
"""
import json,subprocess,time,math,sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
K=int(sys.argv[1]) if len(sys.argv)>1 else 3249
D=int(sys.argv[2]) if len(sys.argv)>2 else 2500
S=10**D;NF=json.load(open(HERE/'nf02_exact.json'));LAM2=997002
sys.set_int_max_str_digits(0)

def cd(x,y):return -(-x//y)
def moments(p):
    # The coefficients are exact rational NF02 coefficients, divided by k0.
    psi=[F(z)/F(NF['k0']) for z in NF[f'psi{p-1}_exact']]
    for j in range(K+1):
        x=sum((a/F(2*j+i+1) for i,a in enumerate(psi)),F(0));v=int(x*S)
        yield f'{v} 1\n'

def transform(p):
    t=time.perf_counter();v=f'{p} {K} {D}\n'+''.join(moments(p))
    r=subprocess.run([str(HERE/'moment_transform')],input=v,text=True,capture_output=True,check=True)
    rows=[tuple(map(int,l.split())) for l in r.stdout.splitlines()];assert [k for k,_,_ in rows]==list(range(K+1));assert rows[-1][2]<S//10**100
    print(f'p={p} K={K}: seconds {time.perf_counter()-t:.3f}; max error < 1e-{D-len(str(rows[-1][2]))}',flush=True)
    # Independent exact rational check at low degree and first 15 k values.
    psi=[F(z)/F(NF['k0']) for z in NF[f'psi{p-1}_exact']];g0=[F(0)];g=[F(1)]
    for k in range(15):
        exact=sum((a*q/F(2*j+i+1) for j,a in enumerate(g) for i,q in enumerate(psi)),F(0));_,mid,e=rows[k];assert abs(exact*S-mid)<=e
        v=[F(0)]*(k+2)
        for j,a in enumerate(g):v[j]+=(2*k*k+2*p*k+p-F(1,4))*a;v[j+1]-=LAM2*a
        for j,a in enumerate(g0):v[j]-=(k+p-1)**2*a
        g0,g=g,[a/(k+1)**2 for a in v]
    return [(v,e) for _,v,e in rows]
A=[(2*v,2*e) for v,e in transform(3)];G=[(LAM2*v,LAM2*e) for v,e in transform(4)]

def physical(n,a,b):
    s=2*b+1;q=b*(b+1);N=n+a;den=N*(N+s);c=F(n*(n+s),den);h=F(1,den);d=1-c
    J=min(n,K-1);v=e=0;sig=1 if b==0 else -1
    for k in range(J,-1,-1):
        # P= (1+qh)A + h D_b A + sigma h G, with 2den as common denominator.
        C0=2*den+2*q-k*k+(8-4*b)*k+9-2*b*(b-1);C1=(k+1)**2
        vv=C0*A[k][0]+C1*A[k+1][0]+2*sig*G[k][0];ee=abs(C0)*A[k][1]+C1*A[k+1][1]+2*G[k][1]
        if k<J:
            j=k;u=n*(n+s)-j*(j+s);w=den-j*(j+s);assert 0<=u<w
            vv=vv*w+2*den*u*v;ee=ee*w+2*den*u*e;v=(abs(vv)//(2*den*w))*(1 if vv>=0 else -1);e=cd(ee,2*den*w)+1
        else:v=(abs(vv)//(2*den))*(1 if vv>=0 else -1);e=cd(ee,2*den)+1
    # Physical angular action is d^3 P/(1+qh)^4; if n>K-1 use low-prefix tail.
    norm=d**3/(1+q*h)**4;mid=norm*F(v,S);err=norm*F(e,S)+(F(1,2**110) if n>J else F(0))
    assert c<=F('.56') and c>=F('.05')
    if n>J:
        assert (c/F('.99'))**(J+1)<=F(1,2**2670)
    print(json.dumps(dict(n=n,b=b,a=a,L=-math.log(float(d)),ideal_mid=float(mid),error_lt=float(err),positive_ideal=mid>err)),flush=True)
for n,a in [(1002,3000),(1200,2400),(1600,1600),(2870,998),(3500,1600)]:
    for b in [0,4]:physical(n,a,b)
print('IDEAL ONLY: signed coefficient transport passed; no realized-weight sign claim.')
