#!/usr/bin/env python3
"""Truncated exact integer tensor Bernstein certificate at all gaps and unbounded J.
Variables q=997002/Lambda^2, V=v/(4/131^2), X=x/64. All in [0,1].
Truncation in q,V ONLY; positive coefficient majorants bound discarded terms.
"""
from check_fast import *
import numpy as np
import argparse

K=60;X=64;MQ=8;MV=16;VD=131**2;v0=F(4,VD);t0=F(1,D);L=D*VD;r=F(1,20)

def majorant(b,A,G,Rq,Rv):
    U=[1+F((2*(b+j)+1)**2,4)*v0*Rv for j in range(K)]
    suf=[F(1)]*(K+1)
    for j in range(K-1,-1,-1):suf[j]=suf[j+1]*(1+t0*Rq*X*U[j])
    pre=F(1);m=F(0)
    for k in range(K+1):
        gamma=F(-k*k+(8-4*b)*k+9-2*b*(b-1),2)
        a=sum(abs(c)*Rq**i for i,c in enumerate(A[k]));aa=sum(abs(c)*Rq**i for i,c in enumerate(A[k+1]));g=sum(abs(c)*Rq**i for i,c in enumerate(G[k]))
        m+=X**k*pre*suf[k]*((1+t0*Rq*X*(1+v0*Rv*(abs(gamma)+F(1,4))))*a+X*v0*Rv*(F((k+1)**2,2)*aa+g))
        if k<K:pre*=U[k]
    return m

def tail_g(p,c,mult):
    gp=gpolys(p,K+1)
    W=[X**k*sum(abs(z)*t0**(k-i) for i,z in enumerate(gp[k])) for k in (K,K+1)]
    C=max(W[0],W[1]/r)
    assert 2*(K+1)**2+2*p*(K+1)+F(p)-F(1,4)<=3*(K+2)**2
    assert (K+p)**2<=2*(K+2)**2
    alpha=X*(3*t0+F(1,(K+2)**2));beta=2*X*X*t0*t0
    assert alpha/r+beta/r**2<1
    bc=bern_1(c);cells=[bc]
    for _ in range(4):cells=[a for z in cells for a in split1(z)]
    norm=mult*max(abs(x) for z in cells for x in z)
    return C*norm

def sqsum(a):return a*a*r/(1-r)+2*a*r/(1-r)**2+r*(1+r)/(1-r)**3

def numerator_tensor(b,A,G):
    sig=1 if b==0 else -1
    # common denominator of F_k coefficients, before multiplying by suffix scale.
    den=lcm(*(x.denominator for a in A+G for x in a))*2*D*VD
    AA=[[int(z*den) for z in a] for a in A];GG=[[int(z*den) for z in a] for a in G]
    # suffix monomials have x-degree = q-degree.
    S=np.zeros((MQ+1,MV+1),dtype=object);S[0,0]=1
    # P[q,V,x], exact coefficients; discarded monomials in q or V cannot affect retained ones.
    P=np.zeros((MQ+1,MV+1,K+2),dtype=object)
    for k in range(K,-1,-1):
        ell=(2*(b+k)+1)**2
        if k<K:
            T=S*L;T[1:,:]+=X*VD*S[:-1,:];T[1:,1:]-=X*ell*S[:-1,:-1];S=T
            T=np.zeros_like(P);T[:,:,1:]+=X*L*P[:,:,:-1];T[:,1:,1:]-=X*D*ell*P[:,:-1,:-1];P=T
        gamma2=-k*k+(8-4*b)*k+9-2*b*(b-1)
        # bracket coefficients, arrays for x^0 V^0, x^1 V^0, x^1 V^1.
        for i in range(MQ+1):
            a=AA[k][i] if i<len(AA[k]) else 0
            aa=AA[k+1][i] if i<len(AA[k+1]) else 0
            g=GG[k][i] if i<len(GG[k]) else 0
            if a:
                for j in range(MQ+1-i):P[i+j,:,j]+=a*S[j,:]
            if i<MQ and a:
                x0,rem=divmod(X*a,D);assert rem==0
                xv,rem=divmod(X*(2*gamma2-1)*a,D*VD);assert rem==0
                for j in range(MQ-i):
                    P[i+1+j,:,j+1]+=x0*S[j,:]
                    P[i+1+j,1:,j+1]+=xv*S[j,:-1]
            xv,rem=divmod(X*(2*(k+1)**2*aa+4*sig*g),VD);assert rem==0
            if xv:
                for j in range(MQ+1-i):P[i+j,1:,j+1]+=xv*S[j,:-1]
    return P,den*L**K

def tensor_bern(P,den):
    for axis in range(3):
        Q=np.moveaxis(P,axis,0);n=Q.shape[0]-1;cs=[comb(n,j) for j in range(n+1)];lc=lcm(*cs)
        w=np.array([[comb(i,j)*(lc//cs[j]) if j<=i else 0 for j in range(n+1)] for i in range(n+1)],dtype=object)
        T=(w@Q.reshape(n+1,-1)).reshape(Q.shape);P=np.moveaxis(T,0,axis);den*=lc
    return P,den

def cert_x(P,den,err):
    # only x subdivisions, full q,V range retained. Exact integer de Casteljau.
    n=P.shape[2]-1;stack=[(P.reshape(-1,n+1).tolist(),0)];leaves=0;minval=None
    while stack:
        A,d=stack.pop();lo=min(min(z) for z in A);dd=den<<(n*d)
        if lo*err.denominator>dd*err.numerator:
            val=F(lo,dd);minval=val if minval is None else min(minval,val);leaves+=1;continue
        if d>=16:return False,leaves,F(lo,dd)
        z=[split_int(row) for row in A];stack.extend([([v[1] for v in z],d+1),([v[0] for v in z],d+1)])
    return True,leaves,minval

def main():
    st=time.time();A=scaled(3,c2,K+1,2);G=scaled(4,c3,K,1);TA=tail_g(3,c2,2);TG=tail_g(4,c3,1)
    # Even if K exceeds a physical termination index, all later products contain the zero numerator factor.
    assert 1+t0*X*(1-F((2*(4+K-1)+1)**2,VD))>0
    assert X*v0<=1
    tail=TA*r/(1-r)+v0*TA*r/2*sqsum(K+1)+X*v0*t0*TA*sqsum(K+10)+X*v0*TG*r/(1-r)
    rows=[]
    for b in (0,4):
        e=majorant(b,A,G,100,1)/100**(MQ+1)+majorant(b,A,G,1,16)/16**(MV+1)
        print('channel',b,'poly truncation',float(e),'ladder tail',float(tail),flush=True)
        P,d=numerator_tensor(b,A,G)
        # Exact check at q=0,V=0, where numerator is just the leading Bessel prefix.
        assert all(F(P[0,0,k],d)==A[k][0]*X**k for k in range(K+1));assert P[0,0,K+1]==0
        print('numerator built',round(time.time()-st,2),flush=True)
        B,den=tensor_bern(P,d);print('Bernstein built',round(time.time()-st,2),flush=True)
        ok,leaves,lo=cert_x(B,den,e+tail*(1+t0*X)**(K+1))
        H=((lo-e)/(1+t0*X)**(K+1)-tail)/(1+t0*X)**3
        print('b',b,'PASS' if ok else 'FAIL','leaves',leaves,'poly lower',float(lo),'H lower',float(H),'seconds',round(time.time()-st,2),flush=True)
        rows.append({'b':b,'pass':ok,'leaves':leaves,'poly_error':str(e),'tail':str(tail),'poly_lower':str(lo),'H_lower':str(H),'H_lower_approx':float(H)})
    Path(__file__).with_name(f'impact_tensor_X{X}_J{(int(VD**0.5)-1)//2}_results.json').write_text(json.dumps({'K':K,'X':X,'MQ':MQ,'MV':MV,'J_min':(int(VD**0.5)-1)//2,'rows':rows},indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--X',type=int,default=64);ap.add_argument('--K',type=int,default=60);ap.add_argument('--mq',type=int,default=8);ap.add_argument('--mv',type=int,default=16);ap.add_argument('--Jmin',type=int,default=65);ap.add_argument('--r',default='1/20');args=ap.parse_args()
    X=args.X;K=args.K;MQ=args.mq;MV=args.mv;VD=(2*args.Jmin+1)**2;v0=F(4,VD);L=D*VD;r=F(args.r);main()
