#!/usr/bin/env python3
"""Faster integer Bernstein sign certificate for ideal angular kernels."""
from check_uniform import *
from math import lcm

def split_int(a):
    n=len(a)-1
    L=[a[0]<<n];R=[a[-1]<<n]
    for i in range(1,n+1):
        a=[x+y for x,y in zip(a,a[1:])];L.append(a[0]<<(n-i));R.append(a[-1]<<(n-i))
    return L,R[::-1]

def cert_int(B,den,maxdepth=20):
    stack=[(B,0)];leaves=0;minval=None;n=len(B)-1
    while stack:
        a,d=stack.pop();lo=min(min(row) for row in a)
        if lo>0:
            l=F(lo,den<<(n*d));leaves+=1;minval=l if minval is None else min(minval,l);continue
        if d>=maxdepth or any(a[i][j]<=0 for i in (0,-1) for j in (0,-1)):
            return False,leaves,F(lo,den<<(n*d))
        vu=max(abs(row[j+1]-row[j]) for row in a for j in range(n))
        vq=max(abs(a[i+1][j]-a[i][j]) for i in range(n) for j in range(n+1))
        trans=vq>vu
        if trans:a=list(map(list,zip(*a)))
        z=[split_int(row) for row in a];left=[v[0] for v in z];right=[v[1] for v in z]
        if trans:left=list(map(list,zip(*left)));right=list(map(list,zip(*right)))
        stack.append((right,d+1));stack.append((left,d+1))
    return True,leaves,minval

def num_int(J,b,A,G):
    n=J-b;C=J*(J+1);sig=1 if b==0 else -1;M=n+1
    # M0 is common denominator for bracket; AA, GG exact integer representations at scale M0.
    den=lcm(*(x.denominator for a in A[:n+2]+G[:n+1] for x in a))*2*D
    AA=[[int(x*den) for x in a] for a in A[:n+2]];GG=[[int(x*den) for x in a] for a in G[:n+1]]
    P=[[0]*(M+1) for _ in range(M+1)]
    B=[C-(b+j)*(b+j+1) for j in range(n)]
    suf=[1];pre=1;dp=[D**(n-j) for j in range(n+1)]
    for v in B:pre*=v
    for k in range(n,-1,-1):
        gamma2=-k*k+(8-4*b)*k+9-2*b*(b-1)
        a=AA[k];v=[0]*(k+2)
        for i,x in enumerate(a):
            xx,r=divmod((2*C+gamma2)*x,2*D);assert r==0;v[i+1]+=xx
        for i,x in enumerate(AA[k+1]):
            xx,r=divmod((k+1)**2*x,2);assert r==0;v[i]+=xx
        for i,x in enumerate(GG[k]):v[i]+=sig*x
        for j,s in enumerate(suf):
            s*=pre*dp[j]
            for i,x in enumerate(a):P[i+j][k+j]+=s*x
            for i,x in enumerate(v):P[i+j][k+j+1]+=s*x
        if k:
            pre//=B[k-1];a1=B[k-1]
            suf=[suf[0]]+[suf[i]+a1*suf[i-1] for i in range(1,len(suf))]+[a1*suf[-1]]
    return P,den*D**n

def bern_int(P,den):
    n=len(P)-1;cs=[comb(n,j) for j in range(n+1)];lc=lcm(*cs);w=[[comb(i,j)*(lc//cs[j]) for j in range(i+1)] for i in range(n+1)]
    Q=[[sum(w[i][k]*P[k][j] for k in range(i+1) if P[k][j]) for j in range(n+1)] for i in range(n+1)]
    R=[[sum(w[j][k]*Q[i][k] for k in range(j+1) if Q[i][k]) for j in range(n+1)] for i in range(n+1)]
    return R,den*lc**2

def main(N,start=0):
    light=exact_light(); A=scaled(3,c2,N+1,2);G=scaled(4,c3,N,1)
    rows=[];t0=time.time()
    for J in range(start,N+1):
        for b in (0,4):
            if (b==0 and J%2) or J<b:continue
            starttime=time.time();P,d=num_int(J,b,A,G);B,e=bern_int(P,d);ok,leaves,mn=cert_int(B,e)
            if J<=12:
                ref=numerator(J,b,A,G);assert all(F(x,d)==y for row,r in zip(P,ref) for x,y in zip(row,r))
                BB=bern_2(ref);assert all(F(x,e)==y for row,r in zip(B,BB) for x,y in zip(row,r))
            print('J',J,'b',b,'PASS' if ok else 'FAIL','leaves',leaves,'min',float(mn),'time',round(time.time()-starttime,2),flush=True)
            rows.append({'J':J,'b':b,'pass':ok,'leaves':leaves,'min_Bernstein':str(mn),'min_approx':float(mn)})
        Path(__file__).with_name('fast_results.json').write_text(json.dumps({'light':light,'rows':rows,'seconds':time.time()-t0},indent=2))
    print('TOTAL',len(rows),'passes',sum(x['pass'] for x in rows),'seconds',time.time()-t0,flush=True)
if __name__=='__main__':main(int(sys.argv[1]) if len(sys.argv)>1 else 64,int(sys.argv[2]) if len(sys.argv)>2 else 0)
