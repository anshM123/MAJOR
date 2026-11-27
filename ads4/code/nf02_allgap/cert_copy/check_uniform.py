#!/usr/bin/env python3
"""Exact rational checks. Ideal shifted spectrum only, no cutoff/interchange claim."""
from fractions import Fraction as F
from math import comb
import json, sys, time
from pathlib import Path
D=997002
k0=F('27.48509817')
c2=[F(x)/k0 for x in ['0','1297957/250','-14066127/1000','2586071/250','36769/250','-2149393/250','42382551/250','-872871009/500','10034471703/1000','-32551961929/1000','62365960687/1000','-9052369139/125','2009389439/40','-19196103417/1000','1557521209/500']]
c3=[F(x)/k0 for x in ['0','1391927/1000','1/1000','1228529/40','-15055687/1000','555153477/1000','-2123008997/250','19839838023/500','-97457517053/1000','28931472071/200','-33422106203/250','73803807381/1000','-5148431833/250','964822933/1000','553835527/1000']]
def mom(c,j): return sum(x/F(i+j+1) for i,x in enumerate(c))
def bern_1(a,n=None):
    if n is None:n=len(a)-1
    return [sum(a[j]*F(comb(i,j),comb(n,j)) for j in range(min(i,len(a)-1)+1)) for i in range(n+1)]
def split1(a):
    L=[a[0]];R=[a[-1]]
    while len(a)>1:
        a=[(x+y)/2 for x,y in zip(a,a[1:])];L.append(a[0]);R.append(a[-1])
    return L,R[::-1]
def exact_light():
    K=[-2*mom(c2,6)+2*mom(c3,4), F(245,2)*mom(c2,4)-35*mom(c3,2), -F(9619,8)*mom(c2,2)+F(569,8)*mom(c3,0), F(51575,32)*mom(c2,0)]
    assert K[1]<0 and K[2]<0 and K[3]>0
    assert K[0]+K[1]/D+K[2]/D**2>F('0.99996')
    assert K[1]+3*K[3]/D**2<0
    g=c2[2:]; coeff=bern_1(g); cells=[coeff]
    for _ in range(4):cells=[p for c in cells for p in split1(c)]
    gmin=min(min(c) for c in cells);gmax=max(max(c) for c in cells)
    assert gmax<0
    gabs=max(abs(gmin),abs(gmax)); finite=mom(g,0)
    print('LIGHT', 'K',list(map(str,K)),'K@1000',float(sum(x/F(D)**i for i,x in enumerate(K))), 'slope',float(c2[1]), 'finite',float(finite),'g range',float(gmin),float(gmax),flush=True)
    eps=F(3,1996)
    err=16*(c2[1]*eps**2/2+F(11,7)*eps*gabs)
    print('16 L <=',float(16*c2[1]),'*log(2 Lambda/3) +',float(16*finite+err),'error bound',float(err),flush=True)
    return {'K_coefficients':list(map(str,K)),'K_min':str(F('0.99996')),'slope':str(c2[1]),'finite_part':str(finite),'g_min':str(gmin),'g_max':str(gmax),'E_constant_upper':str(16*finite+err)}

def gpolys(p,N):
    prev=[];curr=[F(1)];out=[curr]
    for k in range(N):
        nxt=[F(0)]*(k+2)
        for i,v in enumerate(curr):
            nxt[i]+=(8*k*k+8*p*k+4*p-1)*v
            nxt[i+1]-=4*v
        for i,v in enumerate(prev):nxt[i]-=4*(k-1+p)**2*v
        nxt=[v/F(4*(k+1)**2) for v in nxt]
        prev,curr=curr,nxt;out.append(curr)
    return out

def scaled(p,c,N,fac):
    moms=[fac*mom(c,2*i) for i in range(N+1)]
    return [[v*moms[i]/D**(k-i) for i,v in reversed(list(enumerate(g)))] for k,g in enumerate(gpolys(p,N))]

def numerator(J,b,A,G):
    n=J-b;C=J*(J+1);sig=1 if b==0 else -1
    M=n+1
    P=[[F(0)]*(M+1) for _ in range(M+1)]
    B=[C-(b+j)*(b+j+1) for j in range(n)]
    suf=[F(1)];pre=1
    for v in B:pre*=v
    for k in range(n,-1,-1):
        gamma=F(-k*k+(8-4*b)*k+9-2*b*(b-1),2)
        a=A[k];v=[F(0)]*(k+2)
        for i,x in enumerate(a):v[i+1]+=(C+gamma)*x/D
        for i,x in enumerate(A[k+1]):v[i]+=F((k+1)**2,2)*x
        for i,x in enumerate(G[k]):v[i]+=sig*x
        for j,s in enumerate(suf):
            s*=pre
            for i,x in enumerate(a):P[i+j][k+j]+=s*x
            for i,x in enumerate(v):P[i+j][k+j+1]+=s*x
        if k:
            pre//=B[k-1]
            a1=F(B[k-1],D)
            suf=[suf[0]]+[suf[i]+a1*suf[i-1] for i in range(1,len(suf))]+[a1*suf[-1]]
    return P

def bern_2(P):
    n=len(P)-1; w=[[F(comb(i,j),comb(n,j)) for j in range(i+1)] for i in range(n+1)]
    Q=[[sum(w[i][k]*P[k][j] for k in range(i+1) if P[k][j]) for j in range(n+1)] for i in range(n+1)]
    return [[sum(w[j][k]*Q[i][k] for k in range(j+1) if Q[i][k]) for j in range(n+1)] for i in range(n+1)]

def cert(B,maxdepth=20):
    stack=[(B,0)]; leaves=0;minval=None
    while stack:
        a,d=stack.pop();lo=min(min(row) for row in a)
        if lo>0:
            leaves+=1;minval=lo if minval is None else min(minval,lo);continue
        if d>=maxdepth: return False,leaves,lo
        if any(a[i][j]<=0 for i in (0,-1) for j in (0,-1)):
            return False,leaves,lo
        vu=max(abs(row[j+1]-row[j]) for row in a for j in range(len(row)-1))
        vq=max(abs(a[i+1][j]-a[i][j]) for i in range(len(a)-1) for j in range(len(a[0])))
        trans=vq>vu
        if trans:a=list(map(list,zip(*a)))
        z=[split1(row) for row in a]
        left=[v[0] for v in z];right=[v[1] for v in z]
        if trans:left=list(map(list,zip(*left)));right=list(map(list,zip(*right)))
        stack.append((right,d+1));stack.append((left,d+1))
    return True,leaves,minval

def main(N):
    light=exact_light(); A=scaled(3,c2,N+1,2);G=scaled(4,c3,N,1)
    rows=[];t0=time.time()
    for J in range(N+1):
        for b in (0,4):
            if (b==0 and J%2) or J<b:continue
            start=time.time();P=numerator(J,b,A,G);B=bern_2(P)
            ok,leaves,mn=cert(B)
            print('J',J,'b',b,'PASS' if ok else 'FAIL','leaves',leaves,'min',float(mn),'time',round(time.time()-start,2),flush=True)
            rows.append({'J':J,'b':b,'pass':ok,'leaves':leaves,'min_Bernstein':str(mn),'min_approx':float(mn)})
        Path(__file__).with_name('results.json').write_text(json.dumps({'light':light,'rows':rows,'seconds':time.time()-t0},indent=2))
    print('TOTAL',len(rows),'passes',sum(x['pass'] for x in rows),'seconds',time.time()-t0,flush=True)
if __name__=='__main__':main(int(sys.argv[1]) if len(sys.argv)>1 else 24)
