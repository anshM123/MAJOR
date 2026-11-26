#!/usr/bin/env python3
"""Exact physical mass-cap transform. q=998/(Delta_g-2), v=m_gap(J)^2/m^2."""
from check_fast import *
def divide_linear(a,c,d):
    # exact division ascending, by c+d*x
    r=[];prev=0
    for x in a[:-1]:
        v,e=divmod(x-d*prev,c);assert e==0;r.append(v);prev=v
    assert a[-1]==d*prev
    return r

def mul_linear(a,c,d):
    return [c*a[0]]+[c*a[i]+d*a[i-1] for i in range(1,len(a))]+[d*a[-1]]

def transform(P,den,J):
    M=len(P)-1;d=2*J+1;C=998
    S=[comb(M,k)*C**(M-k)*d**k for k in range(M+1)]
    R=[[0]*(M+1) for _ in range(2*M+1)]
    pw=[999**i for i in range(M+1)]
    for j in range(M+1):
        a=S
        for i in range(j+1):
            x=P[i][j]*pw[i]
            if x:
                for k,y in enumerate(a):R[k+2*i][j]+=x*y
            if i<j:a=divide_linear(a,C,1)
        assert all(P[i][j]==0 for i in range(j+1,M+1))
        if j<M:S=mul_linear(divide_linear(S,C,d),C,1)
    return R,den

def bern_rect(P,den):
    n=len(P)-1;m=len(P[0])-1
    def W(d):
        cs=[comb(d,j) for j in range(d+1)];lc=lcm(*cs)
        return [[comb(i,j)*(lc//cs[j]) for j in range(i+1)] for i in range(d+1)],lc
    w,lc=W(n);v,lc2=W(m)
    Q=[[sum(w[i][k]*P[k][j] for k in range(i+1) if P[k][j]) for j in range(m+1)] for i in range(n+1)]
    R=[[sum(v[j][k]*Q[i][k] for k in range(j+1) if Q[i][k]) for j in range(m+1)] for i in range(n+1)]
    return R,den*lc*lc2

def cert_rect(B,den,maxdepth=28):
    n=len(B)-1;m=len(B[0])-1
    stack=[(B,0,0)];leaves=0;minval=None
    while stack:
        a,d,shift=stack.pop();lo=min(min(row) for row in a)
        if lo>0:
            l=F(lo,den<<shift);leaves+=1;minval=l if minval is None else min(minval,l);continue
        if d>=maxdepth or any(a[i][j]<=0 for i in (0,-1) for j in (0,-1)):
            return False,leaves,F(lo,den<<shift)
        vu=m*max(abs(row[j+1]-row[j]) for row in a for j in range(m))
        vq=n*max(abs(a[i+1][j]-a[i][j]) for i in range(n) for j in range(m+1))
        trans=vq>vu
        if trans:a=list(map(list,zip(*a)))
        z=[split_int(row) for row in a];left=[v[0] for v in z];right=[v[1] for v in z]
        if trans:left=list(map(list,zip(*left)));right=list(map(list,zip(*right)))
        s=shift+(n if trans else m)
        stack.append((right,d+1,s));stack.append((left,d+1,s))
    return True,leaves,minval

def ev(P,x,y):return sum(x**i*sum(v*y**j for j,v in enumerate(a) if v) for i,a in enumerate(P))
def main(N,start=0):
    A=scaled(3,c2,N+1,2);G=scaled(4,c3,N,1)
    rows=[];t0=time.time();out=Path(__file__).with_name(f'physical_{start}_{N}.json')
    for J in range(start,N+1):
        for b in (0,4):
            if (b==0 and J%2) or J<b:continue
            st=time.time();P,d=num_int(J,b,A,G);R,e=transform(P,d,J);M=len(P)-1
            # exact cross-check of transformation on 3 independent rational states
            for x,y in [(F(0),F(1)),(F(1),F(1)),(F(3,7),F(4,5))]:
                q=999*x*x/(998+x);u=y*(998+x)/(998+(2*J+1)*x)
                assert ev(R,x,y)==ev(P,q,u)*(998+(2*J+1)*x)**M
            B,z=bern_rect(R,e);ok,leaves,mn=cert_rect(B,z)
            # divide denominator factors by worst-case positive bound (not needed for sign).
            bound=(998+2*J+1)**M * (1+F(J*(J+1),D))**(M+3)
            hn=mn/bound
            print('J',J,'b',b,'PASS' if ok else 'FAIL','leaves',leaves,'lower',float(hn),'time',round(time.time()-st,2),flush=True)
            rows.append({'J':J,'b':b,'pass':ok,'leaves':leaves,'min_Bernstein':str(mn),'normalized_bound':str(hn),'normalized_approx':float(hn)})
        out.write_text(json.dumps({'rows':rows,'seconds':time.time()-t0},indent=2))
    print('TOTAL',len(rows),'passes',sum(x['pass'] for x in rows),'seconds',time.time()-t0,flush=True)
if __name__=='__main__':main(int(sys.argv[1]) if len(sys.argv)>1 else 64,int(sys.argv[2]) if len(sys.argv)>2 else 0)
