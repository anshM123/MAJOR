#!/usr/bin/env python3
"""Rigorous scalar recurrence + anisotropic, sequential real Chebyshev norming.
Uses only Python integer/Fraction arithmetic for proof decisions, and GMP scalar.cpp.
The point of sequential norming is that only one coordinate at a time is complex.
"""
from fractions import Fraction as F
from math import isqrt,prod,lcm
import json,argparse,subprocess,time,concurrent.futures,itertools,os
from pathlib import Path
S=10**80
HERE=Path(__file__).resolve().parent

def cd(x,y):
    assert x>=0 and y>0
    return (x+y-1)//y
def fl(x):x=F(x);return x.numerator*S//x.denominator
def up(x):x=F(x);return -((-x.numerator*S)//x.denominator)
def mu(a,b):return cd(a*b,S)
def dlo(a,b):return a*S//b
def du(a,b):return cd(a*S,b)
def normup(r,i):return isqrt(r*r+i*i)+1
def normlo(r,i):return isqrt(r*r+i*i)
def fmt(x,scale=S):return str(F(x,scale))
def approx(x,scale=S):return float(F(x,scale))
def cmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])

def parameters(c):
    return [F(c[x]) for x in ('d0','dd','h0','dh','nu0','dnu')]+[F(c.get('nuImag','0'))]

def cap_check(c):
    d,dd,h,dh,nu,dnu,ni=parameters(c)
    assert 0<d-dd<=d+dd<1 and h>=dh>=0 and min(dd,dnu)>=0
    assert nu-dnu>0 or c.get('nuAxis','real')=='imag'
    if c.get('fullTwist',False):
        assert h==dh
        root=F(isqrt(fl(1-d-dd)*S),S)
        required=(1 if c.get('coupled',False) else (d+dd)**2)/(998**2*(1+root)**2)
        assert h+dh>=required, ('insufficient twist cap', float(required))
        return float(required)
    return None

def major(c,p,b,radii,J=1500):
    """Bound |P|, |Q_unscaled| over complex affine disks with radii in normalized axes.
    Exact enclosure via center+variation for beta products, h products, d affine term.
    Scales S, no floating comparisons. Returned tail includes j>=K if K supplied.
    """
    d,dd,h,dh,nu,dnu,ni=parameters(c)
    d0=fl(d);h0=fl(h);nr=fl(nu);im=fl(ni)
    # centers exactly represented by scale S (all inputs decimal <=80 digits)
    assert F(d0,S)==d and F(h0,S)==h and F(nr,S)==nu and F(im,S)==ni
    dr,hr,vr=(up(x*y) for x,y in zip([dd,dh,dnu],radii))
    if c.get('coupled',False):
        # h=d^2*x, where h0/hr here initially represent x center/radius.
        xc,xr=h0,hr
        h0=cd(d0*d0*xc,S*S)
        hr=cd((d0+dr)**2*xr+xc*(2*d0*dr+dr*dr),S*S)+1
    dm=d0+dr;hm=h0+hr;s=2*b+1;q=b*(b+1)
    old=cur=S;rp=S;total=tail=0;minD=S
    for j in range(J):
        # beta center real and imag, and product coefficient centers with scale S^2
        re=(2*p-1-2*j)*S//2-im;beta=(re,nr)
        bb=normup(*beta)
        def prod_bound(k1,k2):
            x=(re+k1*S,nr);y=(re+k2*S,nr)
            z=cmul(x,y);dev=vr*(normup(*x)+normup(*y)+vr)
            return z,dev
        BE,BEv=prod_bound(0,-s);B,Bv=prod_bound(0,1)
        A=cmul(beta,((s-1)*S-2*re,-2*nr))
        Av=vr*(normup((s-1)*S-2*re,-2*nr)+2*bb+2*vr)
        # terms d+h A, h B and 1-h BE in scale S^3
        ac=(d0*S*S+h0*A[0],h0*A[1]); av=dr*S*S+hr*normup(*A)+hm*Av
        an=normup(*ac)+av
        bn=hm*(normup(*B)+Bv)
        dc=(S**3-h0*BE[0],-h0*BE[1]);dv=hr*normup(*BE)+hm*BEv
        den=normlo(*dc)-dv
        if den<=0:raise ValueError('D_j pole not excluded')
        minD=min(minD,den//(S*S))
        rn=normup((2*j+1)*S//2+im,nr)+vr
        rd=normlo((j+1)*S+2*im,2*nr)-2*vr
        if rd<=0:raise ValueError('frequency pole not excluded')
        rr=cd(rn*rn,(j+1)*rd)
        val=cur
        if p==3:
            A0,A0v=prod_bound(0,8-4*b);A1,A1v=prod_bound(0,9-4*b)
            C0=((2*q+9-2*b*(b-1))*S*S-A0[0],-A0[1])
            # C0 represents 2(q+A0), C1 represents 2 A1
            ca=normup(2*S**3+h0*C0[0],h0*C0[1])+hr*normup(*C0)+hm*A0v
            cb=hm*(normup(*A1)+A1v)
            val=cd(ca*cur+cd(cb*rp*old,S),2*S**3)
        total+=val
        if j>=int(c.get('K',J)):tail+=val
        if j+1<J:
            nxt=cd(rr*(an*cur+cd(bn*rp*old,S)),den*S)
            old,cur=cur,nxt;rp=rr
    return total,tail,minD

def major_axis(c,p,b,axis,rho):
    r=[F(1)]*3;r[axis]=(rho+1/rho)/2
    # ellipse of each variable lies in disk of radius (rho+rho^-1)/2.
    J=1500 if c.get("ellipseFull",False) else c.get("K",1500)
    from major_rect import major as rect
    vals=[]
    try:vals.append(major(c,p,b,r,J=J)[0])
    except ValueError:pass
    try:vals.append(rect(c,p,b,axis,rho,J=J)[0])
    except ValueError:pass
    if not vals:raise ValueError('no ellipse majorant')
    return min(vals)

def estimate_degree(M,rho,tol):
    # upper bound of Chebyshev truncation tail 2M*rho^-m/(rho-1)
    p,q=rho.numerator,rho.denominator;m=0
    top=2*M*q;bottom=(p-q)*up(tol)
    while top>bottom:top*=q;bottom*=p;m+=1
    err=cd(2*M*q**(m+1),(p-q)*p**m)
    assert err<=up(tol)
    return m,err

def plan(c,p,b):
    cap=cap_check(c);d,dd,h,dh,nu,dnu,ni=parameters(c)
    widths=[dd,dh,dnu];axes=[];tol=F(c.get('epsilon','0.000001'))
    # Each dimension may use a different ellipse. Candidate list selected for cost only.
    rhos=[F(x) for x in c.get('rhos',['1.5','2','3','4','6','8','12','16','24','32','48','64','128','256','512','1024','4096','16384','65536','1048576'])]
    for i,width in enumerate(widths):
        if not width:axes.append(dict(m=0,N=0,rho='0',eps=0,M=0));continue
        best=None
        for rho in rhos:
            if rho<=1:continue
            try:M=major_axis(c,p,b,i,rho)
            except ValueError:continue
            m,err=estimate_degree(M,rho,tol)
            N=max(2, (8*m+4)//5) # ceil(1.6*m), factor 2 incl. node rounding
            if best is None or (N,err)<(best['N'],best['eps']):best=dict(m=m,N=N,rho=str(rho),eps=err,M=M)
        if best is None:raise ValueError(f'axis {i}: no ellipse')
        axes.append(best)
    from major_rect import major as rect
    realvals=[]
    try:realvals.append(major(c,p,b,[F(1)]*3))
    except ValueError:pass
    try:realvals.append(rect(c,p,b))
    except ValueError:pass
    assert realvals, 'no real majorant'
    realM=min(r[0] for r in realvals);tail=min(r[1] for r in realvals);minD=min(r[2] for r in realvals)
    return axes,tail,realM,minD,cap

# Rigorous rational approximation to Lobatto nodes. 12-digit nodes are sufficient;
# Bernstein inequality is applied at the actual nodes, so no huge M*nodeError.
NODE=10**12;CS=10**50

def ataninv(k):
    ss=0;p=k;j=0
    while True:
        x=CS//(p*(2*j+1))
        if not x:break
        ss+=(-1 if j&1 else 1)*x;p*=k*k;j+=1
    assert j<100
    return ss
PI=16*ataninv(5)-4*ataninv(239)

def cosine_nodes(N):
    if N==0:return [F(0)]
    out=[]
    for k in range(N+1):
        x=PI*k//N;x2=x*x//CS;term=ss=CS
        for j in range(1,100):
            # toward zero, sign explicit
            term=(-1 if term>0 else 1)*(abs(term)*x2//(CS*(2*j-1)*(2*j)))
            ss+=term
            if not term:break
        assert j<99
        # error of cos approx <1e-42: Machin series errors <2000/CS,
        # Taylor terms incl. roundoff amplified by <= exp(pi) < 24.
        v=max(-NODE,min(NODE,(ss*NODE+CS//2)//CS));out.append(F(v,NODE))
    out[0]=F(1);out[-1]=F(-1)
    assert all(out[i]>out[i+1] for i in range(N))
    return out

def proof_node_factor(m,N):
    if not N:return
    # |node - exact cosine| <= 1e-12. In angle, <= pi sqrt(err/2) < 3e-6
    # using |cos a-cos b|>=2 sin^2((a-b)/2)>=2 (a-b)^2/pi^2.
    # nearest angle error <= pi/(2N)+3e-6. Require product degree < pi/3.
    assert F(22,7)*m/(2*N)+F(3,10**6)*m<F(314159,100000)/3

def node_rows(c,p,b,axes):
    d,dd,h,dh,nu,dnu,ni=parameters(c)
    xx=[cosine_nodes(a['N']) for a in axes]
    for a in axes:proof_node_factor(a['m'],a['N'])
    # cache frequency-dependent exact coefficients across all (d,h) nodes
    for z in xx[2]:
        nr,im=(nu+dnu*z,ni) if c.get('nuAxis','real')=='real' else (nu,ni+dnu*z)
        w=lcm(nr.denominator,im.denominator);nr0=int(nr*w);ni0=int(im*w)
        for x in xx[0]:
            dn=d+dd*x
            for y in xx[1]:
                hn=h+dh*y
                if c.get('coupled',False):hn*=dn*dn
                yield f'{p} {b} {dn.numerator} {dn.denominator} {hn.numerator} {hn.denominator} {nr0} {ni0} {w}\n'

def eval_chunk(rows,digits,K):
    st=time.perf_counter();proc=subprocess.run([str(HERE/'scalar'),str(digits),str(K)],input=''.join(rows),text=True,capture_output=True,check=True)
    p=q=e=0;count=0
    for line in proc.stdout.splitlines():
        P,Q,Pl,Ql,E=map(int,line.split());p=max(p,P);q=max(q,Q);e=max(e,E);count+=1
    assert count==len(rows)
    return p,q,e,count,time.perf_counter()-st

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--preflight',action='store_true');ap.add_argument('--jobs',type=int,default=1);ap.add_argument('--max-nodes',type=int,default=50000);ap.add_argument('--output',default=None);args=ap.parse_args()
    out=[]
    for c in json.load(open(args.input)):
      for p,b in c.get('channels',[[3,4]]):
        st=time.perf_counter();axes,tail,M,Dgap,cap=plan(c,p,b);K=c.get('K',1500);digits=c.get('digits',245)
        n=prod(a['N']+1 for a in axes);dmin=parameters(c)[0]-parameters(c)[1]
        row=dict(name=c['name'],p=p,b=b,config=c,axes=[dict(m=a['m'],N=a['N'],rho=a['rho'],epsilon=fmt(a['eps']),majorant=fmt(a['M'])) for a in axes],nodes=n,tail_bound=fmt(tail),h_required=cap,preflight_seconds=time.perf_counter()-st)
        if args.preflight or n>args.max_nodes:
            row['status']='PREFLIGHT' if args.preflight else 'SKIP_NODE_LIMIT';print(json.dumps(row),flush=True);out.append(row);continue
        rows=node_rows(c,p,b,axes);chunks=iter(lambda:list(itertools.islice(rows,128)),[]);t=time.perf_counter();P=Q=E=count=0;cpu=0
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            # bounded number of inflight batches
            pending=set()
            for chunk in chunks:
                pending.add(pool.submit(eval_chunk,chunk,digits,K))
                if len(pending)>=2*args.jobs:
                    done,pending=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
                    for f in done:
                        pp,qq,ee,k,secs=f.result();P=max(P,pp);Q=max(Q,qq);E=max(E,ee);count+=k;cpu+=secs
            for f in concurrent.futures.as_completed(pending):
                pp,qq,ee,k,secs=f.result();P=max(P,pp);Q=max(Q,qq);E=max(E,ee);count+=k;cpu+=secs
        assert count==n
        scale=10**digits;pb=F(P,scale);qb=F(Q,scale)
        # scalar returns Q=Q_unscaled/d. For norming use Q_unscaled, avoid extra pole in d.
        # Multiplying max Q by dmax is conservative, then divide by dmin after extension.
        qb*=parameters(c)[0]+parameters(c)[1]
        for a in axes:
            if a['N']:
                ep=F(a['eps'],S);pb=2*pb+3*ep;qb=2*qb+3*ep
        # Ellipse majorants bound the full sum; if K<1500 the sampled polynomial is partial.
        # The majorant is valid for every partial sum as well, and add omitted tail only now.
        pb+=F(tail,S);qb=(qb+F(tail,S))/dmin
        row.update(status='ENCLOSED',P_bound= str(pb),Q_bound=str(qb),P_approx=float(pb),Q_approx=float(qb),sample_P=float(F(P,scale)),sample_Q=float(F(Q,scale)),roundoff_units=str(E),scale_digits=digits,wall_seconds=time.perf_counter()-t,worker_seconds=cpu)
        if c.get('maxTransfer') is not None:row['status']='PASS' if max(pb,qb)<F(c['maxTransfer']) else 'FAIL_TRANSFER'
        print(json.dumps(row),flush=True);out.append(row)
    if args.output:Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
