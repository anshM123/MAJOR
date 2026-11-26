"""Anisotropic rectangular majorant for the real norming theorem.
At most one variable is complex. Other variables stay on their real segments.
All interval endpoints and coefficient operations are outward fixed-point integers.
"""
from check_norm import S,cd,fl,up,mu,parameters,normup,normlo
from fractions import Fraction as F
from math import isqrt

def muliv(a,b):
    z=[x*y for x in a for y in b];return min(z)//S,-(-max(z)//S)
def sqiv(a):
    l,u=a;z=max(l*l,u*u);m=0 if l<=0<=u else min(l*l,u*u)
    return m//S,cd(z,S)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def scale(a,k):return (a[0]*k,a[1]*k) if k>=0 else (a[1]*k,a[0]*k)
def mabs(a):return max(abs(a[0]),abs(a[1]))
def minabs(a):return 0 if a[0]<=0<=a[1] else min(abs(a[0]),abs(a[1]))
def ub(a,b):return normup(mabs(a),mabs(b))
def product_complex(ar,ai,br,bi):return add(muliv(ar,br),neg(muliv(ai,bi))),add(muliv(ar,bi),muliv(ai,br))
def max_affine_prod(ar,ai,br,bi):return ub(ar,ai)*ub(br,bi)//S+1

def major(c,p,b,axis=None,rho=F(1),J=1500):
    d,dd,h,dh,v,dv,ni=parameters(c)
    # Exact rectangle enclosing an axis ellipse; other axes remain real.
    aa=(rho+1/rho)/2 if axis is not None else F(1)
    bb=(rho-1/rho)/2 if axis is not None else F(0)
    w=[dd,dh,dv];ctr=[d,h,v];rad=[w[i]*(aa if i==axis else 1) for i in range(3)];imag=[w[i]*bb if i==axis else F(0) for i in range(3)]
    dr=(fl(d-rad[0]),up(d+rad[0]));di=(-up(imag[0]),up(imag[0]))
    hr=(fl(h-rad[1]),up(h+rad[1]));hi=(-up(imag[1]),up(imag[1]))
    if c.get('coupled',False):
        d2r=add(sqiv(dr),neg(sqiv(di)));d2i=scale(muliv(dr,di),2)
        hr,hi=product_complex(d2r,d2i,hr,hi)
    if c.get('nuAxis','real')=='real':tr=(fl(v-rad[2]),up(v+rad[2]));yr=(fl(ni-imag[2]),up(ni+imag[2]))
    else:tr=(fl(v-imag[2]),up(v+imag[2]));yr=(fl(ni-rad[2]),up(ni+rad[2]))
    ts=sqiv(tr);hm=ub(hr,hi);s=2*b+1;q=b*(b+1);old=cur=S;rp=S;total=tail=0;minD=S
    for j in range(J):
        br=(F(2*p-1-2*j,2)*S-yr[1],F(2*p-1-2*j,2)*S-yr[0]);br=tuple(map(int,br));
        # beta*(beta-s) in a rectangle. Use real quadratic with exact min/max via endpoints + vertex.
        def beta_product(k):
            l,u=br;K=k*S;vs=[l*(l+K),u*(u+K)]
            low=min(vs);high=max(vs)
            if 2*l<=-K<=2*u:low=min(low,-K*K//4)
            re=(low//S-ts[1],-(-high//S)-ts[0]);im=muliv(tr,(2*l+K,2*u+K));return re,im
        BE,BEi=beta_product(-s);B,Bi=beta_product(1)
        # beta(s-1-2 beta)=-2 beta(beta-(s-1)/2)
        A,Ai=beta_product(-(s-1)//2);A=scale(A,-2);Ai=scale(Ai,-2)
        DRe=add((S,S),neg(add(muliv(hr,BE),neg(muliv(hi,BEi)))))
        if DRe[0]<=0:raise ValueError('rect pole not excluded')
        den=DRe[0];minD=min(minD,den)
        Ar,Aii=product_complex(hr,hi,A,Ai);an=ub(add(dr,Ar),add(di,Aii));bn=mu(hm,ub(B,Bi))
        # |r_j| upper envelope quasi-convex as function of t^2, even for y rectangle.
        numx=mabs(((2*j+1)*S//2+yr[0],(2*j+1)*S//2+yr[1]))
        denx=minabs(((j+1)*S+2*yr[0],(j+1)*S+2*yr[1]));rr=0
        for tt in ts:
            denr=isqrt(denx*denx+4*tt*S)
            if denr<=0:raise ValueError('rect frequency pole')
            rr=max(rr,cd(numx*numx+tt*S,(j+1)*denr))
        val=cur
        if p==3:
            c0,c0i=beta_product(8-4*b);c1,c1i=beta_product(9-4*b)
            c0=add(((2*q+9-2*b*(b-1))*S,)*2,neg(c0));c0i=neg(c0i)
            cr,ci=product_complex(hr,hi,c0,c0i);ca=ub(add((2*S,2*S),cr),ci);cb=mu(hm,ub(c1,c1i))
            val=cd(ca*cur+cd(cb*rp*old,S),2*S)
        total+=val
        if j>=int(c.get('K',J)):tail+=val
        if j+1<J:
            nxt=cd(rr*(an*cur+cd(bn*rp*old,S)),den*S)
            old,cur=cur,nxt;rp=rr
    return total,tail,minD
