#!/usr/bin/env python3
"""Directed near-real comparison and direct signed/absolute contour bounds.
Uses rational right edge 998.5, with an explicit extra real segment beyond Lambda.
The ideal-weight comparison is to the *undressed* infinite-twist reference.
"""
from fractions import Fraction as F
from math import ceil, comb
from pathlib import Path
import json,sys,time
from base_bounds import hi,lo,sqrtup,sqrtlo,logu,logl,hcap,base_bound
from check_lowfreq import exp_u
from check_norm import major, S as NS, parameters, cap_check
from assemble_top import psi_sup, B_inf, pow_upper, LAML,LAMU,NF
from assemble import B_real
S=10**40; J=1500; RIGHT=F('998.5'); DELTA=F('0.25')

def cd(n,d): assert n>=0 and d>0; return (n+d-1)//d

def fx(x): return cd(F(x).numerator*S,F(x).denominator)
def mu(a,b): return cd(a*b,S)
def div(a,b): return cd(a*S,b)
def f(x): return F(x,S)
def step_coeff(g,h,dl,du,s,slo,shi,thi):
    # bounds on the complex beta recurrence after multiplication by (1-h beta(beta-s))
    B=max(abs(slo),abs(shi))+thi
    A=1+g*du*B*(s-1+2*B); C=g*B*(B+1)
    den=1-h*B*(B+s)
    force=g*B*(1+s*du+B)
    assert den>0
    return fx(A/den),fx(C/den),fx(force/den),B

def nb_delta(g,du,s,slo,shi,thi):
    assert slo>=F(1,4)
    # Uses telescoping bound with factor 2 for all k. For k<=n/2,
    # j(j+s)h<=.253, so (1-h p_j)^-1<=2. For k>n/2,
    # d n>=497, and the right-hand side is >1 by its cubic term.
    C=exp_u(thi**2*(1/slo**2+1/slo)/2)
    return fx(C*2*g*(shi*(shi+1)*(shi+2)/3+F(s+1,2)*du*shi*(shi+1)))

def near_sum(dl,du,p,b,tlo,thi,ylo,yhi):
    """Upper bound sum |c_j|d^j |dressed X_j - 1|, plus infinite series tail.
    Valid on a>=998, n>=997, d in [dl,du], all indicated frequencies.
    No complex-frequency weight bound needed here.
    """
    s=2*b+1;h=hcap(du,F(998));g=hi(h/du**2);a0=F(2*p-1,2)
    sl=a0-yhi;sh=a0-ylo
    # Start at real exponent >=1/4, descend to alpha and all connection exponents.
    shift=max(0,ceil(F(1,4)-sl))
    x=nb_delta(g,du,s,sl+shift,sh+shift,thi)
    prev=nb_delta(g,du,s,sl+shift+1,sh+shift+1,thi)
    for k in range(shift,0,-1):
        A,C,E,_=step_coeff(g,h,dl,du,s,sl+k,sh+k,thi)
        x,prev=mu(A,x)+mu(C,prev)+E,x
    # e_j bounds |c_j| d^j |Delta_j|. Maintaining weighted values avoids
    # multiplying an artificially rounded-up tiny weight by huge Delta_j.
    weight=S;total=0;e=x;ep=mu(fx(1/dl),prev);rp=fx(dl)
    # At j=0 we use actual W0*Delta_-1 directly. Thereafter rp*ep bounds W_j*Delta_(j-1).
    for j in range(J):
        A,C,E,B=step_coeff(g,h,dl,du,s,sl-j,sh-j,thi)
        wprev=prev if j==0 else mu(rp,ep)
        if p==3:
            DA=8+B*(B+abs(8-4*b))/2; DB=B*(B+abs(9-4*b))/2
            delta=mu(fx(1+h*DA),e)+mu(fx(h*DA),weight)+mu(fx(h*DB/dl),weight+wprev)
        else:delta=e
        total+=delta
        ratio=((F(2*j+1,2)+yhi)**2+thi**2)/((j+1)*(j+1+2*ylo))
        r=fx(du*ratio)
        new=mu(r,mu(A,e)+mu(C,wprev)+mu(E,weight))
        weight=mu(weight,r);ep,e=e,new;rp=r
    # Infinite-twist tail. Uniform upper ratio <=(1+(yhi+.5+thi)/(J+1))^2.
    rr=du*(1+(yhi+F(1,2)+thi)/(J+1))**2;assert rr<1
    total+=cd(weight*fx(1/(1-rr)),S)
    return f(total)

def B_box(tlo,thi,ylo,yhi):
    """Gamma duplication and beta-integral bound, exact variable-y version."""
    assert tlo>=F(1,4) and 0<=ylo<=yhi
    k=int(yhi)+1;siglo=k-yhi;sighi=k-ylo
    prod=F(1)
    for j in range(k):
        a=max(abs(j+F(1,2)-ylo),abs(j+F(1,2)-yhi))
        l=F(j)-yhi;r=F(j)-ylo;dist=F(0) if l<=0<=r else min(abs(l),abs(r))
        # ratio (t^2+a^2)/(t^2+dist^2) monotone; max endpoints.
        ratio=max((t*t+a*a)/(t*t+dist*dist) for t in [tlo,thi]);prod=hi(prod*sqrtup(ratio))
    return hi(pow_upper(F(2),-2*ylo-1)*prod*sqrtup(sighi+1)/(siglo*sqrtlo(F('3.14159'))))

def inf_BH(du,ylo,yhi):
    # Euler: |B H(d)| <= coth(pi Re nu)/(2pi) int_0^1 [u(1-u)]^(y-1/2) (1-du)^(-y-1/2)du.
    # coth<=2 for Re nu>=1/4. For y<.5 use integral<=pi; for y>=.5 use 4^(-y+.5).
    factor=F(1) if ylo<F(1,2) else pow_upper(F(4),F(1,2)-ylo)/F('3.14159')
    return hi(factor*pow_upper(1-du,-yhi-F(1,2)))

def pref(dl,du,p,ylo):
    e=3-F(2*p-1,2)+ylo;fac=F(4)/LAML if p==3 else 2*LAMU*hcap(du,F(998))
    return hi(fac*pow_upper(du if e>=0 else dl,e))

def psidiv0(p,thi):
    poly=[F(x)/F(NF['k0']) for x in NF[f'psi{p-1}_exact']]; assert poly[0]==0
    return hi(sum(abs(poly[j])*(thi/LAML)**(j-1) for j in range(1,len(poly)))/LAML)

def near_contours(dl,du,p,b,Y,step=F('0.25')):
    # bottom [0,delta]; |nu B(nu)|<=sqrt(delta/(4pi)+1/(4pi^2))<=.22.
    assert F('.22')**2>DELTA/(4*F('3.14159'))+1/(4*F('3.14159')**2)
    bottom=hi(DELTA*F('.22')*psidiv0(p,DELTA)*pref(dl,du,p,F(0))*near_sum(dl,du,p,b,F(0),DELTA,F(0),F(0)))
    left=F(0);y=F(0)
    while y<Y:
        y1=min(Y,y+step);W=psi_sup(p,DELTA,DELTA,y1)
        left+=hi((y1-y)*pref(dl,du,p,y)*W*B_box(DELTA,DELTA,y,y1)*near_sum(dl,du,p,b,DELTA,DELTA,y,y1));y=y1
    return bottom,hi(left)

def ds(x):
    # exact fixed decimal; all config points in this module have finite decimal expansions.
    x=F(x);q=x.denominator;n=x.numerator;scale=10**40;assert scale%q==0;v=n*(scale//q);s=str(abs(v)).zfill(41);return ('-' if v<0 else '')+(s[:-40]+'.'+s[-40:]).rstrip('0').rstrip('.')

def config(dl,du,tlo,thi,ylo,yhi,p,b):
    # h cap rounded upward on scale 1e-15. Frequency varies on one axis at a time.
    hc=F(ceil(hcap(du,F(998))*10**15/2),10**15)
    axis='imag' if ylo!=yhi else 'real'
    assert axis=='real' or tlo==thi
    c=dict(d0=ds((dl+du)/2),dd=ds((du-dl)/2),h0=ds(hc),dh=ds(hc),nu0=ds((tlo+thi)/2),dnu=ds((yhi-ylo)/2 if axis=='imag' else (thi-tlo)/2),nuImag=ds((ylo+yhi)/2),nuAxis=axis,fullTwist=True,K=150,digits=80,epsilon='.002',channels=[[p,b]])
    cap_check(c);return c

def abs_transfer(dl,du,tlo,thi,ylo,yhi,p,b):
    c=config(dl,du,tlo,thi,ylo,yhi,p,b)
    from major_rect import major as rect
    vals=[]
    try:vals.append(major(c,p,b,[F(1)]*3)[0])
    except ValueError:pass
    try:vals.append(rect(c,p,b)[0])
    except ValueError:pass
    assert vals
    M=min(vals);return F(M,NS),F(M,NS)/dl

def osc_integrand(dl,du,p,b,tlo,thi,ylo,yhi,P,Q):
    M0,M1=base_bound(p,b,dl,du,tlo,thi,ylo,yhi)
    W=psi_sup(p,tlo,thi,yhi)
    return hi(pref(dl,du,p,ylo)*W*(B_box(tlo,thi,ylo,yhi)*(P*M0+Q*M1)+inf_BH(du,ylo,yhi)))

def top_abs(dl,du,p,b,Y,step=F(5)):
    t=DELTA;total=F(0);worst=(F(0),None)
    while t<RIGHT:
        u=min(t+(F('.25') if t<1 else F(1) if t<10 else step),RIGHT);P,Q=abs_transfer(dl,du,t,u,Y,Y,p,b);v=(u-t)*osc_integrand(dl,du,p,b,t,u,Y,Y,P,Q);total+=v
        if v>worst[0]:worst=(v,(t,u))
        t=u
    return hi(total),worst

def right_abs(dl,du,p,b,Y,step=F('.25')):
    y=F(0);total=F(0);worst=(F(0),None)
    while y<Y:
        u=min(y+step,Y);P,Q=abs_transfer(dl,du,RIGHT,RIGHT,y,u,p,b);v=(u-y)*osc_integrand(dl,du,p,b,RIGHT,RIGHT,y,u,P,Q);total+=v
        if v>worst[0]:worst=(v,(y,u))
        y=u
    return hi(total),worst

def extra(dl,du,p,b):
    # Polynomial extension [Lambda,RIGHT] (width<.000126); bound the difference of the truncated outgoing pieces there.
    # Use enclosing rational left endpoint 998.4998. Right-hand weight known through psi_sup.
    l=F('998.4998');assert l<LAML<LAMU<RIGHT
    P,Q=abs_transfer(dl,du,l,RIGHT,F(0),F(0),p,b)
    return hi((RIGHT-l)*osc_integrand(dl,du,p,b,l,RIGHT,F(0),F(0),P,Q))

if __name__=='__main__':
    dl,du,Y=map(F,sys.argv[1:4]);which=sys.argv[4] if len(sys.argv)>4 else 'all';t=time.perf_counter()
    for p in (3,4):
      y=Y+p-3
      for b in (0,4):
        row=dict(p=p,b=b,Y=str(y),dlo=str(dl),dhi=str(du))
        if which in ('all','near'):
            v,w=near_contours(dl,du,p,b,y);row['bottom'],row['left']=float(v),float(w);row['bottom_exact'],row['left_exact']=str(v),str(w)
        if which in ('all','top'):v=top_abs(dl,du,p,b,y,F(25))[0];row['top']=float(v);row['top_exact']=str(v)
        if which in ('all','right'):v,w=right_abs(dl,du,p,b,y)[0],extra(dl,du,p,b);row['right'],row['extra']=float(v),float(w);row['right_exact'],row['extra_exact']=str(v),str(w)
        print(json.dumps(row),flush=True)
    print('seconds',time.perf_counter()-t)
