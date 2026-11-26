#!/usr/bin/env python3
"""State-local base multiplier estimates; s=1 eigenframe crossing lemma + exact helicity.
All bounds use directed Fraction/integer arithmetic. No spin or twist enumeration.
"""
from fractions import Fraction as F
from math import isqrt,comb
from check_lowfreq import exp_u
import json,argparse
S=10**40

def lo(x):x=F(x);return F(x.numerator*S//x.denominator,S)
def hi(x):x=F(x);return -lo(-x)
def sqrtlo(x):x=F(x);assert x>=0;return F(isqrt((x.numerator*S*S)//x.denominator),S)
def sqrtup(x):return sqrtlo(x)+F(1,S)

def log_bounds(x):
    x=F(x);assert x>0
    k=0
    while x<1:x*=2;k-=1
    while x>2:x/=2;k+=1
    def core(z):
        v=(z-1)/(z+1);vl=lo(v);vu=hi(v);zl=lo(vl*vl);zu=hi(vu*vu);pl=vl;pu=vu;sl=vl;su=vu
        for j in range(1,48):
            pl=lo(pl*zl);pu=hi(pu*zu);sl+=lo(pl/(2*j+1));su+=hi(pu/(2*j+1))
        # remaining terms bounded by geometric: 2 v^97/(97(1-v^2))
        tail=2*hi(pu*zu)/(97*(1-zu));return 2*sl,2*su+tail
    l,u=core(x);l2,u2=core(F(2));return l+k*(l2 if k>=0 else u2),u+k*(u2 if k>=0 else l2)
def logu(x):return log_bounds(x)[1]
def logl(x):return log_bounds(x)[0]
def asinhup(x):
    x=F(x)
    if x>=0:return logu(x+sqrtup(1+x*x))
    return logu(1/(sqrtlo(1+x*x)-x))

def hcap(dhi,a):return hi(dhi*dhi/(a*a*(1+sqrtlo(1-dhi))**2))

def positive_power(slo,shi,tlo,thi,dlo,dhi,a):
    """Uniform |d^alpha S_n(alpha)| with s=1, real sigma>=1/4.
    Exponents sigma in [slo,shi], imaginary magnitude tau in [tlo,thi].
    Returns positive real envelope, not complex phase."""
    assert F(1,4)<=slo<=shi and 0<=tlo<=thi and 0<dlo<=dhi<1 and a>=998
    assert (shi+thi)/(a*(a+1))<F(1,500), 'initial eigenframe estimate'
    # Uniform NB bound. Useful for small frequencies, no crossing singularity.
    # If tau=0 the bound is 1; otherwise for large tau avoid astronomical exponent computation.
    nb=thi*thi*(1/slo**2+1/slo)/2
    if nb<=1:return hi(exp_u(nb))
    if tlo==0:return hi(exp_u(nb)) if nb<500 else (_ for _ in ()).throw(ValueError('split tau interval away from zero'))
    # At endpoint y_n=tau(2n+1)/[(N-n)(N+n+1)]. The midpoint coordinate
    # u=a/(N+1/2) obeys u(2-u)=d(1-1/(2N+1)^2).
    deff=dlo*(1-1/(4*a*a));r=sqrtup(1-deff)
    Y=hi(thi/a*2*r/(1+r))
    kmin=slo/thi;kmax=shi/tlo;qmax=sqrtup(1+kmax*kmax);q2=1+kmin*kmin
    # Exact integral of ||P^{-1}P'||_infty; qmax/qmin <= qmax.
    A=qmax*(asinhup((q2*Y+1)/kmin)+asinhup((q2*Y-1)/kmin))/2
    A=max(F(0),A)
    z=min(F(1),Y);B=shi*logu(2/(1+sqrtlo(1-z*z)))
    if Y>1:
        dy=Y-1
        B+=thi*F(2,3)*dy*sqrtup(2*dy)+shi*sqrtup(2*dy)
    C=sqrtup((1+Y)**2+(kmax*Y)**2)+sqrtup((1-Y)**2+(kmax*Y)**2)
    # c0 <=1.01, verified from |w0| < .002, |r0-1| <= |w0|^2.
    bd=hi(F('1.01')*C*exp_u(hi(A+B)))
    if nb<500:bd=min(bd,hi(exp_u(nb)))
    return bd

def power_s1(slo,shi,tlo,thi,dlo,dhi,a):
    if slo>=F(1,4):return positive_power(slo,shi,tlo,thi,dlo,dhi,a)
    steps=max(0,-((slo-F(1,4)).numerator//(slo-F(1,4)).denominator))
    # floor identity gives ceil(1/4-slo); starting real exponent >=1/4.
    assert slo+steps>=F(1,4)
    h=hcap(dhi,a)
    x=positive_power(slo+steps,shi+steps,tlo,thi,dlo,dhi,a)
    prev=positive_power(slo+steps+1,shi+steps+1,tlo,thi,dlo,dhi,a)
    for k in range(steps,0,-1):
        ss=max(abs(slo+k),abs(shi+k));sp=max(abs(slo+k+1),abs(shi+k+1))
        den=1-h*(ss+1)**2;assert den>0
        A=1+2*h*(ss*ss+thi*thi)/dlo
        B=h*sqrtup((ss*ss+thi*thi)*(sp*sp+thi*thi))/dlo**2
        x,prev=hi((A*x+B*prev)/den),x
    return x

def base_bound(p,b,dlo,dhi,tlo,thi,ylo=F(0),yhi=F(0),a=F(998)):
    """M0,M1 on a real physical d/a domain, frequency rectangle tau+i y.
    Calls retained exact s=9 helicity identity; validates sign/denominator factors.
    """
    assert b in (0,4) and p in (3,4) and 0<=tlo<=thi
    hlo=hcap(dhi,a);q=b*(b+1);d0lo=dlo/(1+q*hlo);d0hi=dhi
    def one(k):
        sl=F(2*p-1,2)+k-yhi;sh=F(2*p-1,2)+k-ylo
        if b==0:return power_s1(sl,sh,tlo,thi,dlo,dhi,a)
        den=(1-dhi-12*hlo)**4;assert 1-dhi-12*hlo>0
        num=sum(comb(4,j)*d0hi**j*power_s1(sl-j,sh-j,tlo,thi,d0lo,d0hi,a) for j in range(5))
        # modulus (d/d0)^beta, <=(1+20h)^max(Re beta,0)
        pref=exp_u(max(F(0),sh)*logu(1+20*hlo))
        return hi(pref*num/den)
    return one(0),one(1)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');args=ap.parse_args()
    for c in json.load(open(args.input)):
        d,dd,n,dn=map(F,[c['d0'],c['dd'],c['nu0'],c['dnu']]);im=F(c.get('nuImag','0'));a=F(c.get('a_min','998'))
        for p,b in c.get('channels',[[3,4]]):
            tlo,thi,ylo,yhi=(n-dn,n+dn,im,im) if c.get('nuAxis','real')=='real' else (n,n,im-dn,im+dn)
            x,y=base_bound(p,b,d-dd,d+dd,tlo,thi,ylo,yhi,a)
            print(json.dumps(dict(name=c['name'],p=p,b=b,M0=str(x),M1=str(y),M0_approx=float(x),M1_approx=float(y))))
if __name__=='__main__':main()
