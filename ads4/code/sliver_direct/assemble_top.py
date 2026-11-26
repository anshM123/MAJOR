#!/usr/bin/env python3
"""Signed top-line contour error assembly from full-twist norming results.
For p=3 use Y=6.25, p=4 Y=7.25 by default. A segment certificate is NOT
coverage of the complete top line. No bottom/left/right contour claim.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import json,argparse
from base_bounds import sqrtup,sqrtlo,logu,logl,hcap,base_bound,hi
from check_lowfreq import exp_u
from check_norm import cap_check,parameters
HERE=Path(__file__).resolve().parent
NF=json.load((HERE/'nf02_exact.json').open());LAML=sqrtlo(F(997002));LAMU=sqrtup(F(997002))
PI_LO=F('3.14159')

def pow_upper(x,e):
    return exp_u(e*(logu(x) if e>=0 else logl(x)))

def B_bound(lo,high,Y):
    """Duplication formula + positive-real beta integral after shift.
    Gautschi Γ(σ)/Γ(σ+1/2) <= sqrt(σ+1)/σ, σ>0.
    Works for x>=1/4. No use of a numerical complex Gamma routine.
    """
    assert lo>=F(1,4) and high>=lo and Y>=0
    k=Y.numerator//Y.denominator+1;sig=k-Y
    product=F(1)
    for j in range(k):
        def rat(t):return (t*t+(j-Y+F(1,2))**2)/(t*t+(j-Y)**2)
        product*=sqrtup(max(rat(lo),rat(high)))
    return hi(pow_upper(F(2),-2*Y-1)*product*sqrtup(sig+1)/(sig*sqrtlo(PI_LO)))

def B_inf(dhi,Y):
    # Euler hypergeometric integral: |B*H(d)| <= coth(pi*x)/(2*pi)
    # times integral of [u(1-u)]^(Y-.5)*(1-du)^(-Y-.5), coth < 2.
    assert Y>=F(1,2)
    return hi(pow_upper(F(4),-(Y-F(1,2)))*pow_upper(1-dhi,-Y-F(1,2))/PI_LO)

def psi_sup(p,lo,high,Y):
    poly=[F(x)/F(NF['k0']) for x in NF[f'psi{p-1}_exact']]
    rl=lo/LAMU;ru=high/LAML;c=(rl+ru)/2;R=sqrtup(((ru-rl)/2)**2+(Y/LAML)**2)
    shifts=[sum(poly[j]*comb(j,k)*c**(j-k) for j in range(k,len(poly))) for k in range(len(poly))]
    return hi(sum(abs(v)*R**k for k,v in enumerate(shifts)))

def bound_segment(row,p,b,dl,dh,lo,high,Y):
    c=row['config'];d,dd,h,hhr,nu,dnu,im=parameters(c)
    assert row['p']==p and row['b']==b and row['status'] in ('PASS','ENCLOSED') and c['fullTwist']
    assert d-dd<=dl and d+dd>=dh and c.get('nuAxis','real')=='real' and im==Y and nu-dnu<=lo<=high<=nu+dnu
    cap_check(c);P,Q=map(F,[row['P_bound'],row['Q_bound']]);M0,M1=base_bound(p,b,dl,dh,lo,high,Y,Y)
    sigma=F(2*p-1,2)-Y;assert 3-sigma>0
    pref=(4/LAML if p==3 else 2*LAMU*hcap(dh,F(998)))
    base=pref*pow_upper(dh,3-sigma)*(B_bound(lo,high,Y)*(P*M0+Q*M1)+B_inf(dh,Y))
    return hi((high-lo)*psi_sup(p,lo,high,Y)*base)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('results');ap.add_argument('--dlo',required=True);ap.add_argument('--dhi',required=True);ap.add_argument('--lo',default='0.25');ap.add_argument('--hi',default='999');args=ap.parse_args()
    rows=json.load(open(args.results));dl,dh,lo,high=map(F,[args.dlo,args.dhi,args.lo,args.hi]);assert F(1,4)<=lo<high
    for p,Y in [(3,F('6.25')),(4,F('7.25'))]:
      for b in [0,4]:
        candidates=[]
        for r in rows:
            c=r['config'];d,dd,h,hhr,nu,dnu,im=parameters(c)
            if r['p']==p and r['b']==b and r['status'] in ('PASS','ENCLOSED') and c['fullTwist'] and d-dd<=dl and d+dd>=dh and im==Y and c.get('nuAxis','real')=='real':candidates.append((nu-dnu,nu+dnu,r))
        if not candidates:continue
        pos=lo;total=F(0);count=0
        while pos<high:
            eligible=[(l,u,r) for l,u,r in candidates if l<=pos<u];assert eligible,('gap',p,b,str(pos))
            l,u,r=max(eligible,key=lambda z:z[1]);end=min(high,u)
            # Subdivide to avoid overestimating the polynomial weight and base bounds.
            while pos<end:
                q=min(end,pos+5);total+=bound_segment(r,p,b,dl,dh,pos,q,Y);count+=1;pos=q
        print(json.dumps(dict(p=p,b=b,Y=str(Y),dlo=str(dl),dhi=str(dh),lo=str(lo),hi=str(high),cells=count,top_error_bound=str(total),approx=float(total),full_top=(lo==F(1,4) and high>=LAMU))))
if __name__=='__main__':main()
