#!/usr/bin/env python3
"""Assemble one L-cell using normed residual kernels and inherited contour certificate.
Does not rerun the contour, spectral-weight or LRINF certificates.
"""
from fractions import Fraction as F
from pathlib import Path
import csv,json,sys
from check_norm import major,cap_check,parameters,S as NS
from base_bounds import base_bound,sqrtup,logl,logu,hi
from check_lowfreq import exp_u
HERE=Path(__file__).resolve().parent
DLO,DHI=F('.02598'),F('.02873')
assert DLO < 1/exp_u(F('3.65')) and DHI > exp_u(F('-3.55'))

def B_real(t):
    assert t>=159
    # coth(pi*t) < 1.001 since exp(2*pi*t)>2001 by its quadratic lower bound.
    assert 1+2*3*t+(2*3*t)**2/2>2001
    return sqrtup(F('1.001')/(4*F('3.14159')*t))
def normalized(P,Q,M0,M1,t):
    return hi(2*sqrtup(DHI)*B_real(t)*(P*M0+Q*M1))

def pow2exp(x):
    """Upper bound exp(x), no underflow; exact rational power of two."""
    x=F(x)
    if x>=0:
        q=x/logl(2);n=-(-q.numerator//q.denominator);return F(2)**n
    q=(-x)/logu(2);n=q.numerator//q.denominator;return F(2)**(-n)

def gaussian_mass(V,w=F(0)):
    """Inherited envelope |f_theta(t)|<= exp(21-(t-1029)^2/(2*.901^2))+exp(-t^2/3).
    Integrates each log-concave tail against exp(w*t), via tangent at V.
    """
    V=F(V);ss=F('.901')**2;assert V>=1040
    e1=21-(V-1029)**2/(2*ss)+w*V;r1=(V-1029)/ss-w;assert r1>0
    e2=-V*V/3+w*V;r2=2*V/3-w;assert r2>0
    return pow2exp(e1)/r1+pow2exp(e2)/r2

def kernel_from_result(row,lo,high):
    c=row['config'];d,dd,h,dh,v,dv,im=parameters(c)
    assert row['p']==3 and row['status'] in ('PASS','ENCLOSED')
    assert c['fullTwist'] and d-dd<=DLO and d+dd>=DHI and v-dv<=lo and v+dv>=high and im==0 and c.get('nuAxis','real')=='real'
    cap_check(c)
    # The full 1500-term sums have been enclosed after adding the real tail.
    M0,M1=base_bound(3,row['b'],DLO,DHI,lo,high)
    P,Q=F(row['P_bound']),F(row['Q_bound'])
    return normalized(P,Q,M0,M1,lo),M0,M1

def tail_kernel(lo,high,b):
    c={'d0':'0.027355','dd':'0.001375','h0':'0.000000000106','dh':'0.000000000106','nu0':str(float((lo+high)/2)),'dnu':str(float((high-lo)/2)),'fullTwist':True}
    cap_check(c);ma,_,_=major(c,3,b,[F(1)]*3)
    P=F(ma,NS);Q=P/DLO;M0,M1=base_bound(3,b,DLO,DHI,lo,high)
    return normalized(P,Q,M0,M1,lo)

def main():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'closure5_results.json'
    rows=json.load(path.open());ct=json.load((HERE/'contour_3p55_3p65.json').open());floor=F(ct['F_lower']);contour=sum(F(v) for p in ['psi2','psi3'] for v in ct[p].values())
    masses=list(csv.DictReader((HERE/'band_masses.csv').open()));out=[]
    for b in [0,4]:
        residual=F(0);maxkernel=F(0);parts=[]
        for r in masses:
            l,u,m=map(F,[r['lo'],r['hi'],r['mass']]);matches=[s for s in rows if s['b']==b and s['p']==3 and F(s['config']['nu0'])-F(s['config']['dnu'])<=l and F(s['config']['nu0'])+F(s['config']['dnu'])>=u]
            assert matches,('missing band',b,l,u)
            K,M0,M1=min((kernel_from_result(s,l,u) for s in matches),key=lambda t:t[0]);E=F(2,998)*K*m;residual+=E;maxkernel=max(maxkernel,K);parts.append(dict(lo=str(l),hi=str(u),mass=str(m),kernel=str(K),M0=float(M0),M1=float(M1),error=float(E)))
        # Full cutoff difference on the finite low and mid frequency range, charged again on overlaps.
        cutoff=F(2,998)*maxkernel*F('1.34e-14')
        # Far tail: ideal psi is zero here, so use the realised spectral envelope directly.
        gauss=F(0);farparts=[]
        for l,u in zip([1040,1050,1070,1100],[1050,1070,1100,1150]):
            K=tail_kernel(F(l),F(u),b);m=gaussian_mass(l);e=F(2,998)*K*m;gauss+=e;farparts.append(dict(lo=l,hi=u,kernel=float(K),mass=float(m),error=float(e)))
        # Above 1150 use the unsplit full kernel and weighted Gaussian envelope.
        # d<=2t, Cauchy radius t/2, K0<1.6 => normalised full kernel <=200 exp(pi nu/2).
        far=F(2,998)*200*gaussian_mass(1150,F(11,7));assert far<F('1e-100')
        low=F('1.46e-7');loc=F('1.1e-20');assert loc>exp_u(-46)
        total=low+residual+cutoff+gauss+far+loc+contour
        r=dict(b=b,low=float(low),residual=float(residual),cutoff=float(cutoff),gaussian=float(gauss),far_lt_1e_minus_100=True,localization=float(loc),contour=float(contour),total=float(total),floor=float(floor),ratio=float(total/floor),lower=float(floor-total),status='PASS_CONDITIONAL' if total<floor else 'FAIL',parts=parts,farparts=farparts)
        print(json.dumps({k:v for k,v in r.items() if k not in ('parts','farparts')}),flush=True);out.append(r)
    (HERE/'assembly_results.json').write_text(json.dumps(out,indent=2)+'\n')
    assert all(r['status']=='PASS_CONDITIONAL' for r in out)
if __name__=='__main__':main()
