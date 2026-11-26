#!/usr/bin/env python3
"""Uniform ideal angular positivity for unbounded spin on a bounded impact region.
All inequalities are exact rational; see RESEARCH_NOTE.md for derivation."""
from check_uniform import *
K=50;t0=F(1,D);v0=F(4,131**2) # J>=65
X=F(225,16);r=F(1,50) # bhat<=7.5, x<=14.0625

def cauchy_tail(p,c,mult):
    # No complex-analytic estimate: direct recurrence majorant at and beyond K.
    g=gpolys(p,K+1)
    W=[X**k*sum(abs(z)*t0**(k-i) for i,z in enumerate(g[k])) for k in (K,K+1)]
    C=max(W[0],W[1]/r)
    # For k>=K+1, |hat g_{k+1}| X^{k+1} <= alpha |hat g_k|X^k + beta |hat g_{k-1}|X^{k-1}.
    # (2k²+2pk+p-1/4)/(k+1)^2 <=3, ((k+p-1)/(k+1))^2 <=2 for p=3,4 and k>=51.
    assert 2*(K+1)**2+2*p*(K+1)+F(p)-F(1,4)<=3*(K+2)**2
    assert (K+p)**2<=2*(K+2)**2
    alpha=X*(3*t0+F(1,(K+2)**2));beta=2*X*X*t0*t0
    assert alpha/r+beta/r**2<1
    b=bern_1(c);cells=[b]
    for _ in range(4):cells=[p for z in cells for p in split1(z)]
    norm=mult*max(abs(x) for z in cells for x in z)
    return C*norm

def series_square(a): # sum_{j>=1} (a+j)^2 r^j
    return a*a*r/(1-r)+2*a*r/(1-r)**2+r*(1+r)/(1-r)**3

def main():
    A=scaled(3,c2,K+1,2);G=scaled(4,c3,K,1)
    MA=[sum(abs(z) for z in a) for a in A];MG=[sum(abs(z) for z in a) for a in G]
    bc=bern_1([a[0]*X**k for k,a in enumerate(A[:K+1])]);cells=[bc]
    for _ in range(4):cells=[p for z in cells for p in split1(z)]
    floor=min(min(z) for z in cells)
    TA=cauchy_tail(3,c2,2);TG=cauchy_tail(4,c3,1)
    tail=TA*r/(1-r)+v0*TA*r/2*series_square(K+1)+X*v0*t0*TA*series_square(K+10)+X*v0*TG*r/(1-r)
    # charge tails of k>=K+1 in physical finite sum, even when termination occurs earlier.
    rows=[]
    for b in (0,4):
        e1=e2=e3=F(0)
        for k in range(K+1):
            e1+=X**k*sum(abs(z) for z in A[k][1:])
            e2+=X**k*abs(A[k][0])*(v0*sum(F((2*(b+j)+1)**2,4) for j in range(k))+k*t0*X)
            gamma=F(-k*k+(8-4*b)*k+9-2*b*(b-1),2)
            assert abs(gamma)<=(k+10)**2
            e3+=X**(k+1)*v0*(F((k+1)**2,2)*MA[k+1]+abs(gamma)*t0*MA[k]+MG[k])
        H=(floor-e1-e2-e3-tail)/(1+t0*X)**3
        assert H>0
        print('b',b,'reference floor',float(floor),'coeff err',float(e1),'ladder err',float(e2),'dress err',float(e3),'tail',float(tail),'H lower',float(H))
        rows.append({'b':b,'reference_floor':str(floor),'coefficient_error':str(e1),'ladder_error':str(e2),'dressing_error':str(e3),'tail':str(tail),'H_lower':str(H),'H_lower_approx':float(H)})
    Path(__file__).with_name('bounded_impact_results.json').write_text(json.dumps({'X':str(X),'J_min':65,'K':K,'rows':rows},indent=2))
if __name__=='__main__':main()
