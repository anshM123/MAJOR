#!/usr/bin/env python3
"""Compare inherited mass/contour/weight inputs to a MAJOR checkout.
This checks exact decimal ceilings and coverage; it does not rerun EPSW or contour scripts.
Usage: python verify_inputs.py /path/to/MAJOR
"""
import json,csv,sys
from fractions import Fraction as F
from decimal import Decimal
from pathlib import Path
H=Path(__file__).resolve().parent

def read(p):return json.load(p.open(),parse_float=Decimal)
def main():
    A=Path(sys.argv[1]);A=A/'ads4' if (A/'ads4').is_dir() else A
    x=read(A/'code/analytic2/cache/rmass3_40_898.json');y=read(A/'code/review_l5/cache/rt_rmass2.json')
    data=x+[r for r in y if 898<=F(r[0])<1059];count=0
    for r in csv.DictReader((H/'band_masses.csv').open()):
        l,u,m=map(F,[r['lo'],r['hi'],r['mass']]);rows=[a for a in data if l<=F(a[0])<u-F(1,10)]
        n=int(u-F(1,10)-l);assert len(rows)==n and all(F(a[0])==l+j and F(a[0])+1<=F(a[1])<=u for j,a in enumerate(rows))
        assert F('1.001')*sum(F(a[2]) for a in rows)<=m
        count+=n
    assert count==881
    C=read(H/'contour_3p5_3p65.json');src=read(A/'code/analytic3/cache/contour3_3.5_3.65_k4.json')
    assert F(C['F_lower'])<=F(src['F'])
    for p in ['psi2','psi3']:
        for k,v in C[p].items():assert F(v)>=F(src[p][k])
    NF=read(H/'nf02_exact.json');S=read(A/'prl/data/nf02_psi.json')
    for k in ['psi2_exact','psi3_exact','k0']:assert NF[k]==S[k]
    # The floor rows tile a strictly larger beta range than each L-cell.
    from check_lowfreq import exp_u
    ref=read(A/'code/lrtight/cache/lrt3_rcert.json')['rows']
    for L,U,bl,bu,floor in [('0.8','3.65','1.9','5.05','1.87e-6'),('2.9','3.2','4.25','4.6','2.5e-6'),('3.5','3.65','4.85','5.025','1.98e-6')]:
        L,U,bl,bu,floor=map(F,(L,U,bl,bu,floor))
        cu=(exp_u(bl/2)+exp_u(-bl/2))/2
        cl=(1/exp_u(-bu/2)+1/exp_u(bu/2))/2
        assert cu**2<1/exp_u(-L) and cl**2>exp_u(U)
        rows=[r for r in ref if F(r['bb'])>bl and F(r['ba'])<bu]
        assert F(rows[0]['ba'])<=bl and F(rows[-1]['bb'])>=bu
        assert all(F(a['bb'])==F(b['ba']) for a,b in zip(rows,rows[1:]))
        assert all(F(r['Fcert'])>=floor for r in rows)
    # Floating-point weights in the source are loaded as exact binary64 rationals.
    W=json.load(open(A/'prl/data/k2/k2_continuum_D_1000_rep.json'))
    sums=[F(0)]*5
    from math import ceil
    for (c,w),x in zip(W['centres_sigma'],W['x_re']):
        c,w,x=map(lambda z:F(float(z)),(c,w,x))
        assert 0<=c<=1029 and F('.449')<w<F('.901')
        sums[max(0,ceil(c)-1025)]+=abs(x)
    assert len(W['x_re'])==len(W['centres_sigma'])==3564
    assert all(a<F(b) for a,b in zip(sums,['1220000','.011','.00053','.000024','.00000026']))
    print('PASS: reference floor rows and repaired spectral inventory.')
    print('PASS: 881 source mass rows, 16 conservative band ceilings, contour ceilings, NF02 exact coefficients.')
if __name__=='__main__':main()
