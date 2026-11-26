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
    C=read(H/'contour_3p55_3p65.json');src=read(A/'code/analytic3/cache/contour3_3.55_3.65_k4.json')
    assert F(C['F_lower'])<=F(src['F'])
    for p in ['psi2','psi3']:
        for k,v in C[p].items():assert F(v)>=F(src[p][k])
    NF=read(H/'nf02_exact.json');S=read(A/'prl/data/nf02_psi.json')
    for k in ['psi2_exact','psi3_exact','k0']:assert NF[k]==S[k]
    print('PASS: 881 source mass rows, 16 conservative band ceilings, contour ceilings, NF02 exact coefficients.')
if __name__=='__main__':main()
