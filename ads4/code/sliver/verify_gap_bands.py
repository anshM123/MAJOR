"""Combine completed near-gap log enclosures. No floating point in sign decision.

The logs contain rounded decimal approximations from mpq_class.get_d(); the
1e-12 relative padding is much larger than both conversion and print errors.
Use the retained EPSW/cutoff bounds described in README.md.
"""
from fractions import Fraction as F
from pathlib import Path
import re

pad=F(1,10**12)

def number(file,key):
    txt=Path(file).read_text()
    hit=re.search(r'(?<![A-Za-z_])'+re.escape(key)+r'=([-+0-9.eE]+)',txt)
    if not hit: raise RuntimeError(f'Missing completed result: {file}: {key}')
    return F(hit.group(1))

def upper(file,key):
    x=number(file,key)
    return x+abs(x)*pad+F(1,10**300)

def lower(file,key):
    x=number(file,key)
    return x-abs(x)*pad-F(1,10**300)

# For the physical interval [998,1002], the mass variable extends slightly
# beyond [-1,1].  Each degree <=300 coefficient l1 bound and analytic tail
# grows by less than 4 (including the geometric Cauchy tail).
r=F(1004,1000)
assert r**301*(1-F(1,100))/(1-r/100)<4
for J,b in [(1002,0),(1005,4)]:
    n,s=J-b,2*b+1
    D0=(n+1000)*(n+1000+s)
    assert max(abs((n+a)*(n+a+s)-D0) for a in (998,1002)) <= 8000*r
    lo=f'gap_low_b{b}.log'; mid=f'gap_mid_b{b}.log'; hi=f'gap_b{b}.log'
    # Low-band ideal weight support is wholly within V=1000; its logged
    # integer error bounds are <1e-140.  Pad by 1e-90 for expansion safety.
    ideal=lower(lo,'ideal_constant')-4*upper(lo,'ideal_variation')-F(1,10**90)
    # EPSW: total residual <3e-6; above 1000 <4e-10 (including cutoff).
    # Above 1040 the direct Gaussian envelope gives total mass <1e-22.
    res=F(8,998)*(F(3,10**6)*upper(lo,'spectral_kernel_norm_upper')+
                   F(4,10**10)*upper(mid,'spectral_kernel_norm_upper')+
                   F(1,10**22)*upper(hi,'spectral_kernel_norm_upper'))+F(1,10**20)
    result=ideal-res
    print(f'J={J} b={b} a=[998,1002] ideal_lower={float(ideal):.12g} residual_upper={float(res):.12g} action_lower={float(result):.12g} {"PASS" if result>0 else "FAIL"}')
    assert result>0, 'Bandwise lower bound insufficient; no positivity conclusion.'
