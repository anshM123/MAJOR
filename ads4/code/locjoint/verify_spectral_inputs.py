#!/usr/bin/env python3
"""Optional inventory check on the actual repository files, no external libraries.
Does not rerun the EPSW or cutoff integral certificates."""
import sys,json,math
from fractions import Fraction as F
from pathlib import Path
root=Path(sys.argv[1]); base=root/'ads4/prl/data'
w=json.loads((base/'k2/k2_continuum_D_1000_rep.json').read_text()); p=json.loads((base/'nf02_psi.json').read_text())
assert len(w['x_re'])==len(w['centres_sigma'])==3564
s=sum(abs(F(float(x))) for x in w['x_re']); assert s<F(1220000)
assert all(0<=F(float(c))<1029 and F(449,1000)<F(float(sigma))<F(901,1000) for c,sigma in w['centres_sigma'])
assert len(w['y_re'])==3 and all(abs(F(float(y)))<1 for y in w['y_re'])
print('PASS: repaired inventory; coefficient absolute sum ceiling =',math.ceil(s))
for k in [2,3]:
    norm=sum(abs(F(x)) for x in p[f'psi{k}_exact'])/F(p['k0']); assert norm<30000000
    print(f'PASS: psi{k}/k0 coefficient l1 ceiling = {math.ceil(norm)}')
print('Inventory only; EPSW and cutoff integral bounds remain retained inputs.')
