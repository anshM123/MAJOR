import sys, os, json
sys.dont_write_bytecode=True
from fractions import Fraction as Fr
import numpy as np, mpmath as mp
sys.path.insert(0, '../largespin')
import ls_core as LS
D=json.load(open('../../prl/data/nf02_psi.json')); k0=Fr(D['k0'])
psi={k:[Fr(c)/k0 for c in D['psi%d_exact'%k]] for k in (2,3)}
def der(c): return [i*c[i] for i in range(1,len(c))]
def chi(c):
    d=der(c); out=[2*x for x in c]
    for i in range(len(d)): out[i+1]-= d[i]  # p*psi' : coefficient i of psi' -> power i+1
    return out
ch2={k:der(der(chi(psi[k]))) for k in (2,3)}
from scipy.special import j0
pg,wg=np.polynomial.legendre.leggauss(400); pg=(pg+1)/2; wg=wg/2
def I(k,b):
    c=np.array([float(x) for x in ch2[k]])
    return (np.polynomial.polynomial.polyval(pg,c)*j0(b*pg)*wg).sum()/24
a1={k:float(psi[k][1]) for k in (2,3)}
print('a1',a1, 'chi2(0)', {k:float(ch2[k][0]) for k in (2,3)})
L=LS.lam(1000)
mp.mp.dps=30
def M1(beta):
    # M1 = (2/pi) int [1/sig^2 - 1/(s^2 cosh(s/2))] dphi
    f=lambda ph: (lambda sig: 1/sig**2 - 1/((2*mp.asinh(sig/2))**2*mp.cosh(mp.asinh(sig/2))))(2*mp.sinh(beta/2)*mp.sin(ph))
    return float(2/mp.pi*mp.quad(f,[0,mp.pi/2]))
for b in [1,5,10,20,30,40,60,80,100,114,150,300]:
    beta=2*np.arcsinh(b/(2*L))
    C={k: L**2*(LS.OmLR(k,[beta])[0]-LS.OmFlat(k,[b])[0]) for k in (2,3)}
    g={k: a1[k]*M1(beta)+I(k,b) for k in (2,3)}
    print(b, ['%.6g/%.6g'%(C[k],g[k]) for k in (2,3)])
