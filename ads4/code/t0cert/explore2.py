import sys, json
sys.dont_write_bytecode=True
from fractions import Fraction as Fr
import numpy as np
from scipy.special import j0
D=json.load(open('../../prl/data/nf02_psi.json')); k0=Fr(D['k0'])
psi={k:[Fr(c)/k0 for c in D['psi%d_exact'%k]] for k in (2,3)}
def der(c): return [i*c[i] for i in range(1,len(c))]
def chi(c):
    d=der(c); out=[2*x for x in c]
    for i in range(len(d)): out[i+1]-= d[i]
    return out
ch2={k:np.array([float(x) for x in der(der(chi(psi[k])))]) for k in (2,3)}
pg,wg=np.polynomial.legendre.leggauss(600); pg=(pg+1)/2; wg=wg/2
a1={k:float(psi[k][1]) for k in (2,3)}
bs=np.arange(0.01,200,0.01)
g={k: a1[k]/24 + (np.polynomial.polynomial.polyval(pg,ch2[k])[None,:]*j0(np.outer(bs,pg))*wg).sum(1)/24 for k in (2,3)}
for c in [0,0.1,0.15,0.2]:
    Kmin=np.full_like(bs,np.inf)
    for s2 in (1,-1):
      for s3 in (1,-1):
        for sig in (1,-1):
          for u in (0,1):
            K=(1-c*s2)*2*g[2]+(1-c*s3)*u*sig*g[3]
            Kmin=np.minimum(Kmin,K)
    neg=bs[Kmin<0]
    print('c=',c,'K<0 at bhat in [%g,%g]'%(neg.min(),neg.max()) if len(neg) else 'never', ' min K on bhat>=b:', [(b, round(Kmin[bs>=b].min(),3)) for b in (5,8,10,12,15,20,30,50,100)])
print('g2,g3 extremes', g[2].min(), g[3].min(), g[3].max())
