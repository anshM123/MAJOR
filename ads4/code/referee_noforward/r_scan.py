# Independent float check of FF/N for NF02 (no forward term): FF = int psi2 (2-u p^2) D + sigma u int psi3 D, x = 1-2up^2.
import json, sys, numpy as np
from fractions import Fraction as Fr
from scipy.special import eval_jacobi, eval_legendre, roots_legendre
D=json.load(open('/home/user/MAJOR/ads4/code/robust/functionals_noforward.json'))['NF02']
X=np.array([float(Fr(c)) for c in D['X']]); Y=np.array([float(Fr(c)) for c in D['Y']])
pv=lambda c,p: np.polyval(c[::-1],p)
_c={}
def GL(n):
    n=int(np.ceil(n/200)*200)
    if n not in _c:
        t,w=roots_legendre(n); _c[n]=((t+1)/2,w/2)
    return _c[n]
def FF(J,u,h,n=None):
    n = n or max(400, J+200)
    p,w=GL(n)
    x=1-2*u*p*p
    Dh = eval_legendre(J,x) if h=='++' else eval_jacobi(J-4,0,8,x)
    fac=(1-p)**5*p
    b2=np.sum(w*fac*pv(X,p)*(2-u*p*p)*Dh); b3=(1 if h=='++' else -1)*u*np.sum(w*fac*pv(Y,p)*Dh)
    return (b2+b3)/(abs(b2)+abs(b3)), b2, b3
if __name__=='__main__':
    # sanity: quadrature convergence and scipy vs mpmath at one point
    import mpmath as mp
    mp.mp.dps=30
    J,u=1001,0.9
    print('n-convergence',FF(J,u,'+-',1400)[0],FF(J,u,'+-',2400)[0])
    xx=-0.7; print('jacobi scipy vs mpmath', eval_jacobi(J-4,0,8,xx), mp.jacobi(J-4,0,8,xx))
    worst=(9,None)
    for J in [4,5,6,7,8,10,20,50,100,200,399,400,401,500,750,1000,1500,2000,3000,4001,6000,10001,20000]:
        for u in np.concatenate([np.linspace(0.375,1,126),[0.999,0.9999,1.0]]):
            for h in (['+-'] if J>=400 else ['++','+-']):
                if h=='++' and J%2: continue
                r=FF(J,u,h)[0]
                if r<worst[0]: worst=(r,(J,u,h))
        print('J=%d running min %.5f at %s'%(J,worst[0],worst[1]),flush=True)
