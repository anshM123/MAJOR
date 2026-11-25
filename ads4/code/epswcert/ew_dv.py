"""epswcert/ew_dv.py (EPSW_CERT.md §2.4) [DV, float, not part of the proof]: (2/pi) int_0^T (W - Phi) dt = (2/pi) int r(mu) sin(mu T)/mu dmu -> r(0)
as T grows (the mechanism of the large-beta failure); the mu > 0.05 part is O(1/T) by one integration by parts and is omitted.
Usage: python3 ew_dv.py DG JSONPATH"""
import json, numpy as np, sys
from scipy.integrate import quad
ADS4='/home/user/MAJOR/ads4'
Dg=int(sys.argv[1]); path=sys.argv[2]
d=json.load(open(path)); Lam=np.sqrt((Dg-2.)*(Dg-1.)); cs=np.array(d['centres_sigma']); xr=np.array([float(s) for s in d['x_re']]); yr=np.array([float(s) for s in d['y_re']])
def Rm(z): return -2*(2*z-1j)*(2*z+1j)**2*(2*z+3j)*(2*z+5j)*(2*z+11j)/(z*(2*z-7j)*(2*z-3j)*(2*z+7j))
def chi(z):
    Q=z*z/(z*z+0.25)**2; c=cs[:,0][None,:]; s=cs[:,1][None,:]; zz=z[:,None]
    g=np.exp(-(zz-c)**2/(2*s*s)) + np.where(c>0, np.exp(-(zz+c)**2/(2*s*s)),0)
    l=(z*z+1)*np.log(z*z+1); e=np.exp(-z*z/2)
    return Q*(g@xr) + sum(yr[j]*z**(2*j)*l*e for j in range(3))
nf=json.load(open(f'{ADS4}/prl/data/nf02_psi.json')); from fractions import Fraction as Fr
co=np.array([float(Fr(x)) for x in nf['psi2_exact']])/27.48509817
def psi(mu): p=mu/Lam; return np.where(p<=1,np.polynomial.polynomial.polyval(p,co),0)

a=0.05
mu=np.concatenate([np.geomspace(1e-12,1e-6,300)[:-1],np.linspace(1e-6,a,500001)])
r=np.concatenate([2*(chi(mu[i:i+5000]-1j)*Rm(mu[i:i+5000]-1j)).real for i in range(0,len(mu),5000)])-psi(mu)
w=np.zeros_like(mu); dm=np.diff(mu); w[:-1]+=dm/2; w[1:]+=dm/2
print('Dg',Dg,'r(0)=',r[0], ' r(a)=',r[-1])
for T in (1e2,1e3,1e4,1e5):
    S=np.sum(w*r*np.sin(mu*T)/mu)
    print('T=%g: (2/pi) int_0^a r sin(mu T)/mu = %.4e   (tail beyond a <= ~ (|r(a)|/a + V)/T)'%(T,2/np.pi*S))

