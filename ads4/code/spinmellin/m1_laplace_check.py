"""M1 check: Mehler-Fock transform L(s) = int [d eta] P_{-1/2+i nu}(eta) e^{-s xi} of a pure exponential (the building block of the Regge moment of a Mellin monomial) is power-law in nu, not ~e^{-pi nu} (numerical quadrature; the Gamma-ratio guess is NOT the closed form)."""
import mpmath as mp
mp.mp.dps=30
def Lnum(s,nu):
    mu=-0.5+1j*nu
    f=lambda x: 2*mp.sinh(x)*mp.exp(-s*x)*mp.legenp(mu,0,mp.cosh(x),type=3)
    return mp.quad(f,[0,1,5,20,mp.inf])
def Lguess(s,nu):
    # guess: Mehler-Fock of e^{-s xi}
    a=(s+mp.mpf(1)/2+1j*nu)/2; b=(s+mp.mpf(1)/2-1j*nu)/2
    return mp.gamma(a-mp.mpf(1)/2)*mp.gamma(b-mp.mpf(1)/2)/(mp.gamma(a+mp.mpf(1)/2)*mp.gamma(b+mp.mpf(1)/2))
for s in (2.5,3.7,6.0):
    for nu in (0.0,1.3,4.0):
        n=Lnum(s,nu); g=Lguess(s,nu); print(s,nu,mp.nstr(n,15),mp.nstr(g,15),mp.nstr(n/g,15))
