"""rv_common.py -- shared definitions for the independent review (written from RESEARCH_NOTE.md statements and the TeX).
NF02 weights: prl/data/nf02_psi.json (read-only)."""
import json, os
from fractions import Fraction as Fr
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..', 'prl', 'data'))
K0 = Fr('27.48509817')
_nf = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
PSI2 = [Fr(c)/K0 for c in _nf['psi2_exact']]
PSI3 = [Fr(c)/K0 for c in _nf['psi3_exact']]

def mom(c, j):
    """int_0^1 psi(p) p^j dp (exact)"""
    return sum(x/Fr(i+j+1) for i, x in enumerate(c))

def psi_mp(c, p):
    return mp.polyval([mp.mpf(x.numerator)/x.denominator for x in reversed(c)], p)

def Rm(z):
    I = mp.mpc(0, 1)
    return -2*(2*z-I)*(2*z+I)**2*(2*z+3*I)*(2*z+5*I)*(2*z+11*I)/(z*(2*z-7*I)*(2*z-3*I)*(2*z+7*I))

def Rp(z):
    return Rm(-z)

def Ppoly(z):
    return (z**2+mp.mpf(1)/4)**2*(z**2+mp.mpf(9)/4)*(z**2+mp.mpf(25)/4)*(z**2+mp.mpf(81)/4)

def Fz(z):
    """F(z) = P(z) cosh(pi z/2) Gamma(5/4+iz/2) Gamma(5/4-iz/2) (complex, mpmath; no overflow in mpmath)."""
    return Ppoly(z)*mp.cosh(mp.pi*z/2)*mp.gamma(mp.mpf(5)/4+1j*z/2)*mp.gamma(mp.mpf(5)/4-1j*z/2)

def Bz(z):
    """B(z) = F(z+i)/R_-(z) (the note's definition)."""
    return Fz(z+1j)/Rm(z)

def Bz_alt(z):
    """B(z) = F(z-i)/R_+(z) (the note's second form)."""
    return Fz(z-1j)/Rp(z)

# factored forms psi_k = (1-p)^5 p X_k(p)/k0 (avoid cancellation of the monomial form near p = 1)
X2 = [Fr(c) for c in _nf['X']]
X3 = [Fr(c) for c in _nf['Y']]
def _expand(Xc):
    # (1-p)^5 p X(p) / k0 as exact monomial list (check vs psi*_exact)
    from math import comb
    base = [Fr(0)] + [Fr((-1)**j*comb(5, j)) for j in range(6)]
    out = [Fr(0)]*(len(base)+len(Xc)-1)
    for i, a in enumerate(base):
        for j, b in enumerate(Xc):
            out[i+j] += a*b
    return [x/K0 for x in out]
assert _expand(X2) == PSI2 and _expand(X3) == PSI3, "factored NF02 form does not reproduce psi*_exact"

def psi2_fact_mp(p):
    return (1-p)**5*p*mp.polyval([mp.mpf(x.numerator)/x.denominator for x in reversed(X2)], p)/(mp.mpf(K0.numerator)/K0.denominator)
