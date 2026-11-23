"""
lsw/lsw_common.py (LARGESPIN_W.md): shared helpers for Route W (weight-shifting / Lemma-I cost model of the spinning
Regge expansion at large J).

Read-only imports (sys.path, no file of theirs is written):
  * n2sub/n2_common.py   : published X, psi2 on the 0.04 grid (prl/data/d_abel_X_1000.txt), P_{-1/2+i nu}(cosh b) for complex nu,
                           the local de-shift operator A (R_-/+ of kernels.json).
  * largespin/ls_core.py : NF02 flat parts + leading-Regge curvature per helicity (hyb_fast), OmLR_k(t), m_gap.

Conventions (as LARGESPIN §3, N2_SUBREGGE §1):
  * k = 2 action of an image correction dI(nu) (relative to the leading raw image, ring weight 2):
        de-shift (R1, 'A-form') : 2 * sum psi2 * f * dnu / Lam           (= 2 int psi2 f dp; Lemma I: = 2 int OmLR sigma)
        stripped (R2)           : 2 * Re sum X * f * dnu / Lam
  * a 'c2-unit' coefficient kappa of a structure S means dI = (kappa/m^2) * S(nu; beta).
  * floor_h = N_flat,h/5 + min(0, curv_h)   (certified flat floor CERT_CLOSE + exact leading-Regge curvature).
"""
import os, sys
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
CODE = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(CODE, 'n2sub'))
sys.path.insert(0, os.path.join(CODE, 'largespin'))
import n2_common as N        # noqa: E402  (read-only)
import ls_core as C          # noqa: E402  (read-only)

Dg = 1000
L = N.lam(Dg)
XD = N.xdata(Dg)
NU, X, PS, DN = XD['nu'], XD['X'], XD['psi2'], XD['dnu']
SEL = NU < L                  # psi2 vanishes beyond Lam; X is ~psi2 beyond nu ~ 25
NU, X, PS = NU[SEL], X[SEL], PS[SEL]


def P_and_dP(beta, nu=None, s=0):
    """P_{-1/2+i(nu - i s)}(cosh beta) and d/dbeta of it (4th-order central differences, h = 1e-3 beta)."""
    nu = NU if nu is None else nu
    nus = nu - 1j*s
    h = 1e-3*max(beta, 1e-3)
    f = lambda b: N.Pnu(b, nus)
    P0 = f(beta)
    d = (-f(beta + 2*h) + 8*f(beta + h) - 8*f(beta - h) + f(beta - 2*h))/(12*h)
    return P0, d


def act(f):
    """(R1 de-shift action, R2 stripped action) of a correction f(nu) on the grid, ring weight 2."""
    return 2*(PS*f).sum()*DN/L, 2*(X*f).real.sum()*DN/L


def floors(J, mfl):
    """per-helicity floor and hybrid pieces at (J, m = mfl*Lam)."""
    out = {}
    for h in ('++', '+-'):
        if (h == '++' and J % 2) or (h == '+-' and J < 4):
            continue
        st = C.hyb_fast(Dg, J, mfl, h)
        Nf = abs(st['B2']) + abs(st['B3']); cv = st['cv2'] + st['cv3']
        out[h] = dict(floor=Nf/5 + min(0.0, cv), Ns=st['Ns'], beta=st['beta'], base=st['base'], Nflat=Nf, curv=cv)
    return out


def structures(beta):
    """CHMRS d = 3 ring basis at O(1/m^2) (eq. 'regge subleading' of 2106.10274 rewritten with P'(eta) = dP/dbeta / sinh beta),
    O(1/m^4) representatives, and index-shifted (non-A-form) harmonics.  Returns dict name -> f(nu) on the grid."""
    nu = NU
    P, dP = P_and_dP(beta)
    P, dP = P.real, dP.real
    ch, sh = np.cosh(beta), np.sinh(beta)
    n2 = nu**2
    S = {
        'A ':   ch*P,                      # eta P
        'B ':   ch*n2*P,                   # eta nu^2 P
        'C1':   dP/sh,                     # P'(eta)
        'C2':   sh*n2*dP,                  # (eta^2-1) nu^2 P'(eta)
        'C3':   sh*dP,                     # (eta^2-1) P'(eta)
        # O(1/m^4) representatives (CHMRS: coefficient of P^(n)/m^{2j} grows <= eta^{j+n}, nu-degree <= 2j)
        '4a':   ch**2*n2**2*P,
        '4b':   ch*sh*n2**2*dP,
    }
    for s in (1, 2):
        Ps, _ = P_and_dP(beta, s=s)
        S['sh%+d' % s] = Ps.real          # Re P_{-1/2 + i(nu - i s)}  (index-shifted, density inside the disc; N3)
        S['sh%+d' % -s] = P_and_dP(beta, s=-s)[0].real
    return S
