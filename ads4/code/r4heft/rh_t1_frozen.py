"""Frozen-angle Mellin asymptotics for a dispersive class (k = 2, class s, D8; k = 4, class s, D8): at fixed angles (dphi0) the light
function is V(y) = y^-1 M(8y) with M(8y) = (1/2pi i) int (8y)^-d Gamma(d)^2 Phi(d) dd exactly (og gauge).  Evaluate M as
Res0(y) [double pole at d = 0, from a circle] + I_{-0.4}(y) [line Re d = -0.4]; check against the line at 2/3 at y0, and show
I_{-0.4}/Res0 -> 0 like y^0.4 while Res0 = A + B log 8y.  Usage: python3 rh_t1_frozen.py <k> <dphi0>"""
import sys
import numpy as np
from scipy.special import loggamma
import rh_common as R
import rh_t1_collinear as T

k, dphi0 = int(sys.argv[1]), float(sys.argv[2])
a, b, p2, p4 = T.cfg(dphi0)
y0 = 2*np.sin(dphi0/2)**2
Sab = R.setup(k, a, b, p2, p4, classes=('s',)); Sba = R.setup(k, b, a, p2, p4, classes=('s',))
f = R.poly_D(8)
def F(d):     # Gamma(d)^2 Phi(d) at y0 (includes (8 y0)^-d and the explicit y0 factors)
    return np.exp(2*loggamma(d))*sum((1.0 if o == 0 else (-1.0)**k)*np.sum(R.residues(S, f, d, 12)['s']) for o, S in enumerate((Sab, Sba)))
th = 2*np.pi*np.arange(64)/64; circ = 0.2*np.exp(1j*th); Fc = np.array([F(d) for d in circ])
h = 0.05; ys_ = np.arange(-10, 10 + h/2, h)
lineL = -0.4 + 1j*ys_; FL = np.array([F(d) for d in lineL])
lineR = 2/3 + 1j*np.arange(-8, 8 + 0.1, 0.2); FR = np.array([F(d) for d in lineR])
def M(y, which):
    s = np.log(y/y0)              # (8y)^-d = (8y0)^-d (y/y0)^-d ; F already carries (8 y0)^-d
    if which == 'res':
        return np.sum(Fc*np.exp(-circ*s)*circ)/64
    if which == 'L':
        return np.sum(FL*np.exp(-lineL*s))*h/(2*np.pi)
    return np.sum(FR*np.exp(-lineR*s))*0.2/(2*np.pi)
L = [f'k = {k}, class s, D8, angles frozen at dphi0 = {dphi0} (y0 = {y0:.4e}); V(y) = (y0/y) M(y), M = Res0 + I_-0.4']
L.append(f'  check at y0: line 2/3 = {M(y0, "R").real:+.8e};  Res0 + I(-0.4) = {(M(y0, "res") + M(y0, "L")).real:+.8e}')
rows = []
for y in [y0] + [10.0**(-e) for e in np.arange(3, 9.01, 1)]:
    r0 = M(y, 'res'); il = M(y, 'L') if y > 1e-6 else np.nan
    rows.append((y, r0))
    L.append(f'  y {y:.1e}  log8y {np.log(8*y):+7.3f}  Res0 = {r0.real:+.6e}   I(-0.4)/Res0 = {(il/r0).real if np.isfinite(il) else float("nan"):+.2e}')
yy = np.array([r[0] for r in rows[1:]]); RR = np.array([r[1].real for r in rows[1:]])
X = np.vstack([np.ones_like(yy), np.log(8*yy)]).T
coef, res, *_ = np.linalg.lstsq(X, RR, rcond=None)
L.append(f'  Res0 = A + B log 8y fitted on y = 1e-3..1e-9: A = {coef[0]:+.6e}, B = {coef[1]:+.6e}, max rel. residual {np.max(np.abs(X@coef/RR - 1)):.1e}'
         f'  (exact for a double pole of Gamma(d)^2 with Phi analytic at 0)')
print('\n'.join(L))
open(f'{R.CACHE}/t1_frozen_k{k}_{dphi0}.txt', 'w').write('\n'.join(L) + '\n')
