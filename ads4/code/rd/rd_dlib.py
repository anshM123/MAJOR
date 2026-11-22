"""rd_dlib: the D' helpers of rd_dprime (log(|Phi| J) at real D points, D' membership), importable without side effects."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import rd_core as C
np.seterr(all='ignore')


def logw(X1, X2, s, Y1, Y2, phi, al, Delta, kind):
    """log(|Phi_kind| x D-Jacobian) at real D points (arrays); kind 'S' or 'EE'. Also returns geometry values."""
    z, y0 = C.raw_from_D(X1, X2, s, Y1, Y2, phi, al, Delta)
    a = 1.0/Delta
    g = C.geometry([z[i] for i in range(7)], y0val=y0.astype(complex), xsqval=(-a*a*X1*X2*(1 - s)*(1 + s)).astype(complex))
    F = C.Flow(0.0, Delta, np.ones(7))
    vals = {k: v for k, v in g.items() if not isinstance(v, list)}
    L = F.logabs_parts(vals)
    psi = al - phi
    L = L + np.log(np.sqrt(X1*X2)*np.sqrt(Y1*Y2)*np.cos(phi)**2*np.abs(np.sin(psi)))
    if kind == 'EE':
        L = L + np.log(np.abs(C.Qratio(g, Delta))) + np.log(a)
    return L, g


def inside(X1, X2, s, Y1, Y2, phi, Delta):
    a = 1.0/Delta
    ux, vx, uy, vy = a*X1, a*X2, a*Y1, a*Y2
    qx = ux*vx*(1 - s*s); qy = uy*vy*np.cos(phi)**2
    return (1 - (ux + vx) + qx > 0) & (1 - (uy + vy) + qy > 0) & (ux < 1) & (vx < 1) & (uy < 1) & (vy < 1)


