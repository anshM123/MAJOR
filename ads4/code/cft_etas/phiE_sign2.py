"""Corrected version of phiE_sign.py: uses the cached table and adds the missing [0, 1e-4] piece (deta ~ 1 there),
which caused the constant -6.37e-5 plateau.  Zero crossings of Phi_E must be spaced by pi/Re(w1) = 1.699 (lowest QNM)."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
ws, vals = np.load('cache_deta_E.npy')
W = ws[-1]; cs = CubicSpline(ws, vals); c4 = -vals[-1]*W**4
fine = np.linspace(ws[0], W, 200001); fv = cs(fine)
def Phi(t):
    main = np.trapezoid(fv*np.cos(fine*t), fine)
    head = vals[0]*(np.sin(ws[0]*t)/t if t > 0 else ws[0])
    tail = quad(lambda w: -c4*w**-4, W, np.inf, weight='cos', wvar=t)[0]
    return 2/np.pi*(main + head + tail)
ts = np.arange(2.0, 9.51, 0.125)
P = np.array([Phi(t) for t in ts])
for t, p in zip(ts[::2], P[::2]):
    print('t=%5.2f  Phi_E=%+.4e   e^{2.664 t} Phi_E = %+.4e' % (t, p, np.exp(2.664*t)*p))
sc = [ts[i] - P[i]*(ts[i+1]-ts[i])/(P[i+1]-P[i]) for i in range(len(ts)-1) if np.sign(P[i]) != np.sign(P[i+1])]
print('zero crossings t =', ['%.3f' % x for x in sc], ' spacings =', ['%.3f' % d for d in np.diff(sc)], ' (QNM: 1.699)')
