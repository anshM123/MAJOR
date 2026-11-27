"""Sign structure of real-time relaxation functions (units 2piT = 1.5, p = 1, eta_E(0) = 1).
Phi(t) := (1/pi) Int_{-inf}^{inf} deta(w) cos(w t) dw  (so that Int_0^inf Phi dt = deta(0); E(s) = Int_0^inf Phi e^{-1.5 s t} dt).
(i) Einstein Phi_E(t), vacuum-subtracted;  (ii) deviation dPhi(t) for the horizon-layer deformation a=b=1+Amp*chi_w (horizon_bump.py).
Checks: Int Phi dt vs 4pi eta/s - 1 (or 1), and the sign of Phi on t in [0, 12]."""
import numpy as np
from scipy.integrate import quad
import einstein_ref as er
import horizon_bump as hb

def deta_E(w):
    return er.deta(w)

def deta_dev(w, Amp, wid):
    g1 = hb.Gdef(w, Amp, wid); g0 = hb.Gdef(w, 0.0, wid)
    return (g1 - g0).imag/w

def table(fun, W=40.0, n=1600):
    ws = np.linspace(1e-4, W, n)
    vals = np.array([fun(x) for x in ws])
    return ws, vals

def Phi(ws, vals, ts, tail=None):
    out = []
    for t in ts:
        integrand = vals*np.cos(ws*t)
        out.append(2/np.pi*np.trapezoid(integrand, ws))
    return np.array(out)

if __name__ == '__main__':
    ts = np.linspace(0.0, 12.0, 49)
    ws, v = table(deta_E, W=40.0, n=2000)
    PE = Phi(ws, v, ts)
    print('# Einstein Phi_E(t) (t in units r+=1; beta = 2pi/1.5 = 4.19):')
    print(' '.join('%.2f:%+.3e' % (t, p) for t, p in zip(ts[::4], PE[::4])))
    print('min Phi_E on t in (0,12] = %.3e at t = %.2f ; sign changes: %d' % (PE[1:].min(), ts[1:][PE[1:].argmin()], int(np.sum(np.diff(np.sign(PE[1:])) != 0))))
    for wid in (0.4, 0.1):
        ws2, v2 = table(lambda x: deta_dev(x, 0.1, wid), W=25.0, n=700)
        P = Phi(ws2, v2, ts)
        print('# layer deformation Amp=0.1 width w=%.2f: dPhi(t)' % wid)
        print(' '.join('%.2f:%+.3e' % (t, p) for t, p in zip(ts[::4], P[::4])))
        print('dPhi: min = %.3e at t=%.2f, max = %.3e at t=%.2f, sign changes on (0,12]: %d ; dPhi -> small at t=12: %.1e' % (
            P[1:].min(), ts[1:][P[1:].argmin()], P.max(), ts[P.argmax()], int(np.sum(np.diff(np.sign(P[1:])) != 0)), P[-1]))
