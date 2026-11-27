"""Sign of the Einstein relaxation function Phi_E(t) = (2/pi) Int_0^inf deta(w) cos(w t) dw at late times (r+=1 units,
2piT=1.5, beta=4.19).  The lowest QNM is complex (w1 = 1.8494 - 2.6639 i), so Phi_E(t) ~ Re[c e^{-i w1 t}] must oscillate.
Numerics: deta tabulated on [0, 60] (cubic spline, 3000 pts), Filon-accurate cos integration on a 2e5 grid,
analytic tail -c4 w^-4 beyond 60 (c4 from the fit at w = 60)."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
import einstein_ref as er
W = 60.0
ws = np.concatenate([np.linspace(1e-4, 10, 1500), np.linspace(10.02, W, 1500)])
vals = np.array([er.deta(x) for x in ws])
np.save('cache_deta_E.npy', np.vstack([ws, vals]))
cs = CubicSpline(ws, vals)
c4 = -vals[-1]*W**4
fine = np.linspace(1e-4, W, 200001)
fv = cs(fine)
def Phi(t):
    main = np.trapezoid(fv*np.cos(fine*t), fine)
    tail = quad(lambda w: -c4*w**-4, W, np.inf, weight='cos', wvar=t)[0] if t > 0 else -c4/(3*W**3)
    return 2/np.pi*(main + tail)
print('check Int_0^inf Phi dt = deta(0) = 1:  (2/pi)*(pi/2)*deta(0) trivially; Phi(0) = (2/pi) Int deta = 4/pi*S*... =', Phi(0.0), '(expect 2*3/4*... eps S_E*2/pi = 0.75*... )')
ts = np.arange(0.0, 9.01, 0.25)
P = np.array([Phi(t) for t in ts])
for t, p in zip(ts, P):
    print('t=%5.2f  Phi_E=%+.6e' % (t, p))
w1 = 1.8494213135830206 - 2.6638526611526463j
print('QNM prediction: zero spacing pi/Re w1 = %.3f ; decay rate %.3f' % (np.pi/w1.real, -w1.imag))
sc = [ (ts[i]+ts[i+1])/2 for i in range(len(ts)-1) if np.sign(P[i]) != np.sign(P[i+1])]
print('sign changes near t =', ['%.2f' % x for x in sc])
