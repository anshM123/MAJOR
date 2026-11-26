"""
tailaudit/ta_contour_bookkeeping.py [DV, double precision]: independent check of the j = 0 contour bookkeeping of TAIL_CERT §1
(quarter circle f(0)/2, segment W(y) = -2B(iy)Im f(iy), the 2Re convention), with the pure power (X_0 = 1) and the realised f of
CONTOUR (theta-cutoff included), at a moderate L where the real-axis integral can be done directly:
   2 Re int_0^oo f(nu) B(nu) w^{1/2 - i nu} dnu  =  sqrt(w) f(0)/2 + sqrt(w) int_0^Y W(y) e^{-Ly} dy + 2 Re int_0^oo f B w^{1/2-i nu}|_{nu=t+iY} dt
(Y = 0.3 < 1/2, no singularity crossed).  Also compares W(y) with CONTOUR's certified cb_w cells, and f(0) with EPSW's r(0).
f evaluated in numpy from the repaired JSON with ew_common's exact kernel polynomials.  Output: cache/ta_contour_bookkeeping.txt
"""
import os, sys, json, math
import numpy as np
from scipy.special import erf, loggamma
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'epswcert'))
sys.dont_write_bytecode = True
import ew_common as EC
out = []
def pr(s):
    print(s, flush=True); out.append(s)

d = json.load(open(os.path.join(ADS4, 'prl', 'data', 'k2', 'k2_continuum_D_1000_rep.json')))
cs = np.array([c for c, s in d['centres_sigma']]); sg = np.array([s for c, s in d['centres_sigma']]); xs = np.array([float(x) for x in d['x_re']])
C = np.concatenate([cs, -cs[cs > 0]]); S = np.concatenate([sg, sg[cs > 0]]); X = np.concatenate([xs, xs[cs > 0]])
order = np.argsort(C); C, S, X = C[order], S[order], X[order]
yv = [float(y) for y in d['y_re']]
LAM = math.sqrt(998*999); NUC = LAM + 100; ACUT = 10.0
cpoly = lambda P: np.array([float(a) + 1j*float(b) for a, b in P])
RN, RD, QN, QD = (cpoly(P) for P in (EC.RM_NUM, EC.RM_DEN, EC.QR_NUM, EC.QR_DEN))
hv = lambda c, z: np.polyval(c[::-1], z)


def G(nu):
    nu = np.asarray(nu, dtype=complex)
    nup = nu - 1j
    th = 0.5*(erf((nup + NUC)/ACUT) - erf((nup - NUC)/ACUT))
    Bg = np.zeros_like(nu)
    for lo in range(0, len(nu), 4000):
        sl = slice(lo, lo + 4000)
        re = nup[sl].real
        sel = (C > re.min() - 12) & (C < re.max() + 12)
        D = nup[sl, None] - C[None, sel]
        Bg[sl] = np.sum(X[None, sel]*np.exp(-D*D/(2*S[None, sel]**2)), axis=1)
    n2 = nup*nup
    Y = (yv[0] + yv[1]*n2 + yv[2]*n2*n2)*np.exp(-n2/2)
    l = (n2 + 1)*(np.log(nu) + np.log(nu - 2j))
    return th*(hv(QN, nup)/hv(QD, nup)*Bg + Y*l*hv(RN, nup)/hv(RD, nup))


def f(nu):
    nu = np.asarray(nu, dtype=complex)
    return G(nu) + np.conj(G(np.conj(nu)))


def B(nu):
    nu = np.asarray(nu, dtype=complex)
    return np.exp(loggamma(2j*nu) - 2*loggamma(0.5 + 1j*nu))


GLX, GLW = np.polynomial.legendre.leggauss(20)
def gl(lo, hi, width):
    k = max(1, int(math.ceil((hi - lo)/width)))
    e = np.linspace(lo, hi, k + 1)
    xs_ = ((e[1:] - e[:-1])[:, None]*(GLX[None, :] + 1)/2 + e[:-1, None]).ravel()
    ws_ = ((e[1:] - e[:-1])[:, None]*GLW[None, :]/2).ravel()
    return xs_, ws_

f0 = f(np.array([1e-12]))[0]
ew = json.load(open(os.path.join(ADS4, 'code', 'epswcert', 'cache', 'ew_cert_1000_re_rep.json')))
pr(f'f(0+) numpy = {f0.real:.6e} (imag {f0.imag:.1e});  EPSW r(0) = {ew["r0_mid"]:.6e}  (double-precision cancellation limits this)')
# W vs CONTOUR cb_w
cw = json.load(open(os.path.join(ADS4, 'code', 'contour', 'cache', 'cb_w.json')))['cells']
bad = 0; mx = 0
for k in range(0, len(cw), 75):
    y0, y1, mn, mxW = cw[k]
    y = 0.5*(y0 + y1)
    Wy = (-2*B(np.array([1j*y]))*f(np.array([1j*y + 1e-13])).imag)[0].real
    mx = max(mx, abs(Wy))
    if not (mn - 1e-9 <= Wy <= mxW + 1e-9):
        bad += 1
pr(f'W(y) = -2B(iy) Im f(iy) (numpy) inside the CONTOUR cb_w cell enclosure at 100 sampled cells: violations = {bad}')
for L in (3.0, 4.0):
    w = math.exp(-L); Yh = 0.3
    # real axis: [eps, 1300]; near 0 the B pole is imaginary, so 2Re is regular
    nus, wts = gl(1e-9, 2.0, 0.01)
    n2_, w2_ = gl(2.0, 1300.0, 0.05)
    nus = np.concatenate([nus, n2_]); wts = np.concatenate([wts, w2_])
    real = np.sum(wts*2*np.real(f(nus)*B(nus)*np.exp((0.5 - 1j*nus)*math.log(w))))
    # contour: arc + segment + horizontal at Y
    ys, wy = gl(1e-9, Yh, 0.002)
    Wseg = (-2*B(1j*ys)*f(1j*ys + 1e-13).imag).real
    seg = math.sqrt(w)*np.sum(wy*Wseg*np.exp(-L*ys))
    arc = math.sqrt(w)*f0.real/2
    ts, wt = gl(0.0, 1300.0, 0.05)
    nu = ts + 1j*Yh
    hor = np.sum(wt*2*np.real(f(nu)*B(nu)*np.exp((0.5 - 1j*nu)*math.log(w))))
    cont = arc + seg + hor
    pr(f'L = {L}: real-axis 2Re int f B w^(1/2-i nu) = {real:.12e};  arc {arc:.3e} + segment {seg:.10e} + horizontal(Y=0.3) {hor:.10e}'
       f' = {cont:.12e};  rel. diff {abs(real-cont)/abs(real):.1e}')
open(os.path.join(HERE, 'cache', 'ta_contour_bookkeeping.txt'), 'w').write('\n'.join(out) + '\n')
