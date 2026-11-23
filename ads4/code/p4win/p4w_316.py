"""p4win/p4w_316.py (P4_WINDOW.md §3.4): the '(3/16) J0 y^2' residual of P2_CONTRACT §2.1 tested against the exact scalar.
If the hybrid's flat O(y^2) term were g_2 + (3/16) J0 while the exact moment's O(1/m^4) flat part is CHMRS's g_2, then
(R_ex - R_ser^(4)) m^4/nu^4 = -(3/16) J0(b) in the window.  Least-squares fit of (R_ex - R_ser^(4)) m^4/nu^4 = a J0(b) + c0 on
window nodes 10 <= nu <= nu_max (b = bhat p); also the direct y-expansion of the scalar hybrid E_2 P_J/2 at fixed b (mpmath).
-> cache/p316.txt"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0, j1, eval_legendre
import mpmath as mp
import p4w_common as C, p4w_exact as X
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))
pr(__doc__.split('->')[0])
for key, r in X.load().items():
    st = r['st']; J = r['J']; m = st['m']; mfl = st['mfl']; beta = st['beta']
    nu = np.asarray(r['nus']); x = np.asarray(r['x'])
    P, dPb = C.Pd(beta, nu); P = P.real
    bb = (2*J + 1)/mfl*nu/C.L
    E2 = (2 - x**2)/(1 - x**2)**2
    f_ex = np.asarray(r['ex'][2])/2
    f_hyb = (E2*eval_legendre(J, 1 - 2*x**2) + 2*(P - j0(bb)))/2
    ser = {k: (C.T(k, beta, nu, P, dPb)[0] - nu**(2*k)*C.gflat(k, bb))/m**(2*k) for k in range(1, 5)}
    R = f_ex - f_hyb - ser[1] - ser[2] - ser[3] - ser[4]
    for numax in (60, 150):
        s = (nu >= 10) & (nu <= numax)
        A = np.vstack([j0(bb[s]), np.ones(s.sum())]).T
        y = R[s]*m**4/nu[s]**4
        co, *_ = np.linalg.lstsq(A, y, rcond=None)
        pr(f"J = {J}: window 10 <= nu <= {numax} ({s.sum()} nodes): (R_ex - R_ser^(4)) m^4/nu^4 = {co[0]:+.4f} J0(b) {co[1]:+.4f};  "
           f"rms resid {np.sqrt(np.mean((A @ co - y)**2)):.3f}, rms data {np.sqrt(np.mean(y**2)):.3f}")
# direct flat expansion of E_2(y) P_J(1 - 2y)/2 at fixed b: c2 coefficient vs g_2 (CHMRS flat part)
mp.mp.dps = 40
for b in (0.8, 1.7, 2.7, 3.6):
    Js = list(range(40, 401, 20)); ys = []; vals = []
    g1 = 11*mp.besselj(0, b)/6 - b*mp.besselj(1, b)/6 - mp.besselj(1, b)/(6*b)
    for Jj in Js:
        yy = (mp.mpf(b)/(2*Jj + 1))**2
        v = (2 - yy)/(1 - yy)**2*mp.legendre(Jj, 1 - 2*yy)/2 - mp.besselj(0, b) - g1*yy
        ys.append(yy); vals.append(v/yy**2)
    y0 = max(ys)
    Am = mp.matrix([[(ys[i]/y0)**k for k in range(6)] for i in range(len(ys))])
    c2 = float(mp.qr_solve(Am, mp.matrix(vals))[0][0])
    g2 = float(C.gflat(2, np.array([b]))[0])
    pr(f"flat E_2 P_J/2 at fixed b = {b}: O(y^2) coefficient {c2:+.6f}, CHMRS g_2 {g2:+.6f}, difference {c2 - g2:+.6f} = {(c2 - g2)/float(mp.besselj(0, b)):+.5f} J0(b)")
open(os.path.join(C.CACHE, 'p316.txt'), 'w').write('\n'.join(OUT) + '\n')
