"""p2con/pc_checks_ksh.py (P2_CONTRACT.md §1.7 checks): (i) the O(r) structure of the exact EE integrand (ee_quad, read-only):
S_nu has no O(r) term and D = D0 (1 + r * 8 [(z2.xh)(z4.n) - (z2.n)(z4.xh)]/(z2.z4)) + O(r^2);  (ii) the B1 monomial algebra vs
direct theta quadrature.  -> cache/checks_ksh.txt"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'eestep'))
import numpy as np
import ee_quad as Q
import pc_ksh as K
out = []
pr = lambda s: (print(s, flush=True), out.append(s))
rng = np.random.default_rng(5)
md = Q.mdot
for trial in range(3):
    B = rng.uniform(0.3, 2.0); pa = rng.uniform(0, 2*np.pi, 4)
    pol = Q.pols(B, *pa)
    th = np.linspace(0.1, 6.1, 7)
    xh = np.array([np.cosh(B), np.sinh(B), 0.0]); n = np.stack([0*th, np.cos(th), np.sin(th)], -1)
    z2, z4 = pol[2], pol[3]
    pred = 8*(md(z2, xh)*md(z4, n) - md(z2, n)*md(z4, xh))/md(z2, z4)
    errD = []
    for r in (1e-4, 2e-4):
        D0 = Q.dressing(th, 1e-9, B, pol)/1.0
        Dr = Q.dressing(th, r, B, pol)
        errD.append(np.max(np.abs((Dr/D0 - 1)/r - pred)))
    # S_nu: logF and logA linear-in-r content (remove explicit powers r^... via comparison at r and 2r)
    R0 = 1e-4; p1 = Q.kernel_parts(R0, B, nth=16, ds=0.5); p2 = Q.kernel_parts(2*R0, B, nth=16, ds=0.5); p4 = Q.kernel_parts(4*R0, B, nth=16, ds=0.5)
    # logF = a log r + c0 + c1 r + c2 r^2: second difference in r removes a log r only if we subtract it; use three radii
    def lin(pa_, pb_, pc_, key):
        f1, f2, f4 = pa_[key], pb_[key], pc_[key]
        # f(r) = a log r + c0 + c1 r + c2 r^2 ; eliminate a and c0 with r, 2r, 4r:
        d1 = f2 - f1; d2 = f4 - f2           # = a log2 + c1 r(1) + c2 r^2 (3) ; a log2 + 2 c1 r + 12 c2 r^2
        return np.max(np.abs((d2 - d1)/R0))   # = c1 + 9 c2 r -> c1
    pr(f'trial {trial}: B = {B:.3f}: max |(D(r)/D0 - 1)/r - pred| = {errD[0]:.2e} (r=1e-4), {errD[1]:.2e} (2e-4); '
       f'O(r) coefficient of log S: logF {lin(p1, p2, p4, "logF"):.2e}, logA {lin(p1, p2, p4, "logA"):.2e}  (|pred| ~ {np.max(np.abs(pred)):.1f})')
for trial in range(3):
    y = np.array([1.3, 0.4, -0.5]); z0 = np.array([1.0, np.cos(1.1 + trial), np.sin(1.1 + trial)])
    z1 = rng.normal(size=3) + 1j*rng.normal(size=3); z1[0] = np.sqrt(z1[1]**2 + z1[2]**2)
    z3 = rng.normal(size=3) + 1j*rng.normal(size=3); z3[0] = np.sqrt(z3[1]**2 + z3[2]**2)
    v = rng.normal(size=3)
    delta = -0.6 + 0.8j
    a = K.B1_direct(y, z0, z1, z3, v, delta, nth=200000); b = K.B1_mono_eval(y, z0, z1, z3, v, delta)
    pr(f'B1 monomials vs direct theta quadrature (delta = {delta}): {a:.10f} vs {b:.10f}, rel {abs(a - b)/abs(a):.1e}')
open(os.path.join(HERE, 'cache', 'checks_ksh.txt'), 'w').write('\n'.join(out) + '\n')
