r"""
m1_interchange/m1i_tw6.py -- per-state actions of Omega^EE on light states of twist > 5, where CLS-D's per-block integral converges
absolutely (m1_direct sec. 4.5), by an independent route: nu-quadrature of the published Delta_g = 1000 continuum weight against the
ladder polynomial (m1_direct/md_ladder, own sympy) for k = 2, exact Q_k (md_k3) for k = 3.  Compared with the moment-ball values of
m1_remainder/out/m1a_hhat_<tw>_b<b>.json.  Also: the share of the k = 2 part coming from nu > Lam (the weight's tail).
Usage: python3 m1i_tw6.py -> out/tw6.txt
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'm1_direct')); sys.path.insert(0, os.path.join(HERE, '..', 'm1_remainder'))
import numpy as np, mpmath as mp, sympy as sp
import md_ladder as ML, md_k3 as K3, m1_weight as MW
mp.mp.dps = 50
Dg = 1000
W = MW.Weight(Dg); Lam = mp.sqrt(W.lam); lam = (Dg - 2)*(Dg - 1)
xg, wg = np.polynomial.legendre.leggauss(24)
edges = np.arange(0.0, float(Lam) + 172, 2.0)
nodes, wts = [], []
for a, b in zip(edges[:-1], edges[1:]):
    nodes += list(0.5*(b - a)*xg + 0.5*(a + b)); wts += list(0.5*(b - a)*wg)
t0 = time.time()
gv = [-(mp.mpf(9216)/Lam)*(4*mp.mpf(v)**2 + 1)/((4*mp.mpf(v)**2 + 9)*(4*mp.mpf(v)**2 + 49))*W.chi_mp(mp.mpf(v))*W.theta_mp(mp.mpf(v))
      for v in nodes]
print(f'weight: {time.time()-t0:.0f}s', flush=True)
Z, Y = ML.zy(40); Q = K3.Qk(40, lam)
out = ['# twist-tau light states, Delta_g = 1000: hhat = hhat2 (quadrature) + hhat3 (exact); [m1a moment-ball log10|hhat|]',
       '# tail = share of |hhat2| from nu > Lam']
for (tw, J, b) in [(6, 2, 0), (6, 4, 0), (6, 8, 0), (6, 8, 4), (6, 20, 0), (6, 20, 4), (8, 8, 0), (10, 8, 0), (22, 8, 0)]:
    D = tw + J; T = (D - 1)*(D - 2)
    S = sp.Poly(ML.h2(b, D, J, Z, Y), ML.s)
    co = [mp.mpf(sp.Rational(c).p)/mp.mpf(sp.Rational(c).q) for c in S.all_coeffs()]
    tot = mp.mpf(0); hi = mp.mpf(0)
    for v, w, g in zip(nodes, wts, gv):
        val = w*g*mp.polyval(co, mp.mpf(v)**2 + mp.mpf(1)/4)
        tot += val
        if v > Lam: hi += val
    h2 = T*tot
    h3 = mp.mpf(sp.Float(K3.hhat3(b, D, J, Q, lam), 40))
    ref = ''
    fn = os.path.join(HERE, '..', 'm1_remainder', 'out', f'm1a_hhat_{tw}_b{b}.json')
    if os.path.exists(fn):
        d = json.load(open(fn))
        if str(J) in d: ref = f'[m1a: sign {d[str(J)][0]:+d}, log10 {d[str(J)][1]:.3f}]'
    tt = h2 + h3
    out.append(f'  twist {tw:3d} J {J:3d} b {b}: hhat2 {mp.nstr(h2, 6):>12s}  hhat3 {mp.nstr(h3, 6):>12s}  total {mp.nstr(tt, 8):>14s}'
               f'  log10|.| {float(mp.log10(abs(tt))):8.3f} {ref}  tail(nu>Lam)/hhat2 = {mp.nstr(T*hi/h2, 3)}')
    print(out[-1], flush=True)
open(os.path.join(HERE, 'out', 'tw6.txt'), 'w').write('\n'.join(out) + '\n')
