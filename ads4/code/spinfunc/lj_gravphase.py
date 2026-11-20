"""
lj_gravphase (HLJ_ANALYTIC.md §4.2): does the SCALAR travelling-wave saddle phase describe the GRAVITON rows?  (Route (c): the
graviton integrand has the same large-parameter exponent (lambda F + spin-J angle) as the scalar, only the O(1) prefactors differ.)
For a graviton row (Q.even_grid(J)[Delta], SM_CURV=none: raw Gamma = m^2 (G/(a env) - D_h)) build two images of the curvature /
J^2/m phase shift, both with the same amplitude A = sqrt(b/sinh b) and Regge term:
   kappa image  (sm_curvfix):  m^2 {A D(x_e) - D(x) + 2[P - A J0(nu b)]/env},   l theta_e = l theta kappa(b)
   saddle image (this file) :  same with l theta_s(nu) = Im[Phi_+(z_+) - Phi_0(z_0)] + const (lj_saddle, scalar exponent; the
                               constant is fixed by matching l theta_e on x in [0.03, 0.08], where kappa is exact)
and report rms of Gamma - image on x in [0.1, 0.45] and [0.45, 0.8] (C2 and C3H, ++ and +-).  No parameter is fitted to the graviton.
Usage: W2_DMIN=1000 [EE1000_DSTATE=1] python3 lj_gravphase.py J:Delta[,J:Delta]
"""
import os, sys
os.environ['SM_CURV'] = 'none'
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
import ee1000_chain as CH
import sm_curvfix as CF
Q = CH.Q


def saddle_phase(Delta, J, nus):
    """l theta_s at the given nu (ascending), continued from the largest nu downwards; uses a separate L1_DEG=2 TPS engine."""
    import subprocess, json
    code = r'''
import os, sys, json
os.environ['L1_DEG'] = '2'
sys.path.insert(0, %r)
import numpy as np, l1_crystal as Cr, lj_saddle as S
Delta, J = %r, %r
nus = np.array(json.loads(sys.stdin.read()))[::-1]
S.SREF[0] = None; S.QREF[0] = None
z, ok, _ = S.track_J(nus[0], Delta, J, 1)
z0, _ = S.newton(Cr.saddle(nus[0], Delta)[0], nus[0], Delta, None, 1)
I = S.info(z, nus[0], Delta, J, 1); S.SREF[0] = I['theta']; S.QREF[0] = S.QLAST[0]
out = []; prev = nus[0]
for nu in nus:
    for nn in np.linspace(prev, nu, max(2, int(abs(prev - nu)/2) + 1))[1:]:
        z, ok = S.newton(z, nn, Delta, J, 1); z0, _ = S.newton(z0, nn, Delta, None, 1)
        I = S.info(z, nn, Delta, J, 1); S.SREF[0] = I['theta']; S.QREF[0] = S.QLAST[0]
    prev = nu
    I = S.info(z, nu, Delta, J, 1); I0 = S.info(z0, nu, Delta, None, 1)
    out.append(float((I['Phi'] - I0['Phi']).imag))
print(json.dumps(out[::-1]))
''' % (HERE, float(Delta), int(J))
    r = subprocess.run([sys.executable, '-c', code], input=__import__('json').dumps(list(map(float, nus))), capture_output=True,
                       text=True, env=dict(os.environ, L1_DEG='2'))
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    return np.array(__import__('json').loads(r.stdout.strip().splitlines()[-1]))


def image(kd, h, J, m, xs, lth, env, SIG):
    b = 2*np.arcsinh((J + 0.5)/m); A = np.sqrt(b/np.sinh(b))
    th = np.clip(lth/(J + 0.5), 0, np.pi - 1e-9); xe = np.sin(th/2)
    osc = A*CF.Dsh(h, J, xe) - CF.Dsh(h, J, xs)
    reg = CF._pdiff(J, m, xs) + j0((2*J + 1)*xs) - A*j0(xs*m*b)
    if kd == 'C2':
        return m*m*(osc + 2*reg/env(h, xs, 2))
    return m*m*SIG[h]*(osc + reg/env(h, xs, 3))


if __name__ == '__main__':
    tag = 'direct' if os.environ.get('EE1000_DSTATE') == '1' else 'K2'
    for spec in sys.argv[1].split(','):
        J, D0 = int(spec.split(':')[0]), float(spec.split(':')[1])
        g = Q.even_grid(J)
        if D0 not in g:
            print(f"J={J} D={D0}: no row ({sorted(g)})"); continue
        r = g[D0]; xs, m = np.asarray(r['xs'], float), r['m']
        b = 2*np.arcsinh((J + 0.5)/m)
        lth_k = (J + 0.5)*2*np.arcsin(xs)*(b/2)/np.sinh(b/2)
        cal = np.array([0.03, 0.04, 0.05, 0.06, 0.07, 0.08])
        nus_all = np.concatenate([cal, xs])*m
        o = np.argsort(nus_all)
        ph = np.empty_like(nus_all); ph[o] = saddle_phase(D0, J, nus_all[o])
        lk_cal = (J + 0.5)*2*np.arcsin(cal)*(b/2)/np.sinh(b/2)
        # constant: phases may differ by multiples of pi only through the constant; fit it on the calibration window
        wrap = lambda a: (a + np.pi) % (2*np.pi) - np.pi          # Im Phi is defined mod 2 pi (log branches)
        c0 = np.angle(np.mean(np.exp(1j*(lk_cal - ph[:len(cal)]))))
        lth_s = lth_k + wrap(ph[len(cal):] + c0 - lth_k)
        ph = np.concatenate([lk_cal + wrap(ph[:len(cal)] + c0 - lk_cal), lth_s]) - c0
        print(f"J={J} Delta={D0:.0f} [{tag}] m={m:.1f} J^2/m={J*J/m:.2f} nodes={len(xs)}; calibration spread "
              f"{np.std(wrap(lk_cal - ph[:len(cal)] - c0)):.2e} rad; l(theta_s - theta_e) at x=0.2/0.4/0.6: " +
              ' '.join(f"{np.interp(xx, xs, lth_s - lth_k):+.3f}" for xx in (0.2, 0.4, 0.6)))
        for kd in ('C2', 'C3H'):
            for h in ('++', '+-'):
                if (kd, h) not in r:
                    continue
                G = np.asarray(r[(kd, h)], float)
                Ik = image(kd, h, J, m, xs, lth_k, Q.A.env, Q.A.SIG)
                Is = image(kd, h, J, m, xs, lth_s, Q.A.env, Q.A.SIG)
                s = []
                for lo, hi in ((0.1, 0.45), (0.45, 0.8)):
                    sel = (xs >= lo - 1e-9) & (xs <= hi + 1e-9)
                    if sel.sum() < 3:
                        continue
                    rms = lambda v: np.sqrt(np.mean(v[sel]**2))
                    s.append(f"x in [{lo},{hi}] (n={sel.sum()}): rms Gamma {rms(G):6.0f} | -kappa {rms(G - Ik):6.0f} | -saddle {rms(G - Is):6.0f}")
                print(f"   {kd:3s} {h}: " + "   ".join(s))
