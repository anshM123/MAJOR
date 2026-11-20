"""
lj_phasecheck (HLJ_ANALYTIC.md §3.2): validation of the travelling-wave saddle against the EXACT scalar skeleton (sm_exact caches,
Arb Laplace/Barnes, independent method).  The saddle predicts the oscillation phase of the action,
    phi(nu) = Im[Phi_+(z_+(nu), nu) - Phi_0(z_0(nu), nu)]       (Phi_0: J = 0 exponent at the same Delta; removes N(nu)'s phase),
so at the zeros nu_k of the exact ex_2(nu) (located by sign changes + cubic interpolation on the cached nodes), phi(nu_k) mod pi
must be constant.  Reports the spread of phi(nu_k) mod pi (circular), and the same for the flat phase (J+1/2) 2 arcsin(nu/m) and the
leading-Regge phase nu beta, for comparison.
Usage: python3 lj_phasecheck.py Dg J [key-tag]
"""
import os, sys
os.environ.setdefault('L1_DEG', '2'); os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import pickle
import numpy as np
from scipy.interpolate import CubicSpline
import l1_crystal as Cr
import lj_saddle as S


def zeros_of(nus, f):
    o = np.argsort(nus); nus, f = nus[o], f[o]
    cs = CubicSpline(nus, f)
    z = []
    for i in range(len(nus) - 1):
        if f[i] == 0 or f[i]*f[i + 1] < 0:
            a, b = nus[i], nus[i + 1]
            for _ in range(60):
                c = 0.5*(a + b)
                if cs(a)*cs(c) <= 0:
                    b = c
                else:
                    a = c
            z.append(0.5*(a + b))
    return np.array(z)


def circ_spread(ph):
    """spread of angles mod pi: map to 2*ph mod 2pi, circular std/2."""
    c = np.exp(2j*ph)
    R = abs(c.mean())
    return np.sqrt(-2*np.log(max(R, 1e-300)))/2, np.angle(c.mean())/2


if __name__ == '__main__':
    Dg = int(sys.argv[1]); J = int(sys.argv[2]); tag = sys.argv[3] if len(sys.argv) > 3 else None
    db = pickle.load(open(os.path.join(HERE, 'cache', 'sm_exact_%d.pkl' % Dg), 'rb'))
    key = [k for k in db if k[0] == J and ((tag is None and len(k) == 3) or (tag is not None and len(k) == 4 and k[3] == tag))][0]
    r = db[key]; st = r['st']; Delta, m, beta = st['Delta'], st['m'], st['beta']
    nus = np.asarray(r['nus']); ex = np.asarray(r['ex'][2])
    sel = (nus > 4*m/(2*J + 1)) & (nus < 0.85*nus.max())
    zs = zeros_of(nus[sel], ex[sel])
    # saddle phase along the zeros (continuation from the largest)
    S.SREF[0] = None; S.QREF[0] = None
    order = np.sort(zs)[::-1]
    z, ok, _ = S.track_J(order[0], Delta, J, 1)
    z0, _ = S.newton(Cr.saddle(order[0], Delta)[0], order[0], Delta, None, 1)
    I = S.info(z, order[0], Delta, J, 1); S.SREF[0] = I['theta']; S.QREF[0] = S.QLAST[0]
    ph = []; prev = order[0]
    for nu in order:
        for nn in np.linspace(prev, nu, max(2, int(abs(prev - nu)/2) + 1))[1:]:
            z, ok = S.newton(z, nn, Delta, J, 1); z0, _ = S.newton(z0, nn, Delta, None, 1)
            I = S.info(z, nn, Delta, J, 1); S.SREF[0] = I['theta']; S.QREF[0] = S.QLAST[0]
        prev = nu
        I = S.info(z, nu, Delta, J, 1); I0 = S.info(z0, nu, Delta, None, 1)
        ph.append((I['Phi'] - I0['Phi']).imag)
    ph = np.array(ph)
    flat = (J + 0.5)*2*np.arcsin(order/m)
    regge = order*beta
    print('Dg=%d J=%d Delta=%g m=%.2f beta=%.4f: %d exact zeros on nu in (%.1f, %.1f)' % (Dg, J, Delta, m, beta, len(zs),
                                                                                            zs.min(), zs.max()))
    for name, p in (('saddle (lj_saddle)', ph), ('flat (J+1/2)theta', flat), ('leading Regge nu*beta', regge)):
        sp, mu = circ_spread(p)
        dev = np.angle(np.exp(2j*(p - mu)))/2
        print('   %-24s phase at zeros mod pi: spread %.4f rad, max |dev| %.4f rad' % (name, sp, np.abs(dev).max()))
    # local frequency check: spacing of consecutive zeros vs pi/Re kappa~
