"""
lj_scan (HLJ_ANALYTIC.md §3): travelling-wave saddles z_+-(nu) of the exact scalar crystal integrand over the H_LJ region.
For each state (Dg, J, Delta) and each wave s = +-1:
  * continuation from the J = 0 crystal (in J at nu = Lam, then down in nu to NU_LO);
  * Re kappa~_s(nu) = Re(1/2 log A(z_s)) - log 2 (local nu-frequency = decay rate in Im nu, lj_saddle docstring), compared with the
    leading-Regge beta (cosh beta = 1 + (2J+1)^2/(2 m^2)) and the flat frequency (2J+1)/sqrt(m^2 - nu^2);
  * Hessian: smallest |eigenvalue| (the fibre/celestial direction x~0 is the soft one at small nu) and the smallest base eigenvalue
    (eigenvalue with the smallest x~0 weight in its eigenvector);
  * complex nu: Re[Phi_s(z_s(nu + i y), nu + i y) - Phi_s(z_s(nu), nu)] / (-y beta) along vertical lines (decay test, >= 1 - eps).
Usage: python3 lj_scan.py Dg J1,J2,... [mfac]   (mfac: m/m_gap via Delta shift; default gap state)  -> cache/lj_scan_<Dg>.txt
"""
import os, sys
os.environ.setdefault('L1_DEG', '2')
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
import lj_saddle as S

NU_LO = float(os.environ.get('LJ_NULO', '3.0'))


def Lam_of(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


def state(Dg, J, mfac=1.0):
    """Delta with m = mfac * m_gap (m^2 = (Delta-J-2)(Delta+J-1))."""
    mg2 = (Dg - 2.0)*(Dg - 1.0 + 2*J)
    m2 = mfac**2*mg2
    # solve (D-J-2)(D+J-1) = m2 for D
    b = -3.0; c = -(J + 2.0)*(J - 1.0) - m2
    D = (-b + np.sqrt(b*b - 4*c))/2
    return D, np.sqrt(m2)


def eig_split(z, nu, Delta, J, sgn):
    zv = [T.var(i, z[i]) for i in range(T.NV)]
    P, _, _, _ = S.phi(zv, nu, Delta, J, sgn)
    _, H = T.grad_hess(P)
    w, V = np.linalg.eig(H)
    fib = np.abs(V[0, :])**2/np.sum(np.abs(V)**2, axis=0)
    i_f = np.argmax(fib)
    base = [abs(w[i]) for i in range(len(w)) if i != i_f]
    return abs(w[i_f]), min(base), fib[i_f]


def scan_state(Dg, J, mfac=1.0, nnu=40, pr=print):
    Lam = Lam_of(Dg)
    Delta, m = state(Dg, J, mfac)
    beta = np.arccosh(1 + (2*J + 1)**2/(2*m*m))
    nus = np.concatenate([np.linspace(Lam, 100, nnu//2), np.geomspace(99.0, NU_LO, nnu//2)])
    res = {}
    for sgn in (1, -1):
        z, ok, Jr = S.track_J(nus[0], Delta, J, sgn)
        rows = []
        if not ok:
            pr('  J-continuation failed at J=%g' % Jr)
        for nu, zz, ok in S.track_nu(z, nus, Delta, J, sgn):
            I = S.info(zz, nu, Delta, J, sgn)
            ef, eb, fw = eig_split(zz, nu, Delta, J, sgn)
            kt = I['kappa'] - np.log(2)
            rows.append((nu, ok, kt.real, kt.imag, ef, eb, fw, zz))
        res[sgn] = rows
    return dict(Dg=Dg, J=J, mfac=mfac, Delta=Delta, m=m, beta=beta, res=res)


def vertical(Dg, J, mfac, nu_r, ys, sgn=1):
    """decay of the wave relative to the J = 0 action at the same Delta (both continued to nu_r + i sgn y, with the branches of
    (Z^2 W^2)^{-1/2} and theta continued along the path): the z-independent normalisation N(nu) and the smooth power-law
    nu-dependence common to all J cancel in the ratio."""
    import l1_crystal as Cr
    Delta, m = state(Dg, J, mfac)
    beta = np.arccosh(1 + (2*J + 1)**2/(2*m*m))
    S.SREF[0] = None; S.QREF[0] = None
    z, ok, _ = S.track_J(nu_r, Delta, J, sgn)
    I = S.info(z, nu_r, Delta, J, sgn); S.SREF[0] = I['theta']; S.QREF[0] = S.QLAST[0]
    z0, _ = S.newton(Cr.saddle(nu_r, Delta)[0], nu_r, Delta, None, sgn)
    R0 = I['Phi'] - S.info(z0, nu_r, Delta, None, sgn)['Phi']
    out = []
    zz, zz0 = z.copy(), z0.copy()
    yprev = 0.0
    for y in ys:
        good = True
        for yy in np.linspace(yprev, y, max(2, int((y - yprev)/0.25) + 1))[1:]:
            zz, ok = S.newton(zz, nu_r + 1j*sgn*yy, Delta, J, sgn); good &= ok
            I = S.info(zz, nu_r + 1j*sgn*yy, Delta, J, sgn); S.SREF[0] = I['theta']; S.QREF[0] = S.QLAST[0]
            zz0, ok = S.newton(zz0, nu_r + 1j*sgn*yy, Delta, None, sgn); good &= ok
        yprev = y
        nu = nu_r + 1j*sgn*y
        R = S.info(zz, nu, Delta, J, sgn)['Phi'] - S.info(zz0, nu, Delta, None, sgn)['Phi']
        out.append((y, good, (R - R0).real/(-y*beta)))
    S.SREF[0] = None; S.QREF[0] = None
    return beta, out


if __name__ == '__main__':
    Dg = float(sys.argv[1]); Js = [int(v) for v in sys.argv[2].split(',')]
    mfac = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    fn = os.path.join(HERE, 'cache', 'lj_scan_%d.txt' % Dg)
    f = open(fn, 'a')
    def pr(*a):
        s = ' '.join(str(x) for x in a); print(s); f.write(s + '\n'); f.flush()
    for J in Js:
        r = scan_state(Dg, J, mfac, pr=pr)
        pr('Dg=%g J=%d mfac=%.3f Delta=%.3f m=%.2f beta=%.5f m/Lam=%.4f' % (Dg, J, mfac, r['Delta'], r['m'], r['beta'],
                                                                             r['m']/Lam_of(Dg)))
        for sgn in (1, -1):
            rows = r['res'][sgn]
            okall = all(x[1] for x in rows)
            nus = np.array([x[0] for x in rows]); kr = np.array([x[2] for x in rows])
            sel = nus >= 2*(2*J + 1 and r['m']/(2*J + 1))       # outside the Bessel (Regge) endpoint zone nu < 2 m/(2J+1)
            kmin = np.min(sgn*kr[sel]) if sel.any() else np.nan
            numin = nus[sel][np.argmin(sgn*kr[sel])] if sel.any() else np.nan
            eb = np.min([x[5] for x in rows]); nub = nus[np.argmin([x[5] for x in rows])]
            pr('  wave %+d: all Newton ok=%s  nu-range [%.1f, %.1f]  min sgn*Re kappa~/beta (nu >= 2m/(2J+1)=%.1f) = %.5f at nu=%.1f;'
               '  min base |eig| = %.3g at nu=%.1f;  fibre |eig| at nu=%.1f: %.3g' % (
                   sgn, okall, nus.min(), nus.max(), 2*r['m']/(2*J + 1), kmin/r['beta'], numin, eb, nub, nus[-1], rows[-1][4]))
        nuB = 5*r['m']/(2*J + 1)
        for nu_r in (0.9*Lam_of(Dg), 0.5*Lam_of(Dg), 0.1*Lam_of(Dg), nuB):
            beta, out = vertical(Dg, J, mfac, nu_r, [1/r['beta'], 3/r['beta'], 10/r['beta'], 30/r['beta']])
            pr('  vertical line nu_r=%.1f (wave +1): ' % nu_r + '  '.join('y=%.2f ok=%d rate/beta=%.4f' % o for o in out))
    f.close()
