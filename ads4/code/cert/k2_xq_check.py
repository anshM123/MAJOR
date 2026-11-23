"""X^Q (Result Q's k = 2 min-norm Abel weight of Q's psi2, RESULT_Q_REVISION §1.3) -- independent re-solve and checks (CERT_CLOSE.md §3.1).
Stored: resultq/cache/rv_k2eng_Y_<Dg>.pkl (u3_chi_EE.SolverAB, LSQR; grid nu_k = (k+1/2) dnu).  Re-solve with DESHIFT_CLOSE §4.1's
independent solver (deshift/d_abel.py: own discretisation dnu' = dnu*5/4, 2N offset rows, FFT type-IV transforms, CGLS; read-only import)
with Q's psi2 in place of NF02's.  Compares X, the graviton-pole identity V3, the off-grid residual of the identity, and
dk0^Q = -2 int (Re X - psi2) p^6 dp.  Usage: python3 k2_xq_check.py [300 1000]"""
import os, sys, json, pickle
import numpy as np
from scipy.interpolate import CubicSpline
HERE = os.path.dirname(os.path.abspath(__file__)); ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'deshift'))
import d_abel  # noqa: E402

q3 = json.load(open(os.path.join(ADS4, 'prl', 'data', 'q3_full_delta0.2.json')))
C2 = np.array([float(x) for x in q3['psi2']])


class AbelQ(d_abel.Abel):
    def __init__(self, Dg, dnu, numax):
        super().__init__(Dg, dnu, numax)
        p = self.nu / self.Lam
        self.phi = np.where(p <= 1, np.polynomial.polynomial.polyval(p, C2), 0.0)
        self.Phi = self.C(self.phi.astype(complex)); self.rhs = self.Phi / self.ch


def main(Dgs):
    for Dg in Dgs:
        st = pickle.load(open(os.path.join(ADS4, 'code', 'resultq', 'cache', 'rv_k2eng_Y_%d.pkl' % Dg), 'rb'))
        nu0 = st['nu']; d0 = nu0[1] - nu0[0]; X0 = st['a'] * st['Y']; phi0 = st['phi']
        numax = float(nu0[-1] + d0 / 2)
        S = AbelQ(Dg, d0 * 1.25, numax)
        Y, res, it = S.solve()
        X = S.a * Y; uX = S.b * Y
        print('## Dg=%d: stored grid dnu=%g N=%d; LSQR residual (stored) %.1e; re-solve dnu=%g, CGLS its %d, residual %.1e'
              % (Dg, d0, len(nu0), st['res'], S.dnu, it, res), flush=True)
        print('   psi2 check: stored phi vs Q3 psi2(nu/Lam): max diff %.1e' % np.max(np.abs(phi0 - np.where(nu0 / S.Lam <= 1, np.polynomial.polynomial.polyval(nu0 / S.Lam, C2), 0))), flush=True)
        Xi = CubicSpline(nu0, X0.real)(S.nu) + 1j * CubicSpline(nu0, X0.imag)(S.nu)
        dX = np.abs(X - Xi); sc = S.phi.max()
        for lo, hi in ((0, 2), (2, 10), (10, 25), (25, S.Lam), (S.Lam, numax)):
            s_ = (S.nu >= lo) & (S.nu < hi)
            print('   max|X_new - X_stored|/max psi2 on nu in [%.0f, %.0f): %.2e' % (lo, hi, dX[s_].max() / sc), flush=True)
        D = lambda z: 1 / (z * z + 2.25)
        V3 = np.sum(S.dnu * Y * (S.hm / S.n * D(S.nu + 1j) + S.hp / S.n * D(S.nu - 1j))); V3r = np.sum(S.dnu * S.phi * D(S.nu))
        print('   V3 graviton pole: new %.12e vs target %.12e (rel %.1e); stored (lhs, rhs) = (%.12e, %.12e) rel %.1e'
              % (V3.real, V3r, abs(V3 - V3r) / V3r, st['V3'][0].real, st['V3'][1], abs(st['V3'][0] - st['V3'][1]) / abs(st['V3'][1])), flush=True)
        p = S.nu / S.Lam; dk = (-2 * np.sum((X - S.phi) * p ** 6) * S.dnu / S.Lam).real
        p0 = nu0 / S.Lam; dk0 = (-2 * np.sum((X0 - phi0) * p0 ** 6) * d0 / S.Lam).real
        print('   dk0^Q = -2 int (Re X - psi2) p^6 dp: new %+.9e, stored %+.9e (rel diff %.1e)' % (dk, dk0, abs(dk - dk0) / abs(dk0)), flush=True)
        ts = np.array([0.3, 1.0, 3.0, 6.0, 10.0, 14.0, 18.0, 22.0])
        import kern
        Rm, Rp = kern.load_D()
        uX0 = (Rm.np(nu0) - Rp.np(nu0)) / st['n'] * st['Y']
        r_new = d_abel.E_direct(S.nu, S.dnu, X, uX, S.phi, ts); r_old = d_abel.E_direct(nu0, d0, X0, uX0, phi0, ts)
        print('   |E[X](t) - Phi(t)| off-grid, new / stored: ' + '; '.join('t=%g: %.1e / %.1e' % (t, v, w) for t, v, w in zip(ts, r_new, r_old)), flush=True)
        np.savez(os.path.join(HERE, 'cache', 'k2_xq_resolve_%d.npz' % Dg), nu=S.nu, X=X, Y=Y, phi=S.phi)


if __name__ == '__main__':
    main([int(a) for a in sys.argv[1:]] or [300, 1000])
