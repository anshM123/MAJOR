"""deshift/d_abel.py -- DESHIFT_CLOSE.md §4: independent re-solve of Result D's k = 2 Abel equation (the NF02 de-shift of C'_EE).

Equation (RESULT_D_STATEMENT §2; EE_FUNCTIONAL §3), in the chi-normalisation:
    C[a Y](t) - i tanh(t) S[b Y](t) = Phi(t)/cosh(t),   Phi = C[psi2],   X = a Y,
    a = (R_- + R_+)/n, b = (R_- - R_+)/n, n = (|R_-| + |R_+|) * nq/|nq|,  nq = R_-^{q2} + R_+^{q2} (phase convention of u3),
    C[f](t) = int_0^oo f cos(nu t) dnu, S likewise;  psi2(nu) = psi2^{NF02}(nu/Lam)/k0 on nu <= Lam (k0 = 27.48509817).
Gauge (the stored 'min-norm' gauge, made explicit): Y = psi2 + dY with dY of minimal L2(dnu) norm among least-squares solutions.
INDEPENDENT discretisation and code (nothing from spinfunc/ or h7i14/ is imported; R_-+ are read as data via kern.load_D):
  * nu-grid nu_m = (m + 1/2) dnu with a DIFFERENT dnu (0.025 at Dg = 300, 0.05 at Dg = 1000) and the same nu_max as stored;
  * t-grid OVERSAMPLED x2 and offset: t_k = (k + 1/2) pi/(2 N dnu), k < 2N (stored: square system, t_k = k pi/(N dnu));
  * transforms by an own FFT-based type-IV cosine/sine transform (numpy.fft), CGLS from zero for the min-norm correction (own code).
Comparisons with the stored u3 weights (spinfunc/cache/u3_chi_EE_<Dg>.pkl): X on common nodes (cubic-spline interpolation of the
stored X), the off-grid de-shift residual |E[X](t) - Phi(t)| at t = 0.3 ... 25 for both (direct double-precision quadrature), the
graviton-pole identity V3, and dk0 = -2 int (Re X - psi2) p^6 dp.
Usage: python3 d_abel.py [300 1000]  -> cache/d_abel.txt, cache/d_abel_<Dg>.npz
"""
import os, sys, time, pickle, json
from fractions import Fraction
import numpy as np
from scipy.interpolate import CubicSpline
import kern

OUT = []


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)


K0 = 27.48509817


def psi2_coeffs():
    d = json.load(open(os.path.join(kern.CODE, 'robust', 'functionals_noforward.json')))['NF02']
    return np.array([float(Fraction(x)) for x in d['psi2_exact']])


class Abel:
    def __init__(self, Dg, dnu, numax):
        self.Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
        N = int(round(numax/dnu)); self.N, self.dnu = N, dnu
        self.nu = (np.arange(N) + 0.5)*dnu
        self.M = 2*N
        self.t = (np.arange(self.M) + 0.5)*np.pi/(2*N*dnu)
        Rm, Rp = kern.load_D()
        q = pickle.load(open(os.path.join(kern.CODE, 'spinfunc', 'cache', 'u2_hclosed.pkl'), 'rb'))['q2']
        import sympy as sp
        ev = lambda s: eval(s, {**sp.__dict__}) if isinstance(s, str) else s
        qm, qp = kern.Rat(ev(q[0])), kern.Rat(ev(q[1]))
        hm, hp = Rm.np(self.nu), Rp.np(self.nu)
        nq = qm.np(self.nu) + qp.np(self.nu)
        self.n = (np.abs(hm) + np.abs(hp))*nq/np.abs(nq)
        self.hm, self.hp = hm, hp
        self.a = (hm + hp)/self.n; self.b = (hm - hp)/self.n
        p = self.nu/self.Lam
        self.phi = np.where(p <= 1, np.polynomial.polynomial.polyval(p, psi2_coeffs())/K0, 0.0)
        self.th = np.tanh(self.t); self.ch = np.cosh(self.t)
        # FFT twiddles for the type-IV kernel exp(+-i pi (2m+1)(2k+1)/(4M)), m, k < M
        L = 2*self.M
        self.L = L
        self.tw_in = np.exp(1j*np.pi*np.arange(self.M)/(2*self.M))
        self.tw_out = np.exp(1j*np.pi*(2*np.arange(self.M) + 1)/(4*self.M))
        self.Phi = self.C(self.phi.astype(complex))
        self.rhs = self.Phi/self.ch

    def _S(self, f, sgn):
        """sum_{m<M} f_m exp(sgn i pi (2m+1)(2k+1)/(4M)), k < M (f zero-padded to length M)."""
        g = np.zeros(self.L, complex); g[:len(f)] = f*(self.tw_in[:len(f)] if sgn > 0 else np.conj(self.tw_in[:len(f)]))
        if sgn > 0:
            v = np.fft.ifft(g)[:self.M]*self.L
            return self.tw_out*v
        v = np.fft.fft(g)[:self.M]
        return np.conj(self.tw_out)*v

    def C(self, f):
        return self.dnu*0.5*(self._S(f, 1) + self._S(f, -1))

    def S(self, f):
        return self.dnu*(self._S(f, 1) - self._S(f, -1))/(2j)

    def A(self, Y):
        return self.C(self.a*Y) - 1j*self.th*self.S(self.b*Y)

    def AH(self, r):
        # kernel is symmetric in (m, k); adjoint of C: dnu sum_k r_k cos(nu_m t_k); of S: dnu sum_k r_k sin
        Ct = self.C(r)[:self.N]; St = self.S(self.th*r)[:self.N]
        return np.conj(self.a)*Ct + 1j*np.conj(self.b)*St

    def cgls(self, b, tol=1e-16, maxit=4000):
        x = np.zeros(self.N, complex); r = b.copy(); s = self.AH(r); p = s.copy(); g = np.vdot(s, s).real
        nb = np.linalg.norm(b); it = 0
        for it in range(maxit):
            q = self.A(p); al = g/np.vdot(q, q).real
            x += al*p; r -= al*q
            s = self.AH(r); gn = np.vdot(s, s).real
            if np.sqrt(gn) < tol*nb:
                break
            p = s + (gn/g)*p; g = gn
        return x, it

    def solve(self):
        dY, it = self.cgls(self.rhs - self.A(self.phi.astype(complex)))
        Y = self.phi + dY
        res = np.linalg.norm(self.A(Y) - self.rhs)/np.linalg.norm(self.rhs)
        return Y, res, it


def E_direct(nu, dnu, X, uX, phi, ts):
    """E[X](t) - Phi(t) = cosh t C[X] - i sinh t S[u X] - C[psi2] by direct midpoint quadrature (double precision)."""
    out = []
    for t in ts:
        c, s = np.cos(nu*t), np.sin(nu*t)
        E = np.cosh(t)*dnu*np.sum(X*c) - 1j*np.sinh(t)*dnu*np.sum(uX*s)
        out.append(abs(E - dnu*np.sum(phi*c)))
    return np.array(out)


def main(Dgs):
    t00 = time.time()
    pr("# d_abel: independent re-solve of Result D's k = 2 Abel equation (NF02 psi2, EE kernel, min-norm gauge)")
    for Dg in Dgs:
        st = pickle.load(open(os.path.join(kern.CODE, 'spinfunc', 'cache', f'u3_chi_EE_{Dg}.pkl'), 'rb'))
        numax = float(st['nu'][-1] + 0.5*(st['nu'][1] - st['nu'][0]))
        dnu = {300: 0.025, 1000: 0.05, 100: 0.0125}[Dg]
        S = Abel(Dg, dnu, numax)
        Y, res, it = S.solve()
        X = S.a*Y; uX = S.b*Y          # u X = chi (h_- - h_+) with u = b/a
        pr(f"\n## Dg = {Dg}: Lam = {S.Lam:.4f}; own grid dnu = {dnu}, N = {S.N} (nu_max {numax:.2f}), rows 2N = {S.M}; CGLS its {it}, "
           f"relative residual {res:.1e}  [{time.time() - t00:.0f}s]")
        # stored on its own grid
        nu0, X0, phi0 = st['nu'], st['X'], st['phi']
        cs_r, cs_i = CubicSpline(nu0, X0.real), CubicSpline(nu0, X0.imag)
        X0i = cs_r(S.nu) + 1j*cs_i(S.nu)
        m = S.nu <= S.Lam
        dX = np.abs(X - X0i)
        pr(f"   phi vs stored phi (interpolated): max diff/max phi = {np.max(np.abs(S.phi - CubicSpline(nu0, phi0)(S.nu)))/phi0.max():.1e}")
        for lo, hi in ((0, 2), (2, 10), (10, 25), (25, S.Lam), (S.Lam, numax)):
            sel = (S.nu >= lo) & (S.nu < hi)
            pr(f"   max|X_new - X_stored|/max psi2 on nu in [{lo:.0f}, {hi:.0f}): {dX[sel].max()/S.phi.max():.2e}")
        # residual of the identity, off-grid, both weights (double precision)
        ts = np.array([0.3, 1.0, 3.0, 6.0, 10.0, 14.0, 18.0, 22.0, 25.0])
        uX0 = st['b']*st['Y']
        r_new = E_direct(S.nu, dnu, X, uX, S.phi, ts); r_old = E_direct(nu0, nu0[1] - nu0[0], X0, uX0, phi0, ts)
        Phi_t = np.array([dnu*np.sum(S.phi*np.cos(S.nu*t)) for t in ts])
        pr("   |E[X](t) - Phi(t)| off-grid, new / stored (|Phi|): " + "; ".join(f"t={t:g}: {a:.1e} / {b:.1e} ({abs(p):.1e})"
                                                                      for t, a, b, p in zip(ts, r_new, r_old, Phi_t)))
        # V3 graviton pole and dk0
        D = lambda z: 1/(z*z + 2.25)
        V3 = np.sum(dnu*Y*(S.hm/S.n*D(S.nu + 1j) + S.hp/S.n*D(S.nu - 1j))); V3r = np.sum(dnu*S.phi*D(S.nu))
        p = S.nu/S.Lam
        dk = -2*np.sum((X - S.phi)*p**6)*dnu/S.Lam
        pr(f"   V3 graviton pole: {V3.real:.12e} vs {V3r:.12e} (rel {abs(V3 - V3r)/V3r:.1e}); stored V3 rel "
           f"{abs(st['V3'][0] - st['V3'][1])/st['V3'][1]:.1e}")
        pr(f"   dk0 = {dk.real:+.6e} (stored {st['dk0'].real:+.6e}; rel diff {abs(dk.real - st['dk0'].real)/abs(st['dk0'].real):.1e}); "
           f"bounds factor (1+dk0)^-1/2: {(1 + dk.real)**-0.5:.7f} vs {(1 + st['dk0'].real)**-0.5:.7f}")
        np.savez(os.path.join(kern.CACHE, f'd_abel_{Dg}.npz'), nu=S.nu, X=X, Y=Y, phi=S.phi, a=S.a, b=S.b, Lam=S.Lam)
    open(os.path.join(kern.CACHE, 'd_abel.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main([int(a) for a in sys.argv[1:]] or [300, 1000])
