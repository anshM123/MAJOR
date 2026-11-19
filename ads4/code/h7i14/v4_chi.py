"""v4 (H7_I14_LIGHTSIDE.md §23.2): the non-local (eta-inversion / Abel-form) weight chi of C^_2' = A_h^{-1}[Psi_2[q^2]].

Problem.  Omega^{A'} k = 2 part = Re int_0^oo psi_2(nu/Lam) C^'_{2,nu} dnu = Re int_0^oo X(nu) Psih_nu dnu, Psih := Psi^T/(h_- + h_+) (forward-
normalised spinning CLS-D functional), X := chi (h_- + h_+) (normalisation-free).  Leading-Regge kernel of Psih (§18):
   K_nu(eta) = 1/2 [(1 + u) P_{l-1} + (1 - u) P_{l+1}](eta),   u = (h_- - h_+)/(h_- + h_+)  (exact rational, v4_hclosed),
and the defining identity is   int X K_nu(eta) dnu = int psi_2 P_l(eta) dnu   for all eta >= 1.
Abel form [D] (generalises h1_etainv_exact): by Mehler-Dirichlet for complex index, P_{-1/2 + i(nu +- i)}(cosh xi) = (sqrt2/pi) int_0^xi
cos((nu +- i) t)/sqrt(cosh xi - cosh t) dt, so the identity is equivalent (Abel transform injective) to
   E[X](t) := cosh t C[X](t) - i sinh t S[u X](t) = Phi(t) := C[psi_2](t)   for all t >= 0,
C, S the cosine/sine transforms on nu in [0, oo).  No analyticity of X, psi_2 or u is used (only real-axis values of u).
For u = i/(2 nu) (eta P_l, the old §R.5 kernel) this is h1's ODE; for u = i c/nu it is solved by H = S[X/nu], H' + c tanh(t) H = Phi/cosh t.
Numerics.  Midpoint nu-grid (step dnu, N points) and the DCT-II/DST-II-paired t-grid (dt = pi/(N dnu)); the collocated equation
   A[X] := C[X] - i tanh(t) S[u X] = Phi/cosh t        (Phi by the same DCT, so u = 0 returns X = psi_2 exactly)
is solved by GMRES, right-preconditioned with the exact ODE solution for the reference u_ref = i c/nu (c = 23/2, the large-nu coefficient of
u_T; for the eta-kernel test c = 1/2 and the preconditioner is the h1 solution).
Checks: (V0) eta-kernel: reproduce h1's psit (cache spinfunc/h1_<Dg>.pkl) and dk0; (V2) residual |A[X] - Phi/cosh|/|Phi/cosh| on the
t-grid (= the eta-identity, by injectivity of the Abel transform); (V3) graviton-pole identity int X 1/2[(1+u) D(nu+i) + (1-u) D(nu-i)] =
int psi_2 D, D = 1/(nu^2 + 9/4) (complex LHS: Re and Im tested); (U) uniqueness probe: GMRES on the homogeneous system from random starts;
(R) refinement (dnu, nu_max).  Outputs: dk0 = -2 int (X - psi_2) p^6 dp (h1 convention), max |X - psi_2|, tail mass, profiles.
Usage: python3 v4_chi.py [Dg ...]  -> cache/v4_chi.txt, cache/v4_chi_<Dg>.pkl
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
sys.dont_write_bytecode = True
import numpy as np
import sympy as sp
from numpy.polynomial import polynomial as Pl
from scipy.fft import dct, dst
from scipy.sparse.linalg import gmres, LinearOperator

_UF = None


def u_spin(nu):
    """u_T = (h_- - h_+)/(h_- + h_+), exact rational (v4_hclosed)."""
    global _UF
    if _UF is None:
        D = pickle.load(open(os.path.join(HERE, 'cache', 'v4_hclosed.pkl'), 'rb'))
        E = {k: eval(v, {**sp.__dict__}) for k, v in D.items()}
        n = sp.Symbol('nu')
        _UF = sp.lambdify(n, sp.cancel(sp.together((E['RmT'] - E['RpT'])/(E['RmT'] + E['RpT']))), 'numpy')
    return _UF(np.asarray(nu, complex))


def u_eta(nu):
    return 1j/(2*np.asarray(nu, float))


class Solver:
    def __init__(self, coeffs, Lam, dnu, N, ufun, c):
        self.c_ = np.asarray(coeffs, float); self.Lam = Lam; self.dnu = dnu; self.N = N
        self.nu = (np.arange(N) + 0.5)*dnu
        self.dt = np.pi/(N*dnu); self.t = np.arange(N)*self.dt
        p = self.nu/Lam
        self.phi = np.where(p <= 1, Pl.polyval(p, self.c_), 0.0)
        self.u = ufun(self.nu); self.c = c
        self.th = np.tanh(self.t); self.ch = np.cosh(self.t)
        self.Phi = self.Cf(self.phi)
        self.b = self.Phi/self.ch
        # exponential-integrator factors for H' + c tanh H = r: q_k = (cosh t_k/cosh t_{k+1})^c
        lc = np.log(np.cosh(self.t))
        self.q = np.exp(-c*np.diff(lc))

    # transforms on the paired grids
    def Cf(self, X):                       # C[X](t_k) = dnu sum_j X_j cos(nu_j t_k)
        return self.dnu*0.5*(dct(X.real, 2) + 1j*dct(X.imag, 2)) if np.iscomplexobj(X) else self.dnu*0.5*dct(X, 2)

    def Sf(self, Y):                       # S[Y](t_k), k = 0..N-1  (DST-II gives t_1..t_N)
        Y = np.asarray(Y, complex)
        s = self.dnu*0.5*(dst(Y.real, 2) + 1j*dst(Y.imag, 2))
        return np.concatenate([[0.0], s[:-1]])

    def Cinv(self, F):                     # X_j = (2/pi) dt [F_0/2 + sum_{k>=1} F_k cos(nu_j t_k)]
        F = np.asarray(F, complex)
        return (2/np.pi)*self.dt*0.5*(dct(F.real, 3) + 1j*dct(F.imag, 3))

    def A(self, X):
        return self.Cf(X) - 1j*self.th*self.Sf(self.u*X)

    def CT(self, r):                       # adjoint of Cf: dnu sum_k r_k cos(nu_j t_k)
        r = np.asarray(r, complex)
        f = lambda v: 0.5*(dct(v, 3) + v[0])
        return self.dnu*(f(r.real) + 1j*f(r.imag))

    def ST(self, r):                       # adjoint of Sf: dnu sum_{k>=1} r_k sin(nu_j t_k)
        r = np.asarray(r, complex)
        x = np.concatenate([r[1:], [0.0]])
        return self.dnu*0.5*(dst(x.real, 3) + 1j*dst(x.imag, 3))

    def AH(self, r):
        return self.CT(r) + 1j*np.conj(self.u)*self.ST(self.th*r)

    def solve_lsqr(self, X0=None, b=None, tol=1e-15, iter_lim=20000):
        """minimum-norm correction: X = X0 + dX, dX = argmin |dX| among the least-squares solutions of A dX = b - A X0 (LSQR from 0
        converges to the minimum-norm solution; it fixes the gauge of the 3-dim null space by dX orthogonal to it in l2(dnu))."""
        from scipy.sparse.linalg import lsqr
        X0 = self.phi.astype(complex) if X0 is None else X0
        b = self.b if b is None else b
        N = self.N
        op = LinearOperator((N, N), matvec=self.A, rmatvec=self.AH, dtype=complex)
        out = lsqr(op, (b - self.A(X0)).astype(complex), atol=tol, btol=tol, iter_lim=iter_lim)
        X = X0 + out[0]
        res = np.linalg.norm(self.A(X) - b)/np.linalg.norm(b)
        return X, res, out[2], out[1]

    def Mref(self, r):
        """exact (continuum) solution of C[X] + c tanh S[X/nu] = r for the reference kernel, discretised: H' + c tanh H = r,
        H(t) = cosh^-c(t) int_0^t r cosh^c (trapezoid; vectorised in blocks of length <= 25 in t to avoid overflow)."""
        r = np.asarray(r, complex); N = len(r); dt = self.dt; c = self.c
        L = np.log(np.cosh(self.t))
        H = np.zeros(N, complex)
        blk = max(2, int(25.0/dt))
        k0 = 0; Hb = 0.0
        while k0 < N - 1:
            k1 = min(N - 1, k0 + blk)
            seg = slice(k0, k1 + 1)
            e = np.exp(c*(L[seg] - L[k0]))                     # cosh^c(t)/cosh^c(t_k0) >= 1
            g = r[seg]*e
            cum = np.concatenate([[0.0], np.cumsum(0.5*dt*(g[1:] + g[:-1]))])
            H[seg] = (Hb + cum)/e
            Hb = H[k1]; k0 = k1
        Fp = r - c*self.th*H
        return self.Cinv(Fp)

    def solve(self, tol=1e-12, maxiter=400, restart=60, x0=None, b=None):
        N = self.N
        b = self.b if b is None else b
        op = LinearOperator((N, N), matvec=lambda y: self.A(self.Mref(y)), dtype=complex)
        it = [0]
        def cb(_):
            it[0] += 1
        y, info = gmres(op, b.astype(complex), x0=x0, rtol=tol, atol=0.0, restart=restart, maxiter=maxiter,
                        callback=cb, callback_type='pr_norm')
        X = self.Mref(y)
        res = np.linalg.norm(self.A(X) - b)/np.linalg.norm(b)
        return X, res, it[0], info

    # diagnostics
    def V3(self, X):
        nu = self.nu; w = self.dnu
        D = lambda z: 1/(z*z + 2.25)
        lhs = np.sum(w*X*0.5*((1 + self.u)*D(nu + 1j) + (1 - self.u)*D(nu - 1j)))
        rhs = np.sum(w*self.phi*D(nu))
        lhs0 = np.sum(w*self.phi*0.5*((1 + self.u)*D(nu + 1j) + (1 - self.u)*D(nu - 1j)))
        return lhs, rhs, lhs0

    def dk0(self, X):
        p = self.nu/self.Lam
        return -2*np.sum((X - self.phi)*p**6)*self.dnu/self.Lam

    def stats(self, X):
        m = self.nu <= self.Lam
        mx = np.max(np.abs(X - self.phi)[m])/np.max(np.abs(self.phi))
        tail = np.sum(np.abs(X)[~m])/np.sum(np.abs(self.phi))
        return mx, tail


def load_weights(Dg):
    from spf_common import NearGap
    NG = NearGap(Dg)
    return np.array([float(v) for v in NG.F.s2]), NG.Lam


SETTINGS = {100: (0.005, 2**16), 300: (0.01, 2**16), 1000: (0.02, 2**17)}     # (dnu, N): nu_max = 328 / 655 / 2621


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    Dgs = [int(a) for a in sys.argv[1:]] or [100, 300, 1000]
    t00 = time.time()
    pr("# v4_chi: Abel-form (eta-inversion) weight X = chi (h_- + h_+) of C^_2' (forward-normalised), GMRES + ODE preconditioner")
    for Dg in Dgs:
        cf, Lam = load_weights(Dg)
        dnu, N = SETTINGS[Dg]
        pr(f"\n## Delta_g = {Dg}: Lam = {Lam:.4f}; dnu = {dnu}, N = {N}, nu_max = {N*dnu:.0f}, dt = {np.pi/(N*dnu):.2e}, t_max = {np.pi/dnu:.0f}")
        # ---- (V0) eta kernel vs h1
        t0 = time.time()
        Se = Solver(cf, Lam, dnu, N, u_eta, 0.5)
        Xe, rese, ite, _ = Se.solve()
        h1 = pickle.load(open(os.path.join(SF, 'cache', 'h1_%d.pkl' % Dg), 'rb'))
        from scipy.interpolate import interp1d
        sel = Se.nu <= Lam
        ph1 = np.interp(Se.nu[sel], h1['nus'], h1['phit'])
        dev = np.max(np.abs(Xe[sel] - ph1))/np.max(Se.phi)
        lhs, rhs, lhs0 = Se.V3(Xe)
        pr(f"  (V0) eta-kernel (u = i/(2nu), old §R.5): GMRES its {ite}, residual {rese:.1e}, max|Im X| {np.max(np.abs(Xe.imag)):.1e}; "
           f"max|X - psit_h1|/max psi = {dev:.1e}; dk0 = {Se.dk0(Xe).real:+.5e} (h1: {h1['dk0']:+.5e}); V3 rel {abs(lhs - rhs)/rhs:.1e}  [{time.time() - t0:.0f}s]")
        # ---- spinning kernel
        t0 = time.time()
        Ss = Solver(cf, Lam, dnu, N, u_spin, 11.5)
        Xs, ress, its, info = Ss.solve()
        lhs, rhs, lhs0 = Ss.V3(Xs)
        mx, tail = Ss.stats(Xs)
        dk = Ss.dk0(Xs)
        pr(f"  (C2') spinning kernel u_T: GMRES its {its} (info {info}), residual (V2) {ress:.1e}  [{time.time() - t0:.0f}s]")
        pr(f"     V3 graviton pole: int X K[D] = {lhs.real:.12e} {lhs.imag:+.1e}i ; int psi D = {rhs:.12e} ; rel diff {abs(lhs - rhs)/rhs:.1e}"
           f"  (without inversion: rel {abs(lhs0 - rhs)/rhs:.1e})")
        pr(f"     dk0 = {dk.real:+.5e} {dk.imag:+.2e}i (eta-kernel/h1: {Se.dk0(Xe).real:+.5e});  max|X - psi|/max psi on [0,Lam] = {mx:.2e}"
           f" (eta: {Se.stats(Xe)[0]:.2e}); tail mass (nu > Lam) {tail:.1e}")
        pr(f"     => A, B_AdS factor 1/(1 + Re dk0) = {1/(1 + dk.real):.6f}; bounds factor {(1 + dk.real)**-0.5:.6f}")
        # profile at small nu and in the interior
        prof = [0.1, 0.3, 0.73, 1, 2, 4, 6.35, 10, 20, 0.25*Lam, 0.5*Lam, 0.75*Lam, 0.95*Lam]
        pr("     profile nu: X/max psi (Re, Im) | psi/max psi | psit_eta/max psi")
        for v in prof:
            j = int(v/dnu)
            pr(f"       {Ss.nu[j]:9.3f}: {Xs[j].real/Ss.phi.max():+.5e} {Xs[j].imag/Ss.phi.max():+.2e}i | {Ss.phi[j]/Ss.phi.max():+.5e} |"
               f" {Xe[j].real/Ss.phi.max():+.5e}")
        # ---- (U) uniqueness probe: homogeneous system from a random start
        rng = np.random.default_rng(1)
        y0 = (rng.standard_normal(N) + 1j*rng.standard_normal(N))*np.exp(-Ss.t/5)
        Xh, resh, ith, _ = Ss.solve(b=np.zeros(N, complex), x0=y0, tol=1e-10, maxiter=60)
        pr(f"  (U) homogeneous system from a random start: |X_h| / |M y0| after {ith} its = "
           f"{np.linalg.norm(Xh)/np.linalg.norm(Ss.Mref(y0)):.1e}")
        # ---- (R) refinement: half dnu at the same nu_max, and double nu_max at the same dnu
        out = dict(nu=Ss.nu, X=Xs, phi=Ss.phi, Xeta=Xe, Lam=Lam, dk0=dk, V3=(lhs, rhs), res=ress)
        for (lab, dn2, N2) in (("dnu/2", dnu/2, 2*N), ("nu_max x2", dnu, 2*N)):
            t0 = time.time()
            S2 = Solver(cf, Lam, dn2, N2, u_spin, 11.5)
            X2, r2, i2, _ = S2.solve()
            Xi = np.interp(Ss.nu, S2.nu, X2.real) + 1j*np.interp(Ss.nu, S2.nu, X2.imag)
            d = np.max(np.abs(Xi - Xs)[Ss.nu <= Lam])/np.max(Ss.phi)
            pr(f"  (R) {lab}: its {i2}, residual {r2:.1e}; max|dX|/max psi = {d:.1e}; dk0 = {S2.dk0(X2).real:+.6e}; "
               f"V3 rel {abs(S2.V3(X2)[0] - S2.V3(X2)[1])/S2.V3(X2)[1]:.1e}  [{time.time() - t0:.0f}s]")
        pickle.dump(out, open(os.path.join(HERE, 'cache', 'v4_chi_%d.pkl' % Dg), 'wb'))
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(HERE, 'cache', 'v4_chi_' + '_'.join(map(str, Dgs)) + '.txt'), 'w').write(buf.getvalue())
