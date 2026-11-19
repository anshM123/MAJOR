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
is solved by LSQR (exact adjoint; minimum-norm correction X - psi_2, which fixes the gauge of the 3-dim null space, §23.2).  For the
eta-kernel check also GMRES right-preconditioned with the exact ODE solution (u_ref = i/(2nu), = h1).  (GMRES with u_ref = 23i/(2nu)
stagnates on u_T because of the null space.)
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


SETTINGS = {100: (0.01, 2**14), 300: (0.02, 2**14), 1000: (0.04, 2**15)}     # (dnu, N): nu_max = 164 / 328 / 1311


def null_space(dnu=0.025, N=4096, ufun=None, k=6):
    """dense SVD of the collocated operator A (depends on u only, not on Lam/psi): smallest singular values and right vectors."""
    S = Solver(np.zeros(10), 1.0, dnu, N, ufun or u_spin, 11.5)
    I = np.eye(N)
    M = np.array([S.A(I[:, j]) for j in range(N)]).T
    U, s, Vh = np.linalg.svd(M)
    return S.nu, s, Vh[-k:].conj()[::-1], s[:2]


def winding(rho_fun):
    nus = np.concatenate([np.linspace(1e-7, 1, 40001), np.linspace(1, 1e5, 400001)[1:]])
    a = np.unwrap(np.angle(rho_fun(nus)))
    return (a[-1] - a[0])/np.pi


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    Dgs = [int(a) for a in sys.argv[1:]] or [100, 300, 1000]
    t00 = time.time()
    pr("# v4_chi: Abel-form (eta-inversion) weight X = chi (h_- + h_+) of C^_2' (forward-normalised); min-norm LSQR solution")
    # ---------------- null space of the Abel operator (u only) and the winding of rho = h_+/h_-
    D = pickle.load(open(os.path.join(HERE, 'cache', 'v4_hclosed.pkl'), 'rb'))
    E = {k: eval(v, {**sp.__dict__}) for k, v in D.items()}; n_ = sp.Symbol('nu')
    F = {k: sp.lambdify(n_, v, 'numpy') for k, v in E.items()}
    w_eta = winding(lambda v: (2*v - 1j)/(2*v + 1j)); w_sc = winding(lambda v: F['RpS'](v)/F['RmS'](v))
    w_T = winding(lambda v: F['RpT'](v)/F['RmT'](v))
    pr(f"(N0) winding Delta arg(rho)/pi on (0, oo), rho = h_+/h_- (rho(0+) = -1, rho(oo) = +1): eta-kernel {w_eta:+.4f}, scalar g {w_sc:+.4f}, "
       f"spinning h {w_T:+.4f}  -> predicted null dimension (w - 1)/2 = {round((w_eta-1)/2)}, {round((w_sc-1)/2)}, {round((w_T-1)/2)}")
    nul = {}
    for lab, (dn, NN) in (("coarse", (0.05, 2048)), ("fine", (0.025, 4096))):
        t0 = time.time()
        nu0, s0, V0, smax = null_space(dn, NN)
        nul[lab] = (nu0, s0, V0)
        pr(f"(N1) dense SVD of the u_T operator, dnu = {dn}, N = {NN}: smallest sv " + " ".join(f"{v:.1e}" for v in s0[::-1][:6]) +
           f"; largest {smax[0]:.2f}  [{time.time() - t0:.0f}s]")
    nu0, s0, V0 = nul["fine"]
    for kk in range(3):
        a = np.abs(V0[kk]); j = np.argmax(a)
        supp = nu0[a > 1e-3*a.max()]
        pr(f"     null vector {kk + 1} (sv {s0[::-1][kk]:.1e}): peak at nu = {nu0[j]:.3f}; |X_h| > 1e-3 max on nu in [{supp.min():.2f}, {supp.max():.2f}]; "
           f"|X_h(20)|/max = {a[int(20/0.025)]/a.max():.1e}, |X_h(50)|/max = {a[int(50/0.025)]/a.max():.1e}")
    for (dn, NN, lab) in ((0.05, 2048, 'scalar g'),):
        S_ = Solver(np.zeros(10), 1.0, dn, NN, u_eta, 0.5)
        hm, hp = F['RmS'](S_.nu), F['RpS'](S_.nu); nn = np.abs(hm) + np.abs(hp)
        Aa, Bb = (hm + hp)/nn, (hm - hp)/nn
        I_ = np.eye(NN)
        M_ = np.array([S_.Cf(Aa*I_[:, j]) - 1j*S_.th*S_.Sf(Bb*I_[:, j]) for j in range(NN)]).T
        sv = np.linalg.svd(M_, compute_uv=False)
        pr(f"(N2) scalar-g kernel (old C^_2 of H0 §24.1, chi-normalisation n = |g| + |g~|; u_sc has a real pole at nu_0): smallest sv "
           + " ".join(f"{v:.1e}" for v in sv[::-1][:4]))
    pickle.dump(dict(nu=nu0, s=s0, V=V0), open(os.path.join(HERE, 'cache', 'v4_chi_null.pkl'), 'wb'))
    for Dg in Dgs:
        cf, Lam = load_weights(Dg)
        dnu, N = SETTINGS[Dg]
        pr(f"\n## Delta_g = {Dg}: Lam = {Lam:.4f}; dnu = {dnu}, N = {N}, nu_max = {N*dnu:.0f}, t_max = {np.pi/dnu:.0f}")
        # ---- (V0) eta kernel vs h1 (GMRES with the exact ODE preconditioner) and with LSQR
        t0 = time.time()
        Se = Solver(cf, Lam, dnu, N, u_eta, 0.5)
        Xe, rese, ite, _ = Se.solve()
        Xe2, rese2, ite2, _ = Se.solve_lsqr()
        h1 = pickle.load(open(os.path.join(SF, 'cache', 'h1_%d.pkl' % Dg), 'rb'))
        sel = Se.nu <= Lam
        ph1 = np.interp(Se.nu[sel], h1['nus'], h1['phit'])
        dev = np.max(np.abs(Xe[sel] - ph1))/np.max(Se.phi)
        lhs, rhs, lhs0 = Se.V3(Xe)
        pr(f"  (V0) eta-kernel (u = i/(2nu), old §R.5): GMRES its {ite}, residual {rese:.1e}; LSQR its {ite2}, residual {rese2:.1e}, "
           f"max|X_gmres - X_lsqr|/max psi {np.max(np.abs(Xe - Xe2))/Se.phi.max():.1e}; max|X - psit_h1|/max psi = {dev:.1e}; "
           f"dk0 = {Se.dk0(Xe).real:+.6e} (h1: {h1['dk0']:+.6e}); V3 rel {abs(lhs - rhs)/rhs:.1e}  [{time.time() - t0:.0f}s]")
        # ---- spinning kernel
        t0 = time.time()
        Ss = Solver(cf, Lam, dnu, N, u_spin, 11.5)
        Xs, ress, its, ist = Ss.solve_lsqr()
        lhs, rhs, lhs0 = Ss.V3(Xs)
        mx, tail = Ss.stats(Xs)
        dk = Ss.dk0(Xs)
        pr(f"  (C2') u_T kernel: LSQR its {its} (istop {ist}), residual (V2) {ress:.1e}  [{time.time() - t0:.0f}s]")
        pr(f"     V3 graviton pole: int X K[D] = {lhs.real:.12e} {lhs.imag:+.1e}i ; int psi D = {rhs:.12e} ; rel diff {abs(lhs - rhs)/rhs:.1e}"
           f"  (without inversion: rel {abs(lhs0 - rhs)/rhs:.1e})")
        jm = np.argmax(np.abs(Xs - Ss.phi)*(Ss.nu <= Lam))
        pr(f"     dk0 = {dk.real:+.5e} {dk.imag:+.2e}i (eta-kernel: {Se.dk0(Xe).real:+.5e});  max|X - psi|/max psi on [0,Lam] = {mx:.2e} at nu = "
           f"{Ss.nu[jm]:.2f} (eta: {Se.stats(Xe)[0]:.2e}); same on nu >= 20: {np.max(np.abs(Xs - Ss.phi)[(Ss.nu >= 20) & (Ss.nu <= Lam)])/Ss.phi.max():.2e};"
           f" tail mass (nu > Lam) {tail:.1e}")
        pr(f"     => A, B_AdS factor 1/(1 + Re dk0) = {1/(1 + dk.real):.6f}; bounds factor {(1 + dk.real)**-0.5:.6f}  (h1 convention, see §23.2)")
        # large-nu symbol: X - psi ~ psi''/2 - (st psi)', st = -i u
        p = Ss.nu/Lam; ph = Ss.phi; dph = Pl.polyval(p, Pl.polyder(cf, 1))/Lam; d2 = Pl.polyval(p, Pl.polyder(cf, 2))/Lam**2
        st = -1j*Ss.u; dst_ = np.gradient(st, Ss.nu)
        sym = 0.5*d2 - (dst_*ph + st*dph)
        pr("     interior p = 0.2/0.4/0.6/0.8: (X - psi)/psi = " + " ".join(
            f"{((Xs[j] - ph[j])/ph[j]).real:+.3e}{((Xs[j] - ph[j])/ph[j]).imag:+.1e}i" for j in [int(q*Lam/dnu) for q in (0.2, 0.4, 0.6, 0.8)])
           + " ; symbol = " + " ".join(f"{(sym[j]/ph[j]).real:+.3e}{(sym[j]/ph[j]).imag:+.1e}i" for j in [int(q*Lam/dnu) for q in (0.2, 0.4, 0.6, 0.8)]))
        msk = (Ss.nu >= 10) & (Ss.nu <= Lam)
        pi1 = np.sum(((0.5*d2 - (dst_*ph + st*dph))*p**6)[msk])*dnu/Lam*(-2)
        pi0 = np.sum(((Xs - ph)*p**6)[msk])*dnu/Lam*(-2)
        pr(f"     dk0 restricted to nu >= 10: exact {pi0.real:+.5e}, first-order symbol {pi1.real:+.5e}")
        prof = [0.1, 0.37, 0.73, 1, 1.5, 2, 4.4, 6.35, 10, 20, 0.25*Lam, 0.5*Lam, 0.75*Lam, 0.95*Lam]
        pr("     profile nu: X/max psi (Re, Im) | psi/max psi | X_eta/max psi")
        for v in prof:
            j = int(v/dnu)
            pr(f"       {Ss.nu[j]:9.3f}: {Xs[j].real/ph.max():+.5e} {Xs[j].imag/ph.max():+.2e}i | {ph[j]/ph.max():+.5e} | {Xe[j].real/ph.max():+.5e}")
        # ---- gauge spread: the 3 null vectors (interpolated), normalised to max |X_h| = max psi
        for kk in range(3):
            vh = np.interp(Ss.nu, nu0, V0[kk].real, right=0) + 1j*np.interp(Ss.nu, nu0, V0[kk].imag, right=0)
            vh = vh/np.max(np.abs(vh))*ph.max()
            r_h = np.linalg.norm(Ss.A(vh))/np.linalg.norm(Ss.Cf(vh))
            lh = Ss.V3(Xs + vh)[0]
            pr(f"     null vector {kk + 1} (max|X_h| = max psi): |A X_h|/|C X_h| = {r_h:.1e}; d(dk0) = {abs(Ss.dk0(Xs + vh) - dk):.1e}; "
               f"V3 shift {abs(lh - lhs)/rhs:.1e}")
        out = dict(nu=Ss.nu, X=Xs, phi=Ss.phi, Xeta=Xe, Lam=Lam, dk0=dk, V3=(lhs, rhs), res=ress)
        # ---- (R) refinement
        for (lab, dn2, N2) in (("dnu/2", dnu/2, 2*N), ("nu_max x2", dnu, 2*N)):
            t0 = time.time()
            S2 = Solver(cf, Lam, dn2, N2, u_spin, 11.5)
            X2, r2, i2, _ = S2.solve_lsqr()
            Xi = np.interp(Ss.nu, S2.nu, X2.real) + 1j*np.interp(Ss.nu, S2.nu, X2.imag)
            m20 = (Ss.nu <= Lam) & (Ss.nu >= 20)
            d = np.max(np.abs(Xi - Xs)[Ss.nu <= Lam])/np.max(Ss.phi); d20 = np.max(np.abs(Xi - Xs)[m20])/np.max(Ss.phi)
            l2, r2_, _ = S2.V3(X2)
            pr(f"  (R) {lab}: its {i2}, residual {r2:.1e}; max|dX|/max psi = {d:.1e} (nu >= 20: {d20:.1e}); dk0 = {S2.dk0(X2).real:+.6e}; "
               f"V3 rel {abs(l2 - r2_)/r2_:.1e}  [{time.time() - t0:.0f}s]")
        pickle.dump(out, open(os.path.join(HERE, 'cache', 'v4_chi_%d.pkl' % Dg), 'wb'))
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(HERE, 'cache', 'v4_chi_' + '_'.join(map(str, Dgs)) + '.txt'), 'w').write(buf.getvalue())
