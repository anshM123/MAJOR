"""
u3 (EE_FUNCTIONAL.md §3): the non-local (Abel-form) weight of C^'_EE = A_h^{-1}[Psi_2[EE]] when the forward normalisation h_- + h_+ has REAL
zeros (nu = 1.17151, 6.64286, exact roots of u2's rational R_-+^EE).

The Abel equation of H7 §23.2,  E[X] = cosh t C[X] - i sinh t S[u X] = C[psi_2],  X = chi (h_- + h_+),  u X = chi (h_- - h_+),
is rewritten in the chi-normalisation (as v4_chi (N2) did for the scalar g):  Y := chi n,  n := -(|R_-| + |R_+|) (> 0 up to the overall sign
that makes X -> psi_2 at large nu),  a := (R_- + R_+)/n,  b := (R_- - R_+)/n  (bounded, |a|, |b| <= 1):
        A[Y] := C[a Y] - i tanh(t) S[b Y] = Phi/cosh t,        X = a Y   (vanishes automatically at the real zeros of h_- + h_+).
Nothing is divided by h_- + h_+; the operator is bounded.  Solved by LSQR (min-norm correction Y - psi_2: chi-normalisation gauge).
Checks: (N) dense SVD of the operator (u only): number of vanishing singular values vs (w-1)/2 = 2 (u2: w = 5) and grid refinement;
(V2) residual; (V3) graviton-pole identity int Y [R_-/n D(nu+i) + R_+/n D(nu-i)] = int psi_2 D, D = 1/(nu^2 + 9/4);
(K) dk0 = -2 int (Re X - psi_2) p^6 dp (h1/v4 forward convention) and the Y-convention value; interior vs the first-order symbol;
(Z) X near the real zeros; refinement dnu/2.  Also q2 through the SAME code (reproduces v4_chi's dk0 up to the gauge convention).
Usage: python3 u3_chi_EE.py [Dg ...]  -> cache/u3_chi_EE.txt, cache/u3_chi_<struct>_<Dg>.pkl
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np
from numpy.polynomial import polynomial as Pl
import v4_chi as C
import u2_hclosed_struct as U2

CACHE = os.path.join(HERE, 'cache')
_F = None


def Rfun(name, lam=0.0, base='q2'):
    """(R_-, R_+) callables for structure `name` (+ lam * base)."""
    global _F
    if _F is None:
        _F = U2.numfuncs()
    a, b = _F[name]; c, d = _F[base]
    if lam == 0.0:
        return a, b
    return (lambda v: a(v) + lam*c(v)), (lambda v: b(v) + lam*d(v))


class SolverAB(C.Solver):
    def __init__(self, coeffs, Lam, dnu, N, Rm, Rp):
        super().__init__(coeffs, Lam, dnu, N, lambda v: np.zeros_like(np.asarray(v, float), dtype=complex), 0.5)
        hm, hp = Rm(self.nu.astype(complex)), Rp(self.nu.astype(complex))
        qm, qp = Rfun('q2')
        nq = qm(self.nu.astype(complex)) + qp(self.nu.astype(complex))
        self.n = (np.abs(hm) + np.abs(hp))*nq/np.abs(nq)      # chi-normalisation with the (smooth, zero-free) phase of q^2's forward norm
        self.hm, self.hp = hm/self.n, hp/self.n
        self.a, self.bb = self.hm + self.hp, self.hm - self.hp

    def A(self, Y):
        return self.Cf(self.a*Y) - 1j*self.th*self.Sf(self.bb*Y)

    def AH(self, r):
        return np.conj(self.a)*self.CT(r) + 1j*np.conj(self.bb)*self.ST(self.th*r)

    def V3Y(self, Y):
        nu = self.nu; w = self.dnu
        D = lambda z: 1/(z*z + 2.25)
        lhs = np.sum(w*Y*(self.hm*D(nu + 1j) + self.hp*D(nu - 1j)))
        rhs = np.sum(w*self.phi*D(nu))
        lhs0 = np.sum(w*self.phi*(self.hm*D(nu + 1j) + self.hp*D(nu - 1j)))
        return lhs, rhs, lhs0


def null_sv(Rm, Rp, dnu, N, k=6):
    S = SolverAB(np.zeros(10), 1.0, dnu, N, Rm, Rp)
    I = np.eye(N)
    M = np.array([S.A(I[:, j]) for j in range(N)]).T
    U, s, Vh = np.linalg.svd(M)
    return S.nu, s, Vh[-k:].conj()[::-1]


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    Dgs = [int(a) for a in sys.argv[1:]] or [100, 300, 1000]
    t00 = time.time()
    pr("# u3: Abel-form weight of C'_EE in the chi-normalisation (bounded operator; X = a Y vanishes at the real zeros of h_- + h_+)")
    for nm in (() if os.environ.get('U3_SKIPN') else ('EE', 'q2')):
        Rm, Rp = Rfun(nm)
        for (dn, NN) in ((0.05, 2048), (0.025, 4096)):
            t0 = time.time()
            nu0, s0, V0 = null_sv(Rm, Rp, dn, NN)
            pr(f"(N) {nm}: dense SVD dnu = {dn}, N = {NN}: smallest sv " + " ".join(f"{v:.1e}" for v in s0[::-1][:6]) +
               f"; largest {s0[0]:.2f}  [{time.time() - t0:.0f}s]")
            if nm == 'EE' and dn == 0.025:
                for kk in range(3):
                    a = np.abs(V0[kk]); j = np.argmax(a); supp = nu0[a > 1e-3*a.max()]
                    pr(f"     vector {kk + 1} (sv {s0[::-1][kk]:.1e}): peak nu = {nu0[j]:.3f}; |Y_h| > 1e-3 max on [{supp.min():.2f}, {supp.max():.2f}]; "
                       f"|Y_h(20)|/max = {a[int(20/dn)]/a.max():.1e}")
                pickle.dump(dict(nu=nu0, s=s0, V=V0), open(os.path.join(CACHE, 'u3_null_EE.pkl'), 'wb'))
    for Dg in Dgs:
        cf, Lam = C.load_weights(Dg)
        dnu, N = C.SETTINGS[Dg]
        pr(f"\n## Delta_g = {Dg}: Lam = {Lam:.4f}; dnu = {dnu}, N = {N}")
        for nm in ('EE', 'q2'):
            t0 = time.time()
            Rm, Rp = Rfun(nm)
            S = SolverAB(cf, Lam, dnu, N, Rm, Rp)
            Y, res, its, ist = S.solve_lsqr()
            X = S.a*Y
            lhs, rhs, lhs0 = S.V3Y(Y)
            p = S.nu/Lam
            dkX = -2*np.sum((X - S.phi)*p**6)*dnu/Lam
            dkY = -2*np.sum((Y - S.phi)*p**6)*dnu/Lam
            pr(f"  {nm}: LSQR its {its} (istop {ist}), residual (V2) {res:.1e}  [{time.time() - t0:.0f}s]")
            pr(f"     V3 graviton pole: {lhs.real:.12e} {lhs.imag:+.1e}i vs {rhs:.12e}: rel {abs(lhs - rhs)/rhs:.1e} (no inversion: {abs(lhs0 - rhs)/rhs:.1e})")
            pr(f"     dk0 (forward/X convention) = {dkX.real:+.5e} {dkX.imag:+.1e}i;  (Y convention) = {dkY.real:+.5e} {dkY.imag:+.1e}i;  "
               f"A, B_AdS factor 1/(1+dk0) = {1/(1 + dkX.real):.6f}, bounds factor {(1 + dkX.real)**-0.5:.6f}")
            m = S.nu <= Lam
            pr(f"     max|X - psi|/max psi on [0, Lam] = {np.max(np.abs(X - S.phi)[m])/S.phi.max():.2e}; on nu >= 20: "
               f"{np.max(np.abs(X - S.phi)[m & (S.nu >= 20)])/S.phi.max():.2e}; tail (nu > Lam) {np.sum(np.abs(X)[~m])/np.sum(np.abs(S.phi)):.1e}")
            # first-order symbol in the interior (u = b/a)
            u = S.bb/S.a
            ph = S.phi; dph = Pl.polyval(p, Pl.polyder(cf, 1))/Lam; d2 = Pl.polyval(p, Pl.polyder(cf, 2))/Lam**2
            st = -1j*u; dst_ = np.gradient(st, S.nu)
            sym = 0.5*d2 - (dst_*ph + st*dph)
            js = [int(q*Lam/dnu) for q in (0.2, 0.4, 0.6, 0.8)]
            pr("     interior p = 0.2/0.4/0.6/0.8: (X - psi)/psi = " + " ".join(f"{((X[j] - ph[j])/ph[j]).real:+.3e}" for j in js) +
               " ; symbol " + " ".join(f"{(sym[j]/ph[j]).real:+.3e}" for j in js))
            if nm == 'EE':
                pr("     near the real zeros of h_- + h_+:  nu: |a|, X/max psi, Y/max psi, psi/max psi")
                for v in (1.0, 1.15, 1.1715, 1.19, 1.4, 6.3, 6.6, 6.6429, 6.7, 7.0, 10.0, 20.0):
                    j = int(v/dnu)
                    pr(f"       {S.nu[j]:8.4f}: {abs(S.a[j]):.2e}  {X[j].real/ph.max():+.4e}{X[j].imag/ph.max():+.1e}i  "
                       f"{Y[j].real/ph.max():+.4e}{Y[j].imag/ph.max():+.1e}i  {ph[j]/ph.max():+.4e}")
                # refinement
                S2 = SolverAB(cf, Lam, dnu/2, 2*N, Rm, Rp)
                Y2, r2, i2, _ = S2.solve_lsqr()
                X2 = S2.a*Y2
                Xi = np.interp(S.nu, S2.nu, X2.real) + 1j*np.interp(S.nu, S2.nu, X2.imag)
                dk2 = -2*np.sum((X2 - S2.phi)*(S2.nu/Lam)**6)*dnu/2/Lam
                l2, r2_, _ = S2.V3Y(Y2)
                pr(f"     (R) dnu/2: residual {r2:.1e}; max|dX|/max psi {np.max(np.abs(Xi - X)[m])/ph.max():.1e} (nu >= 20: "
                   f"{np.max(np.abs(Xi - X)[m & (S.nu >= 20)])/ph.max():.1e}); dk0 {dk2.real:+.5e}; V3 rel {abs(l2 - r2_)/r2_:.1e}")
            pickle.dump(dict(nu=S.nu, X=X, Y=Y, phi=S.phi, Lam=Lam, dk0=dkX, V3=(lhs, rhs), res=res, a=S.a, b=S.bb),
                        open(os.path.join(CACHE, f'u3_chi_{nm}_{Dg}.pkl'), 'wb'))
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(CACHE, 'u3_chi_EE%s.txt' % ('_w' if os.environ.get('U3_SKIPN') else '')), 'w').write(buf.getvalue())
