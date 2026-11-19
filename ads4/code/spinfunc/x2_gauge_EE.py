"""
x2 (EE_FOUNDATIONS.md §3-§4): Lemma N and the canonical (CHMRS-class) gauge for the Abel-form weight of C'_EE.

Set-up (u3 conventions).  chi-normalisation n = (|R_-| + |R_+|) x phase, Y = chi n, hm = R_-/n, hp = R_+/n, a = hm + hp, b = hm - hp,
X = a Y.  Abel equation  E(t) := cosh t C[aY] - i sinh t S[bY] = Phi(t) = C[psi_2](t).
On the real line f(nu) := chi h_- = Y hm (nu > 0), f(-nu) := chi h_+ = Y hp (nu > 0), so that E(t) = 1/2[e^t fh(t) + e^-t fh(-t)],
fh(t) = int_R f e^{-i nu t}.  For t > 0 this gives fh(t) = 2 e^-t Phi(t) - e^-2t fh(-t), and hence the EXACT representation of the
continuation F of f into -1 < Im z < 0 (no exponentially amplified Fourier tail):
   F(z) = (1/2pi) [ int_R f(nu) ( i/(nu - z) - 1/(2 - iz - i nu) ) dnu + int_0^oo psi_2(nu) ( 1/(1 - iz - i nu) + 1/(1 - iz + i nu) ) dnu ],
   F'(z) = (1/2pi) [ int_R f ( i/(nu-z)^2 - i/(2-iz-i nu)^2 ) + int_0^oo psi_2 ( i/(1-iz-i nu)^2 + i/(1-iz+i nu)^2 ) ],
(the psi_2 term is absent for null vectors, Phi = 0).

CHMRS class for EE.  rho(nu) rho(-nu) = 1 (R_+(nu) = R_-(-nu), N_com even), so chi := f/h_- is single-valued and even on R; it is analytic
on the closed strip |Im nu| <= 1 (apart from the nu = 0 term inherited from psi_2 ~ |nu|) iff F/h_- is analytic in -1 <= Im < 0, i.e. iff F
vanishes to the order of the zeros of h_- there.  h_-^EE has exactly one zero in -1 <= Im nu < 0: a DOUBLE zero at -i/2 (and no pole).
So the class is  F(-i/2) = F'(-i/2) = 0:  2 complex conditions for the 2 complex null coefficients.

Checks / outputs
 (S) synthetic validation of the F-representation: f entire, psi := f(nu - i) + f(-nu - i); F(-i/2), F'(-i/2) vs exact.
 (V) reproduction of H7 §24 numbers with the same representation: q^2 null vectors (v5_null.pkl): |F_h(-i/2)|, |F_h(-0.2928i)|; scalar g null
     vector |F_h(-i/2)|/int|f_h| = 1.141327e-3 (H7 §24.6).
 (L) Lemma N for EE: null vectors re-projected on each grid; regularity (f continuous at 0, decay); I_h[g] = int_R f_h g(nu + i) for
     strip-analytic even g (-> 0) and for 1/(nu^2+a^2), a < 1 (-> residue formula -(pi/a) F_h(-i(1-a))); generic test functions.
 (G) canonical gauge at Delta_g = 100, 300, 1000: the 2x2 matrix M = [F_hk(-i/2); F_hk'(-i/2)] (singular values), the coefficients c
     (null vectors normalised to max|X_h| = max psi_2), grid stability (dnu, dnu/2), chi near -i/2 (bounded vs ~ r^-2 in the min-norm gauge),
     gauge-independent quantities (V3, dk0), and the change of the weight X relative to the min-norm one.
Usage: python3 x2_gauge_EE.py  -> cache/x2_gauge_EE.txt, cache/x2_chi_EE_can_<Dg>.pkl   (~5-10 min, one core)
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
H7 = os.path.join(HERE, '..', 'h7i14'); sys.path.insert(0, H7)
import numpy as np
import mpmath as mp
from scipy.sparse.linalg import LinearOperator, lsqr
import v4_chi as C
import u3_chi_EE as U3

buf = io.StringIO()
CACHE = os.path.join(HERE, 'cache')


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


# ------------------------------------------------------------------ the strip continuation
def Fz(z, nu, dnu, fp, fm, phi=None, deriv=0):
    """F(z) (deriv=0) or F'(z) (deriv=1) for -1 < Im z < 0 from f on the midpoint grid (fp at +nu, fm at -nu) and psi_2 (phi) if
    the solution is inhomogeneous."""
    if deriv == 0:
        K1 = lambda v: 1j/(v - z) - 1/(2 - 1j*z - 1j*v)
        K2 = lambda v: 1/(1 - 1j*z - 1j*v) + 1/(1 - 1j*z + 1j*v)
    else:
        K1 = lambda v: 1j/(v - z)**2 - 1j/(2 - 1j*z - 1j*v)**2
        K2 = lambda v: 1j/(1 - 1j*z - 1j*v)**2 + 1j/(1 - 1j*z + 1j*v)**2
    s = dnu*(np.sum(fp*K1(nu)) + np.sum(fm*K1(-nu)))
    if phi is not None:
        s += dnu*np.sum(phi*K2(nu))
    return s/(2*np.pi)


def project_null(S, V, tol=1e-15):
    N = S.N
    op = LinearOperator((N, N), matvec=S.A, rmatvec=S.AH, dtype=complex)
    out = []
    for v in V:
        x = lsqr(op, S.A(v), atol=tol, btol=tol, iter_lim=30000)[0]
        out.append(v - x)
    return np.array(out)


def legP(lv, eta):
    return complex(mp.hyp2f1(-lv, lv + 1, 1, (1 - mp.mpf(eta))/2))


def main():
    t00 = time.time()
    pr("# x2: Lemma N and the CHMRS-class (canonical) gauge for the Abel-form weight of C'_EE (EE_FOUNDATIONS.md §3-§4)")
    # ------------------------------------------------------------ (S) synthetic validation
    pr("\n(S) synthetic check of the F-representation: f(nu) = exp(-(nu - 0.4)^2/6) (1 + 0.3 i nu) (entire), psi := f(nu-i) + f(-nu-i)")
    f0 = lambda z: np.exp(-(z - 0.4)**2/6)*(1 + 0.3j*z)
    df0 = lambda z: np.exp(-(z - 0.4)**2/6)*(0.3j - (z - 0.4)/3*(1 + 0.3j*z))
    for dnu in (0.02, 0.01):
        nu = (np.arange(int(40/dnu)) + 0.5)*dnu
        psi = f0(nu - 1j) + f0(-nu - 1j)
        for z in (-0.5j, -0.3j, 0.7 - 0.5j, -0.85j):
            Fv = Fz(z, nu, dnu, f0(nu), f0(-nu), psi); dF = Fz(z, nu, dnu, f0(nu), f0(-nu), psi, deriv=1)
            pr(f"   dnu = {dnu}: z = {z}: |F - f(z)|/|f(z)| = {abs(Fv - f0(z))/abs(f0(z)):.1e};  |F' - f'(z)|/|f'(z)| = {abs(dF - df0(z))/abs(df0(z)):.1e}")

    # ------------------------------------------------------------ (V) reproduce H7 §24 numbers
    pr("\n(V) the same representation on H7 §24's cases (forward normalisation f = X(1+u)/2 at +nu, X(1-u)/2 at -nu)")
    d5 = pickle.load(open(os.path.join(H7, 'cache', 'v5_null.pkl'), 'rb'))
    nu5, V5, dn5 = d5['nu'], d5['V'], d5['dnu']
    u5 = C.u_spin(nu5)
    for k, v in enumerate(V5):
        fp, fm = v*(1 + u5)/2, v*(1 - u5)/2
        pr(f"   q^2 X_{k + 1} (max|X| = 1): |F_h(-i/2)| = {abs(Fz(-0.5j, nu5, dn5, fp, fm)):.2e}, |F_h(-0.2928i)| = "
           f"{abs(Fz(-0.29280j, nu5, dn5, fp, fm)):.2e}   (H7 §24.2: 3e-6/5e-3/0.11 and 3e-6/3e-3/0.15)")
    # scalar g null vector (chi-normalisation), dense SVD
    F_ = C.pickle.load(open(os.path.join(H7, 'cache', 'v4_hclosed.pkl'), 'rb'))
    import sympy as sp
    E_ = {k: eval(v, {**sp.__dict__}) for k, v in F_.items()}; ns = sp.Symbol('nu')
    RmS, RpS = sp.lambdify(ns, E_['RmS'], 'numpy'), sp.lambdify(ns, E_['RpS'], 'numpy')
    t0 = time.time()
    for (dns, Ns) in ((0.05, 2048), (0.025, 4096)):
        Ssc = U3.SolverAB(np.zeros(10), 1.0, dns, Ns, RmS, RpS)
        I_ = np.eye(Ns)
        M = np.array([Ssc.A(I_[:, j]) for j in range(Ns)]).T
        _, s_, Vh = np.linalg.svd(M)
        y = Vh[-1].conj()
        fp, fm = y*Ssc.hm, y*Ssc.hp
        Fm = Fz(-0.5j, Ssc.nu, dns, fp, fm)
        pr(f"   scalar g (dnu = {dns}): smallest sv {s_[-1]:.1e}, next {s_[-2]:.2f}; |F_h(-i/2)|/int|f_h| = "
           f"{abs(Fm)/(np.sum(np.abs(fp))*dns):.6e}   (H7 §24.6: 1.141327e-3)  [{time.time() - t0:.0f}s]")
        del M, I_

    # ------------------------------------------------------------ EE: null vectors, Lemma N, canonical gauge
    Rm, Rp = U3.Rfun('EE')
    dN = pickle.load(open(os.path.join(CACHE, 'u3_null_EE.pkl'), 'rb'))
    nu0, V0 = dN['nu'], dN['V'][:2]
    pr(f"\n(L) Lemma N for EE (null vectors of u3, dnu = {nu0[1] - nu0[0]:.3f}; u3 SVD singular values {dN['s'][-1]:.1e}, {dN['s'][-2]:.1e})")
    # coarse-grid tests on the SVD grid
    S0 = U3.SolverAB(np.zeros(10), 1.0, nu0[1] - nu0[0], len(nu0), Rm, Rp)
    dn0 = S0.dnu
    Vr = project_null(S0, V0)
    Vn = []
    for v in Vr:
        X = S0.a*v; Vn.append(v/X[np.argmax(np.abs(X))])            # normalise max|X_h| = 1 (X = a Y)
    Vn = np.array(Vn)
    for k, v in enumerate(Vn):
        X = S0.a*v; fp, fm = v*S0.hm, v*S0.hp; aX = np.abs(X)
        pr(f"   Y_{k + 1}: |A Y|/|Y| = {np.linalg.norm(S0.A(v))/np.linalg.norm(v):.1e}; peak of |X| at nu = {S0.nu[np.argmax(aX)]:.3f}; "
           f"X/nu^2 at nu = {S0.nu[0]:.4f}, {S0.nu[1]:.4f}: {abs(X[0])/S0.nu[0]**2:.3f}, {abs(X[1])/S0.nu[1]**2:.3f}; f(0+), f(0-) = "
           f"{fp[0]:.2e}, {fm[0]:.2e}; |X| at nu = 20, 30, 50: {aX[int(20/dn0)]:.1e} {aX[int(30/dn0)]:.1e} {aX[int(50/dn0)]:.1e}")
    nu = S0.nu
    tests = []
    for a in (1.1, 1.5, 3.0):
        tests.append((f"1/(nu^2+{a}^2) [analytic on strip]", lambda z, a=a: 1/(z*z + a*a)))
    tests += [("exp(-nu^2/4) [entire]", lambda z: np.exp(-z*z/4)), ("exp(-nu^2/25) [entire]", lambda z: np.exp(-z*z/25)),
              ("cos(0.5 nu) [entire]", lambda z: np.cos(0.5*z)), ("cos(2 nu) [entire]", lambda z: np.cos(2*z)),
              ("1/(nu^2+9/4) [graviton pole, Psi level]", lambda z: 1/(z*z + 2.25))]
    for lab, g in tests:
        out = []
        for v in Vn:
            fp, fm = v*S0.hm, v*S0.hp
            I = dn0*(np.sum(fp*g(nu + 1j)) + np.sum(fm*g(nu - 1j))); sc = dn0*(np.sum(np.abs(fp*g(nu + 1j))) + np.sum(np.abs(fm*g(nu - 1j))))
            out.append(abs(I)/sc)
        pr(f"   I_h[{lab:42s}]: |I|/scale = " + " ".join(f"{o:.1e}" for o in out))
    # P_l(eta) through h: I = int_R f P_{-1/2+i(nu+i)}(eta)  (Mehler-Fock kernel, entire in nu)
    sel = nu <= 25
    for eta in (1.5, 3.0):
        Pm = np.array([legP(-1.5 + 1j*v, eta) for v in nu[sel]]); Pp = np.array([legP(0.5 + 1j*v, eta) for v in nu[sel]])
        out = []
        for v in Vn:
            fp, fm = (v*S0.hm)[sel], (v*S0.hp)[sel]
            I = dn0*(np.sum(fp*Pm) + np.sum(fm*Pp)); sc = dn0*(np.sum(np.abs(fp*Pm)) + np.sum(np.abs(fm*Pp)))
            out.append(abs(I)/sc)
        pr(f"   I_h[P_l(eta = {eta}) (leading-Regge heavy kernel), nu <= 25]       : |I|/scale = " + " ".join(f"{o:.1e}" for o in out))
    for a in (0.3, 0.6, 0.9):
        g = lambda z, a=a: 1/(z*z + a*a)
        for k, v in enumerate(Vn):
            fp, fm = v*S0.hm, v*S0.hp
            I = dn0*(np.sum(fp*g(nu + 1j)) + np.sum(fm*g(nu - 1j)))
            P = -(np.pi/a)*Fz(-1j*(1 - a), nu, dn0, fp, fm)
            pr(f"   poles in strip, a = {a}, Y_{k + 1}: I_h = {I.real:+.6e}{I.imag:+.6e}i; -(pi/a) F_h(-i(1-a)) = {P.real:+.6e}{P.imag:+.6e}i; "
               f"rel diff {abs(I - P)/abs(I):.1e}")
    for lab, gv in (("nu^2 exp(-nu^2/25) paired with X_h directly [generic]", nu**2*np.exp(-nu**2/25)),
                    ("1/(nu^2+9/4) paired with X_h directly [C-level pole, not A_h of it]", 1/(nu**2 + 2.25))):
        out = [abs(np.sum(S0.a*v*gv))/np.sum(np.abs(S0.a*v*gv)) for v in Vn]
        pr(f"   {lab:62s}: |I|/scale = " + " ".join(f"{o:.2f}" for o in out))
    # F_h at the double zero of h_- (-i/2) and the Lemma N (i) symmetry F_h(nu) = -F_h(-nu-2i) is outside -1<Im<0; report F, F'
    for k, v in enumerate(Vn):
        fp, fm = v*S0.hm, v*S0.hp
        pr(f"   Y_{k + 1}: F_h(-i/2) = {Fz(-0.5j, nu, dn0, fp, fm):.4e}, F_h'(-i/2) = {Fz(-0.5j, nu, dn0, fp, fm, deriv=1):.4e}  (units max|X_h| = 1)")

    # ------------------------------------------------------------ (G) canonical gauge
    pr("\n(G) canonical (CHMRS-class) gauge: F(-i/2) = F'(-i/2) = 0;  Y_can = Y_0 + c1 Y_1 + c2 Y_2  (Y_0: u3's min-norm solution;"
       " Y_k normalised to max|X_k| = max psi_2, the H7 §24 convention for |c|)")
    for Dg in (100, 300, 1000):
        cf, Lam = C.load_weights(Dg)
        dnu, N = C.SETTINGS[Dg]
        res_by_grid = []
        for fac in (1, 2):
            t0 = time.time()
            S = U3.SolverAB(cf, Lam, dnu/fac, N*fac, Rm, Rp)
            Y0, r0, _, _ = S.solve_lsqr()
            Vi = np.array([np.interp(S.nu, nu0, v.real, right=0) + 1j*np.interp(S.nu, nu0, v.imag, right=0) for v in V0])
            Vp = project_null(S, Vi)
            pm = S.phi.max()
            Yk = np.array([v/np.abs(S.a*v).max()*pm for v in Vp])
            nres = [np.linalg.norm(S.A(v))/np.linalg.norm(v) for v in Yk]
            fp0, fm0 = Y0*S.hm, Y0*S.hp
            F0 = np.array([Fz(-0.5j, S.nu, S.dnu, fp0, fm0, S.phi), Fz(-0.5j, S.nu, S.dnu, fp0, fm0, S.phi, deriv=1)])
            M = np.array([[Fz(-0.5j, S.nu, S.dnu, v*S.hm, v*S.hp, deriv=dd) for v in Yk] for dd in (0, 1)])
            sM = np.linalg.svd(M/pm, compute_uv=False)
            c = np.linalg.solve(M, -F0)
            Yc = Y0 + c @ Yk
            Xc, X0 = S.a*Yc, S.a*Y0
            F_after = [abs(Fz(-0.5j, S.nu, S.dnu, Yc*S.hm, Yc*S.hp, S.phi, deriv=dd)) for dd in (0, 1)]
            lhs, rhs, _ = S.V3Y(Yc); l0, _, _ = S.V3Y(Y0)
            p = S.nu/Lam
            dkc = -2*np.sum((Xc - S.phi)*p**6)*S.dnu/Lam; dk0 = -2*np.sum((X0 - S.phi)*p**6)*S.dnu/Lam
            resc = np.linalg.norm(S.A(Yc) - S.b)/np.linalg.norm(S.b)
            res_by_grid.append(dict(nu=S.nu, Xc=Xc, X0=X0, Yc=Yc, c=c, dnu=S.dnu))
            pr(f"  Delta_g = {Dg}, dnu = {S.dnu}: null residuals " + " ".join(f"{v:.1e}" for v in nres) +
               f"; |F_0(-i/2)|, |F_0'(-i/2)| = {abs(F0[0])/pm:.3e}, {abs(F0[1])/pm:.3e} (max psi units); sv(M)/max psi = "
               + " ".join(f"{v:.3e}" for v in sM) + f"  [{time.time() - t0:.0f}s]")
            pr(f"     c = ({c[0]:.5e}, {c[1]:.5e}), |c| = {abs(c[0]):.4f}, {abs(c[1]):.4f}; after: |F(-i/2)|, |F'(-i/2)| = "
               f"{F_after[0]/pm:.1e}, {F_after[1]/pm:.1e}; Abel residual {resc:.1e}")
            pr(f"     gauge-independent: V3 rel (canonical) {abs(lhs - rhs)/rhs:.1e} (min-norm {abs(l0 - rhs)/rhs:.1e}); dk0 canonical "
               f"{dkc.real:+.6e}{dkc.imag:+.1e}i vs min-norm {dk0.real:+.6e}{dk0.imag:+.1e}i")
            m20 = S.nu >= 20
            pr(f"     max|X_can - X_minnorm|/max psi = {np.max(np.abs(Xc - X0))/pm:.3e} (at nu = {S.nu[np.argmax(np.abs(Xc - X0))]:.2f}); on nu >= 20: "
               f"{np.max(np.abs(Xc - X0)[m20])/pm:.1e};  Re X/max psi at nu = 1.0, 3.0, 5.0, 6.0, 6.36, 9.0, 15: canonical " +
               " ".join(f"{Xc[int(v/S.dnu)].real/pm:+.3f}" for v in (1.0, 3.0, 5.0, 6.0, 6.36, 9.0, 15.0)) + " / min-norm " +
               " ".join(f"{X0[int(v/S.dnu)].real/pm:+.3f}" for v in (1.0, 3.0, 5.0, 6.0, 6.36, 9.0, 15.0)))
            # chi near -i/2: chi(z) = F(z)/R_-(z) (N_com zero-free); min-norm vs canonical
            if fac == 1:
                rows = []
                for r in (0.2, 0.1, 0.05):
                    zs = -0.5j + r*np.exp(1j*np.linspace(0, 2*np.pi, 9)[:-1])
                    zs = zs[(zs.imag < -0.02) & (zs.imag > -0.98)]
                    ch0 = max(abs(Fz(z, S.nu, S.dnu, fp0, fm0, S.phi)/Rm(z)) for z in zs)
                    chc = max(abs(Fz(z, S.nu, S.dnu, Yc*S.hm, Yc*S.hp, S.phi)/Rm(z)) for z in zs)
                    rows.append(f"r = {r}: {ch0/pm:.2e} / {chc/pm:.2e}")
                pr("     max|chi| on circles |z + i/2| = r (min-norm / canonical, max psi units): " + "; ".join(rows))
        A, B = res_by_grid
        Xi = np.interp(A['nu'], B['nu'], B['Xc'].real) + 1j*np.interp(A['nu'], B['nu'], B['Xc'].imag)
        pr(f"     grid stability dnu -> dnu/2: |dc| = {abs(A['c'][0] - B['c'][0]):.1e}, {abs(A['c'][1] - B['c'][1]):.1e}; "
           f"max|dX_can|/max psi = {np.max(np.abs(Xi - A['Xc']))/A['Xc'].real.max():.1e}")
        pickle.dump(dict(nu=A['nu'], X=A['Xc'], Y=A['Yc'], X_minnorm=A['X0'], c=A['c'], Lam=Lam),
                    open(os.path.join(CACHE, f'x2_chi_EE_can_{Dg}.pkl'), 'wb'))
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(CACHE, 'x2_gauge_EE.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
