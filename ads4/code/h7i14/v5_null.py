"""v5 (H7_I14_LIGHTSIDE.md §24): the null directions of the Abel-form weight of C^_2' -- theory checks.

Statement tested (§24.1, Lemma N).  Let X_h be a null vector of the Abel operator (E[X_h] = 0), f_h(nu) := X_h (1+u)/2 (nu > 0),
f_h(-nu) := X_h (1-u)/2.  Then (i) f_h continues analytically to the strip -2 < Im nu < 0 and is odd about -i; (ii) for every even Y
analytic on the closed strip |Im nu| <= 1 (with growth slower than the decay of F_h there)
     I_h[Y] := int_0^oo X_h(nu) Ah[Y](nu) dnu = 0,      Ah[Y] := [(1+u) Y(nu+i) + (1-u) Y(nu-i)]/2   (forward-normalised A_h),
and (iii) if Y has simple poles in the strip, I_h[Y] = residue terms; for Y = 1/(nu^2 + a^2), 0 < a < 1:  I_h = -(pi/a) F_h(-i(1-a)).
Checks:
  (T0) null vectors re-projected on a finer grid (LSQR projection), size at nu -> 0, decay, moments int nu^k X_h.
  (T1) I_h[Y] for Y = P_l(eta) (null property), 1/(nu^2+a^2) for a = 0.3 ... 3, exp(-nu^2/s^2), cos(b nu): zero iff Y is strip-analytic.
  (T2) residue formula (iii) with F_h(-i beta) from the Fourier continuation of f_h.
  (T3) the generic (non-range) test functions nu^2 exp(-nu^2/s^2), x-polynomials (the crystal LO shape): I_h != 0.
  (T4) scalar-g kernel (old H0 §24.1 C^_2): its 1-dim null vector, the same tests; and its pairing with the SPINNING leading-Regge
       kernel K^T_nu(eta) (i.e. the old definition applied to <TTTT>), which is not in the range of A_g.
  (T5) index: w = Delta arg rho / pi on (0, oo) = Z_+(rho) - P_+(rho) (zeros minus poles of the rational rho in Im nu > 0), per root.
Usage: python3 v5_null.py  -> cache/v5_null.txt, cache/v5_null.pkl
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.sparse.linalg import LinearOperator, lsqr
import v4_chi as C

buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def proots(expr):
    """roots with multiplicity of a polynomial in nu (square-free factorisation, then nroots)."""
    _, fl = sp.sqf_list(sp.Poly(sp.expand(expr), nsym))
    out = []
    for f, mult in fl:
        out += [complex(r) for r in f.nroots(n=20, maxsteps=200)]*mult
    return out


def funcs():
    D = pickle.load(open(os.path.join(HERE, 'cache', 'v4_hclosed.pkl'), 'rb'))
    E = {k: eval(v, {**sp.__dict__}) for k, v in D.items()}
    n = sp.Symbol('nu')
    F = {k: sp.lambdify(n, v, 'mpmath' if False else 'numpy') for k, v in E.items()}
    return E, F


E_, F_ = funcs()
nsym = sp.Symbol('nu')
uT_sym = sp.cancel(sp.together((E_['RmT'] - E_['RpT'])/(E_['RmT'] + E_['RpT'])))
usc_sym = sp.cancel(sp.together((E_['RmS'] - E_['RpS'])/(E_['RmS'] + E_['RpS'])))
uT = sp.lambdify(nsym, uT_sym, 'numpy')
usc = sp.lambdify(nsym, usc_sym, 'numpy')


def project_null(S, V, tol=1e-14):
    """project vectors V (rows) onto the numerical null space of S.A by removing their row-space component (LSQR min-norm)."""
    N = S.N
    op = LinearOperator((N, N), matvec=S.A, rmatvec=S.AH, dtype=complex)
    out = []
    for v in V:
        x = lsqr(op, S.A(v), atol=tol, btol=tol, iter_lim=20000)[0]
        out.append(v - x)
    Q, _ = np.linalg.qr(np.array(out).T)
    return Q.T


def legP(lv, eta):
    """P_lv(eta), complex degree, eta >= 1."""
    return complex(mp.hyp2f1(-lv, lv + 1, 1, (1 - mp.mpf(eta))/2))


def Kfun(nu, eta, u):
    """forward-normalised leading-Regge kernel 1/2[(1+u) P_{l-1} + (1-u) P_{l+1}], l = -1/2 + i nu."""
    out = np.empty(nu.size, complex)
    for j, v in enumerate(nu):
        l = -0.5 + 1j*v
        out[j] = 0.5*((1 + u[j])*legP(l - 1, eta) + (1 - u[j])*legP(l + 1, eta))
    return out


def Ahat(Y, nu, u):
    return 0.5*((1 + u)*Y(nu + 1j) + (1 - u)*Y(nu - 1j))


def Fcont(nu, X, u, betas, tmax=16.0, nt=4001):
    """F(-i beta) = (1/2pi) int dt fhat(t) e^{beta t}, fhat(t) = int_R f e^{-i nu t}, f = X(1+u)/2 (nu>0), X(1-u)/2 at -nu."""
    dnu = nu[1] - nu[0]
    fp, fm = X*(1 + u)/2, X*(1 - u)/2
    t = np.linspace(-tmax, tmax, nt); dt = t[1] - t[0]
    fh = np.array([dnu*(np.sum(fp*np.exp(-1j*nu*tt)) + np.sum(fm*np.exp(1j*nu*tt))) for tt in t])
    return {b: np.sum(fh*np.exp(b*t))*dt/(2*np.pi) for b in betas}, (t, fh)


def main():
    t00 = time.time()
    pr("# v5_null: null directions of the Abel-form weight of C^_2' (H7_I14_LIGHTSIDE.md §24)")
    # ---------------- (T5) index decomposition
    rho = sp.cancel(E_['RpT']/E_['RmT']); num, den = sp.fraction(sp.together(rho))
    zr = proots(num)
    pl = proots(den)
    cz = [-np.angle(-z)/np.pi for z in zr]; cp = [np.angle(-p)/np.pi for p in pl]
    pr("(T5) index of rho = h_+/h_- (rational): w = Delta_(0,oo) arg rho / pi = sum over roots of -arg(-z)/pi (zeros), +arg(-p)/pi (poles)")
    pr("     zeros: " + ", ".join(f"{z.real:+.4f}{z.imag:+.4f}i [{c:+.3f}]" for z, c in sorted(zip(zr, cz), key=lambda q: q[0].imag)))
    pr("     poles: " + ", ".join(f"{p.real:+.4f}{p.imag:+.4f}i [{c:+.3f}]" for p, c in sorted(zip(pl, cp), key=lambda q: q[0].imag)))
    Zp = sum(1 for z in zr if z.imag > 0); Pp = sum(1 for p in pl if p.imag > 0)
    pr(f"     sum = {sum(cz) + sum(cp):+.6f};  Z_+ - P_+ = {Zp} - {Pp} = {Zp - Pp}  (roots symmetric under z -> -conj z, equal degrees => w = Z_+ - P_+)")
    pr(f"     in the strip 0 < Im nu <= 1: zeros {sum(1 for z in zr if 0 < z.imag <= 1)}; in Im nu > 1: zeros {sum(1 for z in zr if z.imag > 1)}"
       f" -> the index is NOT a strip quantity (the triple zero at 3i/2 and the zero at 5i/2 count)")
    rs = sp.cancel(E_['RpS']/E_['RmS']); ns_, ds_ = sp.fraction(sp.together(rs))
    zs = proots(ns_); ps = proots(ds_)
    pr(f"     scalar g: zeros {[f'{z.imag:+.2f}i' for z in zs]}, poles {[f'{p.imag:+.2f}i' for p in ps]} -> Z_+ - P_+ = "
       f"{sum(z.imag > 0 for z in zs) - sum(p.imag > 0 for p in ps)}; eta P_l: rho = (2nu-i)/(2nu+i) -> 1")
    # ---------------- (T0) null vectors, refined
    d = pickle.load(open(os.path.join(HERE, 'cache', 'v4_chi_null.pkl'), 'rb'))
    nu0, V0 = d['nu'], d['V'][:3]
    dnu, N = 0.0125, 8192
    S = C.Solver(np.zeros(10), 1.0, dnu, N, C.u_spin, 11.5)
    Vi = np.array([np.interp(S.nu, nu0, v.real, right=0) + 1j*np.interp(S.nu, nu0, v.imag, right=0) for v in V0])
    t0 = time.time()
    Vr = project_null(S, Vi)
    nA = np.linalg.norm(S.A(np.exp(-S.nu**2/50)))/np.linalg.norm(np.exp(-S.nu**2/50))
    pr(f"(T0) null vectors projected on dnu = {dnu}, N = {N} (nu_max = {N*dnu:.0f}) [{time.time() - t0:.0f}s]: |A X_h|/|X_h| = " +
       ", ".join(f"{np.linalg.norm(S.A(v)):.1e}" for v in Vr) + f" (a smooth non-null vector: {nA:.2f}); before projection " +
       ", ".join(f"{np.linalg.norm(S.A(v))/np.linalg.norm(v):.1e}" for v in Vi))
    # overlap of refined span with the interpolated one
    Qi, _ = np.linalg.qr(Vi.T)
    sv = np.linalg.svd(Qi.conj().T @ Vr.T, compute_uv=False)
    pr(f"     principal cosines between interpolated and projected spans: {', '.join(f'{v:.8f}' for v in sv)}")
    u = C.u_spin(S.nu)
    # a real-structured basis: the 3 vectors ordered by peak position, normalised to max |X_h| = 1
    pk = [S.nu[np.argmax(np.abs(v))] for v in Vr]
    Vr = Vr[np.argsort(pk)[::-1]]
    Vn = np.array([v/v[np.argmax(np.abs(v))] for v in Vr])
    for k, v in enumerate(Vn):
        a = np.abs(v)
        pr(f"     X_{k + 1}: peak nu = {S.nu[np.argmax(a)]:.3f}; X(nu)/nu^2 at nu = 0.006, 0.019, 0.031: " +
           " ".join(f"{abs(v[j])/S.nu[j]**2:.3f}" for j in range(3)) + f"; |X| at nu = 20, 30, 50: {a[int(20/dnu)]:.1e} {a[int(30/dnu)]:.1e} "
           f"{a[int(50/dnu)]:.1e}; L1 = {np.sum(a)*dnu:.3f}; moments int nu^k X (k = 0, 1, 2, 4, 6): " +
           " ".join(f"{abs(np.sum(v*S.nu**kk)*dnu):.2e}" for kk in (0, 1, 2, 4, 6)))
    # ---------------- (T1)-(T3)
    nu = S.nu
    pr("(T1) I_h[Y] = int X_h Ah[Y] (X_h normalised to max|X_h| = 1); 'scale' = int |X_h| |Ah[Y]|")
    tests = []
    for eta in (1.0, 1.5, 3.0, 10.0):
        sel = nu <= 60
        Kv = np.zeros(nu.size, complex); Kv[sel] = Kfun(nu[sel], eta, u[sel])
        tests.append((f"P_l(eta = {eta})  [MF, entire]", Kv))
    for a in (0.3, 0.6, 0.9, 1.1, 1.5, 3.0):
        tests.append((f"1/(nu^2 + {a}^2)  [{'poles IN strip' if a < 1 else 'analytic on strip'}]", Ahat(lambda z, a=a: 1/(z*z + a*a), nu, u)))
    for s in (2.0, 5.0):
        tests.append((f"exp(-nu^2/{s}^2)  [entire]", Ahat(lambda z, s=s: np.exp(-z*z/(s*s)), nu, u)))
    for b in (0.5, 2.0):
        tests.append((f"cos({b} nu)  [entire]", Ahat(lambda z, b=b: np.cos(b*z), nu, u)))
    tests.append(("Ah[(nu^2+1/4) exp(-nu^2/25)]  [entire]", Ahat(lambda z: (z*z + 0.25)*np.exp(-z*z/25), nu, u)))
    gen = [("nu^2 exp(-nu^2/25)  [generic, not Ah of anything]", nu**2*np.exp(-nu**2/25)),
           ("crystal LO (2 - x^2) P_8(1 - 2x^2), x = nu/100  [generic]", (2 - (nu/100)**2)*np.polynomial.legendre.legval(1 - 2*(nu/100)**2, np.eye(9)[8])),
           ("1/(nu^2 + 9/4)  [graviton pole at the C level, NOT Ah of it]", 1/(nu**2 + 2.25)),
           ("Ah[1/(nu^2 + 9/4)]  [graviton pole, Psi level]", Ahat(lambda z: 1/(z*z + 2.25), nu, u))]
    res = {}
    for lab, Yv in tests + gen:
        I = np.array([np.sum(v*Yv)*dnu for v in Vn]); sc = np.array([np.sum(np.abs(v*Yv))*dnu for v in Vn])
        res[lab] = I
        pr(f"     {lab:62s}: |I_h|/scale = " + " ".join(f"{abs(i)/s_:.1e}" for i, s_ in zip(I, sc)))
    # ---------------- (T2) residue formula
    betas = [1 - a for a in (0.3, 0.6, 0.9)]
    pr("(T2) residue formula for Y = 1/(nu^2 + a^2), a < 1:  I_h = -(pi/a) F_h(-i(1-a)) (F_h by Fourier continuation, t <= 16)")
    for k, v in enumerate(Vn):
        Fc, _ = Fcont(nu, v, u, betas)
        for a in (0.3, 0.6, 0.9):
            I = res[f"1/(nu^2 + {a}^2)  [poles IN strip]"][k]; P = -(np.pi/a)*Fc[1 - a]
            pr(f"     X_{k + 1}, a = {a}: I_h = {I.real:+.6e}{I.imag:+.6e}i ; predicted {P.real:+.6e}{P.imag:+.6e}i ; rel diff {abs(I - P)/abs(I):.1e}")
    # (T2b) values of F_h at the zeros of h_- in the lower strip (poles of chi_h = f_h/h_-) and at -conj of zeros of h_+ in the upper strip
    pts = [-0.5j, -0.29280j]
    pr("     chi_h = f_h/h_- has poles at the zeros of h_- in -1 <= Im < 0 unless F_h vanishes there: |F_h(-i/2)|, |F_h(-0.2928i)| = " +
       "; ".join(" ".join(f"{abs(Fcont(nu, v, u, [-p.imag])[0][-p.imag]):.2e}" for p in pts) for v in Vn))
    pickle.dump(dict(nu=nu, V=Vn, dnu=dnu), open(os.path.join(HERE, 'cache', 'v5_null.pkl'), 'wb'))
    # ---------------- (T6) can 'chi analytic in the strip' fix the gauge?  chi = f/h_- (below the positive axis), f(-nu)/h_+ (above).
    # Poles of chi in {Re nu >= 0, |Im nu| <= 1} sit at the zeros of h_- (-i/2, -0.2928i) and of h_+ (i/2 double, 1.7788+0.6489i), unless
    # F vanishes at  -i/2 (with F' = 0 there), -0.2928i, -1.7788-0.6489i.  Conditions: 4 (complex); free parameters: 3.
    zp = [(-0.5j, 0), (-0.5j, 1), (-0.29280j, 0), (-1.77877 - 0.64888j, 0)]
    def Fpts(nuv, X, uu, pts, tmax=16.0, nt=3201):
        dn = nuv[1] - nuv[0]; fp, fm = X*(1 + uu)/2, X*(1 - uu)/2
        t = np.linspace(-tmax, tmax, nt); dt = t[1] - t[0]
        fh = np.array([dn*(np.sum(fp*np.exp(-1j*nuv*tt)) + np.sum(fm*np.exp(1j*nuv*tt))) for tt in t])
        return np.array([np.sum(fh*(1j*t)**k*np.exp(1j*z*t))*dt/(2*np.pi) for z, k in pts])
    Mh = np.array([Fpts(nu, v, u, zp) for v in Vn]).T            # 4 x 3
    sM = np.linalg.svd(Mh, compute_uv=False)
    pr("(T6) 'chi analytic in the strip' as a gauge condition: values F_h(z) at z = -i/2, F_h'(-i/2), F_h(-0.2928i), F_h(-1.7788-0.6489i):")
    pr("     singular values of the 4x3 matrix [F_h(z_j)]: " + " ".join(f"{v:.2e}" for v in sM))
    for Dg in (100, 300, 1000):
        dd = pickle.load(open(os.path.join(HERE, 'cache', 'v4_chi_%d.pkl' % Dg), 'rb'))
        sel = dd['nu'] <= 60
        F0 = Fpts(dd['nu'][sel], dd['X'][sel], C.u_spin(dd['nu'][sel]), zp, tmax=14.0)/dd['phi'].max()
        c, *_ = np.linalg.lstsq(Mh, -F0, rcond=None)
        r = np.linalg.norm(F0 + Mh @ c)/np.linalg.norm(F0)
        c3, *_ = np.linalg.lstsq(Mh[:3], -F0[:3], rcond=None)
        r3 = abs((F0 + Mh @ c3)[3])/abs(F0[3])
        pr(f"     Delta_g = {Dg}: |F_0(z_j)|/max psi = " + " ".join(f"{abs(v):.2e}" for v in F0) + f"; best gauge: residual/|F_0| = {r:.2e} "
           f"(|c| = " + " ".join(f"{abs(v):.2e}" for v in c) + f"); with the 3 imaginary-axis conditions imposed exactly, the 4th misses by {r3:.2e}")
    # ---------------- (T4) scalar g
    dns, Ns = 0.025, 4096
    Ssc = C.Solver(np.zeros(10), 1.0, dns, Ns, usc, 11/26.)
    nn = np.abs(F_['RmS'](Ssc.nu)) + np.abs(F_['RpS'](Ssc.nu))
    Aa = (F_['RmS'](Ssc.nu) + F_['RpS'](Ssc.nu))/nn; Bb = (F_['RmS'](Ssc.nu) - F_['RpS'](Ssc.nu))/nn
    Msc = lambda Xt: Ssc.Cf(Aa*Xt) - 1j*Ssc.th*Ssc.Sf(Bb*Xt)
    t0 = time.time()
    I_ = np.eye(Ns)
    M = np.array([Msc(I_[:, j]) for j in range(Ns)]).T
    U_, s_, Vh = np.linalg.svd(M)
    xs = Vh[-1].conj()                        # chi-normalised null vector Xt = chi (|g|+|g~|)
    pr(f"(T4) scalar-g kernel, chi-normalisation (dnu = {dns}, N = {Ns}) [{time.time() - t0:.0f}s]: smallest sv " +
       " ".join(f"{v:.1e}" for v in s_[::-1][:3]))
    nus = Ssc.nu; a_ = np.abs(xs); xs = xs/xs[np.argmax(a_)]
    pr(f"     scalar null vector: peak nu = {nus[np.argmax(a_)]:.3f}; |Xt|>1e-3 max on [{nus[a_ > 1e-3*a_.max()].min():.2f}, "
       f"{nus[a_ > 1e-3*a_.max()].max():.2f}]; |Xt(20)|, |Xt(40)| = {a_[int(20/dns)]/a_.max():.1e}, {a_[int(40/dns)]/a_.max():.1e}")
    gm, gp = F_['RmS'](nus), F_['RpS'](nus)
    # CHMRS class for the scalar: chi even and analytic on |Im nu| <= 1 (except the nu = 0 pole term); the null vector's chi_h = f_h/g has a
    # pole at the zero -i/2 of g (and +i/2 of g~) unless F_h(-i/2) = 0.  f = chi g = Xt g/n (g -> R_-^sc: N_com is zero-free in the strip)
    dn_ = nus[1] - nus[0]; fpp, fmm = xs*gm/nn, xs*gp/nn
    tt_ = np.linspace(-16, 16, 3201); dtt = tt_[1] - tt_[0]
    fh_ = np.array([dn_*(np.sum(fpp*np.exp(-1j*nus*t_)) + np.sum(fmm*np.exp(1j*nus*t_))) for t_ in tt_])
    Fm = np.sum(fh_*np.exp(0.5*tt_))*dtt/(2*np.pi); Fref = np.sum(np.abs(fpp))*dn_
    pr(f"     scalar null vector: f_h(0+) = {abs(fpp[0]):.1e}, f_h(0-) = {abs(fmm[0]):.1e} (continuous); |F_h(-i/2)|/int|f_h| = {abs(Fm)/Fref:.3e}"
       f" -> chi_h has poles at -+i/2: the CHMRS class (chi even, analytic in the strip) excludes it")
    def Isc(Ym):      # int chi_h [g Y(nu+i) + g~ Y(nu-i)]  (normalisation of g irrelevant: N_com common)
        Ps = gm*Ym(nus + 1j) + gp*Ym(nus - 1j)
        return np.sum(xs*Ps/nn)*dns, np.sum(np.abs(xs*Ps/nn))*dns
    for lab, Ym in (("1/(nu^2+0.6^2) [in strip]", lambda z: 1/(z*z + 0.36)), ("1/(nu^2+1.5^2)", lambda z: 1/(z*z + 2.25)),
                    ("exp(-nu^2/25)", lambda z: np.exp(-z*z/25)), ("P_l(2)", None)):
        if Ym is None:
            sel = nus <= 60
            Pl = np.zeros(nus.size, complex)
            Pm = np.array([legP(-1.5 + 1j*v, 2.0) for v in nus[sel]]); Pp = np.array([legP(0.5 + 1j*v, 2.0) for v in nus[sel]])
            Pl[sel] = gm[sel]*Pm + gp[sel]*Pp
            I, sc = np.sum(xs*Pl/nn)*dns, np.sum(np.abs(xs*Pl/nn))*dns
        else:
            I, sc = Isc(Ym)
        pr(f"     scalar: I[{lab:28s}] |I|/scale = {abs(I)/sc:.1e}")
    # the old definition on <TTTT>: Psi^T in the scalar forward normalisation; pairing of the scalar null vector with the spinning
    # leading-Regge kernel h_- P_{l-1} + h_+ P_{l+1} (normalised by N_T = (h_- + h_+)/(g + g~), §19.2 (ii))
    hm, hp = F_['RmT'](nus), F_['RpT'](nus)
    for eta in (1.5, 3.0, 10.0):
        sel = nus <= 60
        KT = np.zeros(nus.size, complex); KS = np.zeros(nus.size, complex)
        Pm = np.array([legP(-1.5 + 1j*v, eta) for v in nus[sel]]); Pp = np.array([legP(0.5 + 1j*v, eta) for v in nus[sel]])
        KT[sel] = (hm[sel]*Pm + hp[sel]*Pp)/(hm[sel] + hp[sel])*(gm[sel] + gp[sel])
        KS[sel] = gm[sel]*Pm + gp[sel]*Pp
        I_T = np.sum(xs*KT/nn)*dns; I_S = np.sum(xs*KS/nn)*dns
        sc = np.sum(np.abs(xs*KT/nn))*dns
        pr(f"     old C^_2 on <TTTT> (forward normalisation (ii)): scalar null vector paired with the SPINNING leading-Regge kernel, eta = {eta}:"
           f" |I|/scale = {abs(I_T)/sc:.2e} (scalar kernel: {abs(I_S)/sc:.1e})")
    pr(f"[{time.time() - t00:.0f}s]")
    open(os.path.join(HERE, 'cache', 'v5_null.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
