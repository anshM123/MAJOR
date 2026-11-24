"""
gapscale/gs_weights.py (GAP_SCALING.md §1): how Omega^EE's k = 2 de-shift weight depends on Delta_g.

Solves Result D's k = 2 Abel equation (minimal-norm gauge) at several gaps with the INDEPENDENT solver deshift/d_abel.Abel
(read-only import; its main() is not called, nothing is written outside gapscale/cache), and compares, in the rescaled variable
p = nu/Lam,
    delta_Lam(p) := Lam^2 [Re X(Lam p) - psi2(p)]
with the first-order symbol of H7_I14 §23.2 / DG_SCAN §1.3,   S(p) = 1/2 psi2''(p) - (23/2) (psi2/p)'(p),
i.e. tests  X_Lam(nu) = psi2(nu/Lam) + Lam^-2 S(nu/Lam) + O(Lam^-4)  ("rescaling up to explicit 1/Lam^2"),
and the light-side shift dk0 = -2 int (Re X - psi2) p^6 dp  (DG_SCAN: Lam^2 dk0 -> -168 int psi2 p^4 = -68.08).
Usage: python3 gs_weights.py 300 1000 2000 5000 -> cache/weights.txt, cache/weights_<Dg>.npz
"""
import os, sys, time, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'deshift'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import d_abel as DA

OUT = []


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)


c = DA.psi2_coeffs()/DA.K0
P = np.polynomial.polynomial
psi = lambda p: P.polyval(p, c)
d1 = P.polyder(c); d2 = P.polyder(c, 2)
g = c[1:]                                            # psi2/p
gd = P.polyder(g)
S = lambda p: 0.5*P.polyval(p, d2) - 11.5*P.polyval(p, gd)


def run(Dg, dnu):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    numax = Lam + 40.0
    t0 = time.time()
    A = DA.Abel(Dg, dnu, numax)
    Y, res, it = A.solve()
    X = A.a*Y
    p = A.nu/Lam
    m = p <= 1
    dk = -2*np.sum((X.real - A.phi)[m]*p[m]**6)*dnu/Lam
    np.savez(os.path.join(HERE, 'cache', 'weights_%d.npz' % Dg), nu=A.nu, X=X, Y=Y, phi=A.phi, Lam=Lam)
    pr(f"\n## Dg = {Dg}: Lam = {Lam:.4f}, dnu = {dnu}, N = {A.N}, CGLS its {it}, rel. residual {res:.1e}, {time.time()-t0:.0f}s")
    pr(f"   dk0 = {dk:+.6e};  Lam^2 dk0 = {Lam**2*dk:+.4f}   (first-order symbol: -168 int psi2 p^4 = "
       f"{-168*np.sum(psi(p[m])*p[m]**4)*dnu/Lam:+.4f})")
    return dict(Lam=Lam, p=p, X=X, phi=A.phi, dk=dk)


def main(Dgs):
    R = {}
    for Dg in Dgs:
        R[Dg] = run(Dg, 0.05)
    pr("\n## delta_Lam(p) = Lam^2 [Re X(Lam p) - psi2(p)] vs the symbol S(p)   (max psi2 = %.4f)" % psi(np.linspace(0, 1, 2001)).max())
    ps = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.97]
    pr("   p      S(p)      " + "  ".join(f"Dg={Dg:>5}" for Dg in Dgs))
    for p0 in ps:
        row = []
        for Dg in Dgs:
            r = R[Dg]
            v = np.interp(p0, r['p'], r['Lam']**2*(r['X'].real - r['phi']))
            row.append(f"{v:+10.3f}")
        pr(f"   {p0:<5}  {S(p0):+9.3f}  " + "  ".join(row))
    pr("\n## sup-norm of the residual Lam^2 (Re X - psi2) - S on p in [p_lo, 1] (relative to max|S|), and Im X")
    pg = np.linspace(0, 1, 4001)
    Smax = np.abs(S(pg)).max()
    for Dg in Dgs:
        r = R[Dg]
        line = []
        for plo in (0.02, 0.05, 0.1, 0.2):
            sel = (r['p'] >= plo) & (r['p'] <= 1)
            res = np.abs(r['Lam']**2*(r['X'].real[sel] - r['phi'][sel]) - S(r['p'][sel])).max()/Smax
            line.append(f"p>={plo}: {res:.2e}")
        sel = r['p'] <= 1
        imx = np.abs(r['X'].imag[sel]).max()*r['Lam']**2
        line.append(f"Lam^2 max|Im X| = {imx:.3f}")
        # weight beyond Lam (should decay ~ exp(-(pi/2)(nu - Lam)))
        tail = r['p'] > 1
        nu = r['p']*r['Lam']
        t10 = np.abs(r['X'][tail & (nu > r['Lam'] + 10)]).max() if np.any(tail & (nu > r['Lam'] + 10)) else 0
        line.append(f"max|X| beyond Lam+10 = {t10:.1e}")
        pr(f"   Dg = {Dg}: " + "; ".join(line))
    pr("\n## small nu (fixed nu, not fixed p): Lam*X(nu) and Lam*(X(nu) - psi2(nu/Lam)), Re / Im; a1 = psi2'(0) = %.4f"
       % P.polyval(0, d1))
    pr("   (if X = Lam^-1 Xi(nu) + O(Lam^-3) with a Lam-independent Xi, the columns agree across gaps)")
    for nu0 in (0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0):
        row = []
        for Dg in Dgs:
            r = R[Dg]; nu = r['p']*r['Lam']
            xr = np.interp(nu0, nu, r['X'].real); xi = np.interp(nu0, nu, r['X'].imag); ph = np.interp(nu0, nu, r['phi'])
            row.append(f"{r['Lam']*xr:+9.4f}/{r['Lam']*(xr - ph):+9.4f}/{r['Lam']*xi:+9.4f}")
        pr(f"   nu={nu0:<5} " + "  ".join(row))
    pr("\n## Lam^2 dk0 across gaps (light-side normalisation shift; first-order symbol value -68.08):")
    pr("   " + "; ".join(f"{Dg}: {R[Dg]['Lam']**2*R[Dg]['dk']:+.3f}" for Dg in Dgs))
    open(os.path.join(HERE, 'cache', 'weights.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main([int(a) for a in sys.argv[1:]] or [300, 1000, 2000, 5000])
