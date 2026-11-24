"""prl/data/export_gap_weights.py -- publish the k = 2 weights of Omega^EE used at Delta_g = 2000 and 5000 (Supplement Sec. S8).

The near-gap computations at 2000 / 5000 (ads4/GAP_NEARGAP.md; code/gapnear/gn_bind.py) used the minimal-norm grid solution of
Result D's k = 2 Abel equation produced by code/gapscale/gs_weights.py (independent solver code/deshift/d_abel.Abel, dnu = 0.05,
nu_max = Lam + 40), cached as code/gapscale/cache/weights_<Dg>.npz. This script reads those caches (read-only), re-checks them,
and writes, in this directory only:
    gap/d_abel_X_2000_trunc.txt.gz, gap/d_abel_X_5000_trunc.txt.gz   (columns as d_abel_X_1000.txt; rows nu_k <= Lam + 10)
    gap/gap_weights_checks.txt                                      (the checks below, also for the 1000 cache)
    gap/SHA256SUMS
The grids are written truncated at nu <= Lam + 10, i.e. already in the reading of README.md Sec. 4, 7 (zero beyond). gzip with
mtime 0, so the hashes are reproducible.
Checks per gap: CGLS residual of the discrete equation; X = a Y; the off-grid de-shift residual |E[X](t) - Phi(t)| by direct
quadrature (d_abel.E_direct) at t = 0.3 ... 25; the graviton-pole identity V3; dk0 = -2 int (Re X - psi2) p^6 dp against the
first-order value -68.076/Lam^2; the rescaling Lam^2 (Re X - psi2) - S(p) on p >= 0.1; max |X| beyond Lam + 10; the gauge factor
max |Re X/psi2| on nu <= 25 used by the near-gap budget.
Usage: python3 export_gap_weights.py   (about a minute; numpy, scipy, sympy; imports code/deshift read-only)
"""
import os, sys, gzip, hashlib, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'deshift'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import d_abel as DA

GS = os.path.join(ADS4, 'code', 'gapscale', 'cache')
OUTD = os.path.join(HERE, 'gap')
LOG = []


def pr(s=""):
    print(s, flush=True); LOG.append(s)


P = np.polynomial.polynomial
c = DA.psi2_coeffs()/DA.K0
S = lambda p: 0.5*P.polyval(p, P.polyder(c, 2)) - 11.5*P.polyval(p, P.polyder(c[1:]))


def check(Dg):
    w = np.load(os.path.join(GS, 'weights_%d.npz' % Dg))
    nu, X, Y, phi, Lam = w['nu'], w['X'], w['Y'], w['phi'], float(w['Lam'])
    dnu = nu[1] - nu[0]
    A = DA.Abel(Dg, dnu, Lam + 40.0)
    assert A.N == len(nu) and np.allclose(A.nu, nu, rtol=0, atol=1e-9), "grid mismatch"
    pr(f"## Delta_g = {Dg}: Lam = {Lam:.10f}, dnu = {dnu:g}, N = {len(nu)} (nu_max = Lam + 40), source code/gapscale/cache/weights_{Dg}.npz")
    pr(f"   psi2 column vs fresh psi2: max diff {np.max(np.abs(A.phi - phi)):.1e};  max|X - a Y|/max|X| = "
       f"{np.max(np.abs(X - A.a*Y))/np.max(np.abs(X)):.1e}")
    res = np.linalg.norm(A.A(Y) - A.rhs)/np.linalg.norm(A.rhs)
    pr(f"   discrete equation, relative residual: {res:.1e}")
    ts = np.array([0.3, 1.0, 3.0, 6.0, 10.0, 14.0, 18.0, 22.0, 25.0])
    r = DA.E_direct(nu, dnu, X, A.b*Y, phi, ts)
    Phi = np.array([dnu*np.sum(phi*np.cos(nu*t)) for t in ts])
    pr("   off-grid de-shift residual |E[X](t) - Phi(t)| (|Phi(t)|): " +
       "; ".join(f"t={t:g}: {a:.1e} ({abs(b):.1e})" for t, a, b in zip(ts, r, Phi)))
    D = lambda z: 1/(z*z + 2.25)
    V3 = np.sum(dnu*Y*(A.hm/A.n*D(nu + 1j) + A.hp/A.n*D(nu - 1j))); V3r = np.sum(dnu*phi*D(nu))
    pr(f"   graviton-pole identity V3: {V3.real:.12e} vs {V3r:.12e} (rel {abs(V3 - V3r)/V3r:.1e})")
    p = nu/Lam; m = p <= 1
    dk = -2*np.sum((X.real - phi)[m]*p[m]**6)*dnu/Lam
    pr(f"   dk0 = {dk:+.6e}; Lam^2 dk0 = {Lam**2*dk:+.4f} (first-order value -68.076); (1 + dk0)^(-1/2) = {(1 + dk)**-0.5:.8f}")
    pg = np.linspace(0, 1, 4001); Smax = np.abs(S(pg)).max()
    sel = (p >= 0.1) & m
    pr(f"   rescaling: max |Lam^2 (Re X - psi2) - S(p)| / max|S| on p in [0.1, 1] = "
       f"{np.abs(Lam**2*(X.real[sel] - phi[sel]) - S(p[sel])).max()/Smax:.2e}")
    pr(f"   max |X| on nu > Lam + 10: {np.abs(X[nu > Lam + 10]).max():.1e} (double-precision floor; read as zero)")
    g = nu <= 25
    pr(f"   gauge factor max |Re X/psi2| on nu <= 25: {np.abs(X.real[g]/phi[g]).max():.3f}")
    return nu, X, Y, phi, Lam, dnu


def write_grid(Dg, nu, X, Y, phi, Lam, dnu):
    k = nu <= Lam + 10
    buf = io.StringIO()
    buf.write(f"# Result D k = 2 Abel weight X = chi (h_- + h_+) at Delta_g = {Dg} (code/gapscale/cache/weights_{Dg}.npz, "
              f"gs_weights.py with deshift/d_abel.Abel; NF02 psi2, EE kernel, min-norm gauge = Y - psi2 of minimal l2 norm on this "
              f"grid); Lam = {Lam:.10f}; grid nu_k = (k + 1/2) {dnu:g}; Y = X/a; functional = Re sum_k dnu X_k "
              f"Psi^_(EE,nu_k)/(h_- + h_+)(nu_k)\n")
    buf.write(f"# TRUNCATED at nu_k <= Lam + 10 = {Lam + 10:.4f} (README.md Sec. 4, 7: zero beyond); {int(k.sum())} of {len(nu)} "
              f"rows kept; dropped max|X| = {np.abs(X[~k]).max():.1e}\n")
    buf.write("# nu_k  Re X  Im X  Re Y  Im Y  psi2(nu_k)\n")
    for a in zip(nu[k], X.real[k], X.imag[k], Y.real[k], Y.imag[k], phi[k]):
        buf.write("%.3f %.17e %.17e %.17e %.17e %.17e\n" % a)
    fn = os.path.join(OUTD, f'd_abel_X_{Dg}_trunc.txt.gz')
    with open(fn, 'wb') as f:
        with gzip.GzipFile(filename='', mode='wb', fileobj=f, mtime=0, compresslevel=9) as gz:
            gz.write(buf.getvalue().encode())
    return fn


def main():
    os.makedirs(OUTD, exist_ok=True)
    pr("# export_gap_weights.py: checks of the k = 2 weights of Omega^EE used at Delta_g = 2000, 5000 (1000 cache for comparison)")
    files = []
    for Dg in (1000, 2000, 5000):
        out = check(Dg)
        if Dg != 1000:
            files.append(write_grid(Dg, *out))
        pr()
    pr("Not done here: a continuum (strip-analytic) fit of these weights with an all-t bound eps_W, as published for 300 and 1000")
    pr("(k2/k2_continuum_D_*.json). The k = 3 piece of Omega^EE (psi3, not de-shifted) and the k = 4, 5 weights of Omega^Q_up")
    pr("(functions of nu only) do not depend on the gap; Omega^Q_up was not used at 2000 or 5000.")
    fn = os.path.join(OUTD, 'gap_weights_checks.txt')
    open(fn, 'w').write("\n".join(LOG) + "\n")
    files.append(fn)
    with open(os.path.join(OUTD, 'SHA256SUMS'), 'w') as f:
        for fn in files:
            f.write(f"{hashlib.sha256(open(fn, 'rb').read()).hexdigest()}  {os.path.basename(fn)}\n")


if __name__ == '__main__':
    main()
