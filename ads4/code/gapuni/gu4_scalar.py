"""
gu4_scalar (GAP_UNIFORM.md §4.4): chart-A test of the fixed-phase Laplace structure (Prop. S) on the EXACT scalar heavy moment
(Delta_phi = 3, k = 2; Arb Laplace/Barnes series of agproof/ag_exact, read-only import), over a long range of gaps.
At fixed (J, mu = m/Lam) and Delta_g in DGS it evaluates, exactly as gapscale/gs_edge.py,
    G(p) = ex/2 - E2h P_J(1 - 2x^2) - (P_{-1/2+i nu}(cosh beta) - J0(b p)),   x = nu/m, E2h = (2 - x^2)/(2 (1 - x^2)^2),
    Gamma(x) = m^2 G/E2h (c2), and the whole k = 2 cost S = int_0^1 psi2 G dp (Lam^2 S: chart-A law S ~ 1/Lam^2),
on a fixed x grid (linear interpolation in x of the dense p-grid values), and fits, node by node, Gamma(m) = g0 + g1/m + g2/m^2.
Reported: sup_x |Gamma| per gap (x <= 0.9 and x <= x1), max_x |Gamma_m - Gamma_m0|, and whether sup_{m >= m0} |Gamma_m(x)| is attained
at m0 (the 'monotone envelope' property that makes the 1000 value a bound for every larger gap).
Usage: python3 gu4_scalar.py J:mu[,J:mu...] [Dg,Dg,...]  -> cache/gu4_scalar_<J>_<mu>.json, cache/gu4_scalar.txt (appended)
"""
import os, sys, time, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gu_guard  # noqa: F401
os.environ.setdefault('OMP_NUM_THREADS', '1')
CODE = os.path.join(HERE, '..')
for sub in ('agproof', 'p4win'):
    sys.path.insert(0, os.path.join(CODE, sub))
import numpy as np
from scipy.special import j0, eval_legendre
import ag_exact as AE           # read-only (Arb exact scalar)
import p4w_common as C          # read-only (Pd: Mehler-Dirichlet P_{-1/2+i nu})
from fractions import Fraction

d = json.load(open(os.path.join(CODE, '..', 'prl', 'data', 'nf02_psi.json')))
K0 = float(Fraction(d['k0']))
c2 = np.array([float(Fraction(v)) for v in d['psi2_exact']])/K0
psi2 = lambda p: np.polynomial.polynomial.polyval(p, c2)
XG = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95])
CA = os.path.join(HERE, 'cache')


def one(Dg, J, mu):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    mt = mu*Lam
    al = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*mt*mt))/2
    Delta = round((al + J + 2)*64)/64
    al = Delta - J - 2.0; ab = al + 2*J + 1
    m = np.sqrt(al*ab); eta = (al*al + ab*ab)/(2*al*ab); beta = float(np.arccosh(eta)); mfl = m/Lam
    st = dict(J=J, Delta=float(Delta), m=m, m2=m*m, eta=eta, beta=beta, mfl=mfl)
    t0 = time.time()
    L3 = AE.build(st, numax=Lam*1.02, prec=400, conf=(320, 960, 0.5))
    plo = 0.002
    kmax = (2*J + 1)/mfl/np.sqrt(max(1 - (1/mfl)**2, 1e-12))
    nseg = int(np.ceil((1 - plo)*1.2*kmax/np.pi*2)) + 40
    xs, ws = np.polynomial.legendre.leggauss(16)
    br = np.linspace(plo, 1.0, nseg + 1)
    p = np.concatenate([a + (b - a)*(xs + 1)/2 for a, b in zip(br[:-1], br[1:])])
    wp = np.concatenate([ws*(b - a)/2 for a, b in zip(br[:-1], br[1:])])
    nu = Lam*p; x = nu/m
    ex = np.array([float(L3.value_fast(2, float(v)).mid()) for v in nu])/2
    E2h = (2 - x*x)/(1 - x*x)**2/2
    flat = E2h*eval_legendre(J, 1 - 2*x*x)
    P = C.Pd(beta, nu)[0].real
    G = ex - flat - (P - j0((2*J + 1)/mfl*p))
    Gam = m*m*G/E2h
    S = float(np.sum(wp*psi2(p)*G))
    xg = XG[XG <= x[-1]]
    return dict(Dg=Dg, J=J, mu=float(mfl), m=float(m), Delta=float(Delta), Lam=float(Lam), x1=float(x[-1]), S=S, L2S=float(Lam**2*S),
                xg=xg.tolist(), Gam=np.interp(xg, x, Gam).tolist(), sup09=float(np.abs(Gam[x <= 0.9]).max()),
                sup_all=float(np.abs(Gam).max()), t=time.time() - t0)


def main():
    specs = [tuple(s.split(':')) for s in sys.argv[1].split(',')]
    DGS = [float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1000., 1500., 2000., 3000., 5000., 10000.]
    lines = []
    for J, mu in specs:
        J = int(J); mu = float(mu)
        res = []
        for Dg in DGS:
            r = one(Dg, J, mu); res.append(r)
            print(f"J={J} mu={mu} Dg={Dg:.0f}: m={r['m']:.1f} x1={r['x1']:.4f} Lam^2 S={r['L2S']:+.4e} sup|G|(x<=.9)={r['sup09']:.4e} "
                  f"sup|G|(all)={r['sup_all']:.4e} ({r['t']:.0f}s)", flush=True)
            json.dump(res, open(os.path.join(CA, 'gu4_scalar_%d_%g.json' % (J, mu)), 'w'))
        # analysis
        ms = np.array([r['m'] for r in res]); xg = res[0]['xg']
        n = min(len(r['xg']) for r in res); xg = xg[:n]
        GG = np.array([r['Gam'][:n] for r in res])
        A = np.array([np.ones_like(ms), 1/ms, 1/ms**2]).T
        coef, *_ = np.linalg.lstsq(A, GG, rcond=None)
        resid = np.abs(A @ coef - GG).max(axis=0)
        lines.append(f"## exact scalar, J = {J}, mu = {mu}: gaps {[int(r['Dg']) for r in res]}, m = {[round(v) for v in ms]}")
        lines.append("   Lam^2 S = " + " ".join(f"{r['L2S']:+.4e}" for r in res))
        lines.append("   sup|Gamma| x<=0.9: " + " ".join(f"{r['sup09']:.3e}" for r in res) + "   (all x <= x1: " +
                     " ".join(f"{r['sup_all']:.3e}" for r in res) + ")")
        lines.append(f"   {'x':>5} " + " ".join(f"{'G(' + str(int(r['Dg'])) + ')':>11}" for r in res) +
                     f" | {'g0':>11} {'g1':>11} {'g2':>11} {'fit res':>8} | {'sup_m>=m0/|G_m0|':>16}")
        for i, xx in enumerate(xg):
            col = GG[:, i]
            env = np.abs(col).max()/max(abs(col[0]), 1e-300)
            lines.append(f"   {xx:5.2f} " + " ".join(f"{v:+11.4e}" for v in col) +
                         f" | {coef[0, i]:+11.4e} {coef[1, i]:+11.4e} {coef[2, i]:+11.4e} {resid[i]:8.1e} | {env:16.3f}")
    for l in lines:
        print(l)
    open(os.path.join(CA, 'gu4_scalar.txt'), 'a').write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main()
