"""
gapscale/gs_edge.py (GAP_SCALING.md §3.3): Delta_g dependence of the light-cone-EDGE correction that decides Theorem W's thinnest
premise (EDGE-R), on the EXACT SCALAR (Delta_phi = 3, k = 2; Arb Laplace/Barnes series of agproof/ag_exact, read-only import; the
module's Dg = 1000 constant is not used: each state is built from its own (Delta, J)).
For twist-edge states (Delta = Dg + J) it evaluates, as edger/er_scalar does,
    G(p) = ex/2 - E2h P_J(1 - 2x^2) - (P_{-1/2+i nu}(cosh beta) - J0(b p)),   x = nu/m, p = nu/Lam, E2h = (2 - x^2)/(2 (1 - x^2)^2),
    Gamma(x) = m^2 G / E2h  (c2 units),
and the edge integral with a fixed erfc cutoff chi(p) = erfc((p_e - p)/sigma)/2 (p_e = 0.25, sigma = 0.195: the Theorem W cutoff of
(26, 1.026 Lam)), E = int_{p_lo}^1 psi2 chi G dp, relative to the state's leading-Regge k = 2 action Omega_2 = int psi2 P dp.
Comparisons: chart A (same J, same m/Lam up to the twist-edge shift): (26, 1026) vs (26, 2026);
             chart B (same beta, same m/Lam): (26, 1026) vs (52, 2052).
Usage: [GS_PLO=0.002] [GS_TAG=_full] python3 gs_edge.py "1000:26,2000:26,2000:52,2000:26:1.0257" -> cache/edge<TAG>.txt
With GS_PLO ~ 0 the integral E is (half) the exact scalar's whole k = 2 cost beyond the hybrid, int_0^1 psi2 chi G (chi ~ 1 above p_e).
"""
import os, sys, time, json, builtins
sys.dont_write_bytecode = True
_GSC = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'))
_open = builtins.open


def _gopen(f, mode='r', *a, **k):           # guard: this script writes only under gapscale/cache
    if any(c in mode for c in 'wax+') and not os.path.realpath(os.fspath(f)).startswith(_GSC + os.sep):
        raise PermissionError('gs_edge guard: refusing to write %s' % f)
    return _open(f, mode, *a, **k)


builtins.open = _gopen
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..')
for sub in ('agproof', 'p4win'):
    sys.path.insert(0, os.path.join(CODE, sub))
import numpy as np
from scipy.special import j0, eval_legendre, erfc
import ag_exact as AE           # read-only (Arb exact scalar)
import p4w_common as C          # read-only (Pd: Mehler-Dirichlet P_{-1/2+i nu})
from fractions import Fraction

d = json.load(open(os.path.join(CODE, '..', 'prl', 'data', 'nf02_psi.json')))
K0 = float(Fraction(d['k0']))
c2 = np.array([float(Fraction(v)) for v in d['psi2_exact']])/K0
psi2 = lambda p: np.polynomial.polynomial.polyval(p, c2)
OUT = []
pr = lambda s='': (print(s, flush=True), OUT.append(s))


def run(Dg, J, pe=0.25, sg=0.195, mfl_target=None):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    Delta = Dg + J
    if mfl_target is not None:                    # chart A at fixed u: Delta (in Z/64) with m/Lam ~ mfl_target, above the twist edge
        mt = mfl_target*Lam
        al = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*mt*mt))/2
        Delta = round((al + J + 2)*64)/64
        assert Delta >= Dg + J
    al = Delta - J - 2.0; ab = al + 2*J + 1
    m = np.sqrt(al*ab); eta = (al*al + ab*ab)/(2*al*ab); beta = float(np.arccosh(eta)); mfl = m/Lam
    st = dict(J=J, Delta=float(Delta), m=m, m2=m*m, eta=eta, beta=beta, mfl=mfl)
    t0 = time.time()
    L3 = AE.build(st, numax=Lam*1.02, prec=400, conf=(320, 960, 0.5))
    tb = time.time() - t0
    x1 = 1/mfl
    # Gauss-Legendre on [p_lo, 1] in segments (oscillation ~ (2J+1)/m per unit nu)
    plo = float(os.environ.get('GS_PLO', '0.2'))      # lower end of the p-range (0.2: edge study; ~0: full k = 2 cost)
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
    cur = P - j0((2*J + 1)/mfl*p)
    G = ex - flat - cur
    ch = erfc((pe - p)/sg)/2
    ps = psi2(p)
    E = float(np.sum(wp*ps*ch*G))
    Eabs = float(np.sum(wp*np.abs(ps*ch*G)))
    Sfull = float(np.sum(wp*ps*G))                    # no cutoff: the whole k = 2 cost on [p_lo, 1]
    # leading-Regge k = 2 action of the state (full p range)
    pf, wf = np.polynomial.legendre.leggauss(400)
    pf = (pf + 1)/2; wf = wf/2
    Om2 = float(np.sum(wf*psi2(pf)*C.Pd(beta, Lam*pf)[0].real))
    Gam = m*m*G/E2h
    pr(f"\n## Dg = {Dg}, J = {J}, Delta = {Delta}: m/Lam = {mfl:.5f}, x1 = {x1:.5f}, beta = {beta:.5f}, b = {(2*J+1)/mfl:.2f}; "
       f"build {tb:.0f}s, {len(p)} nodes on p in [{plo}, 1]")
    pr(f"   edge integral E = int psi2 chi G = {E:+.4e} (sign-free {Eabs:.4e}); Omega_2 (LR k=2 action) = {Om2:+.4e}; "
       f"E/Omega_2 = {E/Om2:+.3e}; Lam^2 E = {Lam**2*E:+.4e}; no-cutoff S = int psi2 G = {Sfull:+.4e}, S/Omega_2 = {Sfull/Om2:+.3e}")
    for xx in (0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.96, 0.97):
        if xx < x[0] or xx > x[-1]:
            continue
        i = int(np.argmin(np.abs(x - xx)))
        pr(f"   x = {x[i]:.4f}: Gamma = m^2 G/E2h = {Gam[i]:+.4e} c2;  G/E2h = {G[i]/E2h[i]:+.3e}")
    return dict(Dg=Dg, J=J, mfl=mfl, beta=beta, E=E, Eabs=Eabs, S=Sfull, plo=plo, Om2=Om2, x=x.tolist(), Gam=Gam.tolist())


if __name__ == '__main__':
    res = []
    for spec in sys.argv[1].split(','):
        f = spec.split(':')                       # Dg:J  (twist edge)  or  Dg:J:m/Lam  (fixed u)
        res.append(run(int(f[0]), int(f[1]), mfl_target=float(f[2]) if len(f) > 2 else None))
        open(os.path.join(HERE, 'cache', 'edge%s.txt' % os.environ.get('GS_TAG', '')), 'w').write("\n".join(OUT) + "\n")
    json.dump(res, open(os.path.join(HERE, 'cache', 'edge%s.json' % os.environ.get('GS_TAG', '')), 'w'))
