"""
dc_mstest (DIRCV.md §3.4): pooled raw-level test of the two hypotheses at every per-wave node with a P_J twin on the same CRN.
  ms   = log(G_w/S_w) - log(G_P/S_P)    = log(waves-cv / PJ-cv)   ; CV hypothesis: E[ms] = 0
  md   = log(G_w/G_P) + [log(S_P/ref) - exact] = log(waves-dir / PJ-cv) ; DIR hypothesis: E[md] = 0
(state level; the reference and the exact mean cancel in ms; md uses the exact CV mean, shown exact to 1e-8 in dc_cvmean).
Graviton component: EE (and C3H) helicity-trace of the 2x2 record (all four components move together, dc_raw).  Jackknife errors.
Usage: python3 dc_mstest.py [extra tags in dircv/cache, comma-separated]
"""
import dc_guard, os, sys, pickle
import numpy as np
import z3_nodes as Z3, fw_scalar as FS, pw_gamma as PG

STATES = [(100, 1100.0), (64, 1064.0), (64, 364.0), (100, 400.0)]
extra = [t for t in (sys.argv[1].split(',') if len(sys.argv) > 1 else []) if t]
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

def jk(f, R):
    v0 = f(np.ones(R, bool)); j = np.array([f(np.arange(R) != i) for i in range(R)])
    return v0, np.sqrt((R - 1)/R*np.sum(np.abs(j - j.mean(0))**2, axis=0))

def load(J, D, t):
    for base in (os.path.join(dc_guard.SPIN, 'cache'), dc_guard.CACHE):
        f = os.path.join(base, 'pw_pts_J%d_D%g%s.pkl' % (J, D, t))
        if os.path.exists(f):
            return pickle.load(open(f, 'rb')), pickle.load(open(f.replace(t + '.pkl', t + 'P.pkl'), 'rb'))
    return None, None

rows = []
for (J, D) in STATES:
    FS.install(J, D)
    for t in ['_w'] + extra:
        w, P = load(J, D, t)
        if w is None:
            continue
        for nu in sorted(set(k[1] for k in w if k[0] == D)):
            if (D, nu) not in P or (Z3.DREF, nu) not in w:
                continue
            sw, sp, rf = PG.to_z2(w[(D, nu)], 'EE'), PG.to_z2(P[(D, nu)], 'EE'), PG.to_z2(w[(Z3.DREF, nu)], 'EE')
            R = len(sw['repsS'])
            lxs, _ = Z3.lx_exact(D, J, nu); lxr, _ = Z3.lx_exact(Z3.DREF, J, nu)
            offP = (sp['logNs'] - Z3.logbeta(D, J)) - (rf['logNs'] - Z3.logbeta(Z3.DREF, J))
            offw = (sw['logNs'] - Z3.logbeta(D, J)) - (rf['logNs'] - Z3.logbeta(Z3.DREF, J))
            ex = (lxs - lxr).real
            dw, edw = jk(lambda s: (np.log(sw['repsS'][s].mean()/rf['repsS'][s].mean()) + offw).real - ex, R)
            for kind in ('EE', 'C3H'):
                gw, gp = PG.to_z2(w[(D, nu)], kind), PG.to_z2(P[(D, nu)], kind)
                tr = lambda M: M[..., 0, 0] + M[..., 1, 1]
                def lg(s, G):
                    return np.log(tr(G['reps'][s].mean(0))).real + G['logN'].real
                ms, ems = jk(lambda s: (lg(s, gw) - np.log(sw['repsS'][s].mean()).real - sw['logNs'].real)
                             - (lg(s, gp) - np.log(sp['repsS'][s].mean()).real - sp['logNs'].real), R)
                md, emd = jk(lambda s: lg(s, gw) - lg(s, gp) + (np.log(sp['repsS'][s].mean()/rf['repsS'][s].mean()) + offP).real - ex, R)
                lt = (2*J + 1)*nu/Z3.mphys(D, J)
                rows.append(dict(J=J, D=D, t=t, nu=nu, kind=kind, lt=lt, dw=dw, edw=edw, ms=ms, ems=ems, md=md, emd=emd))
pr("# raw-level hypothesis test (state level; jackknife).  d_w = scalar waves - exact;  ms: CV predicts 0;  md: DIR predicts 0")
pr("  state      tag   nu    l.theta0 kind   d_w (sig)            ms (sig)               md (sig)")
for r in rows:
    pr(f"  ({r['D']:g},{r['J']}) {r['t']:4s} {r['nu']:5.1f}  {r['lt']:4.1f}   {r['kind']:3s}  {r['dw']:+.5f}+-{r['edw']:.5f} ({r['dw']/r['edw']:+.1f})  "
       f"{r['ms']:+.5f}+-{r['ems']:.5f} ({r['ms']/r['ems']:+.1f})  {r['md']:+.5f}+-{r['emd']:.5f} ({r['md']/r['emd']:+.1f})")
for kind in ('EE', 'C3H'):
    for cut, lab in ((lambda r: True, 'all nodes'), (lambda r: r['lt'] <= 3.2 and r['nu'] >= 12, 'window nodes l.theta0 <= 3.2, nu >= 12'),
                     (lambda r: abs(r['dw']) > 2*r['edw'] and r['nu'] >= 12, 'nodes with |d_w| > 2 sigma, nu >= 12')):
        sel = [r for r in rows if r['kind'] == kind and cut(r) and r['nu'] >= 12]
        if not sel:
            continue
        c_cv = sum((r['ms']/r['ems'])**2 for r in sel); c_dir = sum((r['md']/r['emd'])**2 for r in sel)
        # one-parameter interpolation: graviton bias = lam * scalar bias  ->  ms = (lam - 1) d_true ~ (lam-1) d_w; fit lam by ms ~ (lam-1)*dw
        num = sum(r['ms']*r['dw']/r['ems']**2 for r in sel); den = sum(r['dw']**2/r['ems']**2 for r in sel)
        lam = 1 + num/den; elam = 1/np.sqrt(den)
        pr(f"  {kind} [{lab}] n={len(sel)}: chi2_CV = {c_cv:.1f}, chi2_DIR = {c_dir:.1f};  fitted lambda (graviton bias / scalar bias) = {lam:+.2f} +- {elam:.2f}  (CV: 1, DIR: 0)")
pickle.dump(rows, open(os.path.join(dc_guard.CACHE, 'dc_mstest.pkl'), 'wb'))
open(os.path.join(dc_guard.CACHE, 'dc_mstest.txt'), 'w').write('\n'.join(lines) + '\n')
