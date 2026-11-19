"""
t1_band (RESULT_D_CONDITIONAL.md §2-3): map of the per-state failure band of Omega^{A'} in the universal-g model of LARGEJ_OFFDIAG (r5),
for NF02 and for the re-weighted functionals of r5's LP (cache/r5_lp_*.pkl), at Delta_g = 100 / 300 / 1000.

Per state (m/Lam, even J) of a (m/Lam, b-hat) grid (twist gap as r5: alpha >= Delta_g - 2):
    F_hh/N_h, |X|/sqrt(N++ N+-), lambda_min/N_s (r3.lam), kappa_crit = sqrt(F++ F+-)/|X|
kappa_crit is the largest helicity coherence |rho_{++,+-}|/sqrt(rho_{++,++} rho_{+-,+-}) of a PSD heavy density at that state for which
tr(rho F) >= 0 (for kappa <= kappa_crit the 2x2 test passes for EVERY such rho).  Pure model [M]: r5's universal g (validated beta <= 0.11),
F_hh = flat + leading-Regge curvature (no scalar remainder, no graviton diagonal correction).

Also: exact 16 L_A(Lam) = 16 int_0^1 psi2(t) dt/(t^2 + 9/(4 Lam^2)) for each weight vector, and the bound 240 sqrt(16 L_A)/Lam^4 (k(0) = 1),
next to the log form 240 sqrt(A log Lam + B_AdS)/Lam^4.

Usage: T1_DG=300 python3 t1_band.py [nf02|<r5_lp pickle basename> ...]   -> cache/t1_band_<Dg>.txt
Read-only use of r5_universal / r3_redesign / robust.nf_lp.  One core, < 1 GB.
"""
import os, sys, io, time, pickle, glob
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.integrate import quad
import r5_universal as R5
import r3_redesign as R3

CACHE = os.path.join(HERE, 'cache')
BVALID = 0.11
NQ = int(os.environ.get('T1_NQ', '0'))          # 0: r5's 64-node rule; >0: composite GL, NQ geometric panels x 16 nodes (needed for beta >~ 0.3)


class _NG:
    def __init__(self, Dg):
        self.Dg = Dg; self.Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
        x, w = np.polynomial.legendre.leggauss(16)
        ed = np.concatenate([[0.0], np.geomspace(0.2/self.Lam, 1.0, NQ)])
        ed = np.unique(np.concatenate([ed, np.linspace(0, 1, int(os.environ.get('T1_NL', '30')) + 1)]))
        pn, pw = [], []
        for a, b in zip(ed[:-1], ed[1:]):
            pn.append((b - a)*(x + 1)/2 + a); pw.append((b - a)*w/2)
        self.pn, self.pw = np.concatenate(pn), np.concatenate(pw)


if NQ > 0:
    R3.NearGap = _NG


def LA_exact(v, Lm):
    """16 L_A with psi2 from the LP vector v (k(0) = 1)."""
    G = R5.NL
    p = lambda t: float(R3.basis_fn(np.array([t]))[:, 0] @ v[:R5.D + 1])
    a2 = 9.0/(4*Lm*Lm)
    val = quad(lambda t: p(t)/(t*t + a2), 0, 1, points=[3/(2*Lm), 10/Lm], limit=400, epsabs=1e-12, epsrel=1e-11)[0]
    eps = 1e-7
    A = 16*p(eps)/eps
    return 16*val, A


def weights():
    out = [('NF02', R3.nf_vector())]
    for a in sys.argv[1:]:
        if a == 'nf02':
            continue
        d = pickle.load(open(os.path.join(CACHE, a), 'rb'))
        out.append((a, d['v']))
    return out


def main():
    Dg = int(os.environ.get('T1_DG', '300'))
    bmax = float(os.environ.get('T1_BMAX', '2.0'))
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    mfls = [float(v) for v in os.environ.get('T1_MFLS', '1.02,1.05,1.1,1.2,1.3,1.5,1.75,2,2.5,3,4,5').split(',')]
    bhs = np.unique(np.round(np.geomspace(6, 3*Lm, int(os.environ.get('T1_NB', '110'))), 3))
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    sts = [s for s in R5.grid_states(Dg, Lm, bmax, mfls, bhs) if s[1] % 2 == 0]
    blocks, _ = R5.blocks_for(Dg, sts)
    pr(f"# t1_band Delta_g = {Dg}, Lam = {Lm:.3f}; universal-g model (r5), {len(blocks)} even-J states, beta <= {bmax}, "
       f"m/Lam in {mfls}; g validated for beta <= {BVALID} (b-hat <= {BVALID*Lm:.0f})  [{time.time() - t0:.0f}s]")
    for name, v in weights():
        F0 = R5.NL.to_functional(R5.NL.NFGrid(R5.D, R5.R), v); A, B, BA = F0.AB()
        L16, A2 = LA_exact(v, Lm)
        pr(f"\n## weights {name}: A = {A:.2f} (16psi2'(0) = {A2:.2f}), B_AdS = {BA:.2f}; 16 L_A({Lm:.1f}) = {L16:.1f} "
           f"(log form {A*np.log(Lm) + BA:.1f}); bound exact {240*np.sqrt(L16)/Lm**4:.4e}, log form {R3.bound(A, BA, Dg):.4e}")
        rows = []
        for b in blocks:
            Fm = R3.evalF(b, v); Nn = R3.evalN(b, v)
            lam = R3.lam(b, v)[0]
            X = abs(Fm[0, 1]); kc = np.sqrt(max(Fm[0, 0], 0)*max(Fm[1, 1], 0))/X if X > 0 else np.inf
            rows.append((b['mfl'], b['J'], b['bh'], b['beta'], Fm[0, 0]/Nn[0], Fm[1, 1]/Nn[1], X/np.sqrt(Nn.prod()), lam, kc,
                         X*Lm*Lm))
        rows = np.array(rows)
        pr("  per m/Lam: negative b-hat intervals (lambda_min < 0); min diagonal F_hh/N_h over the band; b-hat where |X|/sqrt(NN) = 10, 100")
        for mfl in mfls:
            r = rows[np.isclose(rows[:, 0], mfl)]
            if len(r) == 0:
                continue
            r = r[np.argsort(r[:, 2])]
            neg = r[:, 7] < 0
            ivs = []; i = 0
            while i < len(r):
                if neg[i]:
                    j = i
                    while j + 1 < len(r) and neg[j + 1]:
                        j += 1
                    ivs.append((r[i, 2], r[j, 2], i, j)); i = j + 1
                else:
                    i += 1
            lo_pos = r[~neg][:, 2]
            last_pos_below = max([bb for bb in lo_pos if not ivs or bb < ivs[0][0]], default=None)
            def first(th):
                k = np.where(r[:, 6] >= th)[0]
                return f"{r[k[0], 2]:.0f}" if len(k) else "-"
            pr(f"    m/Lam={mfl:5.2f}: b-hat range [{r[0, 2]:.1f}, {r[-1, 2]:.1f}] (beta <= {r[-1, 3]:.2f}); last positive below band "
               f"{last_pos_below if last_pos_below is None else round(last_pos_below, 1)}; negative intervals "
               + (", ".join(f"[{a:.1f}, {bb:.1f}]" for a, bb, _, _ in ivs) or "none")
               + f"; min F_hh/N in band {np.min(r[neg][:, 4:6]) if neg.any() else np.nan:+.3f}; |X|/sqrtNN >= 10 from {first(10)}, >= 100 from {first(100)}")
        # aggregate over m: band edge b_c = smallest negative b-hat; kappa_crit envelope in b-hat bins
        neg = rows[:, 7] < 0
        bc = rows[neg][:, 2].min() if neg.any() else np.nan
        pos_above = rows[(~neg) & (rows[:, 2] > bc)] if neg.any() else rows[:0]
        pr(f"  band edge (min over m of the first negative b-hat): b-hat_c = {bc:.1f}  (= {bc/Lm**(2/3):.3f} Lam^(2/3); beta_c = "
           f"{rows[neg][np.argmin(rows[neg][:, 2]), 3] if neg.any() else np.nan:.3f}); positive states above b-hat_c: {len(pos_above)}"
           + (f" (max b-hat {pos_above[:, 2].max():.1f}, at m/Lam {pos_above[np.argmax(pos_above[:, 2]), 0]})" if len(pos_above) else ""))
        edges = [0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300, 500, 1000, 2000, 5000]
        pr("  kappa_crit envelope (min over states in the b-hat bin; [M]; '*' = beyond the validated beta <= 0.11):")
        line = []
        for a, bb in zip(edges[:-1], edges[1:]):
            s = rows[(rows[:, 2] >= a) & (rows[:, 2] < bb)]
            if len(s) == 0:
                continue
            k = s[np.argmin(s[:, 8])]
            line.append(f"[{a},{bb}) {k[8]:.3g}{'*' if s[:, 3].max() > BVALID else ''} (lam {s[:, 7].min():+.2f}, |X|L^2 {np.median(s[:, 9]):.0f})")
        for i in range(0, len(line), 3):
            pr("    " + " | ".join(line[i:i + 3]))
        pickle.dump(dict(rows=rows, cols='mfl J bh beta F++/N F+-/N X/sqrtNN lam kappa XLam2'.split()),
                    open(os.path.join(CACHE, f't1_band_{Dg}_{name.replace(".pkl", "")}.pkl'), 'wb'))
    open(os.path.join(CACHE, f't1_band_{Dg}{os.environ.get("T1_TAG", "")}.txt'), 'w').write(buf.getvalue())
    pr(f"[{time.time() - t0:.0f}s]")


if __name__ == '__main__':
    main()
