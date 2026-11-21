"""Result Q, Stage 6 (RESULT_Q_STAGE6.md §1): full-physical-region heavy positivity of the delta = 0.05 / 0.2 functionals at Dg = 100, 300, 1000.

Region A (J <= JCAP = 3000, every integer J with the right parity, twist gap imposed):  HYB model (stage6_model), on
   m/Lam in X = {1 + 10^[-4, 0] (100 pts), geomspace(2, 1000, 360)}  (u from 1 down to 1e-6), plus every gap state u = u_gap(J), J <= 300.
Region B (J > JCAP): LR model (Stage 5's beta_quantities), beta grid as Stage 5 (382 values up to beta = 30; coarse in bhat at small beta --
   superseded by the dense region A'/B' of stage6_scan_ext.py), u in [1e-8, 1] (240 values),
   J continuous, twist gap imposed ('phys'); also without the gap (Stage 5's convention, a superset) for comparison.
EXA: HYB with the exact CHMRS skeleton forward moments on every state where they exist (Stage 3 grid + stage6_fwexact states).
Cross-validation: on all region-A states with 10 <= bhat < 30, HYB vs LR (same states, same normalisation N); EXA vs HYB on the exact states.
Scaling limit u -> 0 at fixed bhat: FF/N -> sign of S_AdS(bhat) = 2 int psi2 P_{-1/2+i Lam p}(1 + bhat^2/(2 Lam^2)) dp: its minimum is listed.
Normalisation: Stage-0 N (the LP's); Stage-5 Nrm also given for the LR rows.
Usage: python3 stage6_scan.py [delta ...]  -> cache/stage6_scan_<delta>.txt, .pkl
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import common_q as CQ
import stage6_model as M
import stage5_positivity as S5P

CACHE = M.CACHE
JCAP = 3000
XS = np.unique(np.concatenate([1 + np.logspace(-4, 0, 100), np.geomspace(2, 1000, 360)]))
BMAX = {100: 800., 300: 1350., 1000: 2300.}
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


REGIONS = [('R1 J<=24, m<=4Lam (Stage-3 box)', lambda J, x, b: (J <= 24) & (x <= 4.0001)),
           ('R2 J>24, m<=4Lam, bhat<20', lambda J, x, b: (J > 24) & (x <= 4.0001) & (b < 20)),
           ('R3 m>4Lam, bhat<20', lambda J, x, b: (x > 4.0001) & (b < 20)),
           ('R4 overlap 10<=bhat<30', lambda J, x, b: (b >= 10) & (b < 30)),
           ('R5 bhat>=20 (J<=3000)', lambda J, x, b: b >= 20),
           ('all region A', lambda J, x, b: b >= 0)]


class Best:
    def __init__(self):
        self.d = {}

    def add(self, key, vals, info):
        if len(vals) == 0:
            return
        i = int(np.argmin(vals))
        if key not in self.d or vals[i] < self.d[key][0]:
            self.d[key] = (float(vals[i]), {k: (v[i] if np.ndim(v) else v) for k, v in info.items()})


def fmt(t):
    v, I = t
    return f"{v:+.4f} ({I['h']} J={I['J']:.0f} m/Lam={I['x']:.3f} bhat={I['b']:.1f})"


def region_A(F, name, Dg, tab, flat, Lam, best, xval):
    for h in ('++', '+-'):
        for (u, Js, B2f, B3f) in flat[h]:
            if len(Js) == 0:
                continue
            sel = Js <= M.jmax_gap(u, Dg)
            if not np.any(sel):
                continue
            J = Js[sel].astype(float); x = 1/np.sqrt(u)
            FF, N, Nrm, P = M.hyb(F, h, u, J, B2f[sel], B3f[sel], tab, Lam)
            FFl, Nl, Nrml = M.lr(F, h, u, J, tab, Lam)
            b = P['bh']; r = FF/N; rl = FFl/Nl
            info = dict(h=h, J=J, x=np.full(len(J), x), b=b)
            for lab, f in REGIONS:
                s = f(J, x, b)
                if np.any(s):
                    best.add(('HYB', lab), r[s], {k: (v[s] if np.ndim(v) else v) for k, v in info.items()})
                    best.add(('LR', lab), rl[s], {k: (v[s] if np.ndim(v) else v) for k, v in info.items()})
                    best.add(('LRnrm', lab), (FFl/Nrml)[s], {k: (v[s] if np.ndim(v) else v) for k, v in info.items()})
            s = (b >= 10) & (b < 30)
            if np.any(s):
                xval.append(np.c_[b[s], r[s], rl[s], np.full(s.sum(), x), J[s], np.full(s.sum(), 1 if h == '++' else -1)])


def gap_rows(F, name, Dg, tab, Lam, best):
    G = M.gap_states(name, Dg, 300)
    rows = []
    for h in ('++', '+-'):
        for (u, J, b2, b3) in G[h]:
            FF, N, Nrm, P = M.hyb(F, h, u, float(J), b2, b3, tab, Lam)
            rows.append((FF/N, h, J, 1/np.sqrt(u), P['bh']))
            info = dict(h=h, J=float(J), x=1/np.sqrt(u), b=float(P['bh']))
            best.add(('HYB', 'gap states J<=300'), np.array([FF/N]), info)
            for lab, f in REGIONS:
                if f(J, 1/np.sqrt(u), P['bh']):
                    best.add(('HYB', lab), np.array([FF/N]), info)
    return rows


def exact_rows(F, name, Dg, tab, Lam, best):
    db = M.load_exact(Dg)
    F_, c2, c3 = M.load(name)
    out = []
    for (Delta, J), rec in db.items():
        x = rec['mfl']; u = 1/x**2
        for h in ('++', '+-'):
            if (h == '++' and J % 2) or (h == '+-' and J < 4):
                continue
            Js, B2, B3 = CQ.parts(h, u, max(J, 4), [c2], [c3]); k = list(Js).index(J)
            fw = (rec[4]['v0'], Lam**2*rec[4]['d0'], rec[5]['v0'])
            FFe, Ne, _, P = M.hyb(F, h, u, float(J), B2[k, 0], B3[k, 0], tab, Lam, fw=fw)
            FFh, Nh, _, _ = M.hyb(F, h, u, float(J), B2[k, 0], B3[k, 0], tab, Lam)
            out.append((FFe/Ne, FFh/Nh, h, J, x, float(P['bh']), 'stage6' if (J > 24 or x > 4.0001) else 'stage3'))
            best.add(('EXA', 'exact states'), np.array([FFe/Ne]), dict(h=h, J=float(J), x=x, b=float(P['bh'])))
    return out


def region_B(F, name, Lam, Dg, best, gap=True, tag='LRB'):
    a, c = F['a'], F['c']
    F_, c2, c3 = M.load(name)
    T = [S5P.CosT(c2), S5P.CosT(np.convolve(c2, [0, 0, 1])), S5P.CosT(c3)]
    betas = np.unique(np.concatenate([np.geomspace(5/Lam, 1.0, 70), np.linspace(1.0, 12.0, 221), np.linspace(12.0, 30.0, 91)]))
    us = np.concatenate([np.geomspace(1e-8, 0.05, 140), np.linspace(0.0505, 1.0, 100)])
    for beta in betas:
        I20, I22, I3, P0, D, H5 = S5P.beta_quantities(beta, Lam, T)
        bh = 2*Lam*np.sinh(beta/2)
        for h in ('++', '+-'):
            s = M.SIG[h]; c4 = M.C4[h]
            J = 0.5*(bh/np.sqrt(us) - 1)
            ok = J > JCAP if gap else J >= (0 if h == '++' else 4)
            if gap:
                ok &= us <= M.u_gap(J, Dg)
            else:
                ok &= bh >= 20
            if not np.any(ok):
                continue
            u = us[ok]; JJ = J[ok]
            B2 = 2*I20 - u*I22; B3 = s*u*I3
            B4 = 2*u*u*P0; dB4 = -2*u*u*Lam**2*D + c4*u**3*P0; B5 = s*u**3*P0
            FW = -a*(dB4 - B5); FF = B2 + B3 + FW - c*B4
            N = np.abs(B2) + np.abs(B3) + np.abs(FW) + c*np.abs(B4)
            Nrm = np.abs(B2) + np.abs(B3) + a*(2*u*u*Lam**2*D + 42*u**3*P0 + u**3*P0) + 2*c*u*u*P0
            info = dict(h=h, J=JJ, x=1/np.sqrt(u), b=np.full(len(u), bh))
            best.add((tag, 'N'), FF/N, info); best.add((tag, 'Nrm'), FF/Nrm, info)
            if beta >= 2:
                best.add((tag, 'N beta>=2'), FF/N, info)


def scaling(tab):
    b = tab.b; S = 2*tab.v[0]; Sf = 2*tab.v[5]
    i = int(np.argmin(S)); k = int(np.argmin(S/np.maximum(np.abs(Sf), 1e-300)))
    return S[i], b[i], (S/np.maximum(np.abs(Sf), 1e-300))[k], b[k]


if __name__ == '__main__':
    deltas = [float(a) for a in sys.argv[1:]] or [0.05, 0.2]
    for delta in deltas:
        name = M.FUNCS[delta]
        F, c2, c3 = M.load(name)
        t0 = time.time()
        flat = M.flat_table(name, XS, JCAP)
        pr(f"\n######## delta = {delta}: '{name}' (A = {F['A']:.4f}, a = {F['a']:.5f}); flat table {len(XS)} u-values x J <= {JCAP} [{time.time() - t0:.0f}s]")
        store = {}
        for Dg in (100, 300, 1000):
            Lam = M.lam(Dg); t1 = time.time()
            tab = M.Tables(name, Lam, BMAX[Dg])
            best = Best()
            xval = []
            region_A(F, name, Dg, tab, flat, Lam, best, xval)
            gr = gap_rows(F, name, Dg, tab, Lam, best)
            ex = exact_rows(F, name, Dg, tab, Lam, best)
            region_B(F, name, Lam, Dg, best, gap=True, tag='LRB')
            region_B(F, name, Lam, Dg, best, gap=False, tag='LR-S5')
            pr(f"\n=== Dg = {Dg} (Lam = {Lam:.2f}) [{time.time() - t1:.0f}s]")
            for (mod, lab), t in sorted(best.d.items(), key=lambda kv: (kv[0][1], kv[0][0])):
                pr(f"   {mod:6s} {lab:34s}: {fmt(t)}")
            X = np.concatenate(xval) if xval else np.zeros((0, 6))
            for lo, hi in ((10, 20), (20, 30)):
                s = (X[:, 0] >= lo) & (X[:, 0] < hi)
                d = X[s, 1] - X[s, 2]
                i = int(np.argmax(np.abs(d)))
                pr(f"   overlap {lo}<=bhat<{hi}: {s.sum()} states; HYB - LR in FF/N: max |d| = {np.abs(d).max():.4f} (at bhat {X[s, 0][i]:.1f},"
                   f" m/Lam {X[s, 3][i]:.2f}, J {X[s, 4][i]:.0f}); mean |d| {np.abs(d).mean():.4f}; min HYB {X[s, 1].min():+.4f}, min LR {X[s, 2].min():+.4f}")
                for xlo, xhi in ((1, 4), (4, 30), (30, 1001)):
                    ss = s & (X[:, 3] >= xlo) & (X[:, 3] < xhi)
                    if ss.sum():
                        dd = X[ss, 1] - X[ss, 2]
                        pr(f"        m/Lam in [{xlo},{xhi}): max|d| {np.abs(dd).max():.4f}, mean d {dd.mean():+.4f}; min HYB {X[ss, 1].min():+.4f}, min LR {X[ss, 2].min():+.4f}")
            if ex:
                e = np.array([(r[0], r[1]) for r in ex])
                for tag in ('stage3', 'stage6'):
                    s = np.array([r[6] == tag for r in ex])
                    if s.any():
                        d = e[s, 0] - e[s, 1]; i = int(np.argmax(np.abs(d)))
                        rr = [r for r, q in zip(ex, s) if q][i]
                        pr(f"   EXA - HYB on {tag} exact states ({s.sum()} cases): max |d| = {np.abs(d).max():.5f} ({rr[2]} J={rr[3]} m/Lam={rr[4]:.2f});"
                           f" min EXA {e[s, 0].min():+.4f}, min HYB {e[s, 1].min():+.4f}")
            S, bS, rS, brS = scaling(tab)
            pr(f"   scaling limit u->0: min_bhat S_AdS = {S:+.3e} at bhat {bS:.1f}; min S_AdS/|S_flat| = {rS:+.4f} at bhat {brS:.1f}  (bhat <= {tab.b[-1]:.0f})")
            store[Dg] = dict(best=best.d, xval=X, exact=ex, gap=gr, scaling=(S, bS, rS, brS))
        pickle.dump(store, open(os.path.join(CACHE, f'stage6_scan_{delta}.pkl'), 'wb'))
        open(os.path.join(CACHE, f'stage6_scan_{delta}.txt'), 'w').write(buf.getvalue())
