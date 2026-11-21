"""Result Q, Stage 3, Task 2 (RESULT_Q_STAGE3.md §2-§3): heavy images of the Stage-0 functionals, lifted to AdS4 via Stage 1, on near-gap
and small-bhat heavy states at Delta_g = 100, 300, 1000.  MODEL [M] ("hybrid-exact skeleton"), per helicity channel h (++: even J; +-: J>=4;
the parity-odd even-J block couples like ++ at this order and is not listed separately):

  smeared k = 2, 3 (as Result D's hybrid):   B2 = B2_flat(h; exact Wigner)  + int psi2 (2 - u p^2) [P_{-1/2+i Lam p}(cosh b) - J0(bh p)] dp
                                             B3 = B3_flat(h)               + sigma u int psi3   [P_{-1/2+i Lam p}(cosh b) - J0(bh p)] dp
  forward k = 4, 5 at nu = 0 (flat graviton normal form + EXACT scalar-skeleton curvature, stage3_fwexact.py, Df = 3):
       B4|0  = 2 u^2 + u^2 (ex4(0) - 2)                         = u^2 ex4(0)
       dB4|0 = (c4_h - 2JJ) u^3 + u^2 [Lam^2 d ex4/d nu^2|_0 - u (5 - 2JJ)]     (scalar flat: ex4 = (2-y)(1-y)^-3 P_J(1-2y), y = nu^2/m^2)
       B5|0  = sigma u^3 ex5(0)/2                                (odd k: Gamma-ratio continuation, a definition [M])
  heavy action FF = B2 + B3 - a (dB4 - B5) - c B4  (c = 1/2),  N = |B2| + |B3| + |a (dB4 - B5)| + |c B4|  (Stage-0 normalisation).
Variants reported side by side: 'flat' (Stage-0 action at the AdS state's u), 'LR' (stage-12 leading-Regge forward pieces:
B4 = 2u^2 P0, dB4 = -2u^2 Lam^2 D + c4 u^3 P0, B5 = sigma u^3 P0, smeared parts as above), 'mult' (dB4 curvature multiplicative).
Lower-bound functional  cL B4 - (dB4 - B5) >= 0:  cL_min = max (dB4 - B5)/B4  =>  g4 Lam^6/(8 pi G) >= |g3|^2 Lam^8/(4 cL_min)  (b0 = 0).
AdS constants: A = 2 psi2'(0); 2 L_AdS = 2 int psi2/(p^2 + 9/(4 Lam^2)); B_eff(Lam) = 2 L_AdS - A log(2 Lam/3).
Usage: python3 stage3_neargap.py [functional name]  -> cache/stage3_neargap.txt (+ .pkl)
"""
import os, sys, io, json, glob, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
from scipy.special import j0
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import common_q as CQ
from stage12_forward_regge import mehler_nodes, P0_D

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()
C4 = {'++': 1.0, '+-': 41.0}; SIG = {'++': 1.0, '+-': -1.0}


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def load_fw(Dg, Df=3):
    db = {}
    for fn in glob.glob(os.path.join(CACHE, 'stage3_fw_%d_Df%g_w*.pkl' % (Dg, Df))):
        db.update(pickle.load(open(fn, 'rb')))
    return db


def LAdS(c2, Lam):
    a = 1.5 / Lam
    edges = np.unique(np.concatenate([[0.0], a * np.logspace(-3, np.log10(1 / a), 80), [1.0]])); edges = edges[edges <= 1.0]
    x, w = np.polynomial.legendre.leggauss(40); tot = 0.
    for lo, hi in zip(edges[:-1], edges[1:]):
        p = lo + (x + 1) * (hi - lo) / 2; ww = w * (hi - lo) / 2
        tot += np.sum(ww * np.polynomial.polynomial.polyval(p, c2) / (p * p + a * a))
    return tot


class Curv:
    """curvature parts int psi (..) [P_{-1/2+i Lam p}(cosh b) - J0(bh p)] dp on 400 Gauss nodes."""

    def __init__(self, c2, c3, Lam, n=400):
        x, w = np.polynomial.legendre.leggauss(n)
        self.p = 0.5 * (x + 1); self.w = 0.5 * w; self.Lam = Lam
        self.ps2 = np.polynomial.polynomial.polyval(self.p, c2); self.ps3 = np.polynomial.polynomial.polyval(self.p, c3)

    def __call__(self, beta, bh, u):
        t, W = mehler_nodes(beta, 1600)
        P = np.cos(np.outer(self.Lam * self.p, t)) @ W
        dP = P - j0(bh * self.p)
        c2 = np.sum(self.w * self.ps2 * (2 - u * self.p ** 2) * dP)
        c3 = u * np.sum(self.w * self.ps3 * dP)
        return c2, c3


def analyse(name, Dgs=(100, 300, 1000), verbose=True):
    F = json.load(open(os.path.join(CACHE, os.environ.get('S3_FUNC', 'functionals_q.json'))))[name]
    c2 = np.array(F['psi2'], float); c3 = np.array(F['psi3'], float) if len(F['psi3']) else np.zeros(1)
    a, c = F['a'], F['c']
    pr(f"\n=== functional '{name}': a = {a:.5f}, c = {c}, flat A = {F['A']:.4f}, B = {F['B']:.4f}, kappa = {F['kappa']:.5f}")
    out = {}
    for Dg in Dgs:
        Lam = np.sqrt((Dg - 2.0) * (Dg - 1.0))
        db = load_fw(Dg)
        if not db:
            pr(f"  Dg={Dg}: no forward records yet"); continue
        CV = Curv(c2, c3, Lam)
        flatcache = {}
        rows = []
        for (Delta, J), rec in sorted(db.items(), key=lambda kv: (kv[0][1], kv[0][0])):
            mfl = rec['mfl']; m = mfl * Lam; u = 1.0 / mfl ** 2
            al = Delta - J - 2.0; m2 = al * (al + 2 * J + 1)
            assert abs(m2 / (m * m) - 1) < 1e-6
            bh = (2 * J + 1) * np.sqrt(u)
            beta = np.arccosh(1 + (2 * J + 1) ** 2 / (2 * m2))
            P0, D = P0_D(beta, 2000)
            cb2, cb3 = CV(beta, bh, u)
            JJ = J * (J + 1.0)
            ex4, d4, ex5 = rec[4]['v0'], rec[4]['d0'], rec[5]['v0']
            for h in ('++', '+-'):
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                key = (h, round(u, 12))
                if key not in flatcache:
                    Js, B2f, B3f = CQ.parts(h, u, 24, [c2], [c3])
                    flatcache[key] = dict(zip(Js, zip(B2f[:, 0], B3f[:, 0])))
                b2f, b3f = flatcache[key][J]
                s = SIG[h]
                B2 = b2f + cb2; B3 = b3f + s * cb3
                fw = {}
                fw['flat'] = (2 * u * u, (C4[h] - 2 * JJ) * u ** 3, s * u ** 3)
                fw['exact'] = (u * u * ex4, (C4[h] - 2 * JJ) * u ** 3 + u * u * (Lam ** 2 * d4 - u * (5 - 2 * JJ)), s * u ** 3 * ex5 / 2)
                rd = Lam ** 2 * d4 / (u * (5 - 2 * JJ))
                fw['mult'] = (u * u * ex4, (C4[h] - 2 * JJ) * u ** 3 * rd, s * u ** 3 * ex5 / 2)
                fw['LR'] = (2 * u * u * P0, -2 * u * u * Lam ** 2 * D + C4[h] * u ** 3 * P0, s * u ** 3 * P0)
                res = dict(Dg=Dg, J=J, mfl=mfl, u=u, bh=bh, beta=beta, h=h, ex4=ex4, d4L=Lam ** 2 * d4, ex5=ex5, P0=P0, DL=Lam ** 2 * D)
                for kv, (b4, db4, b5) in fw.items():
                    if kv == 'flat':
                        bb2, bb3 = b2f, b3f
                    else:
                        bb2, bb3 = B2, B3
                    FF = bb2 + bb3 - a * (db4 - b5) - c * b4
                    N = abs(bb2) + abs(bb3) + abs(a * (db4 - b5)) + abs(c * b4)
                    res[kv] = FF / N
                    res[kv + '_X'] = FF + c * b4; res[kv + '_B4'] = b4
                    res['cL_' + kv] = (db4 - b5) / b4
                res['B2'], res['B3'], res['FWx'] = B2, B3, -a * (fw['exact'][1] - fw['exact'][2]) - c * fw['exact'][0]
                rows.append(res)
        out[Dg] = rows
        L2 = 2 * LAdS(c2, Lam); A = 2 * c2[1]
        pr(f"  Dg = {Dg} (Lam = {Lam:.3f}): {len(rows)} (state, channel) cases, J <= 24, m/Lam in [gap, 4]")
        for kv in ('flat', 'exact', 'mult', 'LR'):
            w = min(rows, key=lambda r: r[kv])
            pr(f"    min FF/N [{kv:5s}] = {w[kv]:+.4f}  at {w['h']} J={w['J']} m/Lam={w['mfl']:.4f} bhat={w['bh']:.2f}")
        for lab, sel in (('bhat < 10', lambda r: r['bh'] < 10), ('10 <= bhat < 20', lambda r: 10 <= r['bh'] < 20),
                         ('bhat >= 20', lambda r: r['bh'] >= 20), ('m <= 1.3 Lam', lambda r: r['mfl'] <= 1.3),
                         ('m > 1.3 Lam', lambda r: r['mfl'] > 1.3)):
            v = [r for r in rows if sel(r)]
            if v:
                w = min(v, key=lambda r: r['exact']); wl = min(v, key=lambda r: r['LR'])
                pr(f"      {lab:16s}: exact min {w['exact']:+.4f} ({w['h']} J={w['J']} m/Lam={w['mfl']:.3f}); LR min {wl['LR']:+.4f}"
                   f" ({wl['h']} J={wl['J']} m/Lam={wl['mfl']:.3f})")
        # forward-piece diagnostics: exact skeleton vs leading Regge
        dev4 = max(abs(r['ex4'] / 2 / r['P0'] - 1) for r in rows)
        devd = max(abs((r['d4L'] - (-2 * r['DL'] + 5 * r['u'] * r['P0'])) / (abs(r['d4L']) + 1e-30)) for r in rows)
        pr(f"    exact skeleton vs leading Regge: max|ex4(0)/(2 P0) - 1| = {dev4:.2e};  max rel |Lam^2 d ex4 - (-2 Lam^2 D + 5u P0)| = {devd:.2e}")
        cLf = max(r['cL_flat'] for r in rows); cLx = max(r['cL_exact'] for r in rows); cLr = max(r['cL_LR'] for r in rows)
        wcl = max(rows, key=lambda r: r['cL_exact'])
        pr(f"    lower bound: cL_min = max (dB4-B5)/B4: flat-at-AdS-states {cLf:.4f}, exact {cLx:.4f} ({wcl['h']} J={wcl['J']} m/Lam={wcl['mfl']:.4f}),"
           f" LR {cLr:.4f}  ->  slope 1/(4 cL) = {1 / (4 * cLf):.4f} / {1 / (4 * cLx):.4f} / {1 / (4 * cLr):.4f}  (flat: cL = 1, slope 1/4)")
        pr(f"    AdS constants (functional unchanged): A = {A:.4f}; 2 L_AdS = {L2:.4f}; B_eff = 2L_AdS - A log(2Lam/3) = {L2 - A * np.log(2 * Lam / 3):.4f}"
           f" (flat B = {F['B']:.4f}); kappa = {F['kappa']:.5f} (+ [I] b0)")
        for kv in ('exact', 'LR'):
            # largest B4|0 coefficient c' keeping every computed image >= 0:  FF(c') = X - c' B4,  X = FF + c B4
            cm = min((r[kv + '_X'] / r[kv + '_B4']) for r in rows)
            w = min(rows, key=lambda r: r[kv + '_X'] / r[kv + '_B4'])
            pr(f"    [{kv}] largest admissible c (all computed states >= 0): c_max = {cm:.4f} ({w['h']} J={w['J']} m/Lam={w['mfl']:.4f});"
               f" upper bound x 1/(2 min(c_max, 1/2)) = {1 / (2 * min(cm, 0.5)):.4f}")
        out[Dg] = dict(rows=rows, L2=L2, A=A)
    return out


if __name__ == '__main__':
    names = sys.argv[1:] or ['Q0-full-robust', 'Q0-min-robust (no B3)', 'Q0-r5-psi3zero-robust (AdS-style weights)']
    res = {}
    for nm in names:
        res[nm] = analyse(nm)
    pickle.dump(res, open(os.path.join(CACHE, 'stage3_neargap.pkl'), 'wb'))
    open(os.path.join(CACHE, 'stage3_neargap.txt'), 'w').write(buf.getvalue())
