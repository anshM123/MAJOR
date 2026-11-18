"""
Opposite-helicity even-kernel moments: diagnosis, kernel normalisation, and the isotropic-mode-free flip-free kernel K* (SPINNING_TTTT_EXT.md sect. 6).
(b) flip response of the one-pair tables per circle mode (M^{q,-q}, laplace.M_table(2)): (0,0) : (+-1,-+1) : (+-2,-+2) = 36 : 16 : 1, independent of nu, a.
    TENSOR-kernel circle weights = pointwise circle weights x n_|q|, n = (1/36, 1/16, 1); fixed by the two Lemma-L flip ratios and then verified on
    18 Lemma-L preserved ratios (REGGE_MOMENT_TTTT.md sect. 4.2). Flip-free iff 36 w_0 + 16 (w_1 + w_-1) + (w_2 + w_-2) = 0 (tensor weights).
    D1^2 = {(0,0): 1/72, (+-2): 1/4}; D1D2 = {(+-1): -1/32}; D2^2 = {(0,0): 1/36};  S_ff = D1^2 + D1D2;  K* = D1^2 + D1D2/2 - D2^2/2 = (D1 + D2)(D1 - D2/2).
    E0, FF0 = the naive pointwise weights used in t16_kernels_ext.py (NOT the project's kernels; kept for the diagnosis).
(a) per-mode convergence of the main modes (k = 4..7, nu = 0, 1, 3, beta = 1, m = 39, 70, 101);
(c) kernel ratios for E0, FF0, D1sq, S_ff, K*, both helicities, with the prediction's cancellation factor;
(d1) k = 2; (d2) mixed wave W[X_s,X_o]; (d3) beta = 0.4, m = 72.
Output: out_t18_kstar.txt
"""
import os, sys, glob
import numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from t13_analysis_l2 import Modes, design, Tpred, CLASSES
from predict4 import Mval, weyl
import t5_analysis_l1 as T5
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
# Circle-mode weights of the TENSOR kernels (index contraction): the pointwise circle weights times n_q = (1/36, 1/16, 1) for |q| = (0, 1, 2).
# n_q is fixed by Lemma L (REGGE_MOMENT_TTTT.md sect. 4.2): with it, the per-mode responses of the one-pair tables reproduce the Lemma-L flip
# responses D1^2 : D1D2 : D2^2 = 1 : -1 : 1 and the preserved ratios D1D2/D1^2 = -0.2709, D2^2/D1^2 = 0.1295 (a = 7, nu = 0) exactly (check in (b)).
# 'E0' and 'FF0' are the NAIVE pointwise circle weights used in t16 (not the project's kernels).
KS = {'E0': {(0, 0): 0.5, (2, -2): 0.25, (-2, 2): 0.25},
      'FF0': {(0, 0): 0.5, (2, -2): 0.25, (-2, 2): 0.25, (1, -1): -0.5, (-1, 1): -0.5},
      'D1sq': {(0, 0): 1/72, (2, -2): 0.25, (-2, 2): 0.25},
      'Sff': {(0, 0): 1/72, (2, -2): 0.25, (-2, 2): 0.25, (1, -1): -1/32, (-1, 1): -1/32},
      'K*': {(2, -2): 0.25, (-2, 2): 0.25, (1, -1): -1/64, (-1, 1): -1/64}}
FILES = ['cache/modes_l2_D45.861045_J20_d.npz', 'cache/modes_l2_D80.484300_J36_d.npz', 'cache/modes_l2_D115.107554_J52_d.npz']

def table_check():
    pr("# (b) per-mode flip response of the one-pair table, ratios to the (2,-2) mode, and kernel flip/preserved responses")
    worst = 0
    for nu in (0.0, 0.5, 1.0, 3.0, 6.0, 10.0):
        for a in (6, 7, 8, 9, 10, 11, 12):
            f = {q: complex(Mval(2, q[0], q[1], 1, 1, nu, a)[0]) for q in ((0, 0), (1, -1), (-1, 1), (2, -2), (-2, 2))}
            r0, r1 = (f[(0, 0)]/f[(2, -2)]).real, (f[(1, -1)]/f[(2, -2)]).real
            worst = max(worst, abs(r0 - 36), abs(r1 - 16), abs((f[(-1, 1)]/f[(1, -1)]) - 1), abs(f[(-2, 2)]/f[(2, -2)] - 1))
            if a in (7, 11) and nu in (0.0, 3.0, 10.0):
                resp = {}
                for name, K in KS.items():
                    fl = sum(w*complex(Mval(2, q[0], q[1], 1, 1, nu, a)[0]) for q, w in K.items())
                    pp = sum(w*complex(Mval(2, q[0], q[1], 1, -1, nu, a)[0]) for q, w in K.items())
                    resp[name] = (pp.real, abs(fl))
                pr(f"  nu={nu} a={a}: f00/f22 = {r0:.6f}, f11/f22 = {r1:.6f};  " + "  ".join(f"{n}: pres {p:+.4f} |flip| {fl:.1e}" for n, (p, fl) in resp.items()))
    pr(f"  max deviation from 36 : 16 : 1 (and +- symmetry) over nu = 0..10, a = 6..12: {worst:.1e}")
    # Lemma-L consistency of the tensor weights: D2^2 = (0,0) only with weight 1/36, D1D2 = (+-1) with -1/32
    D2sq = {(0, 0): 1/36}; D1D2 = {(1, -1): -1/32, (-1, 1): -1/32}
    for nu in (0.0, 1.1, 3.0):
        for a in (7, 9, 11):
            r = {}
            for name, K in (('D1sq', KS['D1sq']), ('D1D2', D1D2), ('D2sq', D2sq)):
                r[name] = (sum(w*complex(Mval(2, q[0], q[1], 1, -1, nu, a)[0]) for q, w in K.items()).real,
                           sum(w*complex(Mval(2, q[0], q[1], 1, 1, nu, a)[0]) for q, w in K.items()))
            nu2 = nu*nu; d5 = nu2 + (a - 2.5)**2; d9 = nu2 + (a - 4.5)**2
            N1 = 64*a**4 - 768*a**3 + 64*a*a*nu2 + 3280*a*a - 320*a*nu2 - 5840*a + 16*nu2**2 + 488*nu2 + 3705
            N2 = 64*a**3 - 32*a*a*nu2 - 584*a*a + 224*a*nu2 + 1720*a - 16*nu2**2 - 424*nu2 - 1641
            N3 = 128*a*a - 128*a*nu2 - 800*a + 16*nu2**2 + 360*nu2 + 1113
            L = {'D1sq': N1, 'D1D2': -N2, 'D2sq': N3}
            pr(f"  Lemma-L check nu={nu} a={a}: pres ratios (tables/LemmaL) D1D2/D1sq {r['D1D2'][0]/r['D1sq'][0]:+.6f}/{L['D1D2']/L['D1sq']:+.6f}, "
               f"D2sq/D1sq {r['D2sq'][0]/r['D1sq'][0]:+.6f}/{L['D2sq']/L['D1sq']:+.6f}; flip ratios D1D2/D1sq {(r['D1D2'][1]/r['D1sq'][1]).real:+.6f} (L: -1), D2sq/D1sq {(r['D2sq'][1]/r['D1sq'][1]).real:+.6f} (L: +1)")

def kernel_ratio(mo, k, nu, wave, K, A, Tp):
    terms = []; pred = 0j
    for (q1, q3), w13 in K.items():
        for (q2, q4), w24 in K.items():
            q = (q1, q3, q2, q4); v, sh = mo.moment(k, nu, wave, q)
            terms.append((w13*w24, v, sh)); pred += w13*w24*sum(A[q][c]*Tp[c] for c in CLASSES)
    shm = max(t[2] for t in terms)
    meas = sum(w*v*np.exp(sh - shm) for w, v, sh in terms)*(-1)**mo.J
    canc = sum(abs(w*sum(A[(q1,q3,q2,q4)][c]*Tp[c] for c in CLASSES)) for (q1,q3),w1 in K.items() for (q2,q4),w2 in K.items() for w in [w1*w2])/abs(pred)
    return complex(meas*np.exp(shm - np.log(abs(pred)))/(pred/abs(pred))).real, canc

if __name__ == "__main__":
    pr("# Opposite-helicity even kernels: diagnosis and the kernel K* (fourpt/tttt/t18_kstar.py)")
    table_check()
    mos = [Modes(f) for f in FILES]
    ms = [float(weyl(mo.Delta, mo.J)[2]) for mo in mos]
    pr("\n# (a) per-mode ratios meas/pred, m^2(R-1) at m = " + ", ".join(f"{m:.1f}" for m in ms) + " (beta = 1)")
    for k in (4, 5, 6, 7):
        for nu in (0.0, 1.0, 3.0):
            As = [design(mo.Delta, mo.J, k, nu) for mo in mos]
            for wave, qs in ((('same', 'same'), [(2, -2, -2, 2), (1, -1, -1, 1), (0, 0, 0, 0), (0, 0, 2, -2)]),
                             (('opp', 'opp'), [(2, -2, 2, -2), (1, -1, 1, -1), (0, 0, 2, -2), (0, 0, 0, 0)])):
                cells = []
                for q in qs:
                    Rs = []
                    for mo, A in zip(mos, As):
                        Tp = Tpred(wave[0], wave[1], mo.Delta, mo.J, 8.0)
                        p = sum(A[q][c]*Tp[c] for c in CLASSES); v, sh = mo.moment(k, nu, wave, q)
                        Rs.append(complex((-1)**mo.J*v*np.exp(sh - np.log(abs(p)))/(p/abs(p))).real)
                    cells.append(f"{q}: " + "/".join(f"{m*m*(R-1):+.0f}" for m, R in zip(ms, Rs)))
                pr(f"  k={k} nu={nu} {wave[0]}: " + " | ".join(cells))
    pr("\n# (c) kernel ratios (R at the three masses; cancellation sum|terms|/|sum| of the prediction; m^2(R-1); c1)")
    for k in (4, 5, 6, 7):
        for nu in (0.0, 1.0, 3.0):
            As = [design(mo.Delta, mo.J, k, nu) for mo in mos]
            for wave in (('same', 'same'), ('opp', 'opp')):
                cells = []
                for name, K in KS.items():
                    res = [kernel_ratio(mo, k, nu, wave, K, A, Tpred(wave[0], wave[1], mo.Delta, mo.J, 8.0)) for mo, A in zip(mos, As)]
                    Rs = [r[0] for r in res]
                    f = T5.fits(ms, Rs)
                    cells.append(f"{name}: " + "/".join(f"{R:.4f}" for R in Rs) + f" (canc {res[-1][1]:.1f}; m2d " + "/".join(f"{m*m*(R-1):+.0f}" for m, R in zip(ms, Rs)) + f"; c1 {f['A'][0]:+.2f})")
                pr(f"  k={k} nu={nu} {wave[0]}: " + " | ".join(cells))
    pr("\n# (d1) k = 2 (C+_2 moment, a = 7), beta = 1, kernels with tensor weights (R at the three masses; m^2(R-1); c1)")
    TK = {n: KS[n] for n in ('D1sq', 'Sff', 'K*')}
    for nu in (0.0, 1.0, 3.0):
        As = [design(mo.Delta, mo.J, 2, nu) for mo in mos]
        for wave in (('same', 'same'), ('opp', 'opp')):
            cells = []
            for name, K in TK.items():
                res = [kernel_ratio(mo, 2, nu, wave, K, A, Tpred(wave[0], wave[1], mo.Delta, mo.J, 8.0)) for mo, A in zip(mos, As)]
                Rs = [r[0] for r in res]; f = T5.fits(ms, Rs)
                cells.append(f"{name}: " + "/".join(f"{R:.4f}" for R in Rs) + f" (m2d " + "/".join(f"{m*m*(R-1):+.0f}" for m, R in zip(ms, Rs)) + f"; c1 {f['A'][0]:+.2f})")
            pr(f"  k=2 nu={nu} {wave[0]}: " + " | ".join(cells))
    pr("\n# (d2) mixed wave |K(W[s,o])|/sqrt|K(W[s,s]) K(W[o,o])|, measured/predicted, at the three masses (beta = 1)")
    def kmom(mo, k, nu, wave, K, A, Tp):
        terms = []; pred = 0j
        for (q1, q3), w13 in K.items():
            for (q2, q4), w24 in K.items():
                q = (q1, q3, q2, q4); v, sh = mo.moment(k, nu, wave, q)
                terms.append((w13*w24, v, sh)); pred += w13*w24*sum(A[q][c]*Tp[c] for c in CLASSES)
        shm = max(t[2] for t in terms)
        return np.log(abs(sum(w*v*np.exp(sh - shm) for w, v, sh in terms)) + 1e-300) + shm, np.log(abs(pred) + 1e-300)
    for k, nu in ((4, 0.0), (6, 0.0), (6, 3.0), (7, 1.0)):
        As = [design(mo.Delta, mo.J, k, nu) for mo in mos]
        cells = []
        for name, K in TK.items():
            vals = []
            for mo, A in zip(mos, As):
                r = {w: kmom(mo, k, nu, w, K, A, Tpred(w[0], w[1], mo.Delta, mo.J, 8.0)) for w in (('same', 'same'), ('opp', 'opp'), ('same', 'opp'))}
                lm = r[('same', 'opp')][0] - 0.5*(r[('same', 'same')][0] + r[('opp', 'opp')][0]); lp = r[('same', 'opp')][1] - 0.5*(r[('same', 'same')][1] + r[('opp', 'opp')][1])
                vals.append(f"{np.exp(lm):.1e}/{np.exp(lp):.1e}")
            cells.append(f"{name}: " + ", ".join(vals))
        pr(f"  k={k} nu={nu}: " + " | ".join(cells))
    pr("\n# (d3) beta = 0.4, m = 72.0 (J = 14), one mass: R for D1sq, Sff, K* and the leading mode")
    mo4 = Modes(os.path.join(HERE, 'cache', 'modes_l2_D74.964099_J14_d.npz'))
    for k in (4, 5, 6, 7):
        for nu in (0.0, 1.0, 3.0):
            A = design(mo4.Delta, mo4.J, k, nu); cells = []
            for wave in (('same', 'same'), ('opp', 'opp')):
                Tp = Tpred(wave[0], wave[1], mo4.Delta, mo4.J, 8.0)
                cells.append(f"{wave[0]}: " + " ".join(f"{n} {kernel_ratio(mo4, k, nu, wave, K, A, Tp)[0]:.4f}" for n, K in TK.items()))
            pr(f"  k={k} nu={nu}: " + " | ".join(cells))
    open(os.path.join(HERE, 'out_t18_kstar.txt'), 'w').write("\n".join(out) + "\n")
