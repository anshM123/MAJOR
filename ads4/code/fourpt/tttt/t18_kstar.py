"""
Diagnosis of the opposite-helicity even-kernel moments and the isotropic-mode-free, flip-free even kernel K* (SPINNING_TTTT_EXT.md sect. 6).
(a) per-mode convergence of the main modes of W[X_s,X_s], W[X_o,X_o] (k = 4..7, nu = 0, 1, 3, beta = 1, m = 39, 70, 101);
(b) the flip response of the one-pair tables per circle mode (M^{q,-q}_{h,h}, laplace.M_table(2)): ratios (0,0) : (+-1,-+1) : (+-2,-+2) = 36 : 16 : 1,
    independent of nu and a (checked on a grid) -> a kernel with circle weights w_q is flip-free iff 36 w_0 + 16 (w_1 + w_-1) + (w_2 + w_-2) = 0;
(c) K* = {(+-2,-+2): 1/4, (+-1,-+1): -1/64, (0,0): 0} per pair (flip-free by (b), no isotropic mode), K* x K* end to end, both helicities;
    its preserved response vs E = D1^2 (circle weights) and the project's S_ff (circle weights ((0,0) 1/2, (+-1) -1/2, (+-2) 1/4)).
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
KS = {'E': {(0, 0): 0.5, (2, -2): 0.25, (-2, 2): 0.25},
      'FF': {(0, 0): 0.5, (2, -2): 0.25, (-2, 2): 0.25, (1, -1): -0.5, (-1, 1): -0.5},
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
    open(os.path.join(HERE, 'out_t18_kstar.txt'), 'w').write("\n".join(out) + "\n")
