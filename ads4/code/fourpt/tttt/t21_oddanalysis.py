"""
Parity-odd exchange end to end (PARITY_ODD_HEAVY.md sect. 4): kernel moments of W[O,O], W[E,E] (and the parity test on W[E,O]) at odd J,
beta = 1, against the factorised heavy formula for the opposite-helicity structure (the identical-T structure at odd J is X^O_opp).

O is the real-coefficient structure S with X^O_opp = i S (tmoments_odd); the physical (CHL-normalised) wave is W[X^O, X^O] = -W[S, S].
Prediction: under flip-free kernels the preserved-helicity content of X^O_opp and X^E_opp is identical (same helicity components), so
W[X^O, X^O] = W[X^E, X^E] = Tpred(opp, opp) (1 + O(1/m^2)).  Ratios R = moment / prediction (t13 design, t18/t19 kernels), for k = 2..7,
nu = 0, 1, 3; raw for k >= 3, t19-cleaned for k = 2.
Parity test (a): the parity reflection of the frame maps q -> -q on all four circles; for a parity-odd wave W_{-q} = -W_q pointwise, so every
q -> -q symmetric kernel (K*, S_ff, A x A) gives exactly zero.  Reported: max_q |W_q + W_{-q}| / max_q |W_q| and the kernel contractions
relative to sqrt(|K(EE)| |K(OO)|) on the same points.
Usage: python3 t21_oddanalysis.py J [J ...]    -> out_t21_oddanalysis.txt
"""
import os, sys, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from t13_analysis_l2 import Modes, design, Tpred, CLASSES
from predict4 import weyl
from cases import case
import t19_k2 as C
import t5_analysis_l1 as T5
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

KS = {'K*': C.KS['K*'], 'A': C.KS['A'], 'Sff': C.KS['Sff']}

def load(J, combo, grid='half'):
    Dv, Jv = case(1.0, J)
    fn = os.path.join(HERE, 'cache', 'modes_odd_D%.6f_J%d_%s_%s.npz' % (Dv, Jv, combo, grid))
    if not os.path.exists(fn): return None
    mo = Modes(fn); mo.fn = fn
    return mo

def ratio(mo, k, nu, K, sign=1.0):
    """raw kernel ratio, wave stored under combo ('same','same') (ci = 0)"""
    wave = ('same', 'same')
    A = design(mo.Delta, mo.J, k, nu); Tp = Tpred('opp', 'opp', mo.Delta, mo.J, 8.0)
    terms = []; pred = 0j
    for (q1, q3), w13 in K.items():
        for (q2, q4), w24 in K.items():
            q = (q1, q3, q2, q4); v, sh = mo.moment(k, nu, wave, q)
            terms.append((w13*w24, v, sh)); pred += w13*w24*sum(A[q][c]*Tp[c] for c in CLASSES)
    shm = max(t[2] for t in terms)
    meas = sum(w*v*np.exp(sh - shm) for w, v, sh in terms)*(-1)**mo.J
    return sign*complex(meas*np.exp(shm - np.log(abs(pred)))/(pred/abs(pred))).real

def cleaned_k2(mo, nu, K, sign=1.0):
    """t19 cleaned extraction (prediction: opp, opp) -- t19 hard-wires LEAD/combo; here the wave is under ('same','same') with the opp
    prediction, so the opp prediction is injected by a thin wrapper"""
    orig = C.Tpred
    C.Tpred = lambda hL, hR, D, J, cT: orig('opp', 'opp', D, J, cT)
    try:
        r = C.cleaned(mo, 2, nu, ('same', 'same'), K)
    finally:
        C.Tpred = orig
    C._SET.clear()
    return sign*r['R']

if __name__ == "__main__":
    Js = [int(a) for a in sys.argv[1:]] or [35]
    pr("# Parity-odd exchange end to end (fourpt/tttt/t21_oddanalysis.py); beta = 1, half grid; R = moment/prediction(opp, opp)")
    data = {}
    for J in Js:
        EE, OO, EO = load(J, 'EE'), load(J, 'OO'), load(J, 'EO', 'half_n100')
        if EE is None or OO is None:
            pr(f"  J={J}: waves missing"); continue
        m = float(weyl(EE.Delta, EE.J)[2])
        pr(f"\n## J = {J}, Delta = {EE.Delta:.6f}, m = {m:.2f}")
        for k in (2, 3, 4, 5, 6, 7):
            for nu in (0.0, 1.0, 3.0):
                cells = []
                for kn, K in KS.items():
                    if (k % 2 == 1) != (kn == 'A'): continue
                    if k == 2:
                        rE = cleaned_k2(EE, nu, K); rO = cleaned_k2(OO, nu, K, -1.0)
                    else:
                        rE = ratio(EE, k, nu, K); rO = ratio(OO, k, nu, K, -1.0)
                    data[(J, k, nu, kn)] = (m, rE, rO)
                    cells.append(f"{kn}: EE {rE:.4f}  -OO {rO:.4f}  (-OO)/EE {rO/rE:.5f}")
                pr(f"  k={k} nu={nu}: " + " | ".join(cells))
        if EO is not None:
            pr("  parity test on W[E,O] (first 100 points): see below")
            mo = EO; M = mo.modes[('same', 'same')]
            qs = list(M.keys())
            num = max(np.max(np.abs(M[q] + M[tuple(-x for x in q)])) for q in qs); den = max(np.max(np.abs(M[q])) for q in qs)
            pr(f"   max_q |W_q + W_-q| / max_q |W_q| = {num/den:.2e}  (parity-odd wave: 0)")
            Mo = OO.modes[('same', 'same')]
            num2 = max(np.max(np.abs(Mo[q] - Mo[tuple(-x for x in q)])) for q in qs); den2 = max(np.max(np.abs(Mo[q])) for q in qs)
            pr(f"   control W[O,O]: max_q |W_q - W_-q| / max_q |W_q| = {num2/den2:.2e}  (parity-even wave: 0)")
    if len(Js) >= 2:
        pr("\n## m-dependence: m^2(R-1) for EE and -OO, and m^2((-OO)/EE - 1)")
        for k in (2, 3, 4, 5, 6, 7):
            for nu in (0.0, 1.0, 3.0):
                for kn in KS:
                    pts = [data[(J, k, nu, kn)] for J in Js if (J, k, nu, kn) in data]
                    if len(pts) < 2: continue
                    pr(f"  k={k} nu={nu} {kn}: EE " + "/".join(f"{m*m*(rE-1):+.0f}" for m, rE, rO in pts) + "  -OO " + "/".join(f"{m*m*(rO-1):+.0f}" for m, rE, rO in pts)
                       + "  ratio " + "/".join(f"{m*m*(rO/rE-1):+.1f}" for m, rE, rO in pts) + "  (m = " + ", ".join(f"{m:.1f}" for m, _, _ in pts) + ")")
    open(os.path.join(HERE, 'out_t21_oddanalysis.txt'), 'w').write("\n".join(out) + "\n")
