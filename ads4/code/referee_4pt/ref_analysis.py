"""Referee scratch analysis (independent of t5's reporting choices).
Reads the authors' cached modes (read-only) or the referee's regenerated ones.
 (1) raw-mode ratios (no class fit) for the dominant modes of each combo, all masses; 1/m fits on raw modes;
 (2) S-term sensitivity: refit with S scaled (x0, x0.5, x2) and with S -> -S at the right pair only;
 (3) own mutation: CHL-type non-Weyl-even helicity factor kappa_same = Delta(Delta-1)/((3-Delta)(2-Delta)) on the right vertex
     of the same-helicity structure (the 'helicity-dependent shadow eigenvalue' scenario), and A_-^2 -> J^2 in V_2(opp);
 (4) os/pf at beta~0.4, nu=0: ansatz comparison with all 5 masses.
Usage: python3 ref_analysis.py DATADIR
"""
import os, sys, glob, itertools
import numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'fourpt'))
sys.path.insert(0, os.path.join(HERE, 'e2e'))
import predict4
from cases import info
from tmoments4 import TensorMoments4
import t5_analysis_l1 as t5
mp.mp.dps = 20
DATADIR = sys.argv[1]
QS, CLASSES = t5.QS, t5.CLASSES
OO, SS, PF = (1, 1, -1, -1), (1, -1, -1, 1), (1, 1, -1, 1)

_orig_Mval = predict4.Mval
def make_Mval(sL=1.0, sR=1.0):
    from predict4 import Mfun, gamma_a
    def Mval(l, q1, q3, h1, h3, nu, a, _right=False):
        R, S, mu = Mfun(l)[(q1, q3, h1, h3)]
        g1 = gamma_a(a - 1, nu)/gamma_a(a, nu)
        return R(nu, a) + g1*S(nu, a), mu
    return Mval

def design(Delta, J, k, nu, sL=1.0, sR=1.0):
    from predict4 import Mfun, gamma_a
    Ap, Am, m, beta = predict4.weyl(Delta, J)
    a = k + 3
    pre = predict4.prefactor(1, Delta, J, k, nu)
    def MV(q1, q3, h1, h3, nuu, s):
        R, S, mu = Mfun(1)[(q1, q3, h1, h3)]
        g1 = gamma_a(a - 1, nuu)/gamma_a(a, nuu)
        return R(nuu, a) + s*g1*S(nuu, a), mu
    A = {}
    for q in QS:
        q1, q3, q2, q4 = q
        row = {}
        for c, hs in CLASSES.items():
            acc = mp.mpc(0)
            for h in hs:
                ML, muL = MV(q1, q3, h[0], h[2], nu, sL)
                if ML == 0: continue
                MR, muR = MV(q2, q4, h[1], h[3], -nu, sR)
                if MR == 0: continue
                acc += ML*MR*predict4.Kcached(nu, muL, muR, beta)
            row[c] = complex(pre*acc)
        A[q] = row
    return A

def Tpred(cL, cR, Delta, J, mut=None):
    Ap, Am, m, beta = predict4.weyl(Delta, J)
    D = mp.mpf(Delta)
    def V(i, h1, h2, side):
        same = (h2 == -h1)
        if mut == 'Jsq' and i == 2 and not same:
            return mp.mpf(J)**2/2
        return predict4.V_l1(i, Ap, Am, m, h1, h2)
    def Vc(c, h1, h2, side):
        return sum(c[i]*V(i+1, h1, h2, side) for i in range(2))
    T = {}
    for cl in CLASSES:
        vl, vr = Vc(cL, cl[0], cl[1], 'L'), Vc(cR, cl[2], cl[3], 'R')
        if mut == 'kappa' and cl[3] == -cl[2]:          # right vertex same helicity gets CHL-type non-Weyl-even factor
            vr *= D*(D - 1)/((3 - D)*(2 - D))
        T[cl] = float(vl*vr)
    return T

def raw_and_fit(mom, k, nu, name, sL=1.0, sR=1.0, mut=None):
    D, J = mom.Delta, mom.J
    cL, cR = t5.combos(D, J)[name]
    A = design(D, J, k, nu, sL, sR)
    Tp = Tpred(cL, cR, D, J, mut)
    pred = {q: sum(A[q][c]*Tp[c] for c in CLASSES) for q in QS}
    big = max(abs(v) for v in pred.values()); logref = np.log(big)
    meas = t5.measured_combo(mom, cL, cR, k, nu, logref)
    R, sol, resid = t5.fit_components(A, meas, Tp, logref)
    # raw ratios of the 4 largest predicted modes
    top = sorted(QS, key=lambda q: -abs(pred[q]))[:6]
    raw = {q: (meas[q]/(pred[q]*np.exp(-logref))) for q in top}
    pp, pm = Tp[OO] + Tp[SS], Tp[OO] - Tp[SS]
    Rm = (sol[OO] - sol[SS])/pm if abs(pm) > 1e-9*abs(pp) else None
    return R, resid, raw, Rm, {q: abs(pred[q])/big for q in top}

def fitc1(ms, Rs, kind='A'):
    ms = np.array(ms); Rs = np.array(Rs)
    if kind == 'A':
        X = np.vstack([1/ms[-3:], 1/ms[-3:]**2, 1/ms[-3:]**4]).T; return np.linalg.solve(X, Rs[-3:] - 1)[0]
    if kind == 'L4':   # least squares c1,c2,c3 on 4 largest
        X = np.vstack([1/ms[-4:], 1/ms[-4:]**2, 1/ms[-4:]**3]).T; return np.linalg.lstsq(X, Rs[-4:] - 1, rcond=None)[0][0]
    if kind == 'L5':   # c1..c4 on 5
        X = np.vstack([1/ms, 1/ms**2, 1/ms**3, 1/ms**4]).T; return np.linalg.lstsq(X, Rs - 1, rcond=None)[0][0]

if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(DATADIR, 'modes4_l1_*.npz')))
    moms = [TensorMoments4.load(f) for f in files]
    byb = {}
    for mom in moms:
        m, beta, _ = info(mom.Delta, mom.J)
        byb.setdefault(round(float(beta), 1), []).append((float(m), mom))
    for b in sorted(byb):
        rows = sorted(byb[b], key=lambda r: r[0]); ms = [r[0] for r in rows]
        print(f"\n##### beta ~ {b}; masses {['%.1f' % x for x in ms]}")
        for nu in (0.0, 1.0):
            for name in ('ss', 'oo', 'os', 'v2v2'):
                res = [raw_and_fit(mom, 6, nu, name) for (_, mom) in rows]
                # raw modes common to all masses
                qs = list(res[-1][2].keys())[:4]
                for q in qs:
                    rr = [r[2].get(q, np.nan) for r in res]
                    if any(np.iscomplexobj(x) and abs(np.imag(x)) > 1e-3*abs(x) for x in rr):
                        pass
                    rr = [np.real(x) for x in rr]
                    if any(np.isnan(rr)): continue
                    print(f"  nu={nu} {name:4s} raw q={q} rel.size {res[-1][4][q]:.2e}: R = " + " ".join(f"{x:.5f}" for x in rr)
                          + f" | c1(A)={fitc1(ms, rr, 'A'):+.3f} c1(L4)={fitc1(ms, rr, 'L4'):+.3f} c1(L5)={fitc1(ms, rr, 'L5'):+.3f}")
        # S-term sensitivity at the largest mass
        mom = rows[-1][1]
        print("  S-sensitivity (largest m, k=6, nu=0): combo v2v2 and v1v2")
        for name in ('v2v2', 'v1v2', 'ss'):
            for (sL, sR, lab) in ((1, 1, 'S as predicted'), (0.5, 0.5, 'S x0.5 both'), (2, 2, 'S x2 both'), (1, -1, 'S -> -S right'), (0, 0, 'S = 0')):
                try:
                    R, resid, raw, Rm, _ = raw_and_fit(mom, 6, 0.0, name, sL, sR)
                    print(f"    {name} {lab:16s}: R_ss={R[SS] if R[SS] is None else round(R[SS],4)} R_oo={R[OO] if R[OO] is None else round(R[OO],4)} R(Sigma-)={Rm if Rm is None else round(Rm,4)} resid={resid:.2e}")
                except Exception as e:
                    print(f"    {name} {lab}: {e}")
        # own mutations, c1 on fitted classes
        for nu in (0.0, 3.0):
            for (name, cl, mut) in (('ss', SS, 'kappa'), ('v2v2', OO, 'Jsq'), ('oo', OO, 'Jsq'), ('ss', SS, None), ('oo', OO, None)):
                Rs = [raw_and_fit(mom_, 6, nu, name, mut=mut)[0][cl] for (_, mom_) in rows]
                print(f"  mutation {mut} nu={nu} {name} {t5.CNAME[cl]}: R = " + " ".join(f"{x:.4f}" for x in Rs) + f" | c1(A)={fitc1(ms, Rs, 'A'):+.3f} c1(L5)={fitc1(ms, Rs, 'L5'):+.3f}")
        # os pf ansatz comparison
        for nu in (0.0, 1.0, 3.0):
            Rs = np.array([raw_and_fit(mom_, 6, nu, 'os')[0][PF] for (_, mom_) in rows])
            msa = np.array(ms)
            for lab, cols in (('c2,c4 (5m)', [2, 4]), ('c2,c4,c6 (5m)', [2, 4, 6]), ('c1,c2,c4 (5m)', [1, 2, 4]), ('c1,c2 (4m)', [1, 2]), ('c2,c3 (4m)', [2, 3])):
                sel = slice(-4, None) if '4m' in lab else slice(None)
                X = np.vstack([1/msa[sel]**p for p in cols]).T; y = Rs[sel] - 1
                c, *_ = np.linalg.lstsq(X, y, rcond=None); r = np.linalg.norm(X@c - y)
                print(f"  os pf nu={nu} fit {lab:14s}: coeffs {np.round(c, 3)}  rms resid {r:.1e}")
