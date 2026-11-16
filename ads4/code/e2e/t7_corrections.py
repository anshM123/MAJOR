"""
Tasks 3, 4: measured/predicted ratios for the helicity-preserving and helicity-flip components and the first
corrections to the heavy-block formula, from the data of t4_spinning_e2e.py (data_t4_l{l}.pkl).

Components.  At leading order every mode satisfies  E_{q1q3} = pref * [ A_{q1q3} X + B_{q1q3} Y ],  with X = P_{-1/2+i nu}(cosh beta)
(helicity preserved at the spinning Regge pair) and Y = F^nu_{2l}(beta) (helicity flipped), A = M_pres, B = M_flip (predict.components).
"Measured" X and Y are obtained from the two purely transverse modes (l,-l) and (l,l) by solving the 2x2 system; the other modes
are then predictions (overdetermination check).   R_pres = X_meas/X, R_flip = Y_meas/Y.
Corrections:  R - 1 = c1/m + c2/m^2 + ... ;  tables of m (R-1) and m^2 (R-1) versus m, and Richardson estimates of c2.
Usage: python3 t7_corrections.py l
"""
import sys, pickle, numpy as np, mpmath as mp
from predict import components, scalar_prefactor, prediction
from formulas import Fnu_mp
from cases import info
mp.mp.dps = 25

def measured_components(l, Delta, J, Dphi, k, nu, res):
    pref, beta = scalar_prefactor(l, Delta, J, Dphi, k, nu)
    X0 = Fnu_mp(nu, 0, beta); Y0 = Fnu_mp(nu, 2*l, beta)
    A1, B1 = components(l, k, nu, l, -l); A2, B2 = components(l, k, nu, l, l)
    E1 = res[(k, nu, l, -l)]/pref; E2 = res[(k, nu, l, l)]/pref
    det = A1*B2 - A2*B1
    X = (E1*B2 - E2*B1)/det; Y = (A1*E2 - A2*E1)/det
    # overdetermination: all other modes
    worst = 0
    for q1 in range(-l, l+1):
        for q3 in range(-l, l+1):
            A, B = components(l, k, nu, q1, q3)
            pred = A*X + B*Y; ex = res[(k, nu, q1, q3)]/pref
            scale = abs(A*X0) + abs(B*Y0)
            if scale > 0:
                worst = max(worst, abs(ex - pred)/scale)
    return X/X0, Y/Y0, worst, beta

if __name__ == "__main__":
    from t4_spinning_e2e import main_cases, get_moments, moments_of, fname, DPHI
    import os
    l = int(sys.argv[1]); Dphi = DPHI
    ks = (2, 4, 6); nus = (0.0, 1.0, 3.0)
    TAILMAX = 1e-4          # points whose estimated numerical error (dropped small-u region) exceeds this are not used in the fits
    store = {}
    tails = {}
    for c in main_cases():
        if os.path.exists(fname(l, *c)):
            mom = get_moments(l, *c)
            store[c] = moments_of(mom, ks, nus)
            for k in ks:
                for nu in nus:
                    tails[c + (k, nu)] = max(mom.tail_estimate(k, nu, l, -l), mom.tail_estimate(k, nu, l, l))
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    cases = sorted(store.keys(), key=lambda c: (round(float(info(*c)[1]), 1), c[0]))
    betas = sorted(set(round(float(info(*c)[1]), 1) for c in cases))
    pr(f"# l = {l} (O = {'J' if l == 1 else 'T'}), Dphi = {Dphi}.  R_pres, R_flip: measured/predicted for the helicity-preserving / helicity-flip component.")
    pr("# 'other modes': largest deviation of the remaining modes from the two-component form, relative to |A X| + |B Y| (it is O(1/m^2)).")
    pr("# 'tail': estimated relative contribution of the dropped small-u region to the modes (l,-l), (l,l) (numerical error of the exact side;")
    pr("#         entries with tail > 1e-3 are not reliable).")
    summary = {}
    for beta in betas:
        for k in ks:
            for nu in nus:
                pr(f"# beta ~ {beta}, k = {k}, nu = {nu}")
                pr("#      m       R_pres        m(R-1)     m^2(R-1)  |   R_flip        m(R-1)     m^2(R-1)  | other modes | tail")
                seq = []
                for c in [c for c in cases if round(float(info(*c)[1]), 1) == beta]:
                    m, b_, eta = info(*c)
                    Rp, Rf, worst, bb = measured_components(l, c[0], c[1], Dphi, k, nu, store[c])
                    Rp, Rf = mp.re(Rp), mp.re(Rf)
                    if tails[c + (k, nu)] < TAILMAX:
                        seq.append((m, Rp, Rf))
                    flag = "" if tails[c + (k, nu)] < TAILMAX else "   (not used in the fits)"
                    pr(f"  {float(m):8.2f}  {float(Rp):11.7f}  {float(m*(Rp-1)):10.4f}  {float(m*m*(Rp-1)):10.3f}  |  {float(Rf):11.7f}  {float(m*(Rf-1)):10.4f}  {float(m*m*(Rf-1)):10.3f}  | {float(worst):.1e}    | {tails[c + (k, nu)]:.0e}{flag}")
                if len(seq) >= 3:
                    # two-point Richardson estimates from the two largest m:  (i) c1/m + c2/m^2  (ii) c2/m^2 + c4/m^4
                    (m1, p1, f1), (m2, p2, f2) = seq[-2], seq[-1]
                    def fit12(r1, r2):
                        # r = c1/m + c2/m^2
                        a11, a12, a21, a22 = 1/m1, 1/m1**2, 1/m2, 1/m2**2
                        det = a11*a22 - a12*a21
                        return ((r1*a22 - r2*a12)/det, (a11*r2 - a21*r1)/det)
                    def fit24(r1, r2):
                        a11, a12, a21, a22 = 1/m1**2, 1/m1**4, 1/m2**2, 1/m2**4
                        det = a11*a22 - a12*a21
                        return ((r1*a22 - r2*a12)/det, (a11*r2 - a21*r1)/det)
                    c1p, c2p = fit12(p1-1, p2-1); c1f, c2f = fit12(f1-1, f2-1)
                    d2p, d4p = fit24(p1-1, p2-1); d2f, d4f = fit24(f1-1, f2-1)
                    # test of each ansatz on the third-largest m
                    m0, p0, f0 = seq[-3]
                    e12p = (c1p/m0 + c2p/m0**2) - (p0-1); e24p = (d2p/m0**2 + d4p/m0**4) - (p0-1)
                    e12f = (c1f/m0 + c2f/m0**2) - (f0-1); e24f = (d2f/m0**2 + d4f/m0**4) - (f0-1)
                    pr(f"#   fits from the two largest m:  preserved: [c1/m + c2/m^2] c1 = {float(c1p):.4f}, c2 = {float(c2p):.2f};   [c2/m^2 + c4/m^4] c2 = {float(d2p):.2f}")
                    pr(f"#                                 flip:      [c1/m + c2/m^2] c1 = {float(c1f):.4f}, c2 = {float(c2f):.2f};   [c2/m^2 + c4/m^4] c2 = {float(d2f):.2f}")
                    pr(f"#   misfit at the third largest m (m = {float(m0):.1f}), relative to R-1:  preserved: odd+even ansatz {float(e12p/(p0-1)):+.3f}, even ansatz {float(e24p/(p0-1)):+.3f};  flip: {float(e12f/(f0-1)):+.3f}, {float(e24f/(f0-1)):+.3f}")
                    c3 = (float('nan'),)*4
                    if len(seq) >= 4:
                        # three-parameter fits with the three largest m, tested on the fourth largest:
                        #   (A) c1/m + c2/m^2 + c4/m^4      (B) c2/m^2 + c4/m^4 + c6/m^6
                        ms = np.array([float(t[0]) for t in seq]); m4 = ms[-4]
                        line = []; c3 = []
                        for name, idx in (("preserved", 1), ("flip", 2)):
                            R = np.array([float(t[idx]) - 1 for t in seq])
                            cA = np.linalg.solve(np.stack([1/ms[-3:], 1/ms[-3:]**2, 1/ms[-3:]**4], 1), R[-3:])
                            cB = np.linalg.solve(np.stack([1/ms[-3:]**2, 1/ms[-3:]**4, 1/ms[-3:]**6], 1), R[-3:])
                            eA = cA @ np.array([1/m4, 1/m4**2, 1/m4**4])/R[-4] - 1; eB = cB @ np.array([1/m4**2, 1/m4**4, 1/m4**6])/R[-4] - 1
                            line.append(f"{name}: (A) c1 = {cA[0]:+.4f}, c2 = {cA[1]:.2f}; (B) c2 = {cB[0]:.2f}; misfit at m = {m4:.1f}: (A) {eA:+.3f}, (B) {eB:+.3f}")
                            c3 += [float(cA[0]), float(cB[0])]
                        pr("#   three-parameter fits (three largest m):  " + "   ||   ".join(line))
                        c3 = tuple(c3)
                    summary[(beta, k, nu)] = (float(d2p), float(d2f), float(c1p), float(c1f)) + c3
    pr("# ---- summary: O(1/m^2) coefficients c2 (even ansatz, two largest m) and the O(1/m) coefficient of the odd+even ansatz")
    pr("#       and of the three-parameter fits: c1 of ansatz (A) [1/m, 1/m^2, 1/m^4], c2 of ansatz (B) [1/m^2, 1/m^4, 1/m^6]")
    pr("#  beta   k   nu     c2_pres     c2_flip   |  c1_pres    c1_flip  ||  3-par: c1_pres   c2_pres    c1_flip   c2_flip")
    for (beta, k, nu), v in sorted(summary.items()):
        pr(f"  {beta:4.1f}   {k}  {nu:4.1f}  {v[0]:10.2f}  {v[1]:10.2f}   | {v[2]:8.4f}  {v[3]:8.4f}  ||       {v[4]:+8.4f}  {v[5]:9.2f}   {v[6]:+8.4f}  {v[7]:9.2f}")
    open(f'out_t7_corrections_l{l}.txt', 'w').write("\n".join(out) + "\n")
    pickle.dump(summary, open(f'data_t7_l{l}.pkl', 'wb'))
