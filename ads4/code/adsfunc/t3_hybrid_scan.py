"""t3: physical heavy states (tau >= Delta_gap, J allowed) of the AdS functional Omega_g3 with smearing scale
Lam = sqrt((Dg-2)(Dg-1)) (the smallest heavy m), in the HYBRID MODEL

   F_hyb = F_flat(J, u; h) + [Regge leading - its flat limit]     (common.hybrid_shift)

which reduces to the flat action for beta -> 0 and to the leading Regge action for u -> 0 at fixed beta.  [model, I]
Unknown AdS 1/m^2 terms are modelled by an envelope C0 cosh(beta)/m^2 x N, with two norms N:
   N_cons = int |Psi_h||D_h| + |Phi| u^3   (conservative),   N_comp = |B2| + |B3| + |FW|   (component-wise).
Criteria reported per (Dg, h, C0):
   (A) F_hyb - C0 cosh(beta)/m^2 N_comp > 0          (model taken at face value)
   (B) F_flat - |shift| - C0 cosh(beta)/m^2 N_cons > 0 (sign of the curvature shift NOT trusted; conservative norm),
       evaluated only on states with beta <= 0.3 (bulk-point side; at beta ~ 1 the flat action is not the reference).
Grid: m_fl = m/Lam in [1, 1.3] step 0.01, 120 log points in [1.3, 60], and [20.5, 23.5] step 0.1; J <= 60 all, then ~4% steps up to
min(3000, 200 m_fl), plus all J in [230, 260] for 20.5 <= m_fl <= 23.5 (the flat +- near-zero).
Floating point; not a certificate.
"""
import sys, json
import numpy as np
from common import FlatGrid, kin, tau_of_m, lam_of_gap, hybrid_shift, SIGMA

DGS = [int(a) for a in sys.argv[1:]] or [20, 50, 100, 300, 1000]
C0S = [0, 100, 300, 1000]
MFL = np.unique(np.concatenate([np.arange(1.0, 1.3, 0.01), np.logspace(np.log10(1.3), np.log10(60), 120),
                               np.arange(20.5, 23.5, 0.1)]))   # + fine box around the flat +- near-zero

def jlist(h, mfl):
    jmax = int(min(3000, 200 * mfl))
    js = list(range(0, 61))
    j = 61
    while j <= jmax:
        js.append(j)
        j = int(j * 1.04) + 1
    if 20.4 < mfl < 23.6:
        js += list(range(230, 261))                                   # all J near the flat +- near-zero (J ~ 244)
    js = sorted(set(x for x in js if x <= jmax))
    if h == '++':
        js = sorted(set(x + (x % 2) for x in js))
    else:
        js = [x for x in js if x >= 4]
    return np.array(js), jmax

G = FlatGrid(1600)
summary = {}
for Dg in DGS:
    Lam = lam_of_gap(Dg)
    for h in ['++', '+-']:
        worstA = {C: (np.inf, None) for C in C0S}
        worstB = {C: (np.inf, None) for C in C0S}
        failA = {C: 0 for C in C0S}
        failB = {C: 0 for C in C0S}
        nst = 0
        maxshift_rel = (0, None)
        for mfl in MFL:
            u = 1 / mfl ** 2
            js, jmax = jlist(h, mfl)
            Jall, B2, B3, FW, ab = G.parts(h, int(js.max()) + 1, u)
            idx = {int(J): i for i, J in enumerate(Jall)}
            m = mfl * Lam
            for J in js:
                tau = tau_of_m(m, J)
                if tau < Dg - 1e-9:
                    continue
                i = idx[int(J)]
                _, beta = kin(tau, J)
                dB2, dB3, dFW = hybrid_shift(h, Lam, m, beta, J)
                Ff = B2[i] + B3[i] + FW[i]
                sh = dB2 + dB3 + dFW
                Fh = Ff + sh
                ncomp = abs(B2[i] + dB2) + abs(B3[i] + dB3) + abs(FW[i] + dFW)
                ncons = ab[i]
                nst += 1
                if abs(sh) / ncons > maxshift_rel[0]:
                    maxshift_rel = (abs(sh) / ncons, (round(mfl, 3), int(J), round(beta, 4)))
                for C in C0S:
                    e = C * np.cosh(beta) / m ** 2
                    a = (Fh - e * ncomp) / ncons
                    b = (Ff - abs(sh) - e * ncons) / ncons
                    if a <= 0: failA[C] += 1
                    if beta > 0.3: b = np.inf
                    if b <= 0: failB[C] += 1
                    if a < worstA[C][0]: worstA[C] = (a, (round(mfl, 4), int(J), round(beta, 5), float(Ff), float(sh)))
                    if b < worstB[C][0]: worstB[C] = (b, (round(mfl, 4), int(J), round(beta, 5), float(Ff), float(sh)))
        summary[(Dg, h)] = (failA, failB, worstA, worstB)
        print('Dg=%5d Lam=%8.3f %s  states=%d   max |shift|/N_cons = %.2e at (m_fl,J,beta)=%s' % (Dg, Lam, h, nst, *maxshift_rel))
        for C in C0S:
            print('    C0=%5d  (A) fails %5d  min rel %.3e at %s' % (C, failA[C], worstA[C][0], worstA[C][1]))
            print('             (B) fails %5d  min rel %.3e at %s' % (failB[C], worstB[C][0], worstB[C][1]))
        sys.stdout.flush()
