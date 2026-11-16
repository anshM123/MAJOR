"""t4: the two near-marginal flat points under the AdS uplift (hybrid model + envelopes).  Floating point.

 (i)  +-, J = 4, at the kinematic threshold tau = Dg:  m^2 = (Dg-2)(Dg+7) > Lam^2 = (Dg-2)(Dg-1), so u = Lam^2/m^2 < 1.
      Flat exact (Lemma P, completed): FF(4,u) = 2591/210 - 41u/18 - 10u^3.  Also with smearing scale Lam' = Dg (CHMRS convention).
 (ii) +-, near-zero J ~ 245, m_fl ~ 21.79: local minimum over J in [230, 260], m_fl in [20.5, 23.5] of the flat action, of the hybrid
      model, and the curvature shift; the envelope sizes C0 cosh(beta)/m^2 x N.
"""
import numpy as np
from common import FlatGrid, kin, tau_of_m, lam_of_gap, hybrid_shift

def FF4(u):
    return 2591 / 210 - 41 * u / 18 - 10 * u ** 3

print('(i) +-, J=4 at threshold')
print('   Dg      Lam    u(Lam)   FF_flat   shift(hyb)  N_cons   | u(Lam=Dg)  FF_flat')
G = FlatGrid(1600)
for Dg in [10, 20, 50, 100, 300, 1000, 10 ** 4]:
    Lam = lam_of_gap(Dg)
    tau = float(Dg)
    m, beta = kin(tau, 4)
    u = (Lam / m) ** 2
    Js, B2, B3, FW, ab = G.parts('+-', 4, u)
    sh = sum(hybrid_shift('+-', Lam, m, beta, 4))
    u2 = (Dg / m) ** 2
    print('  %5d %8.2f  %.5f  %.5f  %+.3e  %.2f   | %.5f  %.5f' % (Dg, Lam, u, FF4(u), sh, ab[0], u2, FF4(u2)))
print('   zero of FF(4,u) at u = %.6f (m_fl = %.6f); dFF/du(1) = %.4f' % (
    float(np.roots([-10, 0, -41 / 18, 2591 / 210])[np.isreal(np.roots([-10, 0, -41 / 18, 2591 / 210]))].real[0]),
    1 / np.sqrt(float(np.roots([-10, 0, -41 / 18, 2591 / 210])[np.isreal(np.roots([-10, 0, -41 / 18, 2591 / 210]))].real[0])),
    -41 / 18 - 30))

print('(ii) +- near-zero region')
mfls = np.linspace(20.5, 23.5, 61)
flat_min = (np.inf, None)
best = {}
for mfl in mfls:
    u = 1 / mfl ** 2
    Js, B2, B3, FW, ab = G.parts('+-', 262, u)
    F = B2 + B3 + FW
    sel = (Js >= 230) & (Js <= 260)
    k = np.argmin(np.where(sel, F, np.inf))
    if F[k] < flat_min[0]:
        flat_min = (F[k], (mfl, int(Js[k]), ab[k], B2[k], B3[k], FW[k]))
print('   flat min %.4e at m_fl=%.3f J=%d (N_cons=%.3f; parts %.5f %.5f %.5f)' % (flat_min[0], *flat_min[1]))
mfl0, J0 = flat_min[1][0], flat_min[1][1]
print('   Dg      Lam      beta     shift      shift*Lam^2   min_hyb(local)   env(C0=100,comp)  env(C0=100,cons)')
for Dg in [10, 20, 50, 100, 300, 1000]:
    Lam = lam_of_gap(Dg)
    mh = (np.inf, None)
    for mfl in np.linspace(20.5, 23.5, 31):
        u = 1 / mfl ** 2
        Js, B2, B3, FW, ab = G.parts('+-', 262, u)
        m = mfl * Lam
        for J in range(236, 256):
            i = J - 4
            tau = tau_of_m(m, J)
            if tau < Dg:
                continue
            _, beta = kin(tau, J)
            sh = sum(hybrid_shift('+-', Lam, m, beta, J))
            v = B2[i] + B3[i] + FW[i] + sh
            if v < mh[0]:
                mh = (v, (mfl, J))
    m = mfl0 * Lam
    tau = tau_of_m(m, J0)
    _, beta = kin(tau, J0)
    sh = sum(hybrid_shift('+-', Lam, m, beta, J0))
    u = 1 / mfl0 ** 2
    Js, B2, B3, FW, ab = G.parts('+-', J0, u)
    i = J0 - 4
    ecomp = 100 * np.cosh(beta) / m ** 2 * (abs(B2[i]) + abs(B3[i]) + abs(FW[i]))
    econs = 100 * np.cosh(beta) / m ** 2 * ab[i]
    print('  %5d %8.2f  %.5f  %+.4e  %+8.3f      %.4e at %s   %.2e        %.2e' % (Dg, Lam, beta, sh, sh * Lam ** 2, mh[0], (round(mh[1][0], 2), mh[1][1]), ecomp, econs))
