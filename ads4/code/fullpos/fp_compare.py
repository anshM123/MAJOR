"""
fullpos/fp_compare.py (FULLPOS_CERT.md §2): finite-mass comparison of the angular action with T0's hybrid action.
  H_ang,b(Delta, J): exact finite ladder sum with the certified moments (fp_poly.Action.H_ang), b = 0 <-> '++' (even J), b = 4 <-> '+-'.
  H_hyb, N_s: T0's hybrid (largespin/ls_core.hyb_fast, read-only import; target psi, i.e. de-shifted weight), at the T0 kinematics
  m^2 = T - C (mfl = m/Lam), sinh(beta/2) = (J+1/2)/m.
Also the fixed-r leading-Regge limit: H_ang(T -> oo, r fixed) vs 2 OmLR_2(beta) (+ sigma u OmLR_3), cosh(beta) = (1+r)/(1-r).   [DV]
"""
import sys, os, time, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
from flint import arb, ctx
import fp_common as FC
import fp_poly as FP


def main():
    Dg = int(sys.argv[1]); K = int(sys.argv[2])
    import ls_core as L
    Z, Y = FP.load_moments(Dg, K)
    lam = (Dg - 2)*(Dg - 1)
    A = FP.Action(Dg, Z, Y, FP.Q_exact(K, lam))
    Lam = math.sqrt(lam)
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    pr('# fixed-r leading-Regge limit: H_ang at T = 1e30*lam vs 2*OmLR2(beta)   (cosh beta = (1+r)/(1-r); J chosen with C = r T)')
    for r in (1e-6, 1e-5, 4e-5, 1e-4):
        # choose T so that J(J+1) = r T for J = 40 (K-limited): T = J(J+1)/r ; then u = Lam^2/(T(1-r)) small
        for J in (20, 40):
            C = J*(J + 1); T = C/r
            Ha = A.H_ang(J, 0, T)
            beta = math.acosh((1 + r)/(1 - r))
            m = math.sqrt(T - C)
            u = lam/(T - C)
            lr = 2*L.OmLR(2, [beta], Dg)[0] + u*L.OmLR(3, [beta], Dg)[0]
            # T0 kinematics: sinh(beta0/2) = (J+1/2)/m
            beta0 = 2*math.asinh((J + 0.5)/m)
            lr0 = 2*L.OmLR(2, [beta0], Dg)[0] + u*L.OmLR(3, [beta0], Dg)[0]
            pr(f'r={r:.0e} J={J} T={T:.3e} bhat={2*Lam*math.sinh(beta/2):8.3f}  H_ang={float(Ha.mid()):.10f}  '
               f'2OmLR2+uOmLR3 (r)={lr:.10f}  (r_beta)={lr0:.10f}  rel.diff={(float(Ha.mid())-lr0)/lr0:.2e}')
    pr('# finite mass: (H_ang - H_hyb)/N_s')
    worst = 0
    for J in (0, 2, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 60):
        for tw in (1000, 1300, 2000):
            for b, h in ((0, '++'), (4, '+-')):
                if (b == 0 and J % 2) or (b == 4 and J < 4):
                    continue
                Delta = J + tw
                T = (Delta - 1)*(Delta - 2); C = J*(J + 1)
                Ha = float(A.H_ang(J, b, T).mid())
                mfl = math.sqrt(T - C)/Lam
                hy = L.hyb_fast(Dg, J, mfl, h)
                d = (Ha - hy['base'])/hy['Ns']
                worst = max(worst, abs(d))
                pr(f'J={J:3d} tw={tw} {h}  H_ang={Ha: .8e}  H_hyb={hy["base"]: .8e}  N_s={hy["Ns"]:.6e}  (H_ang-H_hyb)/N_s={d: .3e}'
                   f'  bhat={hy["bhat"]:.3f} u={hy["u"]:.4f}')
    pr(f'# max |H_ang - H_hyb|/N_s = {worst:.3e}')
    open(os.path.join(FC.CACHE, f'fp_compare_{Dg}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
