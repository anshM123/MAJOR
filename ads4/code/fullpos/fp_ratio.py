"""
fullpos/fp_ratio.py (FULLPOS_CERT.md §2.3): independent [DV] check of the authors' J <= 26 statement
   |H_ang,b - H_hyb,b| <= 0.01 H_hyb,b   for every mass (u = Lam^2/m^2 in (0, umax], umax = 999/(999+2J) at Dg = 1000),
H_ang from the certified moments (exact finite ladder sum, evaluated in Arb), H_hyb from T0's hyb_fast (target psi; the repaired
weight changes cv2 by <= 2 eps_W/Lam ~ 5e-9, recorded but not added).  Sampled in u (NOT a certificate): NU points, geometric
towards u -> 0 plus uniform near umax.
"""
import sys, os, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
import fp_common as FC
import fp_poly as FP


def main():
    Dg = int(sys.argv[1]); K = int(sys.argv[2]); Jmax = int(sys.argv[3]); NU = int(sys.argv[4])
    import ls_core as L
    Z, Y = FP.load_moments(Dg, K)
    lam = (Dg - 2)*(Dg - 1); Lam = math.sqrt(lam)
    A = FP.Action(Dg, Z, Y, FP.Q_exact(K, lam))
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    worst = (0, None)
    for J in range(0, Jmax + 1):
        for b, h in ((0, '++'), (4, '+-')):
            if (b == 0 and J % 2) or (b == 4 and J < 4):
                continue
            C = J*(J + 1)
            Tmin = (J + Dg - 1)*(J + Dg - 2)
            umax = lam/(Tmin - C)
            us = sorted(set([umax*(1 - i/(NU//2)) for i in range(NU//2)] + [umax*10**(-4*i/(NU//2)) for i in range(1, NU//2 + 1)]))
            wJ = (0, None); minrat = 9
            for u in us:
                m2 = lam/u; T = m2 + C
                Ha = float(A.H_ang(J, b, T).mid())
                hy = L.hyb_fast(Dg, J, math.sqrt(m2)/Lam, h)
                rel = (Ha - hy['base'])/hy['base']
                minrat = min(minrat, Ha/hy['Ns'])
                if abs(rel) > wJ[0]:
                    wJ = (abs(rel), u)
            worst = max(worst, (wJ[0], (J, h, wJ[1])))
            pr(f'J={J:3d} {h}  max_u |H_ang/H_hyb - 1| = {wJ[0]:.3e} (at u={wJ[1]:.4f}, umax={umax:.4f})   min_u H_ang/N_s = {minrat:.5f}'
               f'   {"OK(<0.01)" if wJ[0] < 0.01 else "EXCEEDS 0.01"}')
    pr(f'# worst: {worst}')
    open(os.path.join(FC.CACHE, f'fp_ratio_{Dg}_J{Jmax}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
