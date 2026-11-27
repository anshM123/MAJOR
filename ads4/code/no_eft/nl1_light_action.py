r"""
no_eft/nl1_light_action.py  [NE]  --  functional-side input (W*) of NOTE.md, Sec. 4.

Stripped action of the cubic functional Omega^EE (Delta_g = 1000, published endpoint-repaired continuum weight) on LIGHT
double-twist states, from the exact finite ladder formula of sec_fullpos.tex, Eqs. (fp-target), (fp-H):

   Omega^EE_O = (2 pi G Lam^4 / 9) * sum_b |a_{O,b,0}|^2 * hhat_b(Delta, J),
   hhat_b     = T * sum_{k=0}^{J-b} D_k^2 U_{b,k} [ Z_k (1/(k+3) - r/(k+1)) + (2 Y_k + V_{b,k} Z_k + sigma_b Lam^2 Q_k)/T ],

i.e. B * H_ang with the factor B (1-r)^3 = (2 pi G Lam^4/9) T taken out (finite at r = 1, where B and (1-r)^3 are separately
singular).  The detector amplitude |a_{O,b,0}|^2 carries the double zero sin^2 / cos^2 (L1); hhat is the smooth remainder.

CAVEAT (grade): Eq. (fp-target) is derived in the paper for heavy primaries.  Its derivation (ladder ratios from light-transform
conservation, first-moment matrix, finite ladder, nu-moments of polynomials) is algebraic in (Delta, J) and has no step that uses
Delta - J >= Delta_g, so we evaluate it at light (Delta, J); that continuation is [D-sketch], and the numbers below are [NE].

Moments: read-only from ads4/code/fullpos/cache/moments_1000_K<K>.json (arb balls, as used by the FP certificates); Q_k exact.
Output: no_eft/cache/nl1_<K>.txt.   Usage: [NL1_DG=300] python3 nl1_light_action.py [K]   (K = 300 or 1000; Delta_g = 1000 default, 300 with K = 300)
"""
import sys, os, time, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
FPDIR = os.path.join(HERE, '..', 'fullpos')
sys.path.insert(0, FPDIR)
import fp_common as FC          # read-only imports
import fp_poly as FP
from flint import arb, fmpq

Dg = int(os.environ.get('NL1_DG', '1000'))


def hhat(A, J, b, Delta):
    L = FC.ladder(J, b, T=1)
    n, C, u, D, V, sig = L['n'], L['C'], L['u'], L['D'], L['V'], L['sig']
    T = arb(Delta - 1)*arb(Delta - 2)          # Delta may be a Fraction-like float; we pass exact rationals below
    r = arb(C)/T
    U = arb(1); tot = arb(0)
    for k in range(n + 1):
        Vk = arb(fmpq(V[k].numerator, V[k].denominator))
        tot += D[k]**2*U*(A.Z[k]*(1/arb(k + 3) - r/(k + 1)) + (2*A.Y[k] + Vk*A.Z[k] + sig*A.Lam2*A.Q[k])/T)
        if k < n:
            U = U*(C - u[k])/(T - u[k])
    return T*tot, T, r


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    t0 = time.time()
    Z, Y = FP.load_moments(Dg, K)
    lam = (Dg - 2)*(Dg - 1)
    Q = FP.Q_exact(K, lam)
    A = FP.Action(Dg, Z, Y, Q)
    os.makedirs(os.path.join(HERE, 'cache'), exist_ok=True)
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    pr(f'# nl1: Delta_g={Dg}, moments K={K}, setup {time.time()-t0:.1f}s')
    pr('# columns: b  twist  J  sign  log10|hhat_b|  log10(relative ball radius)  local slope d log|hhat| / d log J')
    Js = [J for J in (4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 300, 384, 512, 640, 768, 896, 1000) if J <= K]
    for b in (0, 4):
        for tw in (2, 3, 4, 6, 12, 42, 202, 998):
            prev = None
            for J in Js:
                if b == 0 and J % 2:
                    continue
                if J - b > K:
                    continue
                h, T, r = hhat(A, J, b, J + tw)
                if 0 in h:
                    pr(f'{b} {tw:4d} {J:5d}  ball contains 0 (precision lost): rad={h.rad()}')
                    prev = None; continue
                sgn = '+' if h > 0 else '-'
                lg = float((abs(h)).log()/arb(10).log())
                rel = float((h.rad()/abs(h.mid())).log()/arb(10).log()) if h.rad() != 0 else -999
                sl = ''
                if prev is not None:
                    sl = f'{(lg - prev[1])*math.log(10)/(math.log(J) - math.log(prev[0])):9.3f}'
                pr(f'{b} {tw:4d} {J:5d}  sign {sgn}  log10|hhat| = {lg:10.3f}  log10(rel.rad) = {rel:8.1f}  slope {sl}')
                prev = (J, lg)
    pr(f'# total time {time.time()-t0:.1f}s')
    open(os.path.join(HERE, 'cache', f'nl1_{Dg}_{K}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
