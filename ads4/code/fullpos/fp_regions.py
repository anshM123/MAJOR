"""
fullpos/fp_regions.py (FULLPOS_CERT.md §4): diagnostics for the asymptotic regions, from the certified moment balls.
 (i)  coefficient profile: |Z_k|, Y_k/Z_k, P_k = D_k^2 L_k (L_k = (Z_k - k Z_{k-1}/(k+3))/(k+3)), vs the §29 envelope
      P*_k = ((5/2)_k/k!) (A + B/log(k+2)), A = 9.7e-17, B = 1.21e-4.
 (ii) the comparison constant of ANGULAR_ROUTE §27 (bounded impact parameter), truncated at k <= K (hence a LOWER bound on K_b(r0)):
      K_b(r0) >= sup_{r <= r0} (1-r)^3 sum_{k<=K} D_k^2 [ r^{k-1} S_{b,k} |Z_k| |1/(k+3) - r/(k+1)| + r^k |2Y_k + V_{b,k} Z_k| + r^k Lam^2 |Q_k| ],
      S_{b,k} = sum_{j<k} u_j;  and the leading-Regge value H_LR(r0) = 2 OmLR_2(beta0) (target psi; float) for comparison.
"""
import sys, os, math, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
from flint import arb, fmpq, ctx
import fp_common as FC
import fp_poly as FP


def main():
    Dg = int(sys.argv[1]); K = int(sys.argv[2])
    Z, Y = FP.load_moments(Dg, K)
    lam = (Dg - 2)*(Dg - 1)
    Q = [arb(v) for v in FP.Q_exact(K, lam)]
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    from math import comb, lgamma, log, exp
    pr('# (i) coefficient profile (certified balls; printed to 4 digits)')
    pr('#   k     log10|Z_k|   sign   Y_k/Z_k      log10|P_k|   log10 P*_k   log10 |P_k/P*_k|')
    for k in list(range(0, 50, 5)) + list(range(50, K + 1, 25 if K <= 400 else 50)):
        Lk = (Z[k] - (arb(k)/(k + 3))*Z[k - 1])/(k + 3) if k > 0 else Z[0]/3
        Pk = comb(k + 3, 3)**2*Lk
        Ps = exp(lgamma(k + 2.5) - lgamma(2.5) - lgamma(k + 1))*(9.7e-17 + 1.21e-4/log(k + 2))
        lz = float(abs(Z[k]).log().mid())/log(10); lp = float(abs(Pk).log().mid())/log(10)
        pr(f'  {k:4d}  {lz:11.3f}   {"+" if Z[k] > 0 else "-"}   {float((Y[k]/Z[k]).mid()):11.4e}  {lp:11.3f}  {log(Ps)/log(10):11.3f}  {lp - log(Ps)/log(10):11.3f}')
    pr('# (ii) truncated comparison constant (lower bound on K_b(r0)) and H_LR(r0)')
    import ls_core as L
    for b in (0, 4):
        Lb = FC.ladder(10**9 if b == 0 else 10**9, b, T=1)
        for r0 in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5):
            best = None
            for r in (r0, r0/2, r0/4):
                rr = arb(r); tot = arb(0); S = 0
                for k in range(K + 1):
                    u_prev = (b + k - 1)*(b + k) if k > 0 else 0
                    if k > 0:
                        S += (b + k - 1)*(b + k)
                    V = arb(fmpq(4*k, k + 3)) if b == 0 else 4 + arb(12)/(k + 1)
                    D2 = comb(k + 3, 3)**2
                    t1 = (rr**(k - 1) if k > 0 else arb(0))*S*abs(Z[k])*abs(1/arb(k + 3) - rr/(k + 1))
                    t2 = rr**k*(abs(2*Y[k] + V*Z[k]) + lam*abs(Q[k]))
                    tot += D2*(t1 + t2)
                tot = (1 - rr)**3*tot
                best = tot if best is None or tot > best else best
            beta0 = math.acosh((1 + r0)/(1 - r0)); bh = 2*math.sqrt(lam)*math.sinh(beta0/2)
            hlr = 2*L.OmLR(2, [beta0], Dg)[0]
            Tneed = float(best.mid())/abs(hlr)
            pr(f'  b={b} r0={r0:7.0e} bhat0={bh:9.2f}  K_b(r0) >= {float(best.mid()):.3e}   H_LR(r0) = {hlr:.4e}   '
               f'K/H_LR = {Tneed:.3e}  (needs T > this; T >= Lam^2 C/... )  log10={math.log10(Tneed):.1f}')
    open(os.path.join(FC.CACHE, f'fp_regions_{Dg}_K{K}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
