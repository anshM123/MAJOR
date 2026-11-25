"""
fullpos/fp_tail.py (FULLPOS_CERT.md §4): the ANGULAR_ROUTE §29/§31/§32 coefficient-envelope (tail) certificate, evaluated with the
certified moments.
  P_k = D_k^2 L_k,  L_k = (Z_k - (k/(k+3)) Z_{k-1})/(k+3),  L_0 = Z_0/3;   P*_k = ((5/2)_k/k!)(A + B/log(k+2)), A = 9.7e-17, B = 1.21e-4.
  (i)   envelope test |P_k - P*_k| <= P*_k/4 and |D_k^2 Q_k| <= (k+1) P*_k  for every certified k;
  (ii)  prefix constants D0..D3 over k < K (certified lower/upper balls) and the cutoff w_* of §32 *if* the envelopes held for k >= K;
  (iii) Cauchy majorants of §32.4: |P_k| <= 200 Lam^{2k}, |D_k^2 Q_k| <= 120 Lam^{2k}, for k <= K;
  (iv)  §32 single-sequence representation H^(2) = (1-r)^3 sum_{k<=n+1} V_{b,k} P_k vs the direct k=2 ladder sum (exactness check).
"""
import sys, os, math
sys.dont_write_bytecode = True
from flint import arb, fmpq
from math import comb
import fp_common as FC
import fp_poly as FP


def main():
    Dg = int(sys.argv[1]); K = int(sys.argv[2]); Kgood = int(sys.argv[3]) if len(sys.argv) > 3 else K
    Z, Y = FP.load_moments(Dg, K)
    lam = (Dg - 2)*(Dg - 1)
    Q = [arb(v) for v in FP.Q_exact(K, lam)]
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    A, B = arb('9.7e-17'), arb('1.21e-4')
    P = []
    for k in range(Kgood + 1):
        L = Z[0]/3 if k == 0 else (Z[k] - arb(k)/(k + 3)*Z[k - 1])/(k + 3)
        P.append(comb(k + 3, 3)**2*L)
    Ps = []
    g = arb(1)
    for k in range(Kgood + 1):
        if k > 0:
            g = g*(arb(3)/2 + k)/k                 # (5/2)_k/k!
        Ps.append(g*(A + B/arb(k + 2).log()))
    nenv1 = sum(1 for k in range(Kgood + 1) if abs(P[k] - Ps[k]) <= Ps[k]/4)
    nenv2 = sum(1 for k in range(Kgood + 1) if abs(comb(k + 3, 3)**2*Q[k]) <= (k + 1)*Ps[k])
    pr(f'(i) k <= {Kgood}: envelope |P_k-P*_k| <= P*_k/4 holds for {nenv1} values of k; |D_k^2 Q_k| <= (k+1)P*_k for {nenv2}')
    for k in (0, 1, 2, 5, 10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000):
        if k > Kgood:
            break
        pr(f'   k={k:4d}  log10|P_k| = {float(abs(P[k]).log().mid())/math.log(10):9.2f}   log10 P*_k = '
           f'{float(Ps[k].log().mid())/math.log(10):7.2f}   log10|D_k^2 Q_k| = {float(abs(comb(k+3,3)**2*Q[k]).log().mid())/math.log(10) if Q[k] != 0 else -999:9.2f}'
           f'   rel.rad P_k = {float(P[k].rad()/abs(P[k].mid())):.1e}')
    # (ii) prefix constants
    D0 = D1 = D2 = D3 = arb(0)
    S4 = 0
    for k in range(Kgood):
        d = abs(P[k] - Ps[k])
        D0 += d; D1 += S4*d; D2 += (k + 5)**2*d; D3 += abs(comb(k + 3, 3)**2*Q[k])
        S4 += (4 + k)*(5 + k)
    lg = lambda x: float(x.log().mid())/math.log(10)
    pr(f'(ii) prefix constants over k < {Kgood}: log10 D0 = {lg(D0):.2f}, D1 = {lg(D1):.2f}, D2 = {lg(D2):.2f}, D3 = {lg(D3):.2f}')
    ws = [((A/(4*D0)).log()*arb(2)/5), ((A/(4*D1)).log()*arb(2)/11), ((A/(4*D2)).log()*arb(2)/9), ((A/D3).log()*arb(2)/7)]
    lw = min(float(w.mid()) for w in ws)
    lw = min(lw, math.log(1/math.cosh(50)**2))
    beta_star = 2*(lw/(-2) + math.log(2))                     # 2 arcosh(w^{-1/2}) ~ 2 log(2 w^{-1/2})
    pr(f'     IF the envelopes held for all k >= {Kgood} (they do not, (i)): log10 w_* = {lw/math.log(10):.1f}, beta_* = {beta_star:.0f}')
    # (iii) Cauchy majorants
    Lam2 = arb(lam)
    okP = all(abs(P[k]) <= 200*Lam2**k for k in range(Kgood + 1))
    okQ = all(abs(comb(k + 3, 3)**2*Q[k]) <= 120*Lam2**k for k in range(Kgood + 1))
    mP = max(float((abs(P[k])/Lam2**k).log().mid()) for k in range(Kgood + 1))/math.log(10)
    mQ = max(float((abs(comb(k + 3, 3)**2*Q[k])/Lam2**k).log().mid()) for k in range(Kgood + 1) if Q[k] != 0)/math.log(10)
    pr(f'(iii) Cauchy majorants for k <= {Kgood}: |P_k| <= 200 Lam^2k: {okP} (max log10 |P_k|/Lam^2k = {mP:.2f});'
       f'  |D_k^2 Q_k| <= 120 Lam^2k: {okQ} (max log10 = {mQ:.2f})')
    # (iv) single-sequence representation
    Act = FP.Action(Dg, Z, Y, [fmpq(0)]*(K + 1))       # Q = 0: k = 2 part only
    for (J, b, Delta) in ((10, 0, 1010), (40, 4, 1040), (100, 0, 1500)):
        if J - b + 1 > Kgood:
            continue
        T = (Delta - 1)*(Delta - 2); C = J*(J + 1); n = J - b
        Ta = arb(T); r = arb(C)/Ta
        u = lambda j: (b + j)*(b + j + 1)
        U = [arb(1)]
        for j in range(n + 1):
            U.append(U[-1]*(C - u(j))/(Ta - u(j)) if j < n else arb(0))
        tot = arb(0)
        for k in range(n + 2):
            Uk = U[k]; Ukm = U[k - 1] if k > 0 else arb(0)
            um = u(k - 1)                         # u_{-1} = (b-1)b (= 12 for b = 4) at k = 0
            Vk = Uk*(1 + (arb((k + 3)**2)/2 - um)/Ta) + arb(k*k)/(2*Ta)*Ukm
            tot += Vk*P[k]
        h32 = (1 - r)**3*tot
        hd = Act.H_ang(J, b, T)
        pr(f'(iv) J={J} b={b} Delta={Delta}: §32 form = {h32.str(15)}   direct = {hd.str(15)}   diff = {(h32 - hd).str(3)}')
    open(os.path.join(FC.CACHE, f'fp_tail_{Dg}_K{K}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
