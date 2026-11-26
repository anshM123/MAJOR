"""tailcert/tc_power_dv.py [DV]: exact S(beta) = sum_k U_k (beta)_k/k! at a finite state vs the Theorem-P bound K_osc (sanity check)."""
import mpmath as mp, math, sys
import tc_power as P
a, b = 998, 0; s = 2*b + 1
n = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
N = n + a; T = N*(N + s) + b*(b + 1); C = (n + b)*(n + b + 1)
w = mp.mpf(T - C)/T
L = -mp.log(w)
for tau in (10, 100):
    mp.mp.dps = int(tau*0.7 + 40)
    beta = mp.mpc(0.25, tau)
    U = mp.mpf(1); poch = mp.mpf(1); tot = mp.mpc(0)
    for k in range(n + 1):
        tot += U*poch
        poch *= (beta + k)/(k + 1)
        U *= mp.mpf((n - k)*(n + s + k))/((N - k)*(N + s + k))
    Sh = w**beta*tot
    bnd, _ = P.K_osc(0.25, tau, tau + 0.01, float(L))
    print(f'n={n} a={a} L={float(L):.3f} tau={tau}: |Ŝ| = {mp.nstr(abs(Sh), 8)}   Theorem-P bound (L0 = L) = {bnd:.4g}', flush=True)
