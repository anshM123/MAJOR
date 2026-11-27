"""P1, sector class: F analytic and |F| <= M in S_alpha = {|arg(s+gamma)| < alpha}, pi/2 <= alpha <= pi.
zeta = (s+gamma)^p, p = pi/(2 alpha), maps S_alpha onto the right half plane; Theorem C's Blaschke-Lagrange proof then
applies verbatim in the zeta plane with nodes zeta_n = (n+gamma)^p and target zeta_0 = gamma^p:
   |F(0)| <= eps*Lam_N + (M + eps*Lam_N) * beta_N,
   beta_N = prod (zeta_n - zeta_0)/(zeta_n + zeta_0) ~ exp(-2 gamma^p N^{1-p}/(1-p)),
   Lam_N  = sum_k prod_{j != k} (zeta_k + zeta_j)/|zeta_k - zeta_j| ~ e^{c_p N}.
Hence |F(0)| <= M exp(-C (log M/eps)^{1-p}).  p = 1 (half plane) reproduces Theorem C's log rate.
This script tabulates (1/N) log Lam_N (-> c_p) and -log beta_N / N^{1-p}, and the optimised bound for gamma = 1, M = 1.
Output: out_p1_sector.txt
"""
import mpmath as mp

mp.mp.dps = 60
out = []
g = mp.mpf(1)
for p in (mp.mpf(1), mp.mpf("0.75"), mp.mpf("0.6"), mp.mpf("0.5")):
    rows = {}
    for N in range(1, 121):
        z = [mp.mpf(n + g) ** p for n in range(1, N + 1)]
        z0 = g ** p
        beta = mp.fprod((zz - z0) / (zz + z0) for zz in z)
        Lam = sum(mp.fprod((z[k] + z[j]) / abs(z[k] - z[j]) for j in range(N) if j != k) for k in range(N))
        rows[N] = (beta, Lam)
    out.append(f"p = {mp.nstr(p, 3)} (alpha = {mp.nstr(180 / (2 * p), 4)} deg)")
    for N in (10, 40, 120):
        beta, Lam = rows[N]
        out.append(f"   N={N:3d}  log(Lam)/N = {mp.nstr(mp.log(Lam) / N, 4)}   -log(beta)/N^(1-p) = "
                   f"{mp.nstr(-mp.log(beta) / mp.mpf(N) ** (1 - p), 4)}")
    for k in (6, 12, 24, 48):
        eps = mp.mpf(10) ** (-k)
        W = min(eps * Lam + (1 + eps * Lam) * beta for beta, Lam in rows.values())
        L = mp.log(1 / eps)
        out.append(f"   eps=1e-{k:<2d} bound={mp.nstr(W, 4):>10}   -log(bound)/L^(1-p) = {mp.nstr(-mp.log(W) / L ** (1 - p), 4)}")
txt = "\n".join(out)
print(txt)
open(__file__.replace("p1_sector_rate.py", "out_p1_sector.txt"), "w").write(txt + "\n")
