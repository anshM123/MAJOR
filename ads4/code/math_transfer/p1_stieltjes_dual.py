"""P1, class B (completely monotone relaxation functions = Stieltjes transforms, signed, TV <= 1, support [gamma, inf)).

Dual certificate (upper bound).  For any coefficients c_1..c_N,
   |F(0)| = |int dsigma/lam| <= sup_{lam>=gamma} |1/lam - sum_n c_n/(lam+n)| * ||sigma||_TV + eps * sum_n |c_n|,
because sum_n c_n F(n) = int sum_n c_n/(lam+n) dsigma and |F(n)| <= eps.  We pick c by a 120-digit least-squares fit
in u = gamma/lam (not the optimal minimax, so the bound is conservative) and evaluate the sup on a 4000-point grid in u
(grid sup, not interval arithmetic).  Then W_up(eps) = min_N [e_N + eps*S_N].

Output: out_p1_stieltjes_dual.txt.  Columns: N, e_N, S_N, -ln e_N / sqrt(gamma N), ln S_N / N; then W_up(eps) for
eps = 1e-4 ... 1e-40 with -ln W_up / sqrt(gamma L) and W_up * L^{2 gamma} (the latter would be constant for the
Theorem-C / class-A rate).
"""
import mpmath as mp

mp.mp.dps = 120


def certificate(gamma, N, M=400, G=4000):
    us = [mp.sin(mp.pi * (k + mp.mpf(1) / 2) / (2 * M)) ** 2 for k in range(M)]
    A = mp.matrix(M, N)
    b = mp.matrix(M, 1)
    for i, u in enumerate(us):
        lam = gamma / u
        b[i] = 1
        for n in range(1, N + 1):
            A[i, n - 1] = lam / (lam + n)          # fit lam*r(lam) ~ 1
    c = mp.lu_solve(A.T * A, A.T * b)
    err = mp.mpf(0)
    for k in range(G):
        lam = gamma * G / mp.mpf(k + 1)
        r = sum(c[n - 1] / (lam + n) for n in range(1, N + 1))
        err = max(err, abs(1 / lam - r))
    return err, sum(abs(x) for x in c)


if __name__ == "__main__":
    lines = []
    for gamma in (mp.mpf(1) / 2, mp.mpf(1)):
        lines.append(f"gamma = {gamma}")
        tab = []
        for N in (2, 4, 6, 8, 12, 16, 20, 24, 32, 40):
            e, S = certificate(gamma, N)
            tab.append((N, e, S))
            lines.append(f"  N={N:3d}  e_N={mp.nstr(e, 4):>10}  S_N={mp.nstr(S, 4):>10}  "
                         f"-ln e/sqrt(gN)={mp.nstr(-mp.log(e) / mp.sqrt(gamma * N), 4):>6}  lnS/N={mp.nstr(mp.log(S) / N, 4)}")
        for k in (4, 8, 12, 16, 20, 30, 40):
            eps = mp.mpf(10) ** (-k)
            L = mp.log(1 / eps)
            W, Nb = min((e + eps * S, N) for N, e, S in tab)
            lines.append(f"  eps=1e-{k:<3d} W_up={mp.nstr(W, 4):>10} (N={Nb:2d})  -lnW/sqrt(gL)={mp.nstr(-mp.log(W) / mp.sqrt(gamma * L), 4)}"
                         f"   W*L^(2g)={mp.nstr(W * L ** (2 * gamma), 4)}")
        lines.append("")
    txt = "\n".join(lines)
    print(txt)
    open(__file__.replace("p1_stieltjes_dual.py", "out_p1_stieltjes_dual.txt"), "w").write(txt + "\n")
