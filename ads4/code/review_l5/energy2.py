"""
energy2.py [DV-level evaluation of a rigorous inequality]: Theorem E' (ANALYTIC_SLIVER §6.2), physical-branch (eigen-coordinate) tracking
for the b = 0 dual-Hahn spin recurrence  S_{n+1} = 2 w_n S_n + S_{n-1},  w_n = eps_n alpha/2:
  c_n = P_n^{-1}(S_n, S_{n-1}),  P_n = [[lam+, lam-],[1,1]],  lam± = w ± r, r = sqrt(1+w^2),
  |c_{n+1}|_inf <= |lam+(n)| (1 + (|Δlam+| + |Δlam-|)/(2|r_{n+1}|)) |c_n|_inf,   |S_n| <= (|lam+(n)|+|lam-(n)|)|c_n|_inf  (here n = n_f, using P_{n_f}).
"""
import cmath, math
def log_bound2(nf, N, sig, tau):
    al = complex(sig, tau)
    def lam(n):
        e = 2*(2*n + 1)/((N - n)*(N + n + 1)); w = e*al/2; r = cmath.sqrt(1 + w*w)
        if r.real < 0: r = -r
        return w + r, w - r, r
    lp, lm, r = lam(0)
    c0 = max(abs(1 - lm), abs(lp - 1))/(2*abs(r))
    lg = math.log(c0)
    for n in range(nf):
        lp, lm, r = lam(n); lp2, lm2, r2 = lam(n + 1)
        if n < nf - 1:
            lg += math.log(abs(lp)) + math.log1p((abs(lp2 - lp) + abs(lm2 - lm))/(2*abs(r2)))
        else:
            lg += math.log(abs(lp))
    # final: S_nf = (P_{nf} c_nf)_1 ... we tracked c in the basis P_{nf-1} after the last step: S_nf = (P_{nf-1} Lambda c)_1 -> |S| <= (|lp|+|lm|)|c|
    lp, lm, r = lam(nf - 1)
    lg += math.log(abs(lp) + abs(lm))   # c' = P_{nf-1}^{-1} v_nf = Lambda c (|c'| <= |lp||c|, already counted); |S_nf| <= (|lp|+|lm|)|c'|
    d = 1 - nf*(nf + 1)/(N*(N + 1))
    return lg + sig*math.log(d)
if __name__ == '__main__':
    import sys
    sys.path.insert(0, '.')
    from energy_gen import log_bound
    for nf, a in [(1001, 998), (3000, 998), (100000, 998), (5000, 3000), (20000, 1200), (1000000, 998)]:
        for sig in (0.25, 2.5, 3.5):
            for tau in (10, 300, 700, 900, 990):
                if tau >= a: continue
                b2 = math.exp(log_bound2(nf, nf + a, sig, tau))
                b1 = math.exp(log_bound(nf, nf + a, 1, sig, tau)) if nf <= 100000 else float('nan')
                print(f'n={nf} a={a} alpha={sig}+{tau}i: E-bound {b1:.4g}   E\'-bound {b2:.4g}', flush=True)
