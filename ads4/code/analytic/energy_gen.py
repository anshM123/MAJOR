"""
energy_gen.py: Theorem E (ANALYTIC_SLIVER §6) for general s = 2b+1: quadratic-energy bound for the dual-Hahn spin recurrence
  S_{n+1} = (h_n + alpha e_n) S_n + g_n S_{n-1},   S_0 = 1,
  e_n = (2n+s)(2n+s+1)/((N+n+s)(n+s)(N-n)),  g_n = n(2n+s+1)/((n+s)(2n+s-1)),  h_n = 1 - g_n.
Energy E_n = |S_n|^2 + g_n|S_{n-1}|^2 - 2 y_n Im(S_n conj S_{n-1}),  y_n = tau e_n/2,  u_n = (h_n + sigma e_n)/2,  eta_n = |y_n|/sqrt(g_n).
  E'_{n+1} <= E_n [g_n + (4u_n^2 + 2 sqrt(g_n)|u_n|)/(1-eta_n)],   E_{n+1} <= E'_{n+1}/(1 - delta_{n+1}),
  delta_{n+1} = [(g_{n+1}-g_n)/g_{n+1} + |y_{n+1}-y_n|/sqrt(g_{n+1})]/(1 - eta_{n+1}),   |S_n|^2 <= E_n/(1 - eta_n).
log_bound(n, N, s, sigma, tau): float log of the bound on |d_n^alpha S_n| (d_n = 1 - n(n+s)/(N(N+s))); math in floats [DV-level evaluation;
the certified version (§6.3) uses the closed-form majorant].
"""
import math
def coefs(n, N, s):
    e = (2*n + s)*(2*n + s + 1)/((N + n + s)*(n + s)*(N - n))
    g = 1.0 if s == 1 else n*(2*n + s + 1)/((n + s)*(2*n + s - 1))
    return e, g, 1.0 - g


def log_bound(nf, N, s, sig, tau):
    """log of the Theorem-E bound on |d_nf^alpha S_nf(alpha)|, alpha = sig + i tau (sig >= 0, tau >= 0)."""
    if s == 1:
        logE = math.log(2.0); n0 = 0
    else:
        e0, g0, h0 = coefs(0, N, s); e1, g1, h1 = coefs(1, N, s)
        eta1 = tau*e1/2/math.sqrt(g1); assert eta1 < 1
        logE = math.log((1 + eta1)*((h0 + sig*e0)**2 + (tau*e0)**2 + g1)); n0 = 1
    for n in range(n0, nf):
        e, g, h = coefs(n, N, s); e2, g2, h2 = coefs(n + 1, N, s)
        y = tau*e/2; u = (h + sig*e)/2; sg = math.sqrt(g); eta = y/sg
        assert eta < 1, (n, eta)
        dy = abs(tau*e2/2 - y); dg = g2 - g
        fac = (g + (4*u*u + 2*sg*abs(u))/(1 - eta))*(1 + dy/(sg*(1 - eta))) + dg/(1 - eta)
        logE += math.log(fac)
    ef, gf, hf = coefs(nf, N, s)
    etaf = tau*ef/2/math.sqrt(gf); assert etaf < 1
    d = 1 - nf*(nf + s)/(N*(N + s))
    return 0.5*(logE - math.log(1 - etaf) - math.log(gf if s > 1 else 1.0)*0) + sig*math.log(d)


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '.')
    from flint import arb, acb, ctx
    ctx.prec = 6000
    def exact(nf, N, s, al):
        Nb = arb(N); C0 = nf*(nf+s); D = Nb*(Nb+s); U = arb(1); t = acb(1); tot = acb(1)
        for k in range(1, nf+1):
            U = U*(C0 - (k-1)*(k-1+s))/(D - (k-1)*(k-1+s)); t = t*(al + k - 1)/k; tot += U*t
        d = 1 - arb(C0)/D
        return float(abs(d**al*tot).mid())
    for s in (1, 9):
        for nf, a in [(1001, 998), (1001, 1500), (3000, 998), (5000, 3000), (20000, 1200)]:
            for sig in (0.25, 1.25, 2.5, 3.5):
                for tau in (10, 300, 700, 900):
                    lb = log_bound(nf, nf + a, s, sig, tau)
                    ex = exact(nf, nf + a, s, acb(sig, tau)) if nf <= 5000 else float('nan')
                    print(f's={s} n={nf} a={a} alpha={sig}+{tau}i: exact |S^| = {ex:.5g}   Theorem-E bound = {math.exp(lb):.4g}', flush=True)
