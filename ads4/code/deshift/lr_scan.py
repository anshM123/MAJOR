"""deshift/lr_scan.py -- DESHIFT_CLOSE.md §2.4: Stage 5's leading-Regge (LR) large-beta positivity scan (the HYB region beta >~ 15)
re-done with the REALISED images of the new continuum weights (ab_c4v_h015, ab_c4d_h015, ab_c5_p015), with rigorous worst-sign envelopes.

For each new weight the de-shift residual obeys, for ALL t >= 0 (ab_check.py):
    |W - Phi| <= eps_W,     and  |W(t) - Phi(t)| <= c0 + m1 t^2/2   (c0 = |int E_even|, m1 = int |E_even| mu^2),
so the realised Mehler images are the exact ones with W -> Phi(1 + e(t)), |e| <= e_max(t):
    value, C5 (Phi = 1):       e_max = eps_W
    derivative (Phi = -t^2/2): e_max = min(2 eps_W/t^2, 2 c0/t^2 + m1).
The scan is stage7_s5check.scan (Stage 5's states, normalisation Nrm, both delta), imported READ-ONLY from resultq/; the forward pieces
enter as B4 -> 2u^2 P0v, dB4 -> -2u^2 Lam^2 Dd + c4 u^3 P0d, B5 -> sigma u^3 P05 with P0v = sum W(1+ev), Dd = 1/2 sum W t^2 (1+ed),
P0d = sum W (1+ed), P05 = sum W (1+e5).  All 8 sign patterns (ev, ed, e5) = (+-e_max) are scanned and the minimum kept; the 'exact'
row (e = 0) is the Stage-5 continuum.  beta grid extended from 30 (Stage 5/7) to 45.
Output: cache/lr_scan.txt
"""
import os, sys, re, itertools, time
import numpy as np
import kern
sys.path.insert(0, os.path.join(kern.CODE, 'resultq'))
sys.dont_write_bytecode = True
import stage12_forward_regge as FR       # read-only imports (no files written by these modules at import or in the calls used)
import stage5_positivity as S5P
import stage7_s5check as S7

OUT = []


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def bounds(tag):
    txt = open(os.path.join(kern.CACHE, f'ab_check_{tag}.txt')).read()
    eps = float(re.search(r'eps_W = int_0\^110 \|E_even\| = ([0-9.e+-]+)', txt).group(1))
    mom = [float(x) for x in re.search(r'n = 0\.\.8: (.*)', txt).group(1).split(', ')]
    c0 = abs(float(re.search(r't =   0\.0: .*W - Phi = ([0-9.e+-]+)', txt).group(1)))
    return eps, c0, mom[1]


def main():
    t0 = time.time()
    ev_, c0v, _ = bounds('c4v_h015'); ed_, c0d, m1d = bounds('c4d_h015'); e5_, c05, _ = bounds('c5_p015')
    pr(f"# lr_scan: worst-sign realised images of the continuum weights; eps_W (value, deriv, C5) = {ev_:.2e}, {ed_:.2e}, {e5_:.2e}; "
       f"deriv c0 = {c0d:.1e}, m1 = {m1d:.1e}")
    us = np.concatenate([np.geomspace(1e-6, 0.05, 40), np.linspace(0.06, 1.0, 48)])
    for Lam in (100., 300., 1000.):
        betas = np.unique(np.concatenate([np.geomspace(5/Lam, 1.0, 70), np.linspace(1.0, 12.0, 221), np.linspace(12.0, 45.0, 133)]))
        for name in ('Q0-full-robust', 'Q3-full-delta0.2'):
            Fn, c2n, c3n = S5P.load_func(name)
            Tn = [S5P.CosT(c2n), S5P.CosT(np.convolve(c2n, [0, 0, 1])), S5P.CosT(c3n)]
            base = {}
            for beta in betas:
                n = int(min(3.0e5, 800 + 5*Lam*beta)); n = max(800, (n//20)*20)
                t, W = FR.mehler_nodes(beta, n)
                om = Lam*t
                base[beta] = ((W*Tn[0](om)).sum(), (W*Tn[1](om)).sum(), (W*Tn[2](om)).sum(), t, W)
            res = {}
            for sv, sd, s5 in [(0, 0, 0)] + list(itertools.product((1, -1), repeat=3)):
                BQ = {}
                for beta in betas:
                    I20, I22, I3, t, W = base[beta]
                    tt = np.maximum(t, 1e-300)
                    ev = sv*ev_; e5 = s5*e5_
                    ed = sd*np.minimum(2*ed_/tt**2, 2*c0d/tt**2 + m1d)
                    P0v = (W*(1 + ev)).sum(); Dd = 0.5*(W*t**2*(1 + ed)).sum(); P0d = (W*(1 + ed)).sum(); P05 = (W*(1 + e5)).sum()
                    BQ[(Lam, round(beta, 12))] = (I20, I22, I3, (P0v, Dd, P0d, P05))
                res[(sv, sd, s5)] = S7.scan(Lam, name, betas, us, BQ, None)
            for k in ('b>=20', 'beta>=6', 'beta>=12', 'first_neg'):
                ex = res[(0, 0, 0)][k]
                worst = min((r[k] for kk, r in res.items() if kk != (0, 0, 0)), key=lambda v: v[0])
                pr(f"   Lam = {Lam:5.0f} {name:17s} min {k:9s}: exact {S5P.fmt(ex) if ex[1] else 'n/a'};  worst-sign realised "
                   f"{S5P.fmt(worst) if worst[1] else 'n/a'};  shift {worst[0] - ex[0]:+.1e}" if ex[1] else
                   f"   Lam = {Lam:5.0f} {name:17s} {k}: none (no negative value anywhere)")
            pr(f"   [{time.time() - t0:.0f}s]")
    open(os.path.join(kern.CACHE, 'lr_scan.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
