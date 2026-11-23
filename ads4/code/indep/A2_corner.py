"""Target A (ii): corner power counting.
(1) de-shift survival with the exact (all orders in rho) b-shifted blocks of cblocks.py (own recursion):
    S_b(beta) = eta^n int_0^30 w(nu) Ghat_b(r, eta; nu) dnu  vs  S_LR = eta^n int w P_l, at fixed rho = r e^beta,
    n = 1 (k = 2 weights), n = 3 (k = 4 weights); local exponent d log|S_b - S_LR| / d beta over beta = 8.5 .. 13.5.
(2) edge structure Ghat_b - 2 Re[c_l e^{l beta} F_b(rho)] ~ e^{-2 beta}, F_b = (z/4rho)^h 2F1(h, h+b; 2h; z), h = (1/2 + i nu)/2.
(3) residue formula for the k = 4 weight with a pole at 3i/2:  S_b - S_LR ~ eta^3 Re[2 pi i Res_{3i/2} w c_l e^{l beta} F_b].
(4) the corner count c < a + 1 - J0 re-derived (Stage-3 measure) and the threshold table.
"""
import os, sys, time, math
import numpy as np
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cblocks

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'A2_corner.txt')
lines = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

W = {'w2g': (1, lambda x: x**2 * np.exp(-x**2 / 9)),
     'w2s': (1, lambda x: x**2 / np.cosh(np.pi * x / 2)),
     'w4g': (3, lambda x: x**2 * (x**2 + 1) * (x**2 + 4) * np.exp(-x**2 / 9)),
     'w4p': (3, lambda x: x**2 * (x**2 + 1) * (x**2 + 4) * np.exp(-x**2 / 9) / (x**2 + 9 / 4))}

def c_l(l): return mp.gamma(l + 0.5) / (mp.sqrt(mp.pi) * mp.gamma(l + 1))

def F_b(rho, nu, b):
    h = (mp.mpf(1) / 2 + 1j * nu) / 2
    z = 4 * rho / (1 + rho)**2
    return (z / (4 * rho))**h * mp.hyp2f1(h, h + b, 2 * h, z)

def residue_pred(rho, beta, b, n=3):
    mp.mp.dps = 40
    nu0 = mp.mpc(0, 1.5)
    # w4p = w4g / ((nu - 3i/2)(nu + 3i/2)) -> Res = w4g(3i/2)/(3i); c_l F_b finite at 3i/2 (limit)
    w4g = lambda x: x**2 * (x**2 + 1) * (x**2 + 4) * mp.exp(-x**2 / 9)
    eps = mp.mpf('1e-25')
    vals = []
    for e in (eps, -eps, 1j * eps):
        nu = nu0 + e
        l = -mp.mpf(1) / 2 + 1j * nu
        vals.append(c_l(l) * F_b(rho, nu, b))
    cF = sum(vals) / 3
    l0 = -2
    res = w4g(nu0) / (3j) * cF * mp.exp(l0 * beta)
    return float(mp.cosh(beta)**n * mp.re(2j * mp.pi * res)), complex(cF)

def main():
    t0 = time.time()
    xg, wg = np.polynomial.legendre.leggauss(400)
    nus = 15 * (xg + 1); wq = 15 * wg
    rhos = (0.2, 0.5, 0.85)
    betas = np.arange(8.5, 13.51, 0.5)
    bs = (0, 1, 2, 3, 4)
    NMAX = 200
    Bv = {b: cblocks.recursion_vec(nus, b, NMAX) for b in bs}
    say(f'blocks built ({time.time()-t0:.0f} s); GL 400 nodes on [0,30]; N_max = {NMAX}')

    # (2) edge structure at moderate beta
    say('\n(2) edge structure: |Ghat_b - 2Re[c_l e^{l beta} F_b]| / |Ghat_b| at nu = 0.7, 3.1 (claim: falls like e^{-2 beta})')
    for nu in (0.7, 3.1):
        Bs = {b: cblocks.recursion(nu, b, NMAX) for b in bs}
        for rho in (0.3, 0.85):
            row = []
            for beta in (2.0, 4.0, 6.0, 8.0):
                r = rho * math.exp(-beta); eta = math.cosh(beta)
                Q = cblocks.scaled_legendre_vec([nu], r, eta, NMAX + 2)
                rr = []
                for b in bs:
                    G = cblocks.Ghat_vec(Bs[b], Q, r)[0]
                    mp.mp.dps = 30
                    l = mp.mpc(-0.5, nu)
                    edge = 2 * mp.re(c_l(l) * mp.exp(l * beta) * F_b(rho, nu, b))
                    rr.append(abs(G - float(edge)) / abs(G))
                row.append(max(rr))
            say(f'   nu={nu} rho={rho}: beta = 2,4,6,8 -> max_b rel. dev. ' + '  '.join(f'{x:.1e}' for x in row)
                + '   (e^{-2beta}: ' + '  '.join(f'{math.exp(-2*bb):.1e}' for bb in (2, 4, 6, 8)) + ')')
        # sign convention: F_{-b} instead of F_b
    nu = 0.7; rho = 0.85; beta = 6.0; r = rho * math.exp(-beta)
    Q = cblocks.scaled_legendre_vec([nu], r, math.cosh(beta), NMAX + 2)
    for b in (2, 4):
        G = cblocks.Ghat_vec(cblocks.recursion(nu, b, NMAX), Q, r)[0]
        l = mp.mpc(-0.5, nu)
        e1 = 2 * mp.re(c_l(l) * mp.exp(l * beta) * F_b(rho, nu, b)); e2 = 2 * mp.re(c_l(l) * mp.exp(l * beta) * F_b(rho, nu, -b))
        say(f'   control (rho=0.85, beta=6, nu=0.7, b={b}): with F_b {abs(G-float(e1))/abs(G):.1e}; with F_(-b) {abs(G-float(e2))/abs(G):.1e}')

    # (1) de-shift survival
    say('\n(1) local exponents d log|S_b - S_LR|/d beta (beta 8.5 -> 13.5); rows: weight, rho; columns: b = 0..4, then LR-only rate d log|S_LR|/d beta')
    S = {}
    for rho in rhos:
        for beta in betas:
            r = rho * math.exp(-beta); eta = math.cosh(beta)
            Q = cblocks.scaled_legendre_vec(nus, r, eta, NMAX + 2)
            Pl = Q[0]
            G = {b: cblocks.Ghat_vec(Bv[b], Q, r) for b in bs}
            for wn, (n, wf) in W.items():
                wv = wf(nus) * wq
                SLR = eta**n * np.sum(wv * Pl)
                for b in bs:
                    Sb = eta**n * np.sum(wv * G[b])
                    S[(wn, rho, beta, b)] = (Sb, SLR)
        say(f'   rho={rho} done ({time.time()-t0:.0f} s)')
    for wn in W:
        for rho in rhos:
            cols = []
            for b in bs:
                d = [abs(S[(wn, rho, be, b)][0] - S[(wn, rho, be, b)][1]) for be in betas]
                sl = [math.log(d[k + 1] / d[k]) / (betas[k + 1] - betas[k]) for k in range(len(betas) - 1)]
                cols.append(sl)
            lr = [abs(S[(wn, rho, be, 0)][1]) for be in betas]
            lrs = [math.log(lr[k + 1] / lr[k]) / 0.5 for k in range(len(betas) - 1)]
            # report early window (8.5-10) and late window (12-13.5)
            def fmt(sl): return f'{np.mean(sl[:3]):+.2f}/{np.mean(sl[-3:]):+.2f}'
            say(f'   {wn} rho={rho}: ' + '  '.join(f'b={b}: {fmt(c)}' for b, c in zip(bs, cols)) + f'   | LR {fmt(lrs)}')
    say('   (each entry: mean local exponent over beta in [8.5,10] / [12,13.5])')
    # absolute size at the floor for w4g
    for rho in rhos:
        be = 13.5
        Sb, SLR = S[('w4g', rho, be, 0)]
        say(f'   w4g rho={rho} beta=13.5: |S_0 - S_LR| = {abs(Sb-SLR):.2e}, eta^3 x 1e-16 x int|w| ~ {math.cosh(be)**3*1e-16*np.sum(np.abs(W["w4g"][1](nus)*wq))*math.exp(-be/2):.1e}')

    # (3) residue formula for w4p at beta = 10
    say('\n(3) w4p: (S_b - S_LR) / (eta^3 Re[2 pi i Res_{3i/2} w c_l e^{l beta} F_b]) at beta = 10 (LR has no residue: c_l(3i/2) = 0)')
    for rho in rhos:
        beta = 10.0
        row = []
        for b in bs:
            Sb, SLR = S[('w4p', rho, beta, b)]
            pred, cF = residue_pred(mp.mpf(rho), mp.mpf(beta), b)
            row.append((Sb - SLR).real / pred)
        say(f'   rho={rho}: ' + '  '.join(f'b={b}: {x:.7f}' for b, x in zip(bs, row)))
    say('   same residue for w4g (no pole): identically 0 (w4g regular at 3i/2, c_l F_b finite)')

    # (4) corner count
    say('\n(4) corner count (Stage-3 measure): kernel |K| <~ r^a eta^c on {r <= c0/(2 eta)}, K1\' |dDisc| <= B r^(1-J0):')
    say('    int^oo d eta eta^c int_0^{c0/2eta} r^(a+1-J0) dr ~ int d eta eta^(c - a - 2 + J0)  converges iff  c < a + 1 - J0.')
    rows = [('k=2  C\'_EE (de-shifted)', 0, -0.5, 'LRC*: exact-LR exponent -1/2 (table (1), w2g)'),
            ('k=3  C_3 (not de-shifted)', 0, -0.5, 'P_{-1/2} ~ eta^{-1/2} log eta'),
            ('k=4  C4^ds, f_pos(3i/2)=0', 2, 0.5, 's=1 shell raw eta^{1/2} (bound; measured -1/2, table (1), w4g)'),
            ('k=4  C4^ds, f_pos(3i/2)!=0', 2, 1.0, 'residue eta^{+1} (table (1), w4p; formula (3))'),
            ('k=5  C5 raw', 2, 1.5, 'P_{3/2} ~ eta^{3/2}')]
    for name, a, c, why in rows:
        say(f'    {name:30s} a={a} c={c:+.1f} -> J0 < {a + 1 - c:.1f}   [{why}]')
    # numerical illustration of the count: int_1^E deta eta^c int_0^{1/(2eta)} r^(a+1-J0) dr at J0 = threshold +- 0.1
    say('    numerical check of the count (corner integral up to eta = E, closed form in r):')
    for name, a, c, why in rows:
        J0s = (a + 1 - c - 0.1, a + 1 - c + 0.1)
        out = []
        for J0 in J0s:
            p = a + 2 - J0
            f = lambda E: mp.quad(lambda e: e**c * (1 / (2 * e))**p / p, [1, E])
            out.append((float(f(1e4)), float(f(1e8))))
        say(f'      {name:30s} J0 = thr-0.1: I(1e4) = {out[0][0]:.3e}, I(1e8) = {out[0][1]:.3e} (converges);  '
            f'thr+0.1: {out[1][0]:.3e}, {out[1][1]:.3e} (diverges)')
    say(f'\n(elapsed {time.time()-t0:.0f} s)')
    open(OUT, 'w').write('\n'.join(lines) + '\n')

if __name__ == '__main__':
    main()
