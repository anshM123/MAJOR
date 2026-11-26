"""
xcheck_ladder/xl_extra.py: (a) precision stability / cancellation size of the ladder sum; (b) the (Delta, J) = (1100, 100) window state;
(c) minimum over masses per spin (twist edge .. 5 Lam) at Dg = 1000, against the published minima of Table tab:Dstates.
   python3 xl_extra.py prec | win | scan JMAX
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
from flint import arb, ctx
import xl_common as X
import xl_norm as N

LAM = X.lam_int(1000)


def terms_max(J, b, T, Z, Y, Q, lam):
    """largest |term| of the ladder sum (before the (1-r)^3 factor), as a measure of the cancellation."""
    C = J*(J + 1); r = arb(C)/T; n = J - b; sig = 1 if b == 0 else -1
    U = arb(1); mx = 0.0
    for k in range(n + 1):
        V = arb(Fr(4*k, k + 3).numerator)/(k + 3) if b == 0 else 4 + arb(12)/(k + 1)
        t = math.comb(k + 3, 3)**2*U*(Z[k]*(1/arb(k + 3) - r/(k + 1)) + (2*Y[k] + V*Z[k] + sig*arb(lam)*arb(Q[k]))/T)
        mx = max(mx, abs(float(t.mid())))
        if k < n:
            uk = (b + k)*(b + k + 1); U = U*(C - uk)/(T - uk)
    return mx


def prec():
    out = []
    for J, b, D in [(5, 4, 1005), (18, 4, 1018), (29, 4, 1029), (32, 0, 1032), (32, 4, 1032), (100, 0, 1100), (100, 4, 1100),
                    (300, 4, 1300)]:
        T = Fr((D - 1)*(D - 2))
        res = {}
        for p in (64, 128, 256, 512, 1024, 2048):
            ctx.prec = p
            Z, Y = X.load_cert_moments(1000, 1000, kmax=J + 2)
            Q = X.Q_exact(J + 2, LAM)
            h = X.H_ang(J, b, X.to_arb(T), Z, Y, Q, LAM)
            res[p] = h
        ctx.prec = 2048
        Z, Y = X.load_cert_moments(1000, 1000, kmax=J + 2); Q = X.Q_exact(J + 2, LAM)
        mx = terms_max(J, b, X.to_arb(T), Z, Y, Q, LAM)
        ref = res[2048]
        line = f'(Delta,J,b)=({D},{J},{b}): H = {ref.str(15, radius=False)}  max|term| = {mx:.3e} (ratio {mx/abs(float(ref.mid())):.2e});  ' + \
            '  '.join(f'p={p}: {res[p].str(8)}' for p in (64, 128, 256, 512, 1024))
        print(line, flush=True); out.append(line)
    open(os.path.join(X.OUT, 'precision.txt'), 'w').write('\n'.join(out) + '\n')


def win():
    """(1100, 100): exact action per helicity, hybrid, N_s; quadrature-moment cross-check."""
    ctx.prec = 2048
    J = 100; D = 1100; T = Fr((D - 1)*(D - 2)); Ta = X.to_arb(T)
    Z, Y = X.load_cert_moments(1000, 1000, kmax=J + 2); Q = X.Q_exact(J + 2, LAM)
    d = json.load(open(os.path.join(X.OUT, 'pmom_cont_1000_n24_p1400.json')))
    S = [sum((arb(Sr[j]) for Sr in d['S']), arb(0)) for j in range(len(d['S'][0]))]
    Zq, Yq = X.moments_from_power(S, J + 2)
    mfl = mp.nstr(mp.sqrt(mp.mpf(T.numerator - J*(J + 1)))/mp.sqrt(LAM), 30)
    out = [f'(Delta, J) = (1100, 100): m/Lam = {mfl[:12]}, bhat = {(2*J+1)/float(mfl):.4f}, T = {T}']
    for b, h in ((0, '++'), (4, '+-')):
        H = X.H_ang(J, b, Ta, Z, Y, Q, LAM); Hq = X.H_ang(J, b, Ta, Zq, Yq, Q, LAM)
        hy = N.hybrid(1000, J, mfl, h, npn=320)
        Ns = float(hy['Ns'])
        m2 = float(T - J*(J + 1))
        line = (f'  {h}: H_ang = {H.str(12, radius=False)} (quad-moment route {Hq.str(12, radius=False)}),  N_s = {Ns:.10e},  '
                f'H/N_s = {float(H.mid())/Ns:+.6f},  H_hyb/N_s = {float(hy["Hhyb"])/Ns:+.6f},  (H - H_hyb)/N_s = {(float(H.mid())-float(hy["Hhyb"]))/Ns:+.3e},'
                f'  relative to H_hyb: {(float(H.mid())/float(hy["Hhyb"]) - 1):+.3e}  (x m^2 = {(float(H.mid())/float(hy["Hhyb"]) - 1)*m2:+.1f})')
        print(line, flush=True); out.append(line)
    open(os.path.join(X.OUT, 'window_1100_100.txt'), 'w').write('\n'.join(out) + '\n')


def scan(Jmax, J0=0, step=1):
    sys.path.insert(0, os.path.join(X.ADS4, 'code', 'largespin'))
    import ls_core as hf
    ctx.prec = 1024
    Z, Y = X.load_cert_moments(1000, 1000, kmax=Jmax + 2); Q = X.Q_exact(Jmax + 2, LAM)
    out = ['# Dg=1000: min over m/Lam in [twist edge, 5] of lambda/N_s = min_b H_ang,b/N_s,b  (N_s from ls_core.hyb_fast; '
           'grid: 30 log-spaced points on [tw, 1.3] + 15 on [1.3, 5], then 18 golden-section steps around the minimum)']
    for J in range(J0, Jmax + 1, step):
        if J in (1, 3):
            continue
        bs = [(0, '++'), (4, '+-')] if J % 2 == 0 else [(4, '+-')]
        if J < 4:
            bs = [(0, '++')]
        Ctw = Fr((J + 999)*(J + 998)) - J*(J + 1)
        mtw = math.sqrt(float(Ctw)/LAM)

        def lam_of(mf):
            if mf == mtw:
                T = Fr((J + 999)*(J + 998))
            else:
                T = Fr(mf)**2*LAM + J*(J + 1)
            vals = []
            for b, h in bs:
                H = float(X.H_ang(J, b, X.to_arb(T), Z, Y, Q, LAM).mid())
                vals.append((H/hf.hyb_fast(1000, J, mf, h)['Ns'], h))
            return min(vals)
        grid = [mtw] + [mtw*(1.3/mtw)**(i/30) for i in range(1, 31)] + [1.3*(5/1.3)**(i/15) for i in range(1, 16)]
        vals = [(lam_of(m), m) for m in grid]
        (best, h), mb = min(vals)
        i = [v[1] for v in vals].index(mb)
        if 0 < i < len(grid) - 1:
            a, c = grid[i - 1], grid[i + 1]
            for _ in range(18):
                x1 = c - (c - a)*0.618; x2 = a + (c - a)*0.618
                if lam_of(x1)[0] < lam_of(x2)[0]:
                    c = x2
                else:
                    a = x1
            mb = (a + c)/2; best, h = lam_of(mb)
        line = f'J={J:3d}  min lambda/N_s = {best:+.5f} ({h}) at m/Lam = {mb:.5f}   (twist edge {mtw:.5f}: {lam_of(mtw)[0]:+.5f})'
        print(line, flush=True); out.append(line)
    open(os.path.join(X.OUT, f'minscan_1000_J{J0}-{Jmax}_s{step}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__' and sys.argv[1] != 'kernel':
    {'prec': prec, 'win': win, 'scan': lambda: scan(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))}[sys.argv[1]]()


def kernel():
    """route B literally: per-nu exact finite-T kernel K_T(nu) (Arb), H^(2) = int g K_T; reports int|g K_T|/|H^(2)| (the
    conditioning of the quadrature route) and max|K_T| below and above the state's mass."""
    import numpy as np
    ctx.prec = 1024
    Lam = math.sqrt(LAM)
    W = X.load_cont(1000)
    nu, w = X.cont_nodes(Lam, 24)
    g = X.g_of_chi(X.chi_re(W, nu), nu, Lam)*X.theta(nu, Lam)
    out = []
    for J, b, D in [(5, 4, 1005), (18, 4, 1018), (29, 4, 1029), (32, 0, 1032), (100, 4, 1100)]:
        T = X.to_arb(Fr((D - 1)*(D - 2))); C = J*(J + 1); r = arb(C)/T; n = J - b
        coefZ, coefY = [], []
        U = arb(1)
        for k in range(n + 1):
            V = arb(4*k)/(k + 3) if b == 0 else 4 + arb(12)/(k + 1)
            Dk2 = math.comb(k + 3, 3)**2
            coefZ.append((1 - r)**3*Dk2*U*(1/arb(k + 3) - r/(k + 1) + V/T)); coefY.append((1 - r)**3*Dk2*U*2/T)
            if k < n:
                uk = (b + k)*(b + k + 1); U = U*(C - uk)/(T - uk)
        m = math.sqrt(float(Fr((D - 1)*(D - 2)) - C))
        Kv = np.empty(len(nu))
        for i, v in enumerate(nu):
            s = arb(float(v))**2 + arb(1)/4
            f = [arb(0), arb(1)]
            for mm in range(1, n + 7):
                f.append((-(s - 2*mm*mm)*f[mm] - (mm - 1)**2*f[mm - 1])/(mm + 1)**2)
            tot = arb(0)
            for k in range(n + 1):
                mm = k + 2
                z = 12*f[mm] - (mm - 1)*f[mm - 1] + (mm + 1)*f[mm + 1]
                y = (arb(12 - mm*mm)/2)*f[mm] + (arb(13 - 7*mm)/2)*f[mm - 1] + (arb(13 + 7*mm)/2)*f[mm + 1] \
                    + (arb((mm - 2)**2)/4)*f[mm - 2] + (arb((mm + 2)**2)/4)*f[mm + 2]
                tot += coefZ[k]*z + coefY[k]*y
            Kv[i] = float(tot.mid())
        H2 = float(np.sum(w*g*Kv)); A = float(np.sum(np.abs(w*g*Kv)))
        below = np.abs(Kv[nu <= m]).max(); above = np.abs(Kv[nu > m]).max() if np.any(nu > m) else 0
        tailm = float(np.sum((w*g*Kv)[nu > m]))
        line = (f'(Delta,J,b)=({D},{J},{b}) m={m:.3f}: H2 (route B) = {H2:.12e}; int|gK| = {A:.4e} (ratio {A/abs(H2):.2e}); '
                f'max|K| nu<=m: {below:.3e}, nu>m: {above:.3e} at nu_max={nu.max():.1f}; contribution of nu>m: {tailm:.2e}')
        print(line, flush=True); out.append(line)
    open(os.path.join(X.OUT, 'kernel_routeB.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__' and sys.argv[1] == 'kernel':
    kernel()
