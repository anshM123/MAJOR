"""PROOF.md §2 (c direction): at fixed real y, the Chebyshev coefficients of
c -> F_K(c, y H(c)) on [l,u] (independent evaluator, DCT of N first-kind samples)
versus the Mobius-disk / Bernstein-ellipse majorant  sum_k min_rho 2 M_k(rho) rho^-m/(rho-1),
recomputed here in floating point (test only) for several m, and the certified E_c at the
archived degree.  Also a Monte-Carlo test of the Mobius 'largest modulus at a real endpoint'
claim on random complex points of the disk."""
import sys, os, math, random
import numpy as np
from scipy.fft import dct
from fractions import Fraction as Fr
import rc
os.environ.setdefault('SLIVER_COEFF', rc.COEFF)
sys.path.insert(0, str(rc.ROOT / 'sliver_close'))
from coupled_plan import A as LA, T as LT, Rhos

def majorant(l, u, Hl, Hu, K, ylo, yhi):
    """log M_k(rho) per rho, as in PROOF §2 (floats; test only)."""
    l, u, Hl, Hu = map(float, (l, u, Hl, Hu))
    c0 = (l + u) / 2; dc = (u - l) / 2; H0 = (Hl + Hu) / 2; dH = (Hu - Hl) / 2
    jj = np.arange(K); out = []
    for rho in Rhos:
        rho = float(rho); R = (rho + 1 / rho) / 2
        pairs = [(c0 + zz * dc, (H0 + zz * dH) * y) for zz in (-R, R) for y in (float(ylo), float(yhi))]
        if max(h for c, h in pairs) * (K - 1) * (K + 8) >= 1: continue
        f = np.zeros(K)
        for s in (1, 9):
            pp = jj * (jj + s)
            for c, h in pairs:
                f = np.maximum(f, np.abs(c - pp * h) / (1 - pp * h))
        HH = max(abs(h) for c, h in pairs); nrm = max(abs(1 - c) for c, h in pairs) ** 3
        logM = np.r_[0., np.cumsum(np.log(np.maximum(f, 1e-300)))] + np.logaddexp(LA[:K + 1], LT[:K + 1] + math.log(max(HH, 1e-300))) + math.log(nrm) + math.log(2 / (rho - 1))
        out.append((math.log(rho), logM))
    def bound(m):
        best = np.min(np.array([lm - m * lr for lr, lm in out]), axis=0)
        mx = best.max(); return math.exp(mx) * np.exp(best - mx).sum() if mx < 700 else float('inf')
    return bound

def cheb_coeffs(stem, y, N, b=0):
    z, l, u, Hl, Hu = rc.cell(stem); K = z['K']
    vals = []
    for i in range(N):
        t = math.cos(math.pi * (i + 0.5) / N)
        x = Fr(t).limit_denominator(10 ** 30) / 2 + Fr(1, 2)
        c = l + (u - l) * x
        f = rc.F_K(c, rc.Hsec(c, l, u, Hl, Hu) * y, K, b)
        vals.append(float(f))
    a = dct(np.array(vals), type=2) / N
    a[0] /= 2
    return a

if __name__ == '__main__':
    random.seed(7)
    tests = [('cell_5305_6041', Fr(0), 1024), ('cell_5305_6041', Fr(1, 2), 1024), ('cell_5305_6041', Fr(1), 1024),
             ('cell_8936_9000', Fr(1), 1024)]
    import json
    for stem, y, N in tests:
        z, l, u, Hl, Hu = rc.cell(stem); K = z['K']
        gfile = rc.MERGE / (stem + '_group.json')
        if gfile.exists():
            groups = json.load(open(gfile))['groups']
            g = ([g for g in groups if Fr(g['lo']) <= y <= Fr(g['hi'])] or [min(groups, key=lambda g: abs(Fr(g['lo']) - y))])[0]; y = min(max(y, Fr(g['lo'])), Fr(g['hi'])); ylo, yhi, mcert, Ecert = Fr(g['lo']), Fr(g['hi']), g['deg'], float(Fr(g['error']))
        else:
            ylo, yhi, mcert, Ecert = Fr(0), Fr(1), z['c']['deg'], float(Fr(z['c']['error']))
        bnd = majorant(l, u, Hl, Hu, K, ylo, yhi)
        a = cheb_coeffs(stem, y, N)
        ab = np.abs(a)
        print(f'{stem} y={y} K={K}: |a_n| at n=0,10,50,100,200,300,400,500: ' + ' '.join(f'{ab[n]:.1e}' for n in (0, 10, 50, 100, 200, 300, 400, 500)), flush=True)
        floor = np.median(ab[-200:])
        for m in (20, 50, 100, 150, 200, 300, 400, mcert):
            tail = ab[m + 1:].sum()
            B = bnd(m)
            print(f'   m={m:4d} true tail={tail:.3e} (noise floor {floor:.1e})  majorant={B:.3e}  {"OK" if tail <= B else "VIOLATION"}', flush=True)
        print(f'   certified E_c at m={mcert}: {Ecert:.3e}', flush=True)

    # Mobius disk test: for random complex c in the disk |c-c0|<=R dc, every ladder factor and
    # |A+h T|-type prefactor is bounded by its value at the two real endpoints.
    z, l, u, Hl, Hu = rc.cell('cell_8936_9000'); K = z['K']
    l, u, Hl, Hu = map(float, (l, u, Hl, Hu)); c0 = (l + u) / 2; dc = (u - l) / 2; H0 = (Hl + Hu) / 2; dH = (Hu - Hl) / 2
    worst = 0
    for rho in (1.5, 3.0, 8.0, 20.0):
        R = (rho + 1 / rho) / 2
        for y in (0.0, 0.3, 1.0):
            for s in (1, 9):
                for j in (0, 1, 50, 3000, K // 2, K - 1):
                    p = j * (j + s)
                    fac = lambda c: (c - p * y * (H0 + (c - c0) / dc * dH)) / (1 - p * y * (H0 + (c - c0) / dc * dH))
                    ends = max(abs(fac(c0 - R * dc)), abs(fac(c0 + R * dc)))
                    if (1 - p * y * (H0 - R * dH)) <= 0 or (1 - p * y * (H0 + R * dH)) <= 0: continue
                    for _ in range(300):
                        w = complex(random.uniform(-1, 1), random.uniform(-1, 1))
                        if abs(w) > 1: continue
                        cc = c0 + R * dc * w
                        worst = max(worst, abs(fac(cc)) / ends)
    print('Mobius test: max_{disk}|factor| / max_{real endpoints}|factor| =', worst, 'OK' if worst <= 1 + 1e-12 else 'VIOLATION')
