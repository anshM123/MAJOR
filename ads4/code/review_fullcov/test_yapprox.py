"""PROOF.md §2 (y direction): exact residues of F_K(c, y H(c)) in y, at fixed c.

For a cell (default cell_5305_6041, K=3193) and c in {l, mid, u}:
 1. compute all signed pole coefficients Res_j = (1-c)^3 sum_{k>j} (A_k + T_k/p_j) rho_{j,k}
    in 2600-bit fixed point (independent of the archived code);
 2. verify that F_K(c,yH) - sum_j Res_j/(1-p_j H y) is affine in y (second differences ~ 0);
 3. compare |Res_j| with the per-pole bound used by full_preflight.residue_groups;
 4. compare the exact truncation majorant sum_j |Res_j| 2 r_j^(m+1)/(g_j(1-r_j)) with the
    archived/recomputed Ey(m) at several m, including the certified degree;
 5. compare DCT Chebyshev coefficients of y -> F_K(c,yH(c)) with the residue prediction and
    the true tail sum_{n>m}|a_n| with Ey(m).
"""
import sys, os, time
from fractions import Fraction as Fr
import mpmath as mp
import rc
os.environ.setdefault('SLIVER_COEFF', rc.COEFF)
sys.path.insert(0, str(rc.ROOT / 'sliver_close'))
from full_preflight import residue_degree, residue_groups, lower_exp, maxa, maxt, LOG2_UP, sqrtlo

stem = sys.argv[1] if len(sys.argv) > 1 else 'cell_5305_6041'
z, l, u, Hl, Hu = rc.cell(stem); K = z['K']; my = z['y']['deg']
a, g = rc.coeffs()
W = 2600
mp.mp.dps = 900

def residues(c, b):
    """Floating big-int: rho = r / 2^E with >= 2800 significant bits kept (relative precision
    is required because rho_{j,k} grows by up to e^1000 in k); sum accumulated in fixed 2^-W."""
    s = 2 * b + 1; C, D = c.numerator, c.denominator
    tn = [rc.tcoef_num(k, b) for k in range(K + 1)]
    out = [0] * K
    SIG = 2800
    def mul(r, E, n, d):
        # d > 0 always; keep >= SIG significant bits after the division
        r = r * n
        need = SIG + d.bit_length() - abs(r).bit_length()
        if need > 0:
            r <<= need; E += need
        return r // d, E
    for j in range(1, K):
        pj = j * (j + s)
        r, E = mul(1, 0, C - D, D)
        for i in range(j - 1, -1, -1):          # descending: avoid tiny intermediates
            pi = i * (i + s)
            r, E = mul(r, E, C * pj - D * pi, D * (pj - pi))
        acc = 0
        for k in range(j + 1, K + 1):
            t = r * (2 * pj * a[k] + tn[k])
            acc += (t >> (E - W)) if E >= W else (t << (W - E))
            pk = k * (k + s)
            r, E = mul(r, E, D * pk - C * pj, D * (pk - pj))
        out[j] = mp.mpf(acc) / (2 * pj) / mp.mpf(10) ** 60 / mp.mpf(2) ** W * (mp.mpf((1 - c).numerator) / (1 - c).denominator) ** 3
    return out

def per_pole_bound(j, H):
    # the bound behind residue_groups for the 16-group containing j (without the tail factor)
    a0 = 1 + 16 * ((j - 1) // 16)
    L = min(lower_exp(a0, K, u, s) for s in (1, 9)); e = int(L / LOG2_UP)
    M = Fr(maxa[a0 + 1], 10 ** 60) + Fr(maxt[a0 + 1], 10 ** 60 * a0 * (a0 + 1))
    return Fr(2 * (K - a0), a0) * l ** -8 * M / 2 ** e

report = []
for c in (l, (l + u) / 2, u):
    H = rc.Hsec(c, l, u, Hl, Hu); Hm = mp.mpf(H.numerator) / H.denominator
    for b in (0, 4):
        s = 2 * b + 1
        t0 = time.time()
        import pickle; cf = rc.Path(__file__).resolve().parent / f'.res_{stem}_{c.numerator}_{c.denominator}_{b}.pkl'
        if cf.exists(): R = [mp.mpf(x) for x in pickle.load(open(cf, 'rb'))]
        else:
            R = residues(c, b); pickle.dump([mp.nstr(x, 880) for x in R], open(cf, 'wb'))
        t1 = time.time()
        # 2. affine check
        def pf(y):
            return mp.fsum(R[j] / (1 - j * (j + s) * Hm * y) for j in range(1, K))
        ys = [Fr(0), Fr(1, 3), Fr(2, 3), Fr(1)]
        G = [mp.mpf(rc.F_K(c, H * y, K, b, W=400).numerator) / rc.F_K(c, H * y, K, b, W=400).denominator - pf(mp.mpf(y.numerator) / y.denominator) for y in ys]
        d2a = G[0] - 2 * G[1] + G[2]; d2b = G[1] - 2 * G[2] + G[3]
        maxres = max(abs(x) for x in R[1:])
        # 3. per-pole bound ratio
        ratio = max(abs(R[j]) / (lambda q: mp.mpf(q.numerator) / q.denominator)(per_pole_bound(j, H)) for j in range(1, K) if R[j] != 0)
        # 4. exact truncation majorant vs Ey(m)
        rows = []
        for m in (2, 5, 10, 20, 50, 100, my):
            tot = mp.mpf(0)
            for j in range(1, K):
                aj = j * (j + s) * Hm; gj = mp.sqrt(1 - aj); rj = (1 - gj) / (1 + gj)
                tot += abs(R[j]) * 2 * rj ** (m + 1) / (gj * (1 - rj))
            Ey = mp.mpf(residue_degree(l, u, Hl, K, m)[1])
            rows.append((m, tot, Ey))
            assert tot <= Ey, ('Ey bound violated', c, b, m)
        print(f'{stem} c={float(c):.5f} b={b}: residues {t1-t0:.0f}s  max|Res_j|={mp.nstr(maxres,4)}  2nd-diff of F-PF: {mp.nstr(d2a,3)}, {mp.nstr(d2b,3)}  '
              f'max |Res_j|/bound_j={mp.nstr(ratio,4)}', flush=True)
        for m, tot, Ey in rows:
            print(f'    m={m:4d}  exact majorant={mp.nstr(tot,4):>12s}  Ey(m)={mp.nstr(Ey,4):>12s}', flush=True)
        assert abs(d2a) < mp.mpf(10) ** -60 and abs(d2b) < mp.mpf(10) ** -60
        assert ratio <= 1
        # 5. DCT check (only at mid c, b=0 and b=4)
        if c == (l + u) / 2:
            N = 40; mp.mp.dps = 80
            ts = [mp.cos(mp.pi * (i + mp.mpf(1) / 2) / N) for i in range(N)]
            vals = []
            for t in ts:
                y = Fr(mp.nstr((1 + t) / 2, 70))
                f = rc.F_K(c, H * y, K, b, W=400); vals.append(mp.mpf(f.numerator) / f.denominator)
            coef = [2 * mp.fsum(vals[i] * mp.cos(n * mp.pi * (i + mp.mpf(1) / 2) / N) for i in range(N)) / N for n in range(N)]
            mp.mp.dps = 900
            pred = lambda n: mp.fsum(R[j] * 2 / mp.sqrt(1 - j * (j + s) * Hm) * ((1 - mp.sqrt(1 - j * (j + s) * Hm)) / (1 + mp.sqrt(1 - j * (j + s) * Hm))) ** n for j in range(1, K))
            for n in (2, 3, 5, 8, 12):
                print(f'    DCT a_{n}={mp.nstr(coef[n],8)}  residue prediction={mp.nstr(pred(n),8)}', flush=True)
            for m in (1, 2, 4, 8):
                tail = mp.fsum(abs(x) for x in coef[m + 1:N - 8])
                print(f'    true Chebyshev tail beyond m={m}: ~{mp.nstr(tail,4)}  Ey(m)={mp.nstr(mp.mpf(residue_degree(l, u, Hl, K, m)[1]),4)}', flush=True)
print('test_yapprox PASS', stem)
