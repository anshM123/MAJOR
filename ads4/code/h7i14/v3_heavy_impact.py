"""v3 (H7_I14_LIGHTSIDE.md §20): heavy-side impact of replacing the scalar-g de-shift of the graviton k = 2 functional by the spinning one
(C^_2' of §19), in the existing near-gap chains (read-only imports of spinfunc/l9_analysis, l9_neargap, n1_analysis; nothing written there).

The only change: in l9_analysis.shift for k = 2 the local de-shift coefficients g/(g+g~) -> (1 + u_T)/2, u_T = (h- - h+)/(h- + h+) of the
computed spinning Regge kernel (v3_regge_kernel); after Re only sigma_T = Im u_T matters.  k = 3 (C^_3) keeps the scalar g (its own shift
relation is not derived; unchanged on purpose, so the comparison isolates (HN) for C^_2).
sigma_T(nu): computed on the fly (cached in cache/v3_sigmaT_memo.pkl) for nu <= 300, asymptotic sigma_sc + C3/nu^3 above.
Runs:  (a) l9 Laplace rows, Delta_g = 100, 300, 1000 (variant C23H_K2), baseline and spinning;
       (b) Delta_g = 100 with the direct (n1) rows on x >= 0.2 (as n1_analysis / H0 §24.5), baseline and spinning.
Output: cache/v3_heavy_impact.txt
"""
import os, sys, io, pickle, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
sys.dont_write_bytecode = True
import numpy as np
import v3_regge_kernel as V
import l9_analysis as A
import l9_neargap as NG
import n1_analysis as N1
import l2_scalar_validate as L2

OUT = os.path.join(HERE, 'cache', 'v3_heavy_tmp'); os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True)
NG.HERE = OUT
MEMO_F = os.path.join(HERE, 'cache', 'v3_sigmaT_memo.pkl')
MEMO = pickle.load(open(MEMO_F, 'rb')) if os.path.exists(MEMO_F) else {}
C3 = None


def sigma_T(nu):
    global C3
    key = round(float(nu), 6)
    if key in MEMO:
        return MEMO[key]
    if nu > 300:
        if C3 is None:
            tab = pickle.load(open(os.path.join(HERE, 'cache', 'v3_uT_table.pkl'), 'rb'))
            sel = [c for c in tab if c['nu'] >= 60]
            C3 = np.mean([c['nu']**3*c['ds'] for c in sel[-3:]])
        gp, gm = L2.gfac(nu)
        val = ((gp - gm)/(gp + gm)).imag + C3/nu**3
    else:
        bs = np.linspace(0.12, min(3.0, max(0.6, 30/nu)), 14)
        (hm, h0, hp), res = V.fit3(nu, bs, [V.heavy_image(nu, b, 1, -1) for b in bs])
        val = ((hm - hp)/(hm + hp)).imag
    MEMO[key] = val
    return val


_orig_shift = A.shift


def shift_spin(J, h, k, xs, m):
    if k != 2:
        return _orig_shift(J, h, k, xs, m)
    f = lambda t: A.env(h, t, k)*A.Dc(h, J, t)
    out = []
    for x in xs:
        s = sigma_T(x*m)
        gp, gm = 1 + 1j*s, 1 - 1j*s
        out.append(((gp*f(x + 1j/m) + gm*f(x - 1j/m))/(gp + gm) - f(x)).real)
    return np.array(out)


buf = io.StringIO()
def pr(s=""):
    print(s, flush=True); buf.write(str(s) + "\n")


def run(label, Dgs, direct):
    origgrid = A.grid
    if direct:
        N1.DIRECT = True
        NG.A.grid = lambda J, K, kinds=('C2', 'C3H'), **kw: N1.grid_direct(J, K, kinds, **kw)
    try:
        NG.main('C23H_K2', list(Dgs))
    finally:
        NG.A.grid = origgrid
    txt = open(os.path.join(OUT, 'cache', 'l9_neargap_C23H_K2_%s.txt' % '_'.join(map(str, Dgs)))).read()
    pr(f"## {label}"); pr(txt)


if __name__ == '__main__':
    t0 = time.time()
    for spin in (False, True):
        A.shift = shift_spin if spin else _orig_shift
        lab = 'SPINNING de-shift (C^_2\')' if spin else 'baseline (scalar g, as l9/n1)'
        run(f"{lab}; l9 Laplace rows", (100, 300, 1000), direct=False)
        run(f"{lab}; Delta_g = 100 with direct rows (n1)", (100,), direct=True)
        pickle.dump(MEMO, open(MEMO_F, 'wb'))
    pr("sigma_T vs sigma_sc at the k = 2 nodes used (nu, sigma_T, sigma_sc):")
    for key in sorted(MEMO)[::max(1, len(MEMO)//25)]:
        gp, gm = L2.gfac(key)
        pr(f"   {key:9.3f}  {MEMO[key]:+.5f}  {((gp - gm)/(gp + gm)).imag:+.5f}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'v3_heavy_impact.txt'), 'w').write(buf.getvalue())
