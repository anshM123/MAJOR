"""
Feasibility of the <TTTT> partial waves by numerical Taylor-mode application of the D-words (taylor4) -- coordinator step 2.
At a heavy case (beta = 1, m ~ 40 by default) and a few points of the moment grid of e2e/moments.py:
  * order-8 (P,Q)-jets of the 25 seed blocks g^{a,b} (cblock12: copy of e2e/cblock.c with larger jet arrays; double-double radial sums);
  * all 196 word pairs by Taylor mode, in float64 and in long double (80-bit), contracted with CHL's exact helicity structures
    (word coefficients of X_h from t8_kappa_exact.decompose) on both vertices: W[h_L, h_R] for h in {same, opp};
  * timing per pair and per point; precision: float64 vs long double, and the cancellation factor sum|terms|/|sum|.
Usage: python3 t11_feas_taylor.py [beta m npoints]     Output: out_t11_feas_taylor.txt
"""
import os, sys, time, ctypes
import numpy as np, mpmath as mp, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
import blocks
from blocks import CBlock
from moments import make_grid, frame_data, frame_positions
from cases import cases_fixed_beta
from moments4 import pol_grid4
from spinblock2 import invariants4
from embed2 import shifts
from structures2 import words
from taylor4 import Engine4, seed, pq_jets_M, word_budget, ZN, SN
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def load_lib12(prec):
    key = prec + '12'
    if key not in blocks._libs:
        L = ctypes.CDLL(os.path.join(HERE, 'libcblock12_%s.so' % prec))
        L0 = blocks._lib(prec)
        for name in ('block_size', 'block_build', 'block_free', 'block_get', 'block_eval'):
            f, f0 = getattr(L, name), getattr(L0, name)
            f.restype = f0.restype
            if f0.argtypes is not None: f.argtypes = f0.argtypes
        blocks._libs[key] = L
    return key

def helicity_words(Dv, Jv):
    from t8_kappa_exact import decompose
    res = {}
    for kind in ('same', 'opp'):
        ws, sol, params = decompose(2, sp.Rational(Dv), sp.Integer(Jv), kind)
        sol = sol.subs({p: 0 for p in params})
        v = {w: complex(sp.N(sol[i], 30)) for i, w in enumerate(ws)}
        sc = max(abs(x) for x in v.values())
        res[kind] = {w: x/sc for w, x in v.items()}
    return res

if __name__ == "__main__":
    beta = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    mm = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
    npt = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    (Dv, Jv), = cases_fixed_beta(beta, [mm])
    l = 2; Dx = 3
    pr(f"# Taylor-mode <TTTT> feasibility (fourpt/tttt/t11_feas_taylor.py): beta={beta}, m~{mm}: Delta={Dv:.6f}, J={Jv}")
    t = time.time(); hw = helicity_words(Dv, Jv)
    pr(f"  exact helicity word vectors (CHL X_h decomposed): {time.time() - t:.0f} s")
    N = int(min(8000, max(4000, 400*Dv)))
    key = load_lib12('q')
    ws = words(l)
    blk = {}
    for wL in ws:
        for wR in ws:
            s1, s2 = shifts(wL); s3, s4 = shifts(wR)
            a = -((Dx + s1) - (Dx + s2))/2; b = ((Dx + s3) - (Dx + s4))/2
            if (a, b) not in blk:
                c1, c2 = -a, b
                blk[(a, b)] = (CBlock(Dv, Jv, c1, c2, N, prec=key), 0.0) if c1 + c2 >= 0 else (CBlock(Dv, Jv, -c1, -c2, N, prec=key), c1 + c2)
    pr(f"  {len(blk)} seed blocks g^(a,b), N = {N}")
    g = make_grid(Dv, Jv)
    ok = g['ok']
    U, S = g['U'][ok], g['S'][ok]
    # points spread over the grid (bulk of the Regge region)
    cand = [i for i in range(U.size) if 0.05 < U[i] < 0.6 and 0.2 < S[i] < 0.8]
    sel = [cand[int(k*(len(cand) - 1)/max(1, npt - 1))] for k in range(npt)]
    fd = frame_data(U[sel], S[sel])
    x1, x2, x3, x4 = frame_positions(fd['ru'], fd['xi'])
    z1, z3, z2, z4 = pol_grid4(fd['xi'], l)
    s, z = invariants4({1: x1, 2: x2, 3: x3, 4: x4}, {1: z1, 2: z2, 3: z3, 4: z4})
    mp.mp.dps = 40
    for p in range(len(sel)):
        pr(f"  point {p}: U = {U[sel[p]]:.4f}, S = {S[sel[p]]:.4f} (r = {fd['rs'][p]:.4f}, chi = {fd['chi'][p]:.4f})")
        t = time.time()
        jets = {}
        for ab, (bl, vpow) in blk.items():
            Sq = bl.series(np.array([fd['rs'][p]]), np.array([fd['chi'][p]]), sig=+1, kmax=8, lmax=8, cut=70.0)[0]
            jets[ab] = pq_jets_M(Sq, mp.mpf(fd['rs'][p]), mp.mpf(fd['A0'][p]), mp.mpf(fd['E0'][p]), 8, vpow=vpow, one=mp.mpf(1))
        pr(f"    order-8 jets of {len(blk)} seeds (double-double radial sums, mp composition): {time.time() - t:.1f} s")
        res = {}
        for dtn, dt in (('float64', np.float64), ('longdouble', np.longdouble)):
            sv = {n: dt(s[n][p]) for n in SN}
            t = time.time(); tmax = (0, None)
            Phi = {}
            for wL in ws:
                for wR in ws:
                    s1, s2 = shifts(wL); s3, s4 = shifts(wR)
                    d = (Dx + s1, Dx + s2, Dx + s3, Dx + s4)
                    M = word_budget(wL, wR)
                    a = -(d[0] - d[1])/2; b = (d[2] - d[3])/2
                    G = [[dt(str(jets[(a, b)][i, j])) if i + j <= M else dt(0) for j in range(M + 1)] for i in range(M + 1)]
                    t1 = time.time()
                    E = Engine4(M, sv)
                    Phi[(wL, wR)] = E.apply_words(seed(M, sv, tuple(dt(x) for x in d), G, dtype=dt, one=dt(1)), wL, wR)
                    if time.time() - t1 > tmax[0]: tmax = (time.time() - t1, (wL, wR))
            ttot = time.time() - t
            # contract with polarisations: monomial values at this point, [nmono] x [5^4]
            zp = {n: z[n][p].reshape(-1).astype(dt) for n in ZN}
            Wv = {}
            canc = {}
            for hL in ('same', 'opp'):
                for hR in ('same', 'opp'):
                    tre = np.zeros(zp['z12'].shape, dtype=dt); tim = np.zeros(zp['z12'].shape, dtype=dt); absum = np.zeros(zp['z12'].shape)
                    for (wL, wR), dct in Phi.items():
                        c = hw[hL][wL]*hw[hR][wR]
                        if c == 0: continue
                        val = np.zeros(zp['z12'].shape, dtype=dt)
                        for mono, coef in dct.items():
                            mv = np.ones_like(val)
                            for n, e in zip(ZN, mono):
                                if e: mv = mv*zp[n]**e
                            val = val + coef*mv
                        tre += dt(c.real)*val; tim += dt(c.imag)*val; absum += np.abs(c)*np.abs(val.astype(float))
                    tot = tre.astype(float) + 1j*tim.astype(float)
                    Wv[(hL, hR)] = (tre, tim); canc[(hL, hR)] = float(np.max(absum)/np.max(np.abs(tot)))
            res[dtn] = Wv
            pr(f"    {dtn}: 196 word pairs in {ttot:.1f} s (slowest {tmax[1]} {tmax[0]:.2f} s); cancellation sum|terms|/|W| per helicity pair: "
               + ", ".join(f"{k[0]}-{k[1]}: {v:.1e}" for k, v in canc.items()))
        for k in res['float64']:
            a = res['float64'][k][0] + 1j*res['float64'][k][1]
            b = (res['longdouble'][k][0] - res['float64'][k][0]).astype(float) + 1j*(res['longdouble'][k][1] - res['float64'][k][1]).astype(float)
            dev = np.max(np.abs(b))/np.max(np.abs(a)); b = a
            pr(f"    W[{k[0]},{k[1]}]: max|W| = {np.max(np.abs(b)):.3e}, float64 vs long double rel. deviation {dev:.1e}")
    open(os.path.join(HERE, 'out_t11_feas_taylor.txt'), 'w').write("\n".join(out) + "\n")
