"""
a1p_run (A1PRIME.md section 5): term-wise s-channel functional Omega_k (LIGHTSIDE_DIST section 1) on the new sigma-pieces of a1p_pieces,
with LIGHTSIDE's engine imported read-only (resultq/lsd_core.Setup = H matrices and the separable kernel K on the swapped contour C x C;
resultq/lsd_pos.config = the angle/polarisation sets P1-P3).  Nothing in resultq/ is modified; the Mellin amplitudes come from a1p_pieces.

  profile <k> <tag> [N] : delta13-integrated left residues r[p, class][n] (double-trace poles d12 = -n, n = 0..N) and the s-channel
                          single-trace residue (d12 = 1) of xsN, for both angle orders (a,b), (b,a)  -> cache/a1p_prof_k<k>_<tag>.pkl
  report  <k> <tag>     : objects (XN = xsN + xtN, c6, xu6, xu4 and their top-Regge representatives), partial sums N = 50/100/150,
                          power-law tail, sym = E(a,b) + (-1)^k E(b,a) summed over the six free Wick classes -> cache/a1p_report_k<k>_<tag>.txt
    cd ads4/code/a1p; (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 a1p_run.py profile 4 P1)
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.join(HERE, '..', 'resultq'))
import numpy as np
from scipy.special import loggamma, digamma
import lsd_core as C            # read-only
import lsd_pos as P             # read-only
import a1p_pieces as AP

CACHE = os.path.join(HERE, 'cache')
GRID = dict(s0=1e-9, m=16, q=1.8)          # LIGHTSIDE's production grid
PCS = {2: ('xu4', 'r4u', 'xs3', 'xt3', 'rX3'),
       3: ('xs3', 'xt3', 'rX3'),
       4: ('xs3', 'xt3', 'xs4', 'xt4', 'rX4', 'c6', 'r6', 'xu6', 'r6u'),
       5: ('xs3', 'xt3', 'xs4', 'xt4', 'c6', 'xu6')}
OBJ = {'X3': ('xs3', 'xt3'), 'rX3': ('rX3',), 'X4': ('xs4', 'xt4'), 'rX4': ('rX4',), 'c6': ('c6',), 'r6': ('r6',),
       'xu6': ('xu6',), 'r6u': ('r6u',), 'xu4': ('xu4',), 'r4u': ('r4u',)}
NEFF = {'X3': 2, 'rX3': 2, 'X4': 4, 'rX4': 4, 'c6': 4, 'r6': 4, 'xu6': 4, 'r6u': 4, 'xu4': 2, 'r4u': 2}


def profile(S, pieces, N=150, Y13=9.0, h13=0.1):
    y13 = np.arange(-Y13, Y13 + h13/2, h13)
    nf = np.arange(N + 1).astype(float)
    out = {}
    for p in pieces:
        for c in S.H:
            out[(p, c)] = np.zeros(N + 1, complex); out[(p, c, 'st')] = 0j
    for c2 in sorted(set(AP.CONT[p][1] for p in pieces)):
        ps = [p for p in pieces if AP.CONT[p][1] == c2]
        for yy in y13:
            d13 = c2 + 1j*yy
            w13 = np.exp(2*loggamma(d13))*h13/(2*np.pi)
            dL = -nf.astype(complex)
            Kd = S.K(np.concatenate([dL, [1.0 + 0j]]), d13, deriv=True)
            gL = np.exp(2*loggamma(2 - d13 + nf) - 2*loggamma(nf + 1)); psi14 = digamma(2 - d13 + nf); psin = digamma(nf + 1)
            for p in ps:
                ML = AP.M_amp(p, dL, d13); dML = AP.dM_amp(p, dL, d13)
                st = AP.st_res(p, d13)
                for c, (K, dK) in Kd.items():
                    KL, dKL = K[:N + 1], dK[:N + 1]
                    out[(p, c)] += w13*gL*((-2*psi14*ML + dML)*KL + ML*dKL + 2*psin*ML*KL)
                    if st is not None:
                        out[(p, c, 'st')] += w13*np.exp(2*loggamma(1 - d13))*st*K[N + 1]
    return out


def tail(r, n0=60):
    N = len(r) - 1
    n = np.arange(n0, N + 1)
    pp = -np.polyfit(np.log(n), np.log(np.abs(r[n0:]) + 1e-300), 1)[0]
    A = np.sum(r[n0:]*n**pp)/len(n)
    if pp <= 1.05:
        return pp, np.nan
    return pp, A*(N + 0.5)**(1 - pp)/(pp - 1)


def sec_profile(k, tag, N=150):
    a, b, p2, p4 = P.config(tag)
    out = {'k': k, 'tag': tag, 'grid': GRID, 'N': N}
    for order, (aa, bb) in (('ab', (a, b)), ('ba', (b, a))):
        t0 = time.time()
        S = C.Setup(k, aa, bb, p2, p4, **GRID)
        out[order] = profile(S, PCS[k], N=N)
        print(f'k = {k} {tag} {order}: {len(S.t)} nodes, {time.time() - t0:.0f}s', flush=True)
    with open(os.path.join(CACHE, f'a1p_prof_k{k}_{tag}.pkl'), 'wb') as fh:
        pickle.dump(out, fh)


def sec_report(k, tag):
    L = []
    pr = lambda s: (print(s, flush=True), L.append(s))
    D = pickle.load(open(os.path.join(CACHE, f'a1p_prof_k{k}_{tag}.pkl'), 'rb'))
    sg = 1 if k % 2 == 0 else -1
    pr(f'a1p_run report: k = {k}, tag {tag}, grid {D["grid"]}, N = {D["N"]}; sym = E(a,b) + (-1)^k E(b,a), six free classes; '
       f'"st" = s-channel single-trace residue included term-wise')
    vals = {}
    for o, pcs in OBJ.items():
        if not all(p in PCS[k] for p in pcs):
            continue
        tot = {}
        for order in ('ab', 'ba'):
            R = D[order]
            rs = sum(R[(p, c)] for p in pcs for c in C.CLS)
            st = sum(R[(p, c, 'st')] for p in pcs for c in C.CLS)
            part = {Nc: np.sum(rs[:Nc + 1]) + st for Nc in (50, 100, 150)}
            pp, tl = tail(rs)
            tot[order] = (part, pp, tl, st)
        s = {Nc: tot['ab'][0][Nc] + sg*tot['ba'][0][Nc] for Nc in (50, 100, 150)}
        tl = tot['ab'][2] + sg*tot['ba'][2]
        scale = max(abs(tot['ab'][0][150]), abs(tot['ba'][0][150]), abs(tot['ab'][3]), abs(tot['ba'][3]))
        ext = s[150] + (tl if np.isfinite(tl) else 0)
        vals[o] = (s[150], ext, scale)
        regime = 'n < k (drop-out predicted)' if NEFF[o] < k else ('n = k (borderline)' if NEFF[o] == k else 'n > k')
        pr(f'  {o:4s} [n_eff {NEFF[o]}, {regime}]: (ab) {tot["ab"][0][150].real:+.6e} (ba) {tot["ba"][0][150].real:+.6e}'
           f' | st(ab) {tot["ab"][3].real:+.3e} | sym N=50 {s[50].real:+.5e}, N=100 {s[100].real:+.5e}, N=150 {s[150].real:+.5e},'
           f' tail {tl.real:+.2e} (p {tot["ab"][1]:.2f}/{tot["ba"][1]:.2f}) -> extrap {ext.real:+.6e}; |extrap|/scale {abs(ext)/scale:.1e}')
    # robust extrapolation: S(N) = S_inf + C N^-q fitted to the symmetrised partial sums at N = 60..150 (q free)
    from scipy.optimize import curve_fit
    def sym_partial(pcs):
        acc = {}
        for order in ('ab', 'ba'):
            R = D[order]
            rs = sum(R[(p, c)] for p in pcs for c in C.CLS); st = sum(R[(p, c, 'st')] for p in pcs for c in C.CLS)
            acc[order] = np.cumsum(rs) + st
        return (acc['ab'] + sg*acc['ba']).real
    Ns = np.arange(60, D['N'] + 1)
    def extrap(S):
        y = S[Ns]
        try:
            f = lambda n, a, b, q: a + b*n**(-q)
            popt, _ = curve_fit(f, Ns.astype(float), y, p0=(y[-1], (y[0] - y[-1])*60.0, 1.0), maxfev=20000)
            return popt[0], popt[2]
        except Exception:
            return np.nan, np.nan
    pr('  3-parameter extrapolation of the symmetrised partial sums (S_inf + C N^-q, N = 60..150):')
    for o, pcs in OBJ.items():
        if not all(p in PCS[k] for p in pcs):
            continue
        S = sym_partial(pcs); a, q = extrap(S)
        pr(f'    {o:4s}: S(75) {S[75]:+.6e}, S(125) {S[125]:+.6e}, S(150) {S[150]:+.6e} -> S_inf {a:+.6e} (q {q:.2f})')
    pr('  clean-range estimate (residues of degree-6 pieces are grid-noise dominated beyond n ~ 110): S(100) + power tail fitted on n = 40..100')
    for o, pcs in OBJ.items():
        if not all(p in PCS[k] for p in pcs):
            continue
        rs = sum((D['ab'][(p, c)] + sg*D['ba'][(p, c)]) for p in pcs for c in C.CLS).real
        st = sum((D['ab'][(p, c, 'st')] + sg*D['ba'][(p, c, 'st')]) for p in pcs for c in C.CLS).real
        n = np.arange(40, 101); pp = -np.polyfit(np.log(n), np.log(np.abs(rs[40:101]) + 1e-300), 1)[0]
        A = np.mean(rs[40:101]*n**pp); tl = A*100.5**(1 - pp)/(pp - 1) if pp > 1.05 else np.nan
        pr(f'    {o:4s}: S(100) {np.sum(rs[:101]) + st:+.6e}, residue decay p = {pp:.2f}, tail {tl:+.2e} -> {np.sum(rs[:101]) + st + tl:+.6e}')
    for o, r in (('X3', 'rX3'), ('X4', 'rX4'), ('c6', 'r6'), ('xu6', 'r6u'), ('xu4', 'r4u')):
        if o in vals and r in vals and NEFF[o] == k:
            Sd = sym_partial(OBJ[o]) - sym_partial(OBJ[r])
            pr(f'  difference {o} - {r} (n < k remainder): partial sums N = 50/100/150 {Sd[50]:+.3e} / {Sd[100]:+.3e} / {Sd[150]:+.3e} '
               f'(relative to |{r}| = {abs(sym_partial(OBJ[r])[150]):.4e})')
    for o, r in (('X3', 'rX3'), ('X4', 'rX4'), ('c6', 'r6'), ('xu6', 'r6u'), ('xu4', 'r4u')):
        if o in vals and r in vals and NEFF[o] == k:
            pr(f'  borderline (P): {o} vs its top-Regge term {r}: sym N=150 {vals[o][0].real:+.8e} vs {vals[r][0].real:+.8e}, '
               f'rel. diff {abs(vals[o][0] - vals[r][0])/abs(vals[r][0]):.1e}; extrapolated rel. diff '
               f'{abs(vals[o][1] - vals[r][1])/abs(vals[r][1]):.1e}')
    with open(os.path.join(CACHE, f'a1p_report_k{k}_{tag}.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')


if __name__ == '__main__':
    sec, k, tag = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    if sec == 'profile':
        sec_profile(k, tag, int(sys.argv[4]) if len(sys.argv) > 4 else 150)
    else:
        sec_report(k, tag)
