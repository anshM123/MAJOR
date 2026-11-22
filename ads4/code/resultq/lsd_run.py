"""
lsd_run (LIGHTSIDE_DIST.md): production runs.

  profile <k> <tag>   : Omega_OPE term by term.  For both orders (a,b), (b,a): the delta13-integrated s-channel residue array r[p, c][n]
                        (n = 0..N, double-trace poles d12 = -n) plus the xs single-trace residue, on the fine graded contour.
                        -> cache/lsd_prof_k<k>_<tag>.pkl
  report <k> <tag>    : partial sums, power-law tail extrapolation, 1<->3 (anti)symmetrisation; comparison with EIN's whole-function
                        values (cache/lsd_pos_<tag>_n60[_k3].pkl).  -> cache/lsd_report_k<k>_<tag>.txt

    cd ads4/code/resultq; (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_run.py profile 3 P1)
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import loggamma, digamma
import lsd_core as C
import lsd_mellin as LM
import lsd_pos as P

CACHE = os.path.join(HERE, 'cache')
PCS = ('c0', 'c2', 'c3', 'c4', 'r4', 'xu', 'xu4', 'r4u')
GRID = dict(s0=1e-9, m=16, q=1.8)
if os.environ.get('LSD_GRID') == 'alt':          # robustness grid
    GRID = dict(s0=1e-8, m=14, q=2.2)
SUF = '' if os.environ.get('LSD_GRID') != 'alt' else '_alt'


def profile(S, pieces, N=150, Y13=8.0, h13=0.2, cont=C.CONT):
    """delta13-integrated residues: out[(p, c)] = complex array (N+1,), out[(p, c, 'st')] = single-trace residue (xs)."""
    y13 = np.arange(-Y13, Y13 + h13/2, h13)
    nf = np.arange(N + 1).astype(float)
    out = {}
    for p in pieces:
        for c in S.H:
            out[(p, c)] = np.zeros(N + 1, complex); out[(p, c, 'st')] = 0j
    for c2 in sorted(set(cont[p][1] for p in pieces)):
        ps = [p for p in pieces if cont[p][1] == c2]
        for yy in y13:
            d13 = c2 + 1j*yy
            w13 = np.exp(2*loggamma(d13))*h13/(2*np.pi)
            dL = -nf.astype(complex)
            Kd = S.K(np.concatenate([dL, [1.0 + 0j]]), d13, deriv=True)
            gL = np.exp(2*loggamma(2 - d13 + nf) - 2*loggamma(nf + 1)); psi14 = digamma(2 - d13 + nf); psin = digamma(nf + 1)
            for p in ps:
                ML = LM.M_amp(p, dL, d13); dML = LM.dM_amp(p, dL, d13)
                for c, (K, dK) in Kd.items():
                    KL, dKL = K[:N + 1], dK[:N + 1]
                    out[(p, c)] += w13*gL*((-2*psi14*ML + dML)*KL + ML*dKL + 2*psin*ML*KL)
                    if p == 'xs':
                        out[(p, c, 'st')] += w13*np.exp(2*loggamma(1 - d13))*((d13 - (1 - d13))**2/4 + 1/12)*K[N + 1]
    return out


def tail(r, n0=60):
    """fit |r_n| ~ A n^-p on n >= n0 (log-linear on the complex values via the dominant real/imag part); returns (p, tail estimate
    sum_{n > N} r_n, rel. residual of the fit)."""
    N = len(r) - 1
    n = np.arange(n0, N + 1)
    x = np.log(n)
    yv = r[n0:]
    # fit each of Re, Im as A n^-p with a common p from |r|
    pp = -np.polyfit(x, np.log(np.abs(yv) + 1e-300), 1)[0]
    A = np.sum(yv*n**pp)/len(n)
    model = A*n**-pp
    resid = np.linalg.norm(yv - model)/np.linalg.norm(yv)
    if pp <= 1.05:
        return pp, np.nan, resid
    t = A*(N + 0.5)**(1 - pp)/(pp - 1)
    return pp, t, resid


def sec_profile(k, tag):
    a, b, p2, p4 = P.config(tag)
    out = {'k': k, 'tag': tag, 'grid': GRID}
    for order, (aa, bb) in (('ab', (a, b)), ('ba', (b, a))):
        t0 = time.time()
        S = C.Setup(k, aa, bb, p2, p4, **GRID)
        out[order] = profile(S, PCS)
        print(f'k = {k} {tag} {order}: {len(S.t)} nodes, {time.time() - t0:.0f}s', flush=True)
    with open(os.path.join(CACHE, f'lsd_prof_k{k}_{tag}{SUF}.pkl'), 'wb') as fh:
        pickle.dump(out, fh)


def sec_report(k, tag):
    L = []
    pr = lambda s: (print(s, flush=True), L.append(s))
    D = pickle.load(open(os.path.join(CACHE, f'lsd_prof_k{k}_{tag}{SUF}.pkl'), 'rb'))
    posf = os.path.join(CACHE, f'lsd_pos_{tag}_n60' + ('_k3' if k == 3 else '') + '.pkl')
    pos = pickle.load(open(posf, 'rb')) if os.path.exists(posf) else None
    sg = 1 if k % 2 == 0 else -1
    pr(f'lsd_run report: k = {k}, tag {tag}, grid {D["grid"]}; sym = E(a,b) + (-1)^k E(b,a), summed over the six free classes')
    for p in PCS:
        tot = {}
        for order in ('ab', 'ba'):
            R = D[order]
            part = {Nc: sum(np.sum(R[(p, c)][:Nc + 1]) + R[(p, c, 'st')] for c in C.CLS) for Nc in (50, 100, 150)}
            rsum = sum(R[(p, c)] for c in C.CLS)
            pp, tl, res = tail(rsum)
            tot[order] = (part, pp, tl, res)
        scale = max(abs(tot['ab'][0][150]), abs(tot['ba'][0][150]))
        s50 = tot['ab'][0][50] + sg*tot['ba'][0][50]; s100 = tot['ab'][0][100] + sg*tot['ba'][0][100]
        s150 = tot['ab'][0][150] + sg*tot['ba'][0][150]
        tl = tot['ab'][2] + sg*tot['ba'][2]
        line = (f'  {p:4s}: OPE(ab) {tot["ab"][0][150].real:+.6e}  OPE(ba) {tot["ba"][0][150].real:+.6e}  | sym N=50 {s50.real:+.4e}, '
                f'N=100 {s100.real:+.4e}, N=150 {s150.real:+.4e}, tail est {tl.real:+.2e} (p_ab {tot["ab"][1]:.2f}, p_ba {tot["ba"][1]:.2f})'
                f'  sym/scale {abs(s150)/scale:.2e}')
        if pos is not None and p in ('c0', 'c2', 'c3', 'c4', 'xu'):
            pa = sum(pos[(k, 'ab')][(p, c)] for c in C.CLS); pb = sum(pos[(k, 'ba')][(p, c)] for c in C.CLS)
            line += (f'  | EIN whole-function: (ab) {pa.real:+.6e} (ba) {pb.real:+.6e} sym {(pa + sg*pb).real:+.1e};'
                     f'  OPE - whole: (ab) {(tot["ab"][0][150] - pa).real:+.3e} (ba) {(tot["ba"][0][150] - pb).real:+.3e}')
        pr(line)
    with open(os.path.join(CACHE, f'lsd_report_k{k}_{tag}{SUF}.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')



def sec_exch(k, tag):
    """the s- and t-channel spin-2 exchanges (xs: single-trace s-pole residue included term-wise), term-wise vs EIN whole-function."""
    a, b, p2, p4 = P.config(tag)
    pos = pickle.load(open(os.path.join(CACHE, f'lsd_pos_{tag}_n60' + ('_k3' if k == 3 else '') + '.pkl'), 'rb'))
    sg = 1 if k % 2 == 0 else -1
    res = {}
    for order, (aa, bb) in (('ab', (a, b)), ('ba', (b, a))):
        S = C.Setup(k, aa, bb, p2, p4, **GRID)
        res[order] = profile(S, ('xs', 'xt', 'xu'), h13=0.05, Y13=10.0)
    L = []
    for p in ('xs', 'xt', 'xu'):
        for order in ('ab', 'ba'):
            v = sum(np.sum(res[order][(p, c)]) + res[order][(p, c, 'st')] for c in C.CLS)
            st = sum(res[order][(p, c, 'st')] for c in C.CLS)
            # (xu: no s-channel single-trace pole)
            w = sum(pos[(k, order)][(p, c)] for c in C.CLS)
            L.append(f'  k = {k} {p} ({order}): term-wise {v.real:+.8e} (single-trace residue {st.real:+.3e}), EIN whole {w.real:+.8e}, '
                     f'rel. diff {abs(v - w)/abs(w):.1e}')
    tot = {o: sum(np.sum(res[o][(p, c)]) + res[o][(p, c, 'st')] for p in ('xs', 'xt') for c in C.CLS) for o in ('ab', 'ba')}
    L.append(f'  k = {k} xs + xt: sym {(tot["ab"] + sg*tot["ba"]).real:+.3e} (scale {abs(tot["ab"]):.3e})')
    print('\n'.join(L), flush=True)
    with open(os.path.join(CACHE, f'lsd_exch_k{k}_{tag}.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')


if __name__ == '__main__':
    sec, k, tag = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    {'profile': sec_profile, 'report': sec_report, 'exch': sec_exch}[sec](k, tag)
