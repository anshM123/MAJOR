"""
directgrav/gwaves: exact <TTTT> Taylor-mode partial waves W[X_hL, X_hR] (hL, hR in {same, opp}) for ONE heavy parity-even exchange (Delta, J)
on the moment grid of e2e/moments.py -- a resumable driver around fourpt/tttt/tmoments_l2.TaylorWaves (imported READ-ONLY; nothing under
fourpt/tttt/ is written: the helicity-word cache and the mode caches live in directgrav/cache/).

Identical numerics to tmoments_l2.compute (float64, N = 8000, same grid, same convergence mask, same Fourier modes), but:
  * the state is given by (Delta_g, m/Lam, J) of the near-gap scan (Lam = sqrt((Dg-2)(Dg-1)), Delta = 3/2 + sqrt(m^2 + (J+1/2)^2));
  * the seed radial sums and every finished chunk are checkpointed (container restarts);
  * helicity structures that do not exist at this J (opp for J < 4; same for odd J) are set to zero.
Usage: python3 gwaves.py Dg mfl J [chunk]   -> cache/gw_<Dg>_<mfl>_<J>.npz  (+ cache/gw_<tag>_part/ while running)
"""
import os, sys, time, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
TT = os.path.join(HERE, '..', 'fourpt', 'tttt')
sys.path.insert(0, TT); sys.path.insert(0, os.path.join(TT, '..')); sys.path.insert(0, os.path.join(TT, '..', '..', 'e2e'))
import numpy as np, sympy as sp
import tmoments_l2 as TM
from moments import make_grid, frame_data, frame_positions, series_converged
from moments4 import fourier4
from spinblock2 import invariants4
from structures2 import words
CACHE = os.path.join(HERE, 'cache')
L, DX = 2, 3


def state(Dg, mfl, J):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); m = mfl*Lam
    Delta = 1.5 + np.sqrt(m*m + (J + 0.5)**2)
    return Lam, m, float(Delta)


def helicity_words(Dv, Jv):
    fn = os.path.join(CACHE, 'xwords_%.10f_%d.npz' % (Dv, Jv))
    if os.path.exists(fn):
        d = np.load(fn); return {'same': dict(zip(words(L), d['same'])), 'opp': dict(zip(words(L), d['opp']))}
    from t8_kappa_exact import decompose
    res = {}
    for kind in ('same', 'opp'):
        try:
            if kind == 'opp' and Jv < 4:
                raise ValueError('no opposite-helicity structure for J < 4')
            ws, sol, params = decompose(L, sp.Rational(Dv), sp.Integer(Jv), kind)
            sol = sol.subs({p: 0 for p in params})
            res[kind] = np.array([complex(sp.N(sol[i], 30)) for i in range(len(ws))])
        except Exception as e:
            print(f"  helicity structure {kind} at J={Jv}: {type(e).__name__}: {e} -> set to zero", flush=True)
            res[kind] = np.zeros(len(words(L)), dtype=complex)
    np.savez(fn, **res)
    return {k: dict(zip(words(L), v)) for k, v in res.items()}


TM.helicity_words = helicity_words          # redirect the (write-) cache of the imported class to directgrav/cache


def compute(Dg, mfl, J, chunk=200, N=8000):
    Lam, m, Dv = state(Dg, mfl, J)
    tag = 'gw_%d_%.3f_%d' % (Dg, mfl, J)
    fn = os.path.join(CACHE, tag + '.npz')
    if os.path.exists(fn):
        print('# exists', fn); return fn
    part = os.path.join(CACHE, tag + '_part'); os.makedirs(part, exist_ok=True)
    g = make_grid(Dv, J); ok = g['ok']
    fd = frame_data(g['U'][ok], g['S'][ok])
    conv = series_converged(fd['rs'], fd['chi'], J, N, cut=25.0)
    convfull = np.zeros(g['S'].shape, dtype=bool); convfull[ok] = conv
    fd = {k: v[conv] for k, v in fd.items()}
    npts = fd['rs'].size
    print(f"# <TTTT> waves: Dg={Dg} m/Lam={mfl} m={m:.4f} Delta={Dv:.6f} J={J} N={N}; {npts} grid points", flush=True)
    t0 = time.time()
    TW = TM.TaylorWaves(Dv, J, N, np.float64)
    sfn = os.path.join(part, 'seedS.npz')
    if os.path.exists(sfn):
        d = np.load(sfn, allow_pickle=True); TW.S = d['S'].item()
    else:
        TW.precompute(fd['rs'], fd['chi'])
        np.savez(sfn, S=np.array(TW.S, dtype=object))
    print(f"  setup: {time.time() - t0:.0f} s", flush=True)
    for c0 in range(0, npts, chunk):
        cfn = os.path.join(part, 'c%05d.npz' % c0)
        if os.path.exists(cfn):
            continue
        t = time.time()
        sl = slice(c0, min(npts, c0 + chunk)); n = sl.stop - sl.start
        x1, x2, x3, x4 = frame_positions(fd['ru'][sl], fd['xi'][sl])
        z1, z3, z2, z4 = TM.pol_grid4(fd['xi'][sl], L)
        s, z = invariants4({1: x1, 2: x2, 3: x3, 4: x4}, {1: z1, 2: z2, 3: z3, 4: z4})
        Cs = TW.chunk(s, fd['rs'][sl], fd['chi'][sl], fd['A0'][sl], fd['E0'][sl], sl=sl)
        ls = Dv*np.log(4*fd['rs'][sl]) + DX*(np.log(s['s13']) + np.log(s['s24']))
        keys, arrs = [], []
        for ci, c in enumerate(TM.COMBOS):
            if not Cs[c]:
                continue
            Mv = TM.contract(Cs[c], z, n).reshape(n, -1).reshape((n,) + (5,)*4)
            for q, v in fourier4(Mv, L).items():
                keys.append((ci,) + q); arrs.append(v)
        np.savez(cfn + '.tmp.npz', keys=np.array(keys), modes=np.array(arrs), ls=ls)
        os.replace(cfn + '.tmp.npz', cfn)
        print(f"  points {sl.start}-{sl.stop}: {time.time() - t:.0f} s (total {time.time() - t0:.0f} s)", flush=True)
    # assemble
    modes = {}; logscale = np.zeros(npts)
    for c0 in range(0, npts, chunk):
        d = np.load(os.path.join(part, 'c%05d.npz' % c0))
        sl = slice(c0, min(npts, c0 + chunk)); logscale[sl] = d['ls']
        for key, v in zip(d['keys'], d['modes']):
            modes.setdefault(tuple(int(x) for x in key), []).append((sl, v))
    idx = np.flatnonzero(convfull.ravel())
    keys, arrs = [], []
    for key, lst in modes.items():
        full = np.zeros(npts, dtype=complex)
        for sl, v in lst:
            full[sl] = v
        arr = np.zeros(g['S'].shape, dtype=complex); arr.ravel()[idx] = full
        keys.append(key); arrs.append(arr)
    ls = np.full(g['S'].shape, -np.inf); ls.ravel()[idx] = logscale
    np.savez_compressed(fn, keys=np.array(keys), modes=np.array(arrs), us=g['us'], ts=g['ts'], U=g['U'], S=g['S'], W=g['W'], ok=g['ok'],
                        conv=convfull, logscale=ls, Delta=Dv, J=J, N=N, npts=npts, Dg=Dg, mfl=mfl, m=m, Lam=Lam, time=time.time() - t0)
    print(f"  saved {fn} ({time.time() - t0:.0f} s)", flush=True)
    return fn


if __name__ == "__main__":
    Dg = int(sys.argv[1]); mfl = float(sys.argv[2]); J = int(sys.argv[3])
    chunk = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    compute(Dg, mfl, J, chunk)
