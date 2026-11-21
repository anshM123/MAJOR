"""
bg_eeq (BELOW_FOLD2.md §2.3, §4): production of the analytic-radial EE images.  For a state (Delta, J) and the reference
(Delta_ref = 2e4, same J): EE k = 2 (pref x (v'-u')) and k = 4 (x Lemma-B R) integrands on bg_ser, assembled (bg_geom.assemble),
times the B-correction, integrated in closed form order by order on the nu grid.
  I_kind(nu; Delta) = exp(lconst) sum_j (i nu/2)^j sum_k Delta^{-k} M[kind][(k, j)](nu)      (overall a-power dropped)
For J >= 4 (two conserved structures) the 2 x 2 matrix pref[a][c] is kept (all four entries).
Usage: python3 bg_eeq.py Delta J K depth [Kref]   -> cache/bg_eeq_D{Delta}_J{J}_K{K}.pkl
"""
import os, sys, pickle, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_ser as S, bg_geom as GG, bg_ee as BE

C = os.path.join(HERE, 'cache')
DREF = 2e4
NUS = np.round((np.arange(1025) + 0.5)*0.04, 4)          # the dnu = 0.04 de-shift grid, nu < 41


def fn(Dl, J, K, tag=''):
    return os.path.join(C, 'bg_eeq_D%g_J%d_K%d%s.pkl' % (Dl, J, K, tag))


def images(Dl, J, K, depth, nus=NUS, kinds=('EE4', 'EE2'), verbose=True):
    t = time.time()
    P = BE.ee_pieces(depth, Dl, J, verbose=verbose)
    Bf = GG.bfac(K)
    out = dict(Delta=Dl, J=J, K=K, depth=depth, nus=np.asarray(nus), res={})
    na, nc = len(P['pref']), len(P['pref'][0])
    if verbose:
        print(f'  pieces D={Dl} J={J}: {time.time() - t:.0f}s; pref d/K', [(P['pref'][a][c].d, P['pref'][a][c].K)
                                                                    for a in range(na) for c in range(nc)], flush=True)
    for a in range(na):
        for c in range(nc):
            for kind in kinds:
                polys = [('vmu', P['vmu']), ('pref', P['pref'][a][c])] + ([('R', P['R'])] if kind == 'EE4' else [])
                A = GG.assemble(P['logs'], polys, K)
                M = {}
                for j, cj in A['corr'].items():
                    Pk = S.lmul(cj, Bf.c, K)
                    for k in range(K + 1):
                        M[(k, j)] = S.integrate_fast(Pk[k].E, Pk[k].C, A["E0"], A["E1"], out["nus"]) if len(Pk[k]) else 0*out['nus']
                out['res'][(kind, a, c)] = dict(E0=A['E0'], E1=A['E1'], lconst=A['lconst'], M=M,
                                                nterms=[len(p) for p in A['corr'][0]])
                if verbose:
                    print(f'   {kind} ({a},{c}): terms {out["res"][(kind, a, c)]["nterms"]}, {time.time() - t:.0f}s', flush=True)
    return out


def value(rec, kind, a=0, c=0, Kt=None, Dl=None):
    """I(nu) = exp(lconst) sum_j (i nu/2)^j sum_{k<=Kt} Dl^{-k} M[(k, j)] and per-order pieces (k x nu)."""
    r = rec['res'][(kind, a, c)]
    Kt = rec['K'] if Kt is None else Kt
    Dl = rec['Delta'] if Dl is None else Dl
    nus = rec['nus']
    pre = np.exp(r['lconst'][0] + 1j*nus*r['lconst'][1])
    pieces = np.zeros((Kt + 1, len(nus)), complex)
    for (k, j), v in r['M'].items():
        if k <= Kt:
            pieces[k] += Dl**(-k)*(0.5j*nus)**j*v
    return pre*pieces.sum(0), pre[None, :]*pieces


def get(Dl, J, K, depth, Kref=2, depth_ref=None, tag=''):
    path = fn(Dl, J, K, tag)
    if os.path.exists(path):
        return pickle.load(open(path, 'rb'))
    st = images(Dl, J, K, depth)
    rf = images(DREF, J, Kref, depth_ref or (Kref + (depth - K)))
    db = dict(st=st, ref=rf)
    pickle.dump(db, open(path + '.tmp', 'wb')); os.replace(path + '.tmp', path)
    return db


if __name__ == '__main__':
    Dl, J, K, depth = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    Kref = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    t0 = time.time()
    get(Dl, J, K, depth, Kref)
    print(f'bg_eeq D={Dl} J={J} K={K}: done {time.time() - t0:.0f}s', flush=True)
