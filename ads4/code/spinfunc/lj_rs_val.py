"""
lj_rs_val (HLJ_CLOSE.md §3): calibration-free, entry- and helicity-resolved validation of the closed-form spin-J EE engine (bj_eeq images)
against the independent direct crystal nodes (z2 n3 records, state 336 / reference 2e4 on common random numbers), with the analytic
normalisation rho_an of lj_rs_norm (NO fitted constant).

Per node nu and common-seed replica r (CV estimator of bj_chain/z3: scalar MC ratio replaced by the exact one):
  entry  (a, c):  m^2 Re[ (V_st/V_ref)_ac / (rho_an I_st/I_ref)_ac - 1 ]
  helicity h   :  m^2 Re[ lam_h(V_st V_ref^-1) / lam_h(rho_an I_st I_ref^-1) - 1 ]    (eigenvalues: invariant under any change of the
                  conserved basis common to state and reference; paired by proximity)
  structure    :  entry (a, c) minus entry (0, 0)  (calibration- and normalisation-free), with its own jackknife error
and the engine's own truncation band at the node: the same quantities between consecutive geometric orders (Kr vs Kr-1) and D-hat
orders (N vs N-1).
Usage: python3 lj_rs_val.py J depth N Kr [Nref]   -> cache/lj_rs_val_J<J>_d<depth>_N<N>_K<Kr>.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DST, DREF = 336.0, 20000.0


def eng(Dl, J, depth, N, Kr):
    import bj_ee as BJ, lj_rs_eeq as LE
    f = LE.resolve(Dl, J, depth)
    if f is None and Dl == DREF:              # reference: its a^k terms are ~(5e-5)^k -- any available depth is exact at the c2 level
        f = next(q for q in (LE.resolve(Dl, J, d_) for d_ in (3, 2, 4)) if q)
    db = pickle.load(open(f, 'rb'))
    if (N, Kr) not in db['res'] and Dl == DREF:
        N, Kr = max(k for k in db['res'] if k[0] <= N)
    rec = db['res'][(N, Kr)]
    I = np.array([[BJ.value(rec, 'EE2', a, c)[0] for c in range(2)] for a in range(2)])
    return rec['nus'], np.moveaxis(I, 2, 0)


def rho_an(J, depth, N, Kr):
    f = os.path.join(C, 'lj_rs_norm_J%d_D%g_d%d_N%d_K%d.pkl' % (J, DST, depth, N, Kr))
    if not os.path.exists(f):                 # rho_an is truncation-independent (dropped powers / branch constants): any run of this J
        import glob
        f = sorted(glob.glob(os.path.join(C, 'lj_rs_norm_J%d_D%g_*.pkl' % (J, DST))))[0]
    o = pickle.load(open(f, 'rb'))
    if 'rho_rs0.3' in o:                     # smallest radial scale: truncation-free limit to <= 2e-6 (J = 28) / 2e-5 (J = 32)
        return complex(np.mean(o['rho_rs0.3']))
    return complex(np.mean(o['rho']))


def eig_pair(T, ref=None):
    ev = np.linalg.eigvals(T)
    if ref is None:
        return ev[np.argsort(-np.abs(ev))]
    if abs(ev[0] - ref[0]) + abs(ev[1] - ref[1]) <= abs(ev[1] - ref[0]) + abs(ev[0] - ref[1]):
        return ev
    return ev[::-1]


def jk(vals):
    R = len(vals)
    j = np.array([np.mean(np.delete(vals, i, axis=0), axis=0) for i in range(R)])
    return np.sqrt((R - 1)/R*np.sum(np.abs(j - j.mean(0))**2, axis=0))


def run(J, depth, N, Kr, Nref=4, pr=print, variants=None):
    import bj_chain as BC, z3_nodes as Z3
    rho = rho_an(J, depth, N, Kr)
    nus, Ist = eng(DST, J, depth, N, Kr)
    nur, Iref = eng(DREF, J, depth, Nref, Kr)
    variants = variants or []
    VAR = [(lab, eng(DST, J, d_, n_, k_)) for (lab, d_, n_, k_) in variants]
    P = BC.direct(J)
    m = BC.mphys(DST, J); m2 = m*m
    pr(f'# lj_rs_val J={J}: engine depth {depth} N={N} Kr={Kr} (ref N={Nref}); rho_an = {rho:.8e} (lj_rs_norm, no fit); CV estimator')
    pr('#   columns: entry (0,0) (0,1) (1,0) (1,1) | structure (0,1)-(0,0) (1,0)-(0,0) (1,1)-(0,0) | helicity big/small eigenvalue;'
       ' value +- jackknife [c2]; then engine band rows (same quantities, engine variant - production)')
    out = {}
    for v in sorted(set(k[1] for k in P if k[0] == DST and (DREF, k[1]) in P)):
        if v < float(os.environ.get('LJ_VMIN', '8.0')):
            continue
        st, rf = P[(DST, v)], P[(DREF, v)]
        try:
            Is, Ir = BC.at(nus, Ist, v), BC.at(nur, Iref, v)
        except AssertionError:
            continue
        Vs = st['reps']*np.exp(st['logN']); Vr = rf['reps']*np.exp(rf['logN'])
        o = Z3.scalar_check(J, v, {DST: st, DREF: rf}, DST)
        sr = (st['repsS']*np.exp(st['logNs']))/(rf['repsS']*np.exp(rf['logNs']))
        cvr = o['cv']*sr.mean()/sr
        R = len(cvr)
        Te = rho*Is @ np.linalg.inv(Ir)
        le = eig_pair(Te)
        Ee = rho*Is/Ir
        ent = np.array([m2*((Vs[r]/Vr[r])*cvr[r]/Ee - 1).real.ravel() for r in range(R)])           # (R, 4)
        strc = ent[:, 1:] - ent[:, :1]
        hel = np.array([m2*(eig_pair(Vs[r] @ np.linalg.inv(Vr[r])*cvr[r], le)/le - 1).real for r in range(R)])
        row = dict(ent=(ent.mean(0), jk(ent)), strc=(strc.mean(0), jk(strc)), hel=(hel.mean(0), jk(hel)))
        bands = {}
        for lab, (nv, Iv) in VAR:
            try:
                Isv = BC.at(nv, Iv, v)
            except AssertionError:
                continue
            Tv = rho*Isv @ np.linalg.inv(Ir)
            ev = m2*((rho*Isv/Ir)/Ee - 1).real.ravel()
            bands[lab] = dict(ent=ev, strc=ev[1:] - ev[0], hel=m2*(eig_pair(Tv, le)/le - 1).real)
        row['bands'] = bands
        out[v] = row
        f = lambda x, e: ' '.join(f'{a:+7.0f}+-{b:<5.0f}' for a, b in zip(x, e))
        pr(f'  nu={v:7.3f} entry {f(*row["ent"])} | struct {f(*row["strc"])} | hel {f(*row["hel"])}')
        for lab, b in bands.items():
            pr(f'     band {lab:>10s} entry ' + ' '.join(f'{x:+7.0f}' for x in b['ent']) + '        | struct '
               + ' '.join(f'{x:+7.0f}' for x in b['strc']) + '        | hel ' + ' '.join(f'{x:+7.0f}' for x in b['hel']))
    return out


if __name__ == '__main__':
    J, depth, N, Kr = (int(x) for x in sys.argv[1:5])
    Nref = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    var = []
    if depth == 3:
        var = [('a1(d3)', 3, N, Kr - 1), ('a1d2N5', 2, 5, 1)]
    elif depth == 2:
        var = [('a0', 2, N, 0), ('N-1', 2, N - 1, Kr)]
    lines = []

    def pr(s=''):
        print(s, flush=True); lines.append(s)
    out = run(J, depth, N, Kr, Nref, pr, var)
    tag = 'lj_rs_val_J%d_d%d_N%d_K%d%s' % (J, depth, N, Kr, os.environ.get('LJ_VTAG', ''))
    open(os.path.join(C, tag + '.txt'), 'w').write('\n'.join(lines) + '\n')
    pickle.dump(out, open(os.path.join(C, tag + '.pkl'), 'wb'))
