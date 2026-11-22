"""
lj_rs_gamma (HLJ_CLOSE.md §5): Gamma - kappa in the H_win window nu < nu_B from the closed-form engine, in EXACTLY the fw_gamma
conventions (as bj_gamma / pw_gamma), but with the engine normalised ANALYTICALLY (lj_rs_norm, fw patched basis) instead of being
calibrated on a direct node.  The reference constant per replica is still fixed on the fw direct reference records (the reference
engine equals them to ~1e-5; it carries no engine information).  Then the refined-H_win decomposition (HWIN_THIMBLE §4): cs fit
c0 D_h + s x D_h'/(2J+1) + r on the window, remainder vs T_B (HLJ_ANALYTIC table).
Engine variants (LJ_KR): the production order Kr and Kr - 1, so that the last geometric term is shown.
Usage: [FW_EST=dir] python3 lj_rs_gamma.py J:Delta:depth:N[:Nref] [NUMIN NUMAX]  -> cache/lj_rs_gamma.txt (appended)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('BJ_PATCH', '1')
import numpy as np
import bj_gamma as BG
import fw_gamma as FG
import types
FG.FX.pickle = types.SimpleNamespace(load=pickle.load, dump=lambda *a, **k: None)     # never rewrite fw_sols.pkl from here

C = os.path.join(HERE, 'cache')
KRSEL = os.environ.get('LJ_KR')


def rho_an(J, D):
    import glob
    fs = sorted(glob.glob(os.path.join(C, 'lj_rs_norm_J%d_D%g_*_p.pkl' % (J, D))))
    assert fs, ('no lj_rs_norm (patched) for', J, D)
    return complex(np.mean(pickle.load(open(fs[0], 'rb'))['rho']))


def engine(Dl, J, depth, N, Kr):
    import lj_rs_eeq as LE, bj_ee as BJ
    f = LE.resolve(Dl, J, depth, True)
    db = pickle.load(open(f, 'rb'))
    if (N, Kr) not in db['res']:
        N, Kr = max(k for k in db['res'] if k[1] <= Kr)
    rec = db['res'][(N, Kr)]
    I = np.array([[BJ.value(rec, 'EE2', a, c)[0] for c in range(2)] for a in range(2)])
    return rec['nus'], np.moveaxis(I, 2, 0)


BG.engine = engine
_calib0 = BG.calib


def calib(J, D, spec):
    """bj_gamma.calib with rho_r := rho_an (all replicas); kr_r from the fw reference records as before."""
    depth, N, Nref = spec
    kr_ = int(KRSEL) if KRSEL else depth - 1
    nus, Ist = engine(D, J, depth, N, kr_)
    nur, Iref = engine(BG.DREF, J, depth, Nref, depth - 1)
    BG.engine = lambda Dl, J_, d_, N_, K_: (nus, Ist) if Dl == D else (nur, Iref)
    try:
        out = _calib0(J, D, spec)
    finally:
        BG.engine = engine
    rho = rho_an(J, D)
    # rebuild the state records with rho_an: _calib0 stored reps_s = kr_r rho_r I_s e^{-ls}; replace rho_r by rho_an
    import bj_chain as BC
    m = FG.Z3.mphys(D, J)
    ls = BG.LN[(J, 'st')]
    fixed = {}
    for key, r in out.items():
        if key[0] != D:
            fixed[key] = r
            continue
        v = key[1]
        Is = BC.at(nus, Ist, v)
        reps_ref = out[(BG.DREF, key[1])]['reps']*np.exp(BG.LN[(J, 'ref')])            # kr_r I_ref
        Ir = BC.at(nur, Iref, v)
        kr_r = np.array([np.sum(np.conj(Ir)*x)/np.sum(np.abs(Ir)**2) for x in reps_ref])
        reps = np.array([kr_r[i]*rho*Is for i in range(len(kr_r))])*np.exp(-ls)
        r2 = dict(r); r2['reps'] = reps; r2['E'] = reps.mean(0)
        fixed[key] = r2
    BG.LOG.append(f'J={J} D={D:g}: rho := rho_an = {rho:.8e} (lj_rs_norm, patched basis); engine order Kr = {kr_}; '
                  f'm^2 (rho_fit/rho_an - 1) of the fitted constant above is in the line before')
    # the fitted rho of _calib0 (for the record)
    return fixed


BG.calib = calib

if __name__ == '__main__':
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    lo = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    hi = float(sys.argv[3]) if len(sys.argv) > 3 else 34.0
    spec = sys.argv[1]
    p = spec.split(':')
    J, D = int(p[0]), float(p[1])
    depth = int(p[2]); N = int(p[3]); Nref = int(p[4]) if len(p) > 4 else 4
    BG.SPEC[J] = (depth, N, Nref)
    BG.STATE[J] = D
    pr(f"######## lj_rs_gamma J={J} Delta={D:g} estimator {FG.EST}; engine depth {depth} N={N} Kr={KRSEL or depth - 1}; "
       f"engine nodes nu <= {BG.NUCUT}; rho analytic")
    o = FG.analyse(J, D, 99.0, 99.0, pr, bases=())
    for s in BG.LOG:
        pr('# ' + s)
    if o is not None:
        m = o['m']; nus = o['nus']; xs = nus/m
        eng = nus <= BG.NUCUT + 1e-9
        pr(f'# C2 Gamma - kappa at all nodes (engine nu <= {BG.NUCUT}) and cs fits [c0 D_h + s x D_h\'/(2J+1) + r]:')
        res = {}
        for h in ('++', '+-'):
            y = o[('C2', h)]['y']
            res[h] = y
            pr(f'   C2 {h}: ' + ' '.join(f'{v:.2f}:{yy:+.0f}' for v, yy in zip(nus, y)))
            for (a_, b_) in [(lo, BG.NUCUT), (lo, hi), (2.0, BG.NUCUT)]:
                sel = (nus >= a_ - 1e-9) & (nus <= b_ + 1e-9)
                if sel.sum() < 5:
                    continue
                c, e, ir, rs_, rms = FG.fit(xs[sel], y[sel], h, J, 'cs')
                resid = y[sel] - (np.array([*FG.dshape(h, J, xs[sel])]).T @ c[:2])
                pr(f'      fit nu in [{a_:g}, {b_:g}] ({sel.sum()} nodes, {int((sel & eng).sum())} engine): c0 = {c[0]:+9.1f}, '
                   f's = {c[1]:+9.1f}, r = {c[ir]:+8.1f} +- {e[ir]:6.1f}; rms resid {rs_:7.1f}; max|y - c0 D_h - s..| = '
                   f'{np.max(np.abs(resid)):7.1f}; max|y| {np.max(np.abs(y[sel])):7.1f}')
        pickle.dump(dict(nus=nus, m=m, y=res), open(os.path.join(C, 'lj_rs_gamma_J%d_D%g_K%s.pkl' % (J, D, KRSEL or depth - 1)), 'wb'))
    with open(os.path.join(C, 'lj_rs_gamma.txt'), 'a') as f:
        f.write("\n".join(lines) + "\n")
