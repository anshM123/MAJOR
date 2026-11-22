"""
bj_gamma (BELOW_FOLD_SPINJ.md §4): Gamma - kappa in the window nu < nu_B from the closed-form engine, in EXACTLY the fw_gamma
(FOLD_WINDOW_LARGEJ) conventions, as pw_gamma does for the per-wave points: engine records (z2/fw format, patched conserved basis
BJ_PATCH=1 as the fw records) are injected into fw_gamma's loader at the engine nodes (default nu <= BJ_NUCUT = 9.9, i.e. strictly
below the lowest fw node), the fw direct nodes are kept everywhere else (a_h anchors x = 0.1, 0.2 and the n3/pw window nodes).
Calibration of the engine's nu-independent constants: per replica, state and reference separately, against the fw direct nodes on
nu in [BJ_CALLO, BJ_CALHI] (default 10 ... 17), as bj_chain.  A residual calibration error is a constant x D_h-class term (removed by
the cs fit) to the extent that it is a pure normalisation.
Usage: [FW_EST=dir] BJ_PATCH=1 python3 bj_gamma.py J:Delta[:depth:N] [NUMIN NUMAX]   -> cache/bj_gamma.txt (appended)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('BJ_PATCH', '1')
import numpy as np
import fw_gamma as FG
import bj_eeq as BQ
import bj_ee as BJ

NUCUT = float(os.environ.get('BJ_NUCUT', '9.9'))
CALLO, CALHI = float(os.environ.get('BJ_CALLO', '10')), float(os.environ.get('BJ_CALHI', '17'))
DREF = FG.DREF
_fw0 = FG.fwpts
STATE = {}
LOG = []


def engine(Dl, J, depth, N, Kr):
    db = pickle.load(open(BQ.fn(Dl, J, depth), 'rb'))
    rec = db['res'][(N, Kr)]
    I = np.array([[BJ.value(rec, 'EE2', a, c)[0] for c in range(2)] for a in range(2)])
    return rec['nus'], np.moveaxis(I, 2, 0)


def calib(J, D, spec):
    depth, N, Nref = spec
    nus, Ist = engine(D, J, depth, N, depth - 1)
    nur, Iref = engine(DREF, J, depth, Nref, depth - 1)
    base = _fw0('EE', J, D)
    m = FG.Z3.mphys(D, J)
    per = {}
    for (DD, nu), r in base.items():
        if DD != D or (DREF, nu) not in base or not (CALLO - 1e-6 <= nu <= CALHI + 1e-6):
            continue
        i = int(np.argmin(np.abs(nus - nu))); j = int(np.argmin(np.abs(nur - nu)))
        if abs(nus[i] - nu) > 1e-4 or abs(nur[j] - nu) > 1e-4:
            continue
        st, rf = FG.FX.fix_direct(r, J, D, nu), FG.FX.fix_direct(base[(DREF, nu)], J, DREF, nu)
        Vs = st['reps']*np.exp(st['logN']); Vr = rf['reps']*np.exp(rf['logN'])
        ks = np.einsum('ac,rac->r', np.conj(Ist[i]), Vs)/np.sum(np.abs(Ist[i])**2)
        kr = np.einsum('ac,rac->r', np.conj(Iref[j]), Vr)/np.sum(np.abs(Iref[j])**2)
        per[nu] = (ks, kr)
        LN[(J, 'st')] = st['logN']; LN[(J, 'ref')] = rf['logN']
    assert per, ('no calibration nodes', J, D)
    R = len(next(iter(per.values()))[0])
    rho = {nu: ks/kr for nu, (ks, kr) in per.items()}
    rho_r = np.mean([rho[nu] for nu in rho], axis=0); kr_r = np.mean([per[nu][1] for nu in per], axis=0)
    rm = rho_r.mean()
    LOG.append(f'J={J} D={D:g}: calibration on nu = {sorted(round(v, 2) for v in per)}; per-node m^2 Re(rho_nu/rho - 1) = '
               + ', '.join(f'{nu:.1f}: {m*m*(rho[nu].mean()/rm - 1).real:+.0f}' for nu in sorted(rho))
               + f'; replica spread m^2|d rho|/|rho| = {m*m*np.std(rho_r)/abs(rm)/np.sqrt(R):.0f} c2')
    out = {}
    import bj_chain as BC
    for i, v in enumerate(nus):
        if v > NUCUT:
            continue
        vp = round(float(v)/m, 6)*m                         # l9 record lookup matches x*m (x to 6 decimals) to 1e-7
        Is = BC.at(nus, Ist, vp); Ir = BC.at(nur, Iref, vp)
        # records carry the logN of the calibration node's direct records (q3 strips/uses the phase of logC - logK): E = value/e^logN
        ls, lr = LN[(J, 'st')], LN[(J, 'ref')]
        reps_r = np.array([kr_r[r]*Ir for r in range(R)])*np.exp(-lr)
        reps_s = np.array([kr_r[r]*rho_r[r]*Is for r in range(R)])*np.exp(-ls)
        v = vp
        for DD, reps, ln in ((D, reps_s, ls), (DREF, reps_r, lr)):
            out[(DD, round(float(v), 9))] = dict(ljc=True, E=reps.mean(0), err=0*reps.mean(0), reps=reps, logN=ln, Es=1.0 + 0j,
                                                 errs=0.0, repsS=np.ones(R, complex), logNs=0j, nu=float(v), Delta=DD, J=J,
                                                 m=FG.Z3.mphys(DD, J), method='bj engine', R=R)
    return out


SPEC = {}
MINE = {}
LN = {}


def fwpts(kind, J, D=None):
    D = D if D is not None else FG.CUR.get(J, {}).get('D', STATE.get(J))
    base = _fw0(kind, J, D)
    if J not in SPEC:
        return base
    if J not in MINE:
        MINE[J] = calib(J, D, SPEC[J])
    keep = {k: v for k, v in base.items() if k[1] > NUCUT}
    # C3H at the engine nodes: COPIES of the EE engine records, only so that q3's grid keeps these nodes (it needs both kinds);
    # the C2 values there depend on the EE records alone (a_h: fw anchors nu >= 15; bases: fw_k2 rows).  C3H output at engine
    # nodes is meaningless and is not reported.
    keep.update(MINE[J])
    return keep


FG.fwpts = fwpts
_fixd = FG.FX.fix_direct
FG.FX.fix_direct = lambda r, J, D, nu: r if r.get('method') == 'bj engine' else _fixd(r, J, D, nu)


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    lo = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    hi = float(sys.argv[3]) if len(sys.argv) > 3 else 34.0
    for spec in sys.argv[1].split(','):
        p = spec.split(':')
        J, D = int(p[0]), float(p[1])
        depth = int(p[2]) if len(p) > 2 else 3
        N = int(p[3]) if len(p) > 3 else 12
        SPEC[J] = (depth, N, int(p[4]) if len(p) > 4 else 6)
        STATE[J] = D
        pr(f"######## bj_gamma J={J} Delta={D:g} estimator {FG.EST}; engine depth {depth} N={N}; engine nodes nu <= {NUCUT}")
        o = FG.analyse(J, D, 99.0, 99.0, pr, bases=())
        for s in LOG:
            pr('# ' + s)
        if o is not None:
            m = o['m']; nus = o['nus']; xs = nus/m
            eng = nus <= NUCUT + 1e-9
            pr(f'# C2 at the engine nodes (nu <= {NUCUT}) and cs fits  [c0 D_h + s x D_h\'/(2J+1) + r] (lsq on the listed range):')
            for h in ('++', '+-'):
                r_ = o[('C2', h)]
                y = r_['y']
                pr(f'   C2 {h}: nu: ' + ' '.join(f'{v:.2f}:{yy:+.0f}' for v, yy in zip(nus, y)))
                for (a_, b_) in [(lo, NUCUT), (lo, hi), (2.0, NUCUT), (4.0, hi)]:
                    sel = (nus >= a_ - 1e-9) & (nus <= b_ + 1e-9)
                    if sel.sum() < 5:
                        continue
                    c, e, ir, res, rms = FG.fit(xs[sel], y[sel], h, J, 'cs')
                    pr(f'      fit nu in [{a_:g}, {b_:g}] ({sel.sum()} nodes, {int((sel & eng).sum())} engine): c0 = {c[0]:+9.1f}, s = {c[1]:+9.1f}, '
                       f'r = {c[ir]:+8.1f} +- {e[ir]:6.1f};  rms resid {res:7.1f};  max|y - c0 D_h - s..| = '
                       f'{np.max(np.abs(y[sel] - (np.array([*FG.dshape(h, J, xs[sel])]).T @ c[:2]))):7.1f};  max|y| {np.max(np.abs(y[sel])):7.1f}')
    with open(os.path.join(HERE, 'cache', 'bj_gamma.txt'), 'a') as f:
        f.write("\n".join(lines) + "\n")
