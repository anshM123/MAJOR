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
    assert per, ('no calibration nodes', J, D)
    R = len(next(iter(per.values()))[0])
    rho = {nu: ks/kr for nu, (ks, kr) in per.items()}
    rho_r = np.mean([rho[nu] for nu in rho], axis=0); kr_r = np.mean([per[nu][1] for nu in per], axis=0)
    rm = rho_r.mean()
    LOG.append(f'J={J} D={D:g}: calibration on nu = {sorted(round(v, 2) for v in per)}; per-node m^2 Re(rho_nu/rho - 1) = '
               + ', '.join(f'{nu:.1f}: {m*m*(rho[nu].mean()/rm - 1).real:+.0f}' for nu in sorted(rho))
               + f'; replica spread m^2|d rho|/|rho| = {m*m*np.std(rho_r)/abs(rm)/np.sqrt(R):.0f} c2')
    out = {}
    for i, v in enumerate(nus):
        if v > NUCUT:
            continue
        j = int(np.argmin(np.abs(nur - v)))
        reps_r = np.array([kr_r[r]*Iref[j] for r in range(R)])
        reps_s = np.array([kr_r[r]*rho_r[r]*Ist[i] for r in range(R)])
        for DD, reps in ((D, reps_s), (DREF, reps_r)):
            out[(DD, round(float(v), 6))] = dict(ljc=True, E=reps.mean(0), err=0*reps.mean(0), reps=reps, logN=0j, Es=1.0 + 0j,
                                                 errs=0.0, repsS=np.ones(R, complex), logNs=0j, nu=float(v), Delta=DD, J=J,
                                                 m=FG.Z3.mphys(DD, J), method='bj engine', R=R)
    return out


SPEC = {}
MINE = {}


def fwpts(kind, J, D=None):
    D = D if D is not None else FG.CUR.get(J, {}).get('D', STATE.get(J))
    base = _fw0(kind, J, D)
    if kind != 'EE' or J not in SPEC:
        return base
    if J not in MINE:
        MINE[J] = calib(J, D, SPEC[J])
    keep = {k: v for k, v in base.items() if k[1] > NUCUT}
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
        FG.analyse(J, D, lo, hi, pr, bases=('cs',))
        for s in LOG:
            pr('# ' + s)
    with open(os.path.join(HERE, 'cache', 'bj_gamma.txt'), 'a') as f:
        f.write("\n".join(lines) + "\n")
