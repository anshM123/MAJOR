"""
pw_gamma (HWIN_THIMBLE.md): Gamma - kappa at the pw_points nodes in EXACTLY the fw_gamma (FOLD_WINDOW_LARGEJ) conventions: the
per-wave state records (sum of the modes) and the P_J reference records are converted to the z2/fw direct-record format and injected
into fw_gamma's loader (in-process; no cache of fw_gamma is written except its own appended logs, which are NOT used here: output goes
to cache/pw_gamma.txt).  fw nodes at the same nu are replaced by the pw nodes; fw nodes at other nu (the a_h anchors x = 0.1, 0.2, and
for PW_KEEPFW=1 everything else) are kept.  Helicity bases: fw_gamma.gbasis on the fw_k2 EE rows (state and reference).
Usage: [FW_EST=cv|dir] [PW_TAGS=tag1,tag2] python3 pw_gamma.py J:Delta [NUMIN NUMAX]
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import fw_gamma as FG
import pw_points as PP

TAGS = [t for t in os.environ.get('PW_TAGS', '').split(',')]
KEEPFW = os.environ.get('PW_KEEPFW', '0') == '1'
XANCHOR = float(os.environ.get('PW_XANCHOR', '0.08'))     # fw nodes kept (a_h anchors): x >= this


def to_z2(rec, kind):
    gs = rec['g'][kind]
    ref = gs[0]['logN']
    reps = sum(np.asarray(g['reps'])*np.exp(g['logN'] - ref) for g in gs)
    pm = rec['per_mode']
    refs_ = pm[0]['logN']
    repsS = sum(np.asarray(o['reps'])*np.exp(o['logN'] - refs_) for o in pm)
    R = reps.shape[0]
    return dict(ljc=True, E=reps.mean(0), err=reps.std(0, ddof=1)/np.sqrt(R), reps=reps, logN=ref,
                Es=repsS.mean(), errs=repsS.std(ddof=1)/np.sqrt(R), repsS=repsS, logNs=refs_, nu=rec['nu'], Delta=rec['Delta'],
                J=rec['J'], m=rec['m'], method='pw ' + str(rec['modes']), R=R)


def pw_records(kind, J, D):
    out = {}
    for t in TAGS:
        PP.TAG = t
        db = PP.load(J, D)
        for (DD, nu), rec in db.items():
            if kind in rec.get('g', {}):
                out[(DD, nu)] = to_z2(rec, kind)
    PP.TAG = ''
    return out


_fw0 = FG.fwpts
STATE = {}


def fwpts(kind, J, D=None):
    D = D if D is not None else FG.CUR.get(J, {}).get('D', STATE.get(J))
    base = _fw0(kind, J, D)
    mine = pw_records(kind, J, D)
    mynus = set(nu for (_, nu) in mine)
    m = FG.Z3.mphys(D, J)
    keep = {}
    for (DD, nu), r in base.items():
        if nu in mynus:
            continue
        if not KEEPFW and nu/m < XANCHOR:
            continue
        keep[(DD, nu)] = r
    keep.update(mine)
    # only nodes with both state and reference
    nus_st = set(nu for (DD, nu) in keep if DD == float(D)); nus_rf = set(nu for (DD, nu) in keep if DD == FG.DREF)
    return {k: v for k, v in keep.items() if k[1] in nus_st & nus_rf and k[0] in (float(D), FG.DREF)}


FG.fwpts = fwpts


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    lo = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    hi = float(sys.argv[3]) if len(sys.argv) > 3 else 34.0
    store = {}
    pk = os.path.join(PP.C, 'pw_gamma.pkl')
    if os.path.exists(pk):
        store = pickle.load(open(pk, 'rb'))
    for spec in sys.argv[1].split(','):
        J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
        STATE[J] = D
        pr(f"######## pw_gamma J={J} Delta={D:g} estimator {FG.EST} tags {TAGS} keepfw {KEEPFW}")
        o = FG.analyse(J, D, lo, hi, pr, bases=('cs',))
        if o is not None:
            store[(J, D, lo, hi, FG.EST, tuple(TAGS))] = o
    pickle.dump(store, open(pk, 'wb'))
    with open(os.path.join(PP.C, 'pw_gamma.txt'), 'a') as f:
        f.write("\n".join(lines) + "\n")
