"""
ng_val (NEARGAP2.md §1.2): direct (n1 RQMC) crystal points on an own-mass row against the K2 row, node by node in Gamma (c2) and at the
level of the state's lambda/N_s: the direct records (hc_nodes.to_record, read-only: M0 = all-orders direct matrix in the row's shadow
basis, M2 = M4 = 0) override the K2 records at their nu, exactly as hc_chain's HC_INJ does; jackknife over the 8 RQMC replicas.
Usage: python3 ng_val.py Dg J Delta_own mfl[,mfl...]   -> cache/ng_val_J<J>_D<Delta>.txt
"""
import os, sys, glob, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ng_chain as NC
import numpy as np
E, Q, W = NC.E, NC.Q, NC.W
import hc_nodes as HN
NGC = NC.NGC
HN.C = NGC                                   # direct points of this note (hc_pts_* in neargap2/cache)
_c0 = HN._cands


def _cands(odd):
    save = HN.C; HN.C = NC.SFC                   # the project's shadow caches (read-only) ...
    try:
        d = _c0(odd)
    finally:
        HN.C = save                              # ... while the direct points come from neargap2/cache
    key = 'odd_ng' if odd else 'even_ng'
    if key not in HN._SC:
        dd = {k: list(v) for k, v in d.items()}
        for f in sorted(glob.glob(os.path.join(NGC, 'ng_p2shadow.pkl' if odd else 'ng_shadow_k2.pkl'))) + \
                [os.path.join(NC.SFC, p) for p in (['p2_shadow.pkl'] if odd else ['l9_shadow_k2.pkl'])]:
            if os.path.exists(f):
                for k, v in pickle.load(open(f, 'rb')).items():
                    dd.setdefault(k, []).append((os.path.basename(f), v))
        HN._SC[key] = dd
    return HN._SC[key]


HN._cands = _cands
INJ = {}
JK = [None]


def _wrap(loader, odd):
    def f(kind, J, *a, **k):
        R = loader(kind, J, *a, **k)
        if J not in INJ:
            return R
        kf = 'EE' if kind == 'C2' else kind
        R = dict(R)
        HN.set_resid(R, J, kf)
        for D in INJ[J]:
            add = Q._strip(HN.records(kf, J, D, jk=JK[0]), J)
            for kk, v in add.items():
                for k2 in [k2 for k2 in R if k2[0] == kk[0] and abs(k2[1] - kk[1]) < 1e-6*kk[1]]:
                    R.pop(k2)
                R[kk] = v
        return R
    return f


Q.P4.load = _wrap(Q.P4.load, True)
_wl = _wrap(W._load, False)
W._load = _wl
Q.A._load = _wl

def _fa(a):
    return ', '.join(f'{k}: {complex(v).real:.6g}' for k, v in a.items()) if isinstance(a, dict) else f'{a:.6g}'


if __name__ == '__main__':
    Dg, J = float(sys.argv[1]), int(sys.argv[2])
    Downs = [float(v) for v in sys.argv[3].split(',')]
    mfls = [float(v) for v in sys.argv[4].split(',')]
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    NC.KK = 2
    gK2 = Q.odd_grid(J) if J % 2 else Q.even_grid(J)
    INJ[J] = Downs
    reps = []
    for jk in [None] + list(range(8)):
        JK[0] = jk
        reps.append((Q.odd_grid(J) if J % 2 else Q.even_grid(J)))
    keys = ['C2', 'C3H'] if J % 2 else [(kd, h) for kd in ('C2', 'C3H') for h in (('++', '+-') if J >= 4 else ('++',))]
    pr(f"# ng_val J={J} rows {Downs}: Gamma (c2) at the direct nodes, K2 row vs direct (jackknife over 8 replicas); rows of the model "
       f"{sorted(gK2)}")
    for Down in Downs:
        g0 = gK2[Down]; gd = reps[0][Down]
        P = HN.pts('EE', J, Down); m = HN.mphys(Down, J)
        xdir = sorted(round(k[1]/m, 6) for k in P if k[0] == Down)
        pr(f"#  row Delta={Down:g} (m={m:.3f}); normalisation a: K2 {_fa(g0['a'])}, with direct nodes {_fa(gd['a'])}")
        for x in xdir:
            i0 = int(np.argmin(np.abs(g0['xs'] - x))); i1 = int(np.argmin(np.abs(gd['xs'] - x)))
            if abs(gd['xs'][i1] - x) > 1e-5:
                pr(f"   x={x:.4f}: node dropped (despike / incomplete)"); continue
            k2 = abs(g0['xs'][i0] - x) <= 1e-5
            s = f"   x={x:.4f} nu={x*m:8.2f}:"
            for k in keys:
                vj = np.array([r[Down][k][i1] for r in reps[1:]])
                err = np.sqrt(7/8*((vj - vj.mean())**2).sum())
                kv = f"{g0[k][i0]:+9.1f}" if k2 else "    (none)"
                dv = f"{gd[k][i1] - g0[k][i0]:+8.1f}" if k2 else "      -"
                s += f"  {k}: K2 {kv} dir {gd[k][i1]:+9.1f} +- {err:6.1f} (diff {dv})"
            pr(s)
    for mfl in mfls:
        lam = []
        for r in reps:
            grids = {'even': {}, 'odd': {}}; grids['odd' if J % 2 else 'even'][J] = r
            em, om = NC.models(grids)
            f = (om if J % 2 else em)[J][0]
            lam.append(NC.lam_of(NC.dens(Dg, J, mfl, f, 256))[0])
        grids = {'even': {}, 'odd': {}}; grids['odd' if J % 2 else 'even'][J] = gK2
        em, om = NC.models(grids); f = (om if J % 2 else em)[J][0]
        lk = NC.lam_of(NC.dens(Dg, J, mfl, f, 256))[0]
        lj = np.array(lam[1:]); e = np.sqrt(7/8*((lj - lj.mean())**2).sum())
        pr(f"   state m/Lam={mfl}: lam/N_s K2 rows {lk:+.6f} | with direct nodes {lam[0]:+.6f} +- {e:.6f} (jackknife) | shift {lam[0] - lk:+.6f}")
    open(os.path.join(NGC, 'ng_val_J%d_%s.txt' % (J, '_'.join('%g' % D for D in Downs))), 'w').write("\n".join(lines) + "\n")
