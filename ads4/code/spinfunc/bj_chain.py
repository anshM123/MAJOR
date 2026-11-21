"""
bj_chain (BELOW_FOLD_SPINJ.md §3): the closed-form EE engine's values (bj_eeq, kind EE2 = the C2/C^'_EE records) fed into Result D's
Delta_g = 300 chain (hc_chain z4 mode, HYPOTHESES_CLOSE §5 protocol, unchanged) as direct small-nu nodes below the fold.

1. calibrate(): the engine images carry an unknown nu-independent constant per Delta (dropped overall powers / normalisations).  The
   constants are fixed against the z2 direct records (n3 contour, state 336 and reference 2e4 on common random numbers) at the overlap
   nodes nu = 8.33 ... 16.7 (default; Z2NODES of bj_eeq), separately for the state and the reference and per common-seed replica, so that
   the CRN-correlated errors largely cancel in the state/reference ratio kappa_st/kappa_ref; jackknife over the 8 replicas.
   Shape check: per node, m^2 Re(rho_nu/rho - 1) (rho = kappa_st,nu/kappa_ref,nu), i.e. the engine-vs-direct difference in c2 units of
   the state/reference ratio, at every overlap node 8.33 ... 33.3.
2. records(): z2-format records (E 2x2 in the conserved-structure basis of l5_spin.conserved, reps = 8 replicas carrying the replica-
   wise calibration, logN = 0) at the engine nodes nu < NUCUT, state AND reference; written to cache/z2bj_pts_EE_J<J>.pkl.
3. chain(): runs hc_chain.py z4 300 in-process with z3_nodes.zpts patched: EE = z2 (+x) records with the engine nodes added below
   NUCUT (replacing the z2 nodes at 5.0 / 6.67); C3H below NUCUT: 'keep' (z2 x records, as the main run) | 'drop'.
Usage: python3 bj_chain.py calib J [N Kr]
       BJ_C3H=keep|drop BJ_NUCUT=8.3 HC_TAG=bj28 [hc/z4 env] python3 bj_chain.py chain J [N Kr]
"""
import os, sys, pickle, runpy
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bj_eeq as BQ
import bj_ee as BJ

C = os.path.join(HERE, 'cache')
DREF = 20000.0
DST = 336.0
CAL = [float(v) for v in os.environ.get('BJ_CAL', '8.332091574,9.998509889,11.664928204,13.331346519,16.664183148').split(',')]
NUCUT = float(os.environ.get('BJ_NUCUT', '8.3'))


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def engine(Dl, J, depth, N, Kr):
    db = pickle.load(open(BQ.fn(Dl, J, depth), 'rb'))
    rec = db['res'][(N, Kr)]
    I = np.array([[BJ.value(rec, 'EE2', a, c)[0] for c in range(2)] for a in range(2)])     # (a, c, nu)
    return rec['nus'], np.moveaxis(I, 2, 0)                                                 # (nu, a, c)


def direct(J):
    import z3_nodes as Z3
    return Z3.zpts('EE', J)


def _near(nus, v):
    i = int(np.argmin(np.abs(nus - v)))
    assert abs(nus[i] - v) < 1e-5, (v, nus[i])
    return i


def at(nus, I, v):
    """engine matrix at nu = v: the grid value if v is a grid node, else a cubic spline of the unwrapped log (grid step 0.25 below 12)."""
    i = int(np.argmin(np.abs(nus - v)))
    if abs(nus[i] - v) < 1e-5:
        return I[i]
    from scipy.interpolate import CubicSpline
    sel = (nus > v - 1.6) & (nus < v + 1.6)
    assert sel.sum() >= 6 and np.max(np.diff(nus[sel])) < 0.3, ('no dense grid around', v)
    out = np.zeros((2, 2), complex)
    for a in range(2):
        for c in range(2):
            L = np.log(I[sel, a, c]); L = L.real + 1j*np.unwrap(L.imag)
            out[a, c] = np.exp(CubicSpline(nus[sel], L.real)(v) + 1j*CubicSpline(nus[sel], L.imag)(v))
    return out


def calibrate(J, depth=3, N=8, Kr=2, depth_ref=None, Nref=4, pr=print):
    depth_ref = depth if depth_ref is None else depth_ref
    nus, Ist = engine(DST, J, depth, N, Kr)
    nur, Iref = engine(DREF, J, depth_ref, Nref, Kr)
    P = direct(J)
    m = mphys(DST, J)
    allnodes = sorted(set(k[1] for k in P if k[0] == DST and (DREF, k[1]) in P))
    R = len(P[(DST, allnodes[0])]['reps'])

    def kap(I, V):          # per replica: sum conj(I) V / sum |I|^2 over the 2x2 entries
        return np.einsum('ac,rac->r', np.conj(I), V)/np.sum(np.abs(I)**2)
    per = {}
    for v in allnodes:
        if v < 8.0:
            continue
        st, rf = P[(DST, v)], P[(DREF, v)]
        Vs = st['reps']*np.exp(st['logN']); Vr = rf['reps']*np.exp(rf['logN'])
        try:
            Is, Ir = at(nus, Ist, v), at(nur, Iref, v)
        except AssertionError:
            continue
        ks = kap(Is, Vs); kr = kap(Ir, Vr)
        # entrywise residual of the shape (state and reference separately), relative
        es = np.max(np.abs(Vs.mean(0)/(ks.mean()*Is) - 1)); er = np.max(np.abs(Vr.mean(0)/(kr.mean()*Ir) - 1))
        per[v] = dict(ks=ks, kr=kr, rho=ks/kr, es=es, er=er)
    cal = [v for v in per if any(abs(v - c) < 0.05 for c in CAL)]
    # replica-wise calibration: mean over the calibration nodes of rho_r (weights 1/var of the jackknife)
    w = {}
    for v in cal:
        rr = per[v]['rho']
        jk = np.array([np.mean(np.delete(rr, i)) for i in range(R)])
        var = (R - 1)/R*np.sum(np.abs(jk - jk.mean())**2)
        w[v] = 1.0/max(var, 1e-30)
    W = sum(w.values())
    rho_r = sum(w[v]*per[v]['rho'] for v in cal)/W
    kr_r = sum(w[v]*per[v]['kr'] for v in cal)/W
    rho = rho_r.mean()
    jk = np.array([np.mean(np.delete(rho_r, i)) for i in range(R)])
    err = np.sqrt((R - 1)/R*np.sum(np.abs(jk - jk.mean())**2))
    pr(f'# bj_chain calibration J={J} (state {DST:g}, ref {DREF:g}; engine depth {depth} N={N} Kr={Kr}, ref depth {depth_ref} N={Nref}): '
       f'rho = kappa_st/kappa_ref fitted on nu = {[round(v, 2) for v in cal]}; m^2 |d rho/rho| (jackknife) = {m*m*err/abs(rho):.1f} c2')
    pr('#   node nu: shape check m^2 Re(rho_nu/rho - 1) [c2] +- jk | max entrywise |V/(kappa I) - 1| state / reference (MC + shape)')
    out = {}
    for v in sorted(per):
        rr = per[v]['rho']
        jkv = np.array([np.mean(np.delete(rr, i)) for i in range(R)])
        e = np.sqrt((R - 1)/R*np.sum(np.abs(jkv - jkv.mean())**2))
        d = m*m*(rr.mean()/rho - 1)
        out[v] = (d, m*m*e/abs(rho))
        pr(f'    {v:7.3f}: {d.real:+9.1f} {d.imag:+9.1f}i +- {m*m*e/abs(rho):7.1f}   | {per[v]["es"]:.2e} / {per[v]["er"]:.2e}'
           + ('   [cal]' if v in cal else ''))
    return dict(rho_r=rho_r, kr_r=kr_r, rho=rho, err=err, nodes=out, nus=nus, Ist=Ist, nur=nur, Iref=Iref, R=R, m=m)


def k2check(J, cal, pr=print, sfx='b'):
    """engine vs the K2 rows (w1_pts_EE_J<J>_k2_<sfx>, l9_analysis.pred form, K = 2) above the fold, entrywise in the row's
    structure/shadow basis: m^2 Re[(D(nu) - D(nu0))_engine - (same)_K2], D = log(state/reference) of the M = I b^-T entries."""
    import l9_analysis as L9, hc_nodes as HN
    d = pickle.load(open(os.path.join(C, 'w1_pts_EE_J%d_k2_%s.pkl' % (J, sfx)), 'rb'))
    nus, Ist, nur, Iref = cal['nus'], cal['Ist'], cal['nur'], cal['Iref']
    m = cal['m']
    rows = []
    for nu in sorted(set(k[1] for k in d if k[0] == DST)):
        if (DREF, nu) not in d:
            continue
        i = int(np.argmin(np.abs(nus - nu))); j = int(np.argmin(np.abs(nur - nu)))
        if abs(nus[i] - nu) > 2e-3 or abs(nur[j] - nu) > 2e-3:
            continue
        rs, rr = d[(DST, nu)], d[(DREF, nu)]
        (bs, _, _), _ = HN.shadow(DST, J, rs.get('shadow_resid')); (br, _, _), _ = HN.shadow(DREF, J, rr.get('shadow_resid'))
        Ms = Ist[i] @ np.linalg.inv(bs).T; Mr = Iref[j] @ np.linalg.inv(br).T
        De = np.log(Ms/Mr)
        a = [rs[k] for k in ('M0', 'M2', 'M4', 'M6')]; r = [rr[k] for k in ('M0', 'M2', 'M4', 'M6')]
        Dk = np.log(L9.pred(a, r, 2)) + (rs['logC'] - rs['logK']) - (rr['logC'] - rr['logK'])
        rows.append((nu, De, Dk))
    if len(rows) < 2:
        pr(f'# K2 check J={J}: fewer than 2 common nodes'); return
    nu0, De0, Dk0 = rows[0]
    pr(f'# K2 check J={J} (w1 k2_{sfx} rows, pred K=2) vs engine: m^2 Re[(D(nu) - D({nu0:.2f}))_eng - (same)_K2] per entry (a,c), c2')
    for nu, De, Dk in rows[1:]:
        dd = m*m*((De - De0) - (Dk - Dk0)).real
        pr(f'    nu={nu:7.2f} (x={nu/m:.3f}): ' + '  '.join(f'({a},{c}) {dd[a, c]:+8.1f}' for a in range(2) for c in range(2))
           + f'   | engine m^2 Re dD(++-ish 00) {m*m*(De - De0)[0, 0].real:+9.1f}')


def records(J, cal, nucut=NUCUT):
    nus, Ist, nur, Iref = cal['nus'], cal['Ist'], cal['nur'], cal['Iref']
    R = cal['R']; m = cal['m']
    out = {}
    for i, v in enumerate(nus):
        if v >= nucut:
            continue
        j = _near(nur, v)
        repr_ = np.array([cal['kr_r'][r]*Iref[j] for r in range(R)])
        reps_ = np.array([cal['kr_r'][r]*cal['rho_r'][r]*Ist[i] for r in range(R)])
        for D, reps, mm in ((DST, reps_, m), (DREF, repr_, mphys(DREF, J))):
            out[(D, round(float(v), 9))] = dict(E=reps.mean(0), err=0*reps.mean(0), reps=reps, logN=0j, Es=1.0 + 0j, errs=0.0,
                                                 repsS=np.ones(R, complex), logNs=0j, nu=float(v), Delta=D, J=J, m=mm,
                                                 method='bj closed-form engine (calibrated)', R=R)
    return out


def chain(J, recs, c3h='keep', nucut=NUCUT):
    import z3_nodes as Z3
    orig = Z3.zpts

    def zpts(kind, J_):
        P = dict(orig(kind, J_))
        if J_ != J:
            return P
        drop = [k for k in P if k[1] < nucut and k[0] in (DST, DREF)]
        if kind == 'EE' or c3h == 'drop':
            for k in drop:
                P.pop(k)
        if kind == 'EE':
            P.update(recs)
        return P
    Z3.zpts = zpts
    sys.argv = ['hc_chain.py', 'z4', '300']
    runpy.run_path(os.path.join(HERE, 'hc_chain.py'), run_name='__main__')


if __name__ == '__main__':
    mode, J = sys.argv[1], int(sys.argv[2])
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    Kr = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    depth = int(sys.argv[5]) if len(sys.argv) > 5 else Kr + 1
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    cal = calibrate(J, depth=depth, N=N, Kr=Kr, pr=pr)
    if mode == 'calib':
        try:
            k2check(J, cal, pr)
        except Exception as ex:
            pr(f'# K2 check failed: {ex}')
        open(os.path.join(C, 'bj_calib_J%d_N%d_K%d_d%d.txt' % (J, N, Kr, depth)), 'w').write('\n'.join(lines) + '\n')
    else:
        recs = records(J, cal)
        pickle.dump(recs, open(os.path.join(C, 'z2bj_pts_EE_J%d.pkl' % J), 'wb'))
        chain(J, recs, os.environ.get('BJ_C3H', 'keep'))
