"""
z2 (EE_DG300_SMALLNU.md): EE graviton (C^'_EE, kind 'EE' = w0 integrand) and C^_3 (kind 'C3H', o_H q) crystal integrals, with the
scalar on the same samples, at SMALL nu (5-35) for the Delta_g = 300 near-gap rows at large J (24, 28, 32) and their heavy references,
STATE AND REFERENCE ON COMMON RANDOM NUMBERS (proposal adapted on the state and re-used for every reference; same Sobol seeds):
  * nu >= NU_SPLIT (above the fold nu* ~ 11.3-11.5, located with z1_fold): the n3 hybrid wall-free contour C(rho1 = 3.5, nu_f)
    (n3_run rule: nu_f = max(nu, 15) for nu >= 13, nu_f = 13 below), graviton-adapted Gaussian proposal (n3_run.adapt_g);
  * nu <  NU_SPLIT (below the fold): the A+ + A- Lefschetz-thimble sum of n4 (n4_run.run_pair: already CRN; saddles tracked from
    nu_c = 10 or lower).
Project modules are imported read-only; only n1_pw.P9.build_spin is wrapped in-process so that kind 'EE' evaluates the w0 integrand.
Stored per point (key (Delta, round(nu, 9))) in n1_run format (E, err, reps, logN; Es, errs, repsS, logNs; ljc = True):
  cache/z2_pts_<kind>_J<J>.pkl ; raw per call: cache/z2_raw.pkl (key (method, Delta, J, round(nu, 9))).
Usage: OMP_NUM_THREADS=1 python3 z2_points.py J nu1,nu2,... Dstate[,Dref1,Dref2...]  [env Z2_M (log2 pts/replica), Z2_R]
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
KINDS = tuple(os.environ.get('Z2_KINDS', 'EE,C3H').split(','))
NU_SPLIT = float(os.environ.get('Z2_NUSPLIT', '11.0'))
MQ = int(os.environ.get('Z2_M', '10'))
RQ = int(os.environ.get('Z2_R', '8'))
SEED = int(os.environ.get('Z2_SEED', '4242'))


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


_P = [False]


QSC_SRC_OLD = "        ZW = dot(Z, W); ZZWW = T.mul(dot(Z, Z), dot(W, W))\n"
QSC_SRC_NEW = "        ZW = dot(Z, W)/QSC[0]; ZZWW = T.mul(dot(Z, Z), dot(W, W))/(QSC[0]*QSC[0])\n"


def patch_qscale():
    """The scalar prefactor Q_J = |Z|^J |W|^J P_J(cos) of l1_crystal.build overflows (inf) at Delta_ref = 2e4 for J >~ 20 (|Z||W| ~ 1e13
    there), which masks EVERY contour point (the scalar ok-mask also gates the graviton kinds).  Fix (point-wise copy only, in-process):
    Q is evaluated with ZW -> ZW/s, ZZWW -> ZZWW/s^2 for a per-block constant s (set by set_qsc at the block's saddle), i.e. the scalar
    integrand is multiplied by the constant s^-J; the stored scalar logNs then gets + J log s ('lqsc').  Graviton kinds untouched."""
    import n1_pw as PW
    if 'QSC' in PW.Cr.__dict__:
        return
    src = open(os.path.join(HERE, 'l1_crystal.py')).read()
    assert src.count(QSC_SRC_OLD) == 1
    src = src.replace(QSC_SRC_OLD, QSC_SRC_NEW)
    saved = sys.modules.get('l1_tps')
    sys.modules['l1_tps'] = PW
    try:
        PW.Cr.__dict__['QSC'] = [1.0]
        exec(compile(src, 'l1_crystal_pw_qsc', 'exec'), PW.Cr.__dict__)
    finally:
        sys.modules['l1_tps'] = saved


def set_qsc(Delta, J, c, wc, L):
    """s := sqrt|Z^2 W^2| at the block's saddle (xi = 0); returns J log s."""
    import n1_pw as PW
    import n1_direct as D
    PW.Cr.QSC[0] = 1.0
    if J == 0:
        return 0.0
    Wp = wc[None, :].astype(complex)
    PW.set_batch(1)
    zl, lj = D.zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
    zz = D._zzww(zl, M=PW.Cr)
    s = float(np.sqrt(np.abs(np.asarray(zz).ravel()[0])))
    PW.Cr.QSC[0] = s
    return J*np.log(s)


def patch():
    if _P[0]:
        return
    patch_qscale()
    import n1_pw as PW
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(PW, PW.S, PW.Cr, PW.A6)
    orig = PW.P9.build_spin

    def bs(zl, nu, Delta, J, mons, sol, solt, kind):
        if kind == 'EE':
            return bEE(zl, nu, Delta, J, mons, sol, solt)
        return orig(zl, nu, Delta, J, mons, sol, solt, kind)
    PW.P9.build_spin = bs
    _P[0] = True


def _save(path, key, val):
    import fcntl
    with open(path + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        db = pickle.load(open(path, 'rb')) if os.path.exists(path) else {}
        db[key] = val
        tmp = path + '.tmp%d' % os.getpid()
        pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, path)
        fcntl.flock(lk, fcntl.LOCK_UN)


def fn(kind, J):
    return os.path.join(C, 'z2_pts_%s_J%d.pkl' % (kind, J))


def load(kind, J):
    f = fn(kind, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


# ------------------------------------------------------------------------------------------------ n3 contour, CRN
def n3_block(J, Delta, nu, prop, m, R, seed0=SEED, rho1=3.5):
    import n3_hybrid as Hy
    import n3_flow as F
    import n3_run as NR
    import l1_crystal as Cr
    nuf = max(nu, F.NUF) if nu >= NR.NUF_LOW else NR.NUF_LOW
    c, wc, L, _ = F.saddle(nuf, Delta, J)
    lqsc = set_qsc(Delta, J, c, wc, L)
    evS = F.Evaluator(nu, Delta, J, c, wc, L, 'S', nuf=nuf); F.set_branch(evS); Hy.attach_nuflow(evS)
    evs = {k: F.Evaluator(nu, Delta, J, c, wc, L, k, nuf=nuf) for k in KINDS}
    for e in evs.values():
        F.set_branch(e)
    t0 = time.time()
    if prop is None:
        prop = NR.adapt_g(evS, evs, rho1, m=min(m, 10))
    tad = time.time() - t0
    mu, Ch, ldC = prop
    valsS = []; valsG = {k: [] for k in KINDS}
    diag = dict(nbad=0, nbad_g=0, wind=0.0, imdrift=0.0, sdrift=0.0, n=0, nouter=0)
    for r in range(R):
        eta = Hy._etas(m, seed0 + r); eta = np.concatenate([eta, -eta])
        accS = 0.0; accG = {k: None for k in KINDS}
        for c0 in range(0, eta.shape[0], Hy.CHUNK):
            et = eta[c0:c0 + Hy.CHUNK]; Z = mu[None, :] + et @ Ch.T
            wgt = np.exp(0.5*ldC + 0.5*np.sum(et*et, axis=1))
            ct = Hy.contour(evS, Z, rho1, keep_path=True)
            accS = accS + np.sum(wgt*ct['Fval'])
            d = ct['diag']; diag['nbad'] += int((~ct['ok']).sum()); diag['n'] += Z.shape[0]; diag['nouter'] += d['nouter']
            for kk in ('wind', 'imdrift', 'sdrift'):
                diag[kk] = max(diag[kk], d[kk])
            for k in KINDS:
                Fg, Pgc, pref0, okg = Hy.grav_on_paths(evs[k], ct, Z.shape[0])
                v = np.sum(wgt*Fg, axis=-1)
                accG[k] = v if accG[k] is None else accG[k] + v
                diag['nbad_g'] += int((~okg).sum())
        valsS.append(accS/eta.shape[0])
        for k in KINDS:
            valsG[k].append(accG[k]/eta.shape[0])
    valsS = np.array(valsS)
    (Pc, _), _, _ = evS.run(evS.pint, np.zeros((1, 7), complex), None)
    logNs = Pc[0] + np.log(np.linalg.det(L))
    out = {}
    for k in KINDS:
        vg = np.array(valsG[k])
        (Pg, _), _, _ = evs[k].run(evs[k].grav, np.zeros((1, 7), complex), None)
        out[k] = dict(ljc=True, E=vg.mean(axis=0), err=vg.std(axis=0, ddof=1)/np.sqrt(R), reps=vg, logN=Pg[0] + np.log(np.linalg.det(L)),
                      Es=valsS.mean(), errs=valsS.std(ddof=1)/np.sqrt(R), repsS=valsS, logNs=logNs + lqsc, lqsc=lqsc, nu=nu, Delta=Delta, J=J,
                      m=mphys(Delta, J), nuf=nuf, rho1=rho1, mq=m, R=R, diag=dict(diag), method='n3_hybrid CRN', seed0=seed0)
    raw = dict(method='n3', nu=nu, Delta=Delta, J=J, nuf=nuf, prop=prop, logbeta=Cr.log_beta(Delta, J), logNs=logNs + lqsc, lqsc=lqsc,
               repsS=valsS, diag=diag, tadapt=tad, kinds=KINDS, mq=m, R=R)
    return out, raw, prop


def n3_pair(J, nu, Ds, m=MQ, R=RQ):
    prop = None
    res = {}
    for D in Ds:
        t = time.time()
        out, raw, prop = n3_block(J, D, nu, prop, m, R)
        res[D] = (out, raw)
        o = out[KINDS[0]]
        rel = {k: float(np.max(np.abs(out[k]['err']))/np.max(np.abs(out[k]['E']))) for k in KINDS}
        print(f"  n3 J={J} Delta={D:g} nu={nu:.4f} (nu_f {raw['nuf']}): scalar rel err {abs(o['errs']/o['Es']):.1e}; max rel err {rel}; "
              f"nbad {o['diag']['nbad']}/{o['diag']['nbad_g']} of {o['diag']['n']}; {time.time() - t:.0f}s", flush=True)
    return res


# ------------------------------------------------------------------------------------------------ n4 thimbles (below the fold)
def _patch_n4():
    import n4_thimble as N4
    if getattr(N4, '_z2', False):
        return
    orig = N4.thimble

    def th(rec, Delta, J, kinds=(), **kw):
        lq = set_qsc(Delta, J, rec['c'], rec['wc'], rec['L'])
        o = orig(rec, Delta, J, kinds, **kw)
        o['S']['logN'] = o['S']['logN'] + lq; o['S']['logNs'] = o['S']['logNs'] + lq; o['S']['lqsc'] = lq
        return o
    N4.thimble = th
    N4._z2 = True


def n4_pair(J, nu, Ds, m=MQ, R=RQ, more=()):
    import n4_run as NRun
    _patch_n4()
    NRun.SEED = SEED
    res = None
    for nuc in (10.5, 10.0, 9.5, 9.0, 8.0):
        if nuc > 11.0 or nuc < nu - 1e-9:
            continue
        if nuc != NRun.NUC:
            NRun._TR.clear()
        NRun.NUC = nuc
        try:
            res = NRun.run_pair(J, nu, Ds, KINDS, m=m, R=R, more=more)
            break
        except RuntimeError as ex:
            print(f"  n4 J={J} nu={nu}: nu_c = {nuc}: {ex}; retrying lower", flush=True)
    if res is None:
        raise RuntimeError('n4 failed')
    out = {}
    for D, raw in res.items():
        raw['nuc_used'] = NRun.NUC
        o = {k: NRun.combine_kind(raw, k) for k in KINDS}
        for k in KINDS:
            o[k]['method'] = 'n4_thimble CRN (A+ + A-)'
        out[D] = (o, dict(method='n4', **raw))
    return out


def run(J, nus, Ds, m=MQ, R=RQ):
    patch()
    for nu in nus:
        keys = [(D, round(nu, 9)) for D in Ds]
        if all(all(k in load(kd, J) for kd in KINDS) for k in keys):
            print(f"# J={J} nu={nu}: done", flush=True); continue
        t = time.time()
        if nu >= NU_SPLIT:
            res = n3_pair(J, nu, Ds, m, R)
        else:
            res = n4_pair(J, nu, Ds, m, R, more=[v for v in nus if v < NU_SPLIT])
        for D, (o, raw) in res.items():
            _save(os.path.join(C, 'z2_raw.pkl'), (raw['method'], round(D, 9), J, round(nu, 9)), raw)
            for k in KINDS:
                _save(fn(k, J), (D, round(nu, 9)), o[k])
        print(f"## J={J} nu={nu:.4f} Ds={Ds}: done in {time.time() - t:.0f}s", flush=True)


if __name__ == '__main__':
    J = int(sys.argv[1])
    nus = [float(v) for v in sys.argv[2].split(',')]
    Ds = [float(v) for v in sys.argv[3].split(',')]
    if len(Ds) == 1:
        Ds.append(20000.0)
    run(J, nus, Ds)
