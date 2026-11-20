"""
z3 (EE_DG300_SMALLNU.md): the z2 crystal points (state Delta = 336 and reference 2e4 on common random numbers) turned into
  (a) a SCALAR VALIDATION against the exact Mellin object (d_D - d_ref must vanish; n4_exact.record, own cache cache/z3_exact.pkl --
      the project's n4_exact.pkl is not written), with jackknife errors over the 8 common-seed replicas;
  (b) chain records in the K2-record format of the l9/q3/w2 chains (M0 = direct all-orders matrix E binv^T, M2 = M4 = M6 = 0,
      logC = logN, logK = logK of the l9/w1 shadow cache), for the state AND the reference (both direct), with three estimators:
        'dir'  : plain state/reference ratio (CRN);
        'cv'   : scalar control variate  G -> G x [(s_st/s_ref)_exact / (s_st/s_ref)_MC]   (the exact scalar ratio is known for even J);
        replica-jackknife variants ('dir', i) / ('cv', i), i = 0..R-1, for error propagation through the chain;
  (c) per-node relative errors of the helicity-diagonal ratios.
Usage: python3 z3_nodes.py validate [J ...]   -> cache/z3_validate.txt
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
EXPK = os.path.join(C, 'z3_exact.pkl')


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def zpts(kind, J):
    f = os.path.join(C, 'z2_pts_%s_J%d.pkl' % (kind, J))
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def exact_rec(D, J):
    db = pickle.load(open(EXPK, 'rb')) if os.path.exists(EXPK) else {}
    key = (round(float(D), 9), int(J))
    if key not in db:
        import n4_exact as E4
        db = pickle.load(open(EXPK, 'rb')) if os.path.exists(EXPK) else {}
        db[key] = E4.record(D, J)
        tmp = EXPK + '.tmp%d' % os.getpid()
        pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, EXPK)
    return db[key]


def log_q(D, J):
    sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
    from exactmom import log_q as lq
    return lq(D, J, 3.0)


def lx_exact(D, J, nu):
    import n3_exact as NX
    r = exact_rec(D, J)
    return log_q(D, J) + np.log(complex(NX.X(nu, r)[0])), r


def logbeta(D, J):
    import l1_crystal as Cr
    return Cr.log_beta(D, J)


def nodes(J, kind='EE'):
    P = zpts(kind, J)
    out = {}
    for (D, nu), r in P.items():
        out.setdefault(nu, {})[D] = r
    return {nu: d for nu, d in sorted(out.items()) if DREF in d and len(d) >= 2}


def _jk(f, R):
    v0 = f(np.ones(R, bool))
    jk = np.array([f(np.arange(R) != i) for i in range(R)])
    return v0, jk, np.sqrt((R - 1)/R*np.sum(np.abs(jk - jk.mean(0))**2, axis=0))


def scalar_check(J, nu, dd, Dst):
    """Re/Im of (d_st - d_ref) with jackknife errors and Delta c2; returns dict incl. the complex CV factor exact/MC (all reps and JK)."""
    st, rf = dd[Dst], dd[DREF]
    R = len(st['repsS'])
    lxs, rs = lx_exact(Dst, J, nu); lxr, _ = lx_exact(DREF, J, nu)
    off = (st['logNs'] - logbeta(Dst, J)) - (rf['logNs'] - logbeta(DREF, J))
    def lr(sel):
        return np.log(st['repsS'][sel].mean()/rf['repsS'][sel].mean()) + off
    v0, jk, err = _jk(lr, R)
    dm = v0 - (lxs - lxr)
    dm = complex(dm.real, (dm.imag + np.pi) % (2*np.pi) - np.pi)
    m = mphys(Dst, J); x = nu/m; E2 = (2 - x*x)/(1 - x*x)**2
    ex2 = float(np.interp(nu, rs['nus'], rs['ex'][2]))
    fac = m*m*ex2/E2
    jkre = jk.real; jkim = (jk.imag - v0.imag + np.pi) % (2*np.pi) - np.pi
    e = complex(np.sqrt((R - 1)/R*np.sum((jkre - jkre.mean())**2)), np.sqrt((R - 1)/R*np.sum((jkim - jkim.mean())**2)))
    # CV factor: G_cv = G_dir x exp((lxs - lxr) - lr): the MC scalar ratio replaced by the exact one (per jackknife sample too)
    cvf = np.exp((lxs - lxr) - v0); cvf_jk = np.exp((lxs - lxr) - jk)
    return dict(d=dm, err=e, dc2=fac*dm.real, ec2=abs(fac)*e.real, m=m, cv=cvf, cv_jk=cvf_jk, R=R, mth=st.get('method', ''))


def validate(Js, pr):
    pr("# z3 scalar validation: d_D - d_ref = [log Y_D - log(q_D X_D)] - [same at Delta_ref = 2e4]  (must be 0; X exact, n3_exact/n4_exact);"
       " jackknife over the common-seed replicas;  Delta c2 = m^2 Re(d - d_ref) ex2/E2")
    res = {}
    for J in Js:
        for nu, dd in nodes(J).items():
            for Dst in sorted(D for D in dd if D != DREF):
                try:
                    o = scalar_check(J, nu, dd, Dst)
                except Exception as ex:
                    pr(f"  J={J} nu={nu:.4f} Delta={Dst:g}: exact object failed ({ex})"); continue
                res[(J, nu, Dst)] = o
                pr(f"  J={J:2d} nu={nu:7.4f} x={nu/o['m']:.4f} Delta={Dst:g} [{o['mth'][:10]}]: Re {o['d'].real:+.5f} +- {o['err'].real:.5f}, "
                   f"Im {o['d'].imag:+.5f} +- {o['err'].imag:.5f};  Delta c2 = {o['dc2']:+8.1f} +- {o['ec2']:6.1f}  "
                   f"({abs(o['dc2'])/max(o['ec2'], 1e-9):.1f} sigma)")
    return res


# ------------------------------------------------------------------------------------------------ chain records
_SH = {}


def shadow(D, J):
    if not _SH:
        for f in sorted(glob.glob(os.path.join(C, '*shadow_k2*.pkl'))):
            _SH.update(pickle.load(open(f, 'rb')))
    return _SH[(float(D), int(J))]


def to_record(r, D, J, fac=1.0, sel=None):
    """n1-format point -> K2-record format (M0 = mean over replicas sel of E binv^T, times fac (complex scalar))."""
    b, logK = shadow(D, J)[0], shadow(D, J)[1]
    binv = np.linalg.inv(b)
    reps = r['reps'] if sel is None else r['reps'][sel]
    na = r['E'].shape[0]
    H = np.array([rep[:na] @ binv.T for rep in reps]).mean(0)
    Z = 0*H
    return dict(M0=H*fac, M2=Z, M4=Z, M6=Z, logC=complex(r['logN']), logK=complex(logK), direct=True)


def records(kind, J, est='dir', jk=None, numin=0.0, numax=1e9, Ds=None):
    """chain records {(D, nu): rec} for all z2 nodes of (kind, J) with numin <= nu <= numax; est 'dir' | 'cv'; jk: replica left out."""
    out = {}
    P = zpts(kind, J)
    by = {}
    for (D, nu), r in P.items():
        by.setdefault(nu, {})[D] = r
    for nu, dd in by.items():
        if DREF not in dd or not (numin <= nu <= numax):
            continue
        R = len(dd[DREF]['reps'])
        sel = None if jk is None else (np.arange(R) != jk)
        for D, r in dd.items():
            if Ds is not None and D not in Ds and D != DREF:
                continue
            fac = 1.0
            if est == 'cv' and D != DREF:
                ee = zpts('EE', J)          # the scalar on the common samples is stored with every kind; use this kind's own
                o = scalar_check(J, nu, dd, D)
                fac = o['cv'] if jk is None else o['cv_jk'][jk]
            out[(float(D), round(nu, 9))] = to_record(r, D, J, fac, sel)
    return out


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    if sys.argv[1] == 'validate':
        Js = [int(a) for a in sys.argv[2:]] or [24, 28]
        validate(Js, pr)
        open(os.path.join(C, 'z3_validate.txt'), 'w').write("\n".join(lines) + "\n")
