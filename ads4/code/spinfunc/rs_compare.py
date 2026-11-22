"""
rs_compare (HRS_CHECK.md): crystal-cycle (rs_points) vs analytic-radial engine (bg_eeq caches, read-only) for the Result Q EE k = 4
image, node by node, plus the scalar validation against the exact Mellin skeleton (records from bf_exact.load(), read-only; missing
records are built with n4_exact.record into the own cache cache/rs_exact.pkl).

Per node nu (state Delta, reference 2e4, common random numbers, R replicas, jackknife errors):
  scalar:  d_s = log rho_s(MC) - log rho_s(exact)   (rho = state/reference ratio incl. logN - log beta)    -> must vanish
           Delta c2 = m^2 Re(d_s) ex2/E2  (z3_nodes convention)
  EE4:     d = log rho(crystal) - log rho(engine)   -> must be a nu-independent constant kappa (the engine drops Delta-only constants)
           per-harmonic Gamma difference  dG(nu) = Re[m^2 T_st(nu) (d(nu) - kappa)] / sum_j |R_j(nu)|   (c2 per harmonic, as BELOW_FOLD2 §4.2)
           estimators: 'dir' (plain CRN ratio) and 'cv' (scalar control variate: x exact/MC scalar ratio).
  kappa: fitted (weighted, complex) over the nodes given, or fixed at a node (RS_KAPPA_NU).
Usage: python3 rs_compare.py J Dstate tag1[,tag2,...] [kind]     -> stdout (+ cache/rs_compare_D<D>_J<J>.txt)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
EXPK = os.path.join(C, 'rs_exact.pkl')


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def exact_rec(D, J):
    import bf_exact as BE
    key = (round(float(D), 9), int(J))
    db = BE.load()
    if key in db and np.all(np.isfinite(db[key]['ex'][2])):
        return db[key]
    own = pickle.load(open(EXPK, 'rb')) if os.path.exists(EXPK) else {}
    if key not in own:
        import n4_exact as N4
        own[key] = N4.record(float(D), int(J))
        pickle.dump(own, open(EXPK + '.tmp', 'wb')); os.replace(EXPK + '.tmp', EXPK)
    return own[key]


def lx_exact(D, J, nu):
    import n3_exact as NX
    sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
    from exactmom import log_q
    r = exact_rec(D, J)
    return log_q(D, J, 3.0) + np.log(complex(NX.X(nu, r)[0])), r


def logbeta(D, J):
    import l1_crystal as Cr
    return Cr.log_beta(D, J)


def jk(f, R):
    v0 = f(np.ones(R, bool))
    s = np.array([f(np.arange(R) != i) for i in range(R)])
    e = np.sqrt((R - 1)/R*np.sum(np.abs(s - s.mean(0))**2, axis=0))
    ere = np.sqrt((R - 1)/R*np.sum((s.real - s.real.mean(0))**2, axis=0))
    eim = np.sqrt((R - 1)/R*np.sum((s.imag - s.imag.mean(0))**2, axis=0))
    return v0, s, ere + 1j*eim


def wrap(z):
    return complex(z.real, (z.imag + np.pi) % (2*np.pi) - np.pi)


def load_pts(kind, J, tags, Dst):
    P = {}
    for t in tags:
        f = os.path.join(C, 'rs_pts_%s_J%d_D%g%s.pkl' % (kind, J, Dst, t))
        if os.path.exists(f):
            P.update(pickle.load(open(f, 'rb')))
    return P


EKIND = {'EE4': 'EE4', 'EE': 'EE2'}           # crystal kind -> engine kind (k = 4 / k = 2)


def engine(Dst, J, kind='EE4'):
    import rs_engine as RE
    db = RE.get(Dst, J)                       # bg_eeq cache if present (read-only), else rs_eng_*
    r = RE.gamma_K(db, '++', kind=EKIND[kind])
    r3 = RE.gamma_K(db, '++', kind=EKIND[kind], Kt=3)
    return db, r, r3


def node_compare(J, nu, st, rf, Dst, eng):
    """returns dict with scalar check and EE4 log-ratio difference (dir and cv) with jackknife errors."""
    R = len(st['repsS'])
    m = mphys(Dst, J)
    lxs, rs = lx_exact(Dst, J, nu); lxr, _ = lx_exact(DREF, J, nu)
    offs = (st['logNs'] - logbeta(Dst, J)) - (rf['logNs'] - logbeta(DREF, J))
    fs = lambda sel: np.log(st['repsS'][sel].mean()/rf['repsS'][sel].mean()) + offs
    vs, sjk, es = jk(fs, R)
    ds = wrap(vs - (lxs - lxr))
    x = nu/m; E2 = (2 - x*x)/(1 - x*x)**2
    ex2 = float(np.interp(nu, rs['nus'], rs['ex'][2]))
    fac = m*m*ex2/E2
    # engine ratio at nu (grid point)
    db, r, r3 = eng
    i = int(np.argmin(np.abs(r['nus'] - nu)))
    lr_ = np.log(r['ratio'])           # off-grid nodes: linear interpolation of the (smooth) engine log-ratio
    le = complex(np.interp(nu, r['nus'], lr_.real), np.interp(nu, r['nus'], np.unwrap(lr_.imag)))
    offg = st['logN'] - rf['logN']
    fg = lambda sel: np.log(st['reps'][sel].mean(0)[0, 0]/rf['reps'][sel].mean(0)[0, 0]) + offg
    vg, gjk, eg = jk(fg, R)
    # control variate: x (exact scalar ratio / MC scalar ratio), replica-consistent
    fcv = lambda sel: fg(sel) - (fs(sel) - (lxs - lxr))
    vc, cjk, ec = jk(fcv, R)
    out = dict(nu=nu, m=m, R=R, ds=ds, es=es, dc2s=fac*ds.real, ec2s=abs(fac)*es.real, fac=fac,
               d_dir=wrap(vg - le), e_dir=eg, d_cv=wrap(vc - le), e_cv=ec,
               jkb_dir=(R - 1)*(np.mean(gjk) - vg), jkb_cv=(R - 1)*(np.mean(cjk) - vc), jk_dir=gjk - le, jk_cv=cjk - le,
               Tst=r['Tst'][i], SR=float(sum(abs(v[i]) for v in __import__('bf_qan').Rj(r['nus']).values())),
               GK=r['GK'][i].real, Gex=r['Gex'][i].real, trunc=(r['GK'][i] - r3['GK'][i]).real,
               relS=float(abs(st['errs']/st['Es'])), relG=float(np.max(np.abs(st['err']/st['E']))),
               mth=st.get('method', ''))
    return out


def fit_kappa(rows, key, nus=None):
    sel = [o for o in rows if (nus is None or any(abs(o['nu'] - v) < 1e-6 for v in nus))]
    w = np.array([1.0/max(abs(o['e_' + key])**2, 1e-30) for o in sel])
    d = np.array([o['d_' + key] for o in sel])
    ref = d[int(np.argmax(w))]
    d = np.array([ref + wrap(x - ref) for x in d])
    k = np.sum(w*d)/np.sum(w)
    return k, 1.0/np.sqrt(np.sum(w))


def main(J, Dst, tags, kind='EE4'):
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    P = load_pts(kind, J, tags, Dst)
    eng = engine(Dst, J, kind)
    nus = sorted(nu for (D, nu) in P if abs(D - Dst) < 1e-6 and (DREF, nu) in P)
    rows = [node_compare(J, nu, P[(Dst, nu)], P[(DREF, nu)], Dst, eng) for nu in nus]
    pr(f'# rs_compare: J={J} Delta={Dst} kind={kind} tags={tags}; crystal (n3/n4 thimble, CRN vs Delta_ref=2e4) vs engine (bg_eeq K=4)')
    pr('# scalar: d_s = log rho_MC - log rho_exact (must vanish); Delta c2 = m^2 Re d_s ex2/E2')
    for o in rows:
        pr(f"  nu={o['nu']:6.2f} [{o['mth'][:9]}] R={o['R']} relerr S {o['relS']:.1e} G {o['relG']:.1e}: Re d_s = {o['ds'].real:+.2e} +- {o['es'].real:.1e}, "
           f"Im {o['ds'].imag:+.2e} +- {o['es'].imag:.1e};  Delta c2 = {o['dc2s']:+9.1f} +- {o['ec2s']:7.1f} ({abs(o['dc2s'])/max(o['ec2s'],1e-12):.1f} sig)")
    for key in ('dir', 'cv'):
        kap, ek = fit_kappa(rows, key, [float(v) for v in os.environ['RS_KAPPA_NUS'].split(',')] if os.environ.get('RS_KAPPA_NUS') else None)
        pr(f'# EE4 [{key}]: kappa (fitted) = {kap.real:+.3e}{kap.imag:+.3e}i +- {ek:.1e};  per node: d - kappa (log ratio), and the per-harmonic'
           ' Gamma difference dG = Re[m^2 T_st (d - kappa)]/sum|R_j| vs the engine Gamma_K/sum|R_j|, Gamma_exact/sum|R_j|')
        for o in rows:
            dd = wrap(o['d_' + key] - kap); e = o['e_' + key]
            g = o['m']**2*o['Tst']/o['SR']
            dG = (g*dd).real; edG = np.hypot(abs(g.real)*e.real, abs(g.imag)*e.imag)
            o['dG_' + key] = dG; o['edG_' + key] = edG
            o['sys_' + key] = abs(g)*abs(o['ds'])          # scalar log-ratio residual in the same units (method-bias indicator)
            o['bias_' + key] = (g*o['jkb_' + key]).real     # jackknife estimate of the ratio-of-means bias, same units
            pr(f"  nu={o['nu']:6.2f}: Re(d-k) = {dd.real:+.2e} +- {e.real:.1e}, Im {dd.imag:+.2e} +- {e.imag:.1e};  dG = {dG:+10.2f} +- {edG:8.2f} [scal-sys {o['sys_' + key]:.1f}; jk-bias {o['bias_' + key]:+.1f}] "
               f"({abs(dG)/max(edG, 1e-12):.1f} sig) | engine GK {o['GK']/o['SR']:+7.2f} [a4 {o['trunc']/o['SR']:+.2f}] Gex {o['Gex']/o['SR']:+8.2f}; "
               f"m2|T|/SR = {abs(g):.2e}")
    def chi2(sel, key):
        c = 0.0; n = 0
        kap, _ = fit_kappa(rows, key)
        for o in rows:
            if sel(o['nu']):
                dd = wrap(o['d_' + key] - kap); e = o['e_' + key]
                c += (dd.real/e.real)**2 + (dd.imag/e.imag)**2; n += 2
        cs = 0.0; ns = 0
        for o in rows:
            if sel(o['nu']):
                cs += (o['ds'].real/o['es'].real)**2 + (o['ds'].imag/o['es'].imag)**2; ns += 2
        return c, n, cs, ns
    for lab, sel in (('nu < 9 (n4 thimbles)', lambda v: v < 9), ('nu >= 9 (n3)', lambda v: v >= 9), ('all', lambda v: True)):
        c, n, cs, ns = chi2(sel, 'dir')
        pr(f'# chi2 [{lab}]: EE4 crystal-engine (Re, Im of d - kappa) {c:.1f} / {n} values (kappa fitted on all: 2 dof); scalar check {cs:.1f} / {ns}')
    open(os.path.join(C, 'rs_compare_D%g_J%d_%s_%s.txt' % (Dst, J, kind, '_'.join(t.strip('_') or 'main' for t in tags))), 'w').write('\n'.join(lines) + '\n')
    return rows


if __name__ == '__main__':
    J, Dst = int(sys.argv[1]), float(sys.argv[2])
    tags = sys.argv[3].split(',') if len(sys.argv) > 3 else ['']
    kind = sys.argv[4] if len(sys.argv) > 4 else 'EE4'
    main(J, Dst, tags, kind)
