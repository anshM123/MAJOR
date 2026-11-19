"""
p8 (PARITY_ODD_ROUTEA.md sect. 11): NON-PERTURBATIVE (direct randomised-QMC, H0 doc 24.2) values of the parity-odd (odd J, +-) route-A'
crystal action, kinds C2 (q^2/V^4) and C3H (o_H q/V^4), n0 = +1, for the near-gap rows at Delta_g = 100.

Point-wise copy of p3_odd_action (recipe sect. 10): p3's source is exec'd through n1_pw._load with the TPS engine swapped for point arrays
  {'l1_tps': PW, 'l1_crystal': PW.Cr, 'l5_spin': PW.S, 'l6_spin_action': PW.A6}
and three textual edits (asserted to apply exactly once): the `assert T.DEG in (6, 8)` and the `T.DEG >= 8` K3 import branch are removed
(the point-wise backend has no expansion; EXP is never used here).  p1_tto_odd / p2_shadow_odd are the REAL modules (imported before the
swap, so they keep the TPS engine): structure coefficient vectors from P1.conserved_odd (reduced TPS), shadow beta_OO from
P2.shadow_odd(Delta, J) (cache/p2_shadow.pkl, DEG 6 key).  Only P8.build (the integrand) is used.

Driver (copy of n1_run.point): contour, masks and RQMC of n1_direct (re-centred log contour, path-continued logs, antithetic scrambled
Sobol, R = 8 replicas); the structure axis is [odd] only (1x1 at odd J); the scalar Q_J is integrated on the same samples (n1_run._both).
Stored per point (resumable; cache/p8_pts_<kind>_J<J>.pkl, key (Delta, nu rounded 1e-9)):
  E (1x1), err, reps (replicas), logN = log[exp(lamF + G + logJac)(w_c) det L_w]  (so that the full integral = exp(logN) E, in the same
  normalisation as p3's exp(logC) (M0 + M2 + ...) b), Es/errs (scalar), nmask, m.
Usage: OMP_NUM_THREADS=1 L1_DEG=6 python3 p8_odd_pw.py run kind J [J ...]     env P8_ROWS (explicit state rows Delta, default 100),
       P8_MFL (extra rows with m(Delta, J) = mfl Lam(100), e.g. 1.05,1.1,1.2,1.3), P8_XS (default 0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8),
       P8_M (log2 points, state, 14), P8_MHI (state rows other than Delta = 100 at x >= 0.3; default P8_M), P8_MREF (12), P8_R (8).
       Every state point (row, nu = x m(row, J)) is followed by its reference (Delta_ref = 2e4, same nu, same kind).
       python3 p8_odd_pw.py validate   -> direct vs p3 K2 (and K3 where cached) at Delta = 400 / 1600, x in [0.3, 0.7]
"""
import os, sys, pickle, time, warnings
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
import numpy as np
import p1_tto_odd as P1                      # REAL modules first (TPS engine), cached in sys.modules before any swap
import p2_shadow_odd as P2
import n1_direct as D
import n1_pw as PW
from n1_run import _both

assert P2.T.DEG == 6, 'p8 expects the DEG-6 shadow cache (L1_DEG=6)'


def _p8_source():
    src = open(os.path.join(HERE, 'p3_odd_action.py')).read()
    edits = [("assert T.DEG in (6, 8)\n", ""),
             ("if T.DEG >= 8:                                     # K3 (doc 23.1 engine, imported read-only); shadow at K3 as in l9_points\n"
              "    import l9_tps as K9\n"
              "    EXP = K9.expand_pref6\n"
              "    T.expand_pref = lambda lamF, G, q=None: (lambda s_: (s_[0], s_[1], s_[2] + s_[3]))(K9.expand_pref6(lamF, G, q))\n"
              "else:\n"
              "    EXP = T.expand_pref\n", "EXP = None\n"),
             ("os.environ.setdefault('L1_DEG', '6')\n", "")]
    for a, b in edits:
        assert src.count(a) == 1, a[:60]
        src = src.replace(a, b)
    return src


def _load_p8():
    import types
    swaps = {'l1_tps': PW, 'l1_crystal': PW.Cr, 'l5_spin': PW.S, 'l6_spin_action': PW.A6}
    saved = {k: sys.modules.get(k) for k in swaps}
    sys.modules.update(swaps)
    try:
        M = types.ModuleType('p3_odd_action_pw'); M.__dict__['__file__'] = os.path.join(HERE, 'p3_odd_action.py')
        exec(compile(_p8_source(), 'p3_odd_action_pw', 'exec'), M.__dict__)
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    assert M.P1 is P1 and M.P2 is P2 and M.T is PW and M.S is PW.S and M.Cr is PW.Cr and M.A6 is PW.A6
    return M


P8 = _load_p8()

XS = [float(v) for v in os.environ.get('P8_XS', '0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8').split(',')]
MQ = int(os.environ.get('P8_M', '14'))
MQREF = int(os.environ.get('P8_MREF', '12'))
RQ = int(os.environ.get('P8_R', '8'))
ROWS = [float(v) for v in os.environ.get('P8_ROWS', '100').split(',') if v]
MFL = [float(v) for v in os.environ.get('P8_MFL', '').split(',') if v]     # extra rows: m(Delta, J) = mfl * Lam(Delta_g = 100)
MQHI = int(os.environ.get('P8_MHI', str(MQ)))                              # log2 points for state rows at x >= 0.3
LAM100 = 98.49873095629202                                                 # spf_common.NearGap(100).Lam
DREF = 20000.0


def delta_of(mfl, J):
    """row Delta (rounded 1e-4) with m(Delta, J) = mfl * Lam(100)."""
    M = mfl*LAM100; b = 2*J + 1
    return round((-b + np.sqrt(b*b + 4*M*M))/2 + J + 2, 4)


def rows_for(J):
    return ROWS + [delta_of(r, J) for r in MFL]


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def fn(kind, J):
    return os.path.join(HERE, 'cache', 'p8_pts_%s_J%d.pkl' % (kind, J))


def load(kind, J):
    f = fn(kind, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def point(kind, J, Delta, nu, mq=None):
    assert J % 2 == 1 and J >= 5
    c, wc, L = D.full_saddle(nu, Delta, J)
    gfun = lambda zl: P8.build(zl, nu, Delta, J, kind, 1.0)[:3]
    both = lambda zl: _both(gfun, zl, nu, Delta, J)
    mq = mq or (MQ if Delta < 1e4 else MQREF)
    E, err, nm = D.integrate_qmc(both, c, wc, L, nu, Delta, J, mq, RQ, 1.2)
    PW.set_batch(1)
    zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
    lamF, G, _ = gfun(zl)
    logN = lamF[0] + G[0] + lj[0] + np.log(np.linalg.det(L))
    na = 1
    return dict(E=E[:na], err=err[:na], Es=E[na, 0], errs=err[na, 0], logN=logN, nmask=nm, wc=wc, m=mphys(Delta, J),
                reps=D.LAST_VALS[0], mq=mq, R=RQ)


def Mphys(r, Delta, J):
    """direct full action in p3's normalisation: exp(logN - logC_p3) ... returned as (log prefactor, value/replicas) with beta^{-1}:
    full = exp(lp) * v,  v = E / beta_OO  (1x1)."""
    b, logK, _ = P2.shadow_odd(Delta, J)
    bi = 1.0/b[-1, -1]
    return r['logN'], r['E'][0, 0]*bi, np.array([rep[0, 0]*bi for rep in r['reps']]), logK


def run(kind, Js):
    for J in Js:
        for row in rows_for(J):
            for x in XS:
                for Dl in [row, DREF]:
                    nu = x*mphys(row, J)
                    key = (Dl, round(nu, 9))
                    if key in load(kind, J):
                        continue
                    t = time.time()
                    mq = MQREF if Dl >= 1e4 else (MQ if (x < 0.3 or row <= 100.0) else MQHI)
                    r = point(kind, J, Dl, nu, mq=mq)
                    r.update(row=row, x=x)
                    db = load(kind, J); db[key] = r
                    tmp = fn(kind, J) + '.tmp%d' % os.getpid()
                    pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, fn(kind, J))
                    rel = np.max(np.abs(r['err']))/np.max(np.abs(r['E']))
                    print(f"{kind} J={J} Delta={Dl:.2f} (row {row:.2f}) x={x} nu={nu:.4f}: rel err {rel:.1e}, scalar rel err "
                          f"{abs(r['errs']/r['Es']):.1e}, nmask {r['nmask']}, {time.time() - t:.0f}s", flush=True)


def validate(pts=((5, 1600.0, 0.5, 'C2'), (5, 1600.0, 0.5, 'C3H'), (5, 400.0, 0.3, 'C2'), (5, 400.0, 0.7, 'C2'), (7, 400.0, 0.4, 'C2'),
                  (9, 1600.0, 0.6, 'C3H'), (9, 400.0, 0.4, 'C3H')), mq=13):
    """direct vs p3 Laplace (K0..K2 cumulative; K3 where a _k3 cache exists) for the 1x1 odd element at large Delta (point-wise)."""
    import p4_odd_analysis as A4
    out = []
    for (J, Dl, x, kind) in pts:
        nu = x*mphys(Dl, J)
        R2 = A4.load(kind, J); st = A4.getrec(R2, Dl, nu)
        if st is None:
            out.append(f"J={J} Delta={Dl} x={x} {kind}: no p3 record"); continue
        t = time.time()
        r = point(kind, J, Dl, nu, mq=mq)
        lp, v, vr, logK = Mphys(r, Dl, J)
        dv = v*np.exp(lp - st['logC'])                      # direct, in p3 units (M0 + M2 + M4 at K2)
        sd = np.std(vr, ddof=1)/np.sqrt(len(vr))*abs(np.exp(lp - st['logC']))
        o = [st[k][-1, -1] for k in ('M0', 'M2', 'M4')]
        K0, K1, K2 = o[0], o[0] + o[1], o[0] + o[1] + o[2]
        line = (f"J={J} Delta={Dl:.0f} x={x} {kind} (nu={nu:.1f}): direct {dv:.8e} (rel RQMC err {sd/abs(dv):.1e}); "
                f"|direct/K0-1| {abs(dv/K0 - 1):.1e}, |direct/K1-1| {abs(dv/K1 - 1):.1e}, |direct/K2-1| {abs(dv/K2 - 1):.1e}; "
                f"|M4/M0| {abs(o[2]/o[0]):.1e}; nmask {r['nmask']}; {time.time() - t:.0f}s")
        print(line, flush=True); out.append(line)
    open(os.path.join(HERE, 'cache', 'p8_validate.txt'), 'w').write("\n".join(out) + "\n")


if __name__ == '__main__':
    if sys.argv[1] == 'validate':
        validate()
    else:
        run(sys.argv[2], [int(a) for a in sys.argv[3:]])
