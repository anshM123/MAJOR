"""
w3 (EE_NEARGAP.md §2): NON-PERTURBATIVE (direct RQMC, n1_direct contour; n1_run.point unchanged) crystal values of C^_EE, as q4_direct did
for q^2: the point-wise backend n1_pw gets the EE integrand of w0 (make_build_even_EE on the point-wise module copies), kind 'EE'.
Cache: cache/w3_n1_EE_J<J>.pkl keyed (Delta, nu rounded 1e-9) (state and reference Delta_ref = 2e4 at the same nu).
Usage:  OMP_NUM_THREADS=1 python3 w3_direct.py run J Delta x1,x2,...
        W2_BASIS=q2 python3 w3_direct.py analyse "J:Delta,J:Delta,..."    -> cache/w3_direct_compare.txt, cache/w3_direct_grids.pkl
"""
import os, sys, glob, pickle, time
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
C = os.path.join(HERE, 'cache')


def wfn(J):
    return os.path.join(C, 'w3_n1_EE_J%d.pkl' % J)


def wload(J):
    f = wfn(J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def setup_pw():
    import n1_pw as PW
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(PW, PW.S, PW.Cr, PW.A6)
    orig = PW.P9.build_spin
    PW.P9.build_spin = lambda zl, nu, Delta, J, mons, sol, solt, kind: (bEE(zl, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                       orig(zl, nu, Delta, J, mons, sol, solt, kind))
    return PW


def run(J, Delta, xs):
    setup_pw()
    import n1_run as R
    m = R.mphys(Delta, J)
    for D in (Delta, 20000.0):
        for x in xs:
            nu = x*m
            key = (D, round(nu, 9))
            if key in wload(J):
                continue
            t = time.time()
            r = R.point('EE', J, D, nu)
            db = wload(J); db[key] = r
            pickle.dump(db, open(wfn(J) + '.tmp', 'wb')); os.replace(wfn(J) + '.tmp', wfn(J))
            rel = np.max(np.abs(r['err']))/np.max(np.abs(r['E']))
            print(f"EE J={J} Delta={D:.0f} x={x} nu={nu:.3f}: max rel err {rel:.1e}, nmask {r['nmask']}, {time.time() - t:.0f}s", flush=True)


def selfcheck(J=18, D=400.0, x=0.15):
    """point-wise EE integrand at the saddle vs the TPS one (constant term): same function?"""
    PW = setup_pw()
    import n1_direct as Dd
    import l5_spin as S
    import w1_rows as W1
    P = W1.setup_even()
    m = W1.mphys(D, J); nu = x*m
    mons, sol, _ = S.conserved(D, J); _, solt, _ = S.conserved(3.0 - D, J)
    zs, Lc = Dd.saddle_L(nu, D)
    PW.set_batch(1)
    lamF, G, pref = PW.P9.build_spin([np.array([zs[i]]) for i in range(7)], nu, D, J, mons, sol, solt, 'EE')
    import l1_tps as T
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF2, G2, pref2 = P.build_spin(zv, nu, D, J, mons, sol, solt, 'EE')
    a = pref[..., 0]; b = pref2[..., 0]
    print("pointwise vs TPS at the saddle: max |diff|/max|.| =", np.max(np.abs(a - b))/np.max(np.abs(b)), " lamF", lamF[0], lamF2[0],
          " G", G[0], G2[0])


def analyse(specs):
    import w2_chain as W2                  # EE-mode loaders/basis (W2_BASIS), phase strip, u_EE de-shift
    import n1_analysis as N1
    W2.set_mode(True)
    Q3 = W2.Q
    for fn in (sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl'))) + sorted(glob.glob(os.path.join(C, 'w1_shadow_k2*.pkl'))) +
               [os.path.join(C, 'u5_shadow_k2.pkl')]):
        if os.path.exists(fn):
            N1.SHK2.update(pickle.load(open(fn, 'rb')))
    N1.dload = lambda kind, J: wload(J) if kind == 'C2' else {}
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# w3: direct (RQMC) vs K2 Laplace for C^_EE (C^_3 = K2 throughout), l9 units; basis mode {W2.BASIS}")
    store = {}
    for J, D0 in specs:
        N1.DIRECT = True
        gd = N1.grid_direct(J, 2, ('C2', 'C3H'), Ds=[D0])
        N1.DIRECT = False
        g2 = N1.grid_direct(J, 2, ('C2', 'C3H'), Ds=[D0])
        N1.DIRECT = True
        if D0 not in gd:
            pr(f"  J={J} Delta={D0}: no row"); continue
        r, r2 = gd[D0], g2[D0]
        store[(J, D0)] = dict(direct=r, k2=r2)
        pr(f"  J={J} Delta={D0:.0f} (m={r['m']:.1f}, J^2/m={J*J/r['m']:.2f}): direct points used {r['ndirect']}; a_h direct/K2: " +
           ", ".join(f"{h} {r['a'][h]/r2['a'][h]:.4f}" for h in r['a']))
        pr("    x       : " + " ".join(f"{x:7.3f}" for x in r['xs']))
        m = r['m']
        for e in ('++', '+-'):
            Dh = Q3.A.Dsh(e, J, r['xs'])
            ra = m*m*((r[('C2', e)]/m**2 + Dh)*r['a'][e]/r2['a'][e] - Dh)
            pr(f"    EE  {e:3s} direct (a_K2): " + " ".join(f"{v:7.0f}" for v in ra))
            pr(f"    EE  {e:3s} K2           : " + " ".join(f"{v:7.0f}" for v in r2[('C2', e)]))
            pr(f"    EE  {e:3s} RQMC err     : " + " ".join(f"{v:7.0f}" for v in r[('C2', e, 'err')]))
        pr(f"    EE  off direct       : " + " ".join(f"{v:7.3f}" for v in r[('C2', 'off')]*np.sqrt(r['a']['++']*r['a']['+-']/(r2['a']['++']*r2['a']['+-']))))
        pr(f"    EE  off K2           : " + " ".join(f"{v:7.3f}" for v in r2[('C2', 'off')]))
    pickle.dump(store, open(os.path.join(C, 'w3_direct_grids.pkl'), 'wb'))
    open(os.path.join(C, 'w3_direct_compare.txt'), 'w').write("\n".join(lines) + "\n")


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(int(sys.argv[2]), float(sys.argv[3]), [float(v) for v in sys.argv[4].split(',')])
    elif sys.argv[1] == 'check':
        selfcheck()
    else:
        analyse([(int(a.split(':')[0]), float(a.split(':')[1])) for a in sys.argv[2].split(',')])
