"""
y4 (EE_DG100.md): non-perturbative crystal values of C^_EE on the Delta = 100 row (Delta_g = 100), state + reference (Delta_ref = 2e4, same
nu), with the project's direct integrators imported read-only and only the graviton builder patched in-process (w0 integrand):
  n3even J..  : n3_run.point (hybrid Picard-Lefschetz contour), x = 0.15 (2^11 x 2 x 8 pts), 0.1 (2^12 x 2 x 8), kind 'EE'
  n3odd  J..  : n3_run_odd.point (odd +-, p8 P8.build), x = 0.15, 0.1 (2^12 x 2 x 8), kind 'EE'
  n1     J..  : n1_run.point (re-centred real contour, RQMC), x = 0.2, 0.25, 0.3, 0.4 ... 0.8, kind 'EE'
  p8     J..  : p8_odd_pw.point (odd), x = 0.2, 0.25, 0.3, 0.4 ... 0.8, kind 'EE' (Delta = 100 row)
Files (n1_run / p8 formats, key (Delta, round(nu, 9))): cache/y4_n3_pts_EE_J<J>.pkl, y4_n3_odd_pts_EE_J<J>.pkl, y4_n1_pts_EE_J<J>.pkl,
cache/y4_p8_pts_EE_J<J>.pkl.
Usage: OMP_NUM_THREADS=1 python3 y4_direct.py n3even|n3odd|n1|p8 J [J ...]
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
XS_HI = [float(v) for v in os.environ.get('Y4_XS', '0.2,0.25,0.3').split(',')]   # x >= 0.4: K2 (accurate there, H0 §24.4)


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def fn(tag, J):
    return os.path.join(C, 'y4_%s_pts_EE_J%d.pkl' % (tag, J))


def load(tag, J):
    f = fn(tag, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def store(tag, J, key, r):
    db = load(tag, J); db[key] = r
    tmp = fn(tag, J) + '.tmp%d' % os.getpid()
    pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, fn(tag, J))


def patch():
    import y5_n4 as Y5
    Y5.patch()


if __name__ == '__main__':
    mode = sys.argv[1]; Js = [int(a) for a in sys.argv[2:]]
    patch()
    if mode in ('n3even', 'n3odd'):
        import n3_run as N3
        import n3_run_odd as N3O
        mod = N3 if mode == 'n3even' else N3O
        mod.KINDS = ('EE',)
        tag = 'n3' if mode == 'n3even' else 'n3_odd'
        for J in Js:
            for x, mq in ((0.15, 11 if mode == 'n3even' else 12), (0.1, 12)):
                for Dl in (100.0, DREF):
                    nu = x*mphys(100.0, J); key = (Dl, round(nu, 9))
                    if key in load(tag, J):
                        continue
                    t = time.time()
                    res = mod.point(J, Dl, nu, m=mq)['EE']
                    store(tag, J, key, res)
                    rel = float(np.max(np.abs(res['err']))/np.max(np.abs(res['E'])))
                    print(f"{mode} EE J={J} Delta={Dl:.0f} x={x} nu={nu:.4f}: max rel err {rel:.2e}, scalar rel err "
                          f"{abs(res['errs']/res['Es']):.1e}, {time.time() - t:.0f}s", flush=True)
    elif mode in ('n1', 'p8'):
        import n1_run as N1
        import p8_odd_pw as P8M
        tag = mode
        for J in Js:
            for x in XS_HI:
                for Dl in (100.0, DREF):
                    nu = x*mphys(100.0, J); key = (Dl, round(nu, 9))
                    if key in load(tag, J):
                        continue
                    t = time.time()
                    if mode == 'n1':
                        r = N1.point('EE', J, Dl, nu)
                    else:
                        r = P8M.point('EE', J, Dl, nu, mq=(P8M.MQREF if Dl >= 1e4 else P8M.MQ))
                        r.update(row=100.0, x=x)
                    store(tag, J, key, r)
                    rel = float(np.max(np.abs(r['err']))/np.max(np.abs(r['E'])))
                    print(f"{mode} EE J={J} Delta={Dl:.0f} x={x} nu={nu:.4f}: max rel err {rel:.2e}, {time.time() - t:.0f}s", flush=True)
