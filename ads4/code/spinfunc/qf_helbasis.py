"""
qf_helbasis (RESULT_Q_FWD_CLOSE.md): copy of bh_helbasis (BELOW_FOLD3.md §1.2) with its own cache cache/qf_helbasis.pkl (bh_helbasis.pkl
is read, never written).  T0(Delta): the project's order-0 helicity basis (l8_analysis.basis rule: rank-1 decomposition over x of the
order-0 crystal matrix of the C^_2 action, ++ = flattest d/W_h), from fresh l6_spin_action.action points, x = 0.15 ... 0.8.
Usage: L1_DEG=6 python3 qf_helbasis.py J Delta [Delta ...]
"""
import os, sys, pickle, time
os.environ.setdefault('L1_DEG', '6')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

PK = os.path.join(HERE, 'cache', 'qf_helbasis.pkl')
PK_BH = os.path.join(HERE, 'cache', 'bh_helbasis.pkl')
XS = (0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def load():
    db = pickle.load(open(PK_BH, 'rb')) if os.path.exists(PK_BH) else {}
    if os.path.exists(PK):
        db.update(pickle.load(open(PK, 'rb')))
    return db


def records(J, D, xs=XS):
    import l6_spin_action as L6
    m = mphys(D, J)
    for x in xs:
        key = (J, round(float(D), 6), x)
        if key in load():
            continue
        r = L6.action(x*m, D, J)
        own = pickle.load(open(PK, 'rb')) if os.path.exists(PK) else {}
        own[key] = dict(M0=r['M'][:, :, 0], M=r['M'], logC=r['logC'], logK=r['logK'], slog=r['scal_log'], sS=r['scal_S'], m=m, x=x,
                        b=r['b'], resid=r['shadow_resid'])
        pickle.dump(own, open(PK + '.tmp', 'wb')); os.replace(PK + '.tmp', PK)
    return {k[2]: v for k, v in load().items() if k[0] == J and k[1] == round(float(D), 6)}


def basis(J, D, xs_min=0.15):
    import l7_dc2_analysis as An
    from l8_analysis import realvec
    R = records(J, D)
    pts = sorted((x, R[x]) for x in R if x >= xs_min)
    xs = np.array([p[0] for p in pts])
    sv, vs = An.rank1([p[1]['M0'] for p in pts])
    vv = [realvec(v) for v in vs]
    T = np.array([v[0] for v in vv]).T
    Ti = np.linalg.inv(T)
    flat = []
    for j in range(2):
        q = np.array([(Ti @ p[1]['M0'] @ Ti.T)[j, j]*np.exp(p[1]['logC'] - p[1]['logK'] - p[1]['slog'])/p[1]['sS'][0] for p in pts])
        q = q/An.W('++', J, xs)
        flat.append(np.std(np.abs(q/q[0])))
    if flat[0] > flat[1]:
        T = T[:, ::-1]; Ti = np.linalg.inv(T)
    offd = max(abs((Ti @ p[1]['M0'] @ Ti.T)[0, 1])/np.sqrt(abs((Ti @ p[1]['M0'] @ Ti.T)[0, 0]*(Ti @ p[1]['M0'] @ Ti.T)[1, 1]))
               for p in pts)
    return dict(T=T, Ti=Ti, sv3=sv[2], imag=max(v[1] for v in vv), flat=flat, offdiag0=offd, xs=xs)


if __name__ == '__main__':
    J = int(sys.argv[1])
    for D in [float(v) for v in sys.argv[2:]]:
        t = time.time()
        B = basis(J, D)
        print(f'J={J} Delta={D}: T (cols ++, +-) = {B["T"].round(6).tolist()}  sv3/sv1 = {B["sv3"]:.2e}  imag = {B["imag"]:.1e}  '
              f'flatness (col0, col1) = {B["flat"][0]:.3e}, {B["flat"][1]:.3e}  max order-0 offdiag = {B["offdiag0"]:.1e}  [{time.time() - t:.0f}s]',
              flush=True)
