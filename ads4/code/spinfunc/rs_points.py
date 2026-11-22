"""
rs_points (HRS_CHECK.md): independent crystal-cycle (thimble) evaluation of the Result Q EE k = 4 image for low-J states, to test H_RS
(analytic-radial engine on the real section D == crystal-cycle integral beyond leading order).

Protocol = EE_DG300_SMALLNU's z2 direct-node protocol, unchanged (z2_points imported read-only, nothing patched on disk):
  * nu >= RS_NUSPLIT: n3 hybrid contour C(rho1 = 3.5, nu_f) (nu_f = 13 below 13), CRN state / Delta_ref = 2e4;
  * nu <  RS_NUSPLIT: n4 A+ + A- Lefschetz-thimble sum (saddles tracked from nu_c <= 10.5), CRN.
New here: kind 'EE4' = the w0 EE (k = 2) integrand times the Lemma-B factor R = ((v'-u')^2 - 1)/(4u'v') (rational in z, so no new
branch), i.e. exactly the integrand of bg_eeq's 'EE4' image (bg_ee was checked point-wise against the same w0 builder) and of bf_q.
The scalar is computed on the same samples (validation against the exact Mellin skeleton).
Output: cache/rs_pts_<kind>_J<J><tag>.pkl (key (Delta, round(nu, 9)); z2/n1 record format), cache/rs_raw<tag>.pkl.
Usage: OMP_NUM_THREADS=1 python3 rs_points.py J nu1,nu2,... Dstate[,Dref]   [env RS_TAG, RS_KINDS, Z2_M, Z2_R, RS_NUSPLIT, Z2_SEED]
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
os.environ['Z2_KINDS'] = os.environ.get('RS_KINDS', 'EE4')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import z2_points as Z

C = os.path.join(HERE, 'cache')
TAG = os.environ.get('RS_TAG', '')
NU_SPLIT = float(os.environ.get('RS_NUSPLIT', '11.0'))
_P = [False]


def lemma_b_R(zl):
    """R = ((v'-u')^2 - 1)/(4 u'v') from raw z (list of 7 arrays), l1_crystal gauge (as bf_core.cross_ratios)."""
    def dot(p, q):
        return -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]
    ux, vx, uy, vy, xp, yp = (np.asarray(zl[i]) for i in range(1, 7))
    x = [ux, vx, xp]; y = [uy, vy, yp]
    xsq = dot(x, x); x4 = [c/(-xsq) for c in x]; x2 = [-c for c in x4]; x1 = [-c for c in y]; x3 = y
    sub = lambda p, q: [p[i] - q[i] for i in range(3)]
    msq = lambda p, q: -dot(sub(p, q), sub(p, q))
    m12, m34, m13, m24, m14, m23 = msq(x1, x2), msq(x3, x4), msq(x1, x3), msq(x2, x4), msq(x1, x4), msq(x2, x3)
    up = m12*m34/(m13*m24); vp = m14*m23/(m13*m24)
    return ((vp - up)**2 - 1)/(4*up*vp)


def patch():
    if _P[0]:
        return
    Z.patch()                                   # z2: Q_J rescaling + kind 'EE' -> w0 integrand
    import n1_pw as PW
    inner = PW.P9.build_spin

    def bs(zl, nu, Delta, J, mons, sol, solt, kind):
        if kind == 'EE4':
            lamF, G, pref = inner(zl, nu, Delta, J, mons, sol, solt, 'EE')
            return lamF, G, pref*lemma_b_R(zl)[None, None, :]
        return inner(zl, nu, Delta, J, mons, sol, solt, kind)
    PW.P9.build_spin = bs
    _P[0] = True


def fn(kind, J):
    return os.path.join(C, 'rs_pts_%s_J%d%s.pkl' % (kind, J, TAG))


def load(kind, J, tag=None):
    f = os.path.join(C, 'rs_pts_%s_J%d%s.pkl' % (kind, J, TAG if tag is None else tag))
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def run(J, nus, Ds, m=Z.MQ, R=Z.RQ):
    patch()
    for nu in nus:
        keys = [(D, round(nu, 9)) for D in Ds]
        if all(all(k in load(kd, J) for kd in Z.KINDS) for k in keys):
            print(f"# J={J} nu={nu}: done", flush=True); continue
        t = time.time()
        if nu >= NU_SPLIT:
            res = Z.n3_pair(J, nu, Ds, m, R)
        else:
            res = Z.n4_pair(J, nu, Ds, m, R, more=[v for v in nus if v < NU_SPLIT])
        for D, (o, raw) in res.items():
            raw = dict(raw); raw.pop('prop', None)
            Z._save(os.path.join(C, 'rs_raw%s.pkl' % TAG), (raw['method'], round(D, 9), J, round(nu, 9)), raw)
            for k in Z.KINDS:
                Z._save(fn(k, J), (D, round(nu, 9)), o[k])
        print(f"## J={J} nu={nu:.4f} Ds={Ds}: done in {time.time() - t:.0f}s", flush=True)


if __name__ == '__main__':
    J = int(sys.argv[1])
    nus = [float(v) for v in sys.argv[2].split(',')]
    Ds = [float(v) for v in sys.argv[3].split(',')]
    if len(Ds) == 1:
        Ds.append(20000.0)
    run(J, nus, Ds)
