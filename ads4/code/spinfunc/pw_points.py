"""
pw_points (HWIN_THIMBLE.md): per-wave direct crystal points, state (Delta) and heavy reference (Delta_ref = 2e4) on COMMON RANDOM
NUMBERS, scalar + graviton kinds on the same contour points.

  state : modes PW_MODES (default "1,-1": the two travelling waves, each on the n3_hybrid contour of its OWN wave saddle;
          "0": P_J on the standard n3 contour -- the fw_points/z2 method, for comparison);
  ref   : mode 0 (P_J; at Delta_ref = 2e4 the spin phase is negligible: (2J+1) nu/m_ref <~ 1).
  contour nu_f: PW_NUF if set, else pw_core.nuf_of(nu) (= the n3/z2 rule: max(nu, 15) above 13, 13 below).
  proposal: Gaussian N(0, PW_S^2) (default 1), or (PW_ADAPT=1) n3_hybrid.adapt on the state's first mode, re-used everywhere (CRN).
Stored per node (cache/pw_pts_J<J>_D<Dstate><PW_TAG>.pkl, key (Delta, round(nu, 9))): per mode the scalar replicas 'reps' and
'logN', and per graviton kind the replica matrices; meta (m, R, nu_f, diag, time).  Resumable (finished nodes skipped).
Usage: OMP_NUM_THREADS=1 python3 pw_points.py J Dstate nu1,nu2,...   [env PW_KINDS (default 'S'; e.g. 'EE,C3H'), PW_M, PW_R, ...]
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
M = int(os.environ.get('PW_M', '9'))
R = int(os.environ.get('PW_R', '8'))
SEED = int(os.environ.get('PW_SEED', '4242'))
MODES = [int(v) for v in os.environ.get('PW_MODES', '1,-1').split(',')]
KINDS = [k for k in os.environ.get('PW_KINDS', 'S').split(',') if k and k != 'S']
NUFX = float(os.environ['PW_NUF']) if os.environ.get('PW_NUF') else None
ADAPT = os.environ.get('PW_ADAPT', '0') == '1'
SCALE = float(os.environ.get('PW_S', '1.0'))
TAG = os.environ.get('PW_TAG', '')


def fn(J, D):
    return os.path.join(C, 'pw_pts_J%d_D%g%s.pkl' % (J, D, TAG))


def load(J, D):
    f = fn(J, D)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def save(J, D, key, val):
    import fcntl
    path = fn(J, D)
    with open(path + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        db = load(J, D); db[key] = val
        tmp = path + '.tmp%d' % os.getpid()
        pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, path)
        fcntl.flock(lk, fcntl.LOCK_UN)


def setup():
    import fw_frame as FW
    import n1_pw as PW
    import l5_spin as S
    FW.patch_conserved_clean(S)
    FW.patch_conserved_clean(PW.S)


def node(J, Ds, nu):
    import pw_core as Co
    import n3_hybrid as Hy
    import z3_nodes as Z3
    out = {}
    prop = None
    for D, modes in ((Ds, MODES), (DREF, [0])):
        t0 = time.time()
        evs = [Co.make_ev(nu, D, J, md, nuf=NUFX) for md in modes]
        gk = {}
        for ev in evs:
            gk[id(ev)] = {k: Co.make_gev(ev, k) for k in KINDS}
        if prop is None:
            prop = Hy.adapt(evs[0], 3.5) if ADAPT else (np.zeros(7), SCALE*np.eye(7), 7*np.log(SCALE**2))
        res = Co.block(evs, prop, M, R, SEED, gkinds=gk)
        rec = dict(J=J, Delta=D, nu=nu, m=Z3.mphys(D, J), modes=modes, R=R, mq=M, nuf=evs[0].nuf, adapt=ADAPT, scale=SCALE,
                   per_mode=[{k: v for k, v in o.items() if k != 'g'} for o in res],
                   g={k: [o['g'][k] for o in res] for k in KINDS}, time=time.time() - t0,
                   wbranch=[getattr(e, 'wbranch', None) for e in evs])
        out[D] = rec
        rel = [abs(o['err']/o['E']) for o in res]
        print(f"  J={J} Delta={D:g} nu={nu:.4f} modes {modes} nu_f {evs[0].nuf:g}: scalar rel err {np.round(rel, 4).tolist()}; "
              f"diag {[o['diag'] for o in res]}; {time.time() - t0:.0f}s", flush=True)
    return out


def run(J, Ds, nus):
    setup()
    for nu in nus:
        key = (Ds, round(nu, 9))
        db = load(J, Ds)
        if key in db and (DREF, round(nu, 9)) in db:
            print(f"# J={J} Delta={Ds:g} nu={nu}: done", flush=True); continue
        t = time.time()
        o = node(J, Ds, nu)
        for D, rec in o.items():
            save(J, Ds, (D, round(nu, 9)), rec)
        print(f"## J={J} Delta={Ds:g} nu={nu:.4f}: done in {time.time() - t:.0f}s", flush=True)


if __name__ == '__main__':
    J = int(sys.argv[1]); Ds = float(sys.argv[2])
    nus = [float(v) for v in sys.argv[3].split(',')]
    run(J, Ds, nus)
