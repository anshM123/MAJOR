"""n3fix/d2_logs.py: log-by-log comparison (n1_pw record/replay sequence) of the C3H integrand at a few real inner points: n1's
two-level continuation from the saddle (principal logs at the saddle) vs n3's replay along the stored ray (set_branch saddle references).
Usage: python3 d2_logs.py J Delta x [kind]"""
import nf_guard  # noqa: F401
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc')); sys.path.insert(0, SF)
os.environ['HC_CACHE'] = os.path.join(HERE, 'cache')
os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import hc_points as HC
import z2_points as Z2

J = int(sys.argv[1]); Dl = float(sys.argv[2]); x = float(sys.argv[3]); kind = sys.argv[4] if len(sys.argv) > 4 else 'C3H'
HC.setup(odd=(J % 2 == 1))
import n1_pw as PW, n1_direct as D, n3_flow as F, n3_hybrid as Hy
nu = x*HC.mphys(Dl, J)
c, wc, L, _ = F.saddle(nu, Dl, J)
Z2.set_qsc(Dl, J, c, wc, L)
ev = F.Evaluator(nu, Dl, J, c, wc, L, kind, nuf=nu); F.set_branch(ev)
rng = np.random.default_rng(1)
Z = 0.8*rng.normal(size=(3, 7))
# n3 route: saddle refs (set_branch), replay along 14 ray steps
(Pg0, p0), ref0, _ = ev.run(ev.grav, np.zeros((1, 7), complex), None)
nlog = len(ref0)
print(f"{kind}: number of recorded logs per evaluation: {nlog}; saddle refs from set_branch: {len(ev.saddle_refs['grav'])}")
refs = [np.repeat(r, Z.shape[0], axis=-1) for r in ref0]
for k in range(1, 15):
    _, refs, _ = ev.run(ev.grav, (k/14)*Z.astype(complex), refs, rep=1)
# n1 route: record from the saddle (index 0), plain principal logs there, then replay at the full point
PW.path_mode('record')
Wp = wc[None, :] + np.vstack([np.zeros((1, 7)), 0.5*Z]) @ L.T
PW.set_batch(Wp.shape[0]); zl, lj = D.zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
ev_ = PW.P9.build_spin(zl, nu, Dl, J, ev.mons, ev.sol, ev.solt, kind) if J % 2 == 0 else None
PW.path_mode('replay_record')
Wp = wc[None, :] + np.vstack([np.zeros((1, 7)), Z]) @ L.T
PW.set_batch(Wp.shape[0]); zl, lj = D.zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
PW.P9.build_spin(zl, nu, Dl, J, ev.mons, ev.sol, ev.solt, kind)
n1refs = PW.REC[1]
PW.path_mode('plain')
print(f"n1 recorded {len(n1refs)} logs")
# saddle comparison: n3 saddle ref vs n1 principal at saddle (index 0)
for i, (a, b) in enumerate(zip(ref0, n1refs)):
    a = np.asarray(a); b = np.asarray(b)
    ds = (a[..., 0] - b[..., 0])                        # saddle
    dp = (refs[i] - b[..., 1:])                         # points
    dsi = np.round(ds.imag/(2*np.pi), 3); dpi = np.round(dp.imag/(2*np.pi), 3)
    if np.abs(ds).max() > 1e-8 or np.abs(dp).max() > 1e-8:
        print(f" log {i:3d} shape {a.shape}: saddle diff/2pi i: {np.unique(dsi)[:6]} (re {np.abs(ds.real).max():.1e}); "
              f"points diff/2pi i: {np.unique(dpi)[:6]} (re {np.abs(dp.real).max():.1e})")
# prefactors: n1 route (replay at the full point, batch [saddle, pts]) vs n3 route (replay with refs, batch pts only)
PW.path_mode('record')
Wp = wc[None, :] + np.vstack([np.zeros((1, 7)), 0.5*Z]) @ L.T
PW.set_batch(Wp.shape[0]); zl, lj = D.zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
PW.P9.build_spin(zl, nu, Dl, J, ev.mons, ev.sol, ev.solt, kind)
PW.path_mode('replay_record')
Wp = wc[None, :] + np.vstack([np.zeros((1, 7)), Z]) @ L.T
PW.set_batch(Wp.shape[0]); zl, lj = D.zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
_, _, p1 = PW.P9.build_spin(zl, nu, Dl, J, ev.mons, ev.sol, ev.solt, kind)
PW.path_mode('plain')
_, refs_b = refs, [r.copy() for r in refs]
(Pg3, p3), _, _ = ev.run(ev.grav, Z.astype(complex), [r for r in refs_b], rep=1)
(Pg3b, p3b), _, _ = ev.run(ev.grav, np.vstack([np.zeros((1, 7)), Z]).astype(complex), [np.concatenate([np.asarray(a)[..., :1], b], axis=-1) for a, b in zip(ref0, refs_b)], rep=1)
print('pref n1 (point 0)  ', np.round(p1[..., 1].ravel(), 12))
print('pref n3 (point 0)  ', np.round(p3[..., 0].ravel(), 12))
print('pref n3 w/ saddle  ', np.round(p3b[..., 1].ravel(), 12))
