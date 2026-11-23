"""
dc_ptwin (DIRCV.md §7): extra statistics for the UNBIASED P_J twin only: state P_J (mode 0, n3 contour) + reference, on the same
contours and the same proposal as production (n3_hybrid.adapt on the W+ evaluator, seed 31), new Sobol seeds; no wave modes.
Saved in pw_points format to dircv/cache/pw_pts_J<J>_D<D><tag>P.pkl (state P_J and its own reference), for dc_pool.
Usage: python3 dc_ptwin.py J Delta nu R seed tag
"""
import dc_guard, os, sys, pickle, time
os.environ.setdefault('L1_DEG', '6'); os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
import pw_points as PP
import pw_core as Co
import n3_hybrid as Hy
import z3_nodes as Z3

J, Ds, nu, R, seed, tag = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
M = 9
PP.setup()
t0 = time.time()
prop = Hy.adapt(Co.make_ev(nu, Ds, J, 1), 3.5)
out = {}
for D in (Ds, Z3.DREF):
    ev = Co.make_ev(nu, D, J, 0)
    gk = {id(ev): {k: Co.make_gev(ev, k) for k in ('EE', 'C3H')}}
    res = Co.block([ev], prop, M, R, seed, gkinds=gk)
    out[(D, round(nu, 9))] = dict(J=J, Delta=D, nu=nu, m=Z3.mphys(D, J), modes=[0], R=R, mq=M, nuf=ev.nuf, adapt=True, scale=1.0,
                                  per_mode=[{k: v for k, v in o.items() if k != 'g'} for o in res],
                                  g={k: [o['g'][k] for o in res] for k in ('EE', 'C3H')}, time=time.time() - t0, wbranch=[None])
    print(f"  Delta={D:g} nu={nu}: scalar rel err {abs(res[0]['err']/res[0]['E']):.4f}; diag {res[0]['diag']}; {time.time()-t0:.0f}s", flush=True)
pickle.dump(out, open(os.path.join(dc_guard.CACHE, 'pw_pts_J%d_D%g%sP.pkl' % (J, Ds, tag)), 'wb'))
print(f"## done {tag} in {time.time()-t0:.0f}s", flush=True)
