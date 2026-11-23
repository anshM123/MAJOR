"""
dc_tail (DIRCV.md §3.7): per-POINT samples (the production block stores only replica means) of the per-wave scalar and graviton
integrands at one node, same contour/proposal as pw_points (PW_ADAPT=1), to test heavy tails (Hill tail index, top-k share),
graviton-only point drops (branch/winding filter), and the point-wise graviton/scalar ratio that the control variate assumes ~ const.
Also mode 0 (P_J, n3 contour) for comparison.  Usage: python3 dc_tail.py J Delta nu m R seed -> cache/tail_J_D_nu.pkl
"""
import dc_guard, os, sys, pickle, time
os.environ.setdefault('L1_DEG', '6')
import numpy as np
import pw_points as PP
import pw_core as Co
import n3_hybrid as Hy

J, D, nu = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
m, R, seed = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
PP.setup()
t0 = time.time()
modes = [1, -1, 0]
evs = [Co.make_ev(nu, D, J, md) for md in modes]
gk = {id(ev): {k: Co.make_gev(ev, k) for k in ('EE', 'C3H')} for ev in evs}
prop = Hy.adapt(evs[0], 3.5)
mu, Ch, ldC = prop
out = dict(J=J, D=D, nu=nu, m=m, R=R, seed=seed, modes=modes, data={})
for ev in evs:
    S, G, OKs, OKg = [], {k: [] for k in ('EE', 'C3H')}, [], {k: [] for k in ('EE', 'C3H')}
    for r in range(R):
        eta = Hy._etas(m, seed + r); eta = np.concatenate([eta, -eta])
        for c0 in range(0, eta.shape[0], Hy.CHUNK):
            et = eta[c0:c0 + Hy.CHUNK]; Z = mu[None, :] + et @ Ch.T
            wgt = np.exp(0.5*ldC + 0.5*np.sum(et*et, axis=1))
            ct = Hy.contour(ev, Z, 3.5, keep_path=True)
            S.append(wgt*ct['Fval']); OKs.append(ct['ok'])
            for k, ge in gk[id(ev)].items():
                Fg, Pgc, pref0, okg = Co.grav_on_paths(ge, ct, Z.shape[0])
                G[k].append(wgt*Fg); OKg[k].append(okg)
    (Pc, _), _, _ = ev.run(ev.pint, np.zeros((1, 7), complex), None)
    lg = {}
    for k, ge in gk[id(ev)].items():
        (Pg, _), _, _ = ge.run(ge.grav, np.zeros((1, 7), complex), None)
        lg[k] = Pg[0] + np.log(np.linalg.det(ge.L))
    out['data'][ev.mode] = dict(S=np.concatenate(S), okS=np.concatenate(OKs), logN=Pc[0] + np.log(np.linalg.det(ev.L)),
                                G={k: np.concatenate(v, axis=-1) for k, v in G.items()}, okG={k: np.concatenate(v) for k, v in OKg.items()},
                                logNg=lg)
    print(f"mode {ev.mode}: {len(out['data'][ev.mode]['S'])} points, scalar drops {(~out['data'][ev.mode]['okS']).sum()}, "
          f"graviton drops " + ", ".join(f"{k}: {(~v).sum()}" for k, v in out['data'][ev.mode]['okG'].items()) + f" ({time.time()-t0:.0f}s)", flush=True)
pickle.dump(out, open(os.path.join(dc_guard.CACHE, 'tail_J%d_D%g_nu%g_s%d.pkl' % (J, D, nu, seed)), 'wb'))
