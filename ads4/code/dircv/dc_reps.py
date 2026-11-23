"""dc_reps (DIRCV.md §3.3): per-replica view at one node: scalar waves vs exact, P_J vs exact, graviton waves/P_J (state), ms."""
import dc_guard, os, sys, pickle
import numpy as np
import z3_nodes as Z3, fw_scalar as FS, pw_gamma as PG
J, D = 100, 1100.0
FS.install(J, D)
tagfiles = sys.argv[2].split(',') if len(sys.argv) > 2 else ['_w']
nu = float(sys.argv[1]) if len(sys.argv) > 1 else 13.3
for t in tagfiles:
    base = dc_guard.SPIN + '/cache/' if os.path.exists(dc_guard.SPIN + '/cache/pw_pts_J100_D1100%s.pkl' % t) else dc_guard.CACHE + '/'
    w = pickle.load(open(base + 'pw_pts_J100_D1100%s.pkl' % t, 'rb')); P = pickle.load(open(base + 'pw_pts_J100_D1100%sP.pkl' % t, 'rb'))
    if (D, nu) not in w: continue
    sw, sp, rf = PG.to_z2(w[(D, nu)], 'EE'), PG.to_z2(P[(D, nu)], 'EE'), PG.to_z2(w[(Z3.DREF, nu)], 'EE')
    lxs, _ = Z3.lx_exact(D, J, nu); lxr, _ = Z3.lx_exact(Z3.DREF, J, nu)
    off = (sw['logNs'] - Z3.logbeta(D, J)) - (rf['logNs'] - Z3.logbeta(Z3.DREF, J))
    offP = (sp['logNs'] - Z3.logbeta(D, J)) - (rf['logNs'] - Z3.logbeta(Z3.DREF, J))
    dw = np.log(sw['repsS']/rf['repsS']) + off - (lxs - lxr)
    dP = np.log(sp['repsS']/rf['repsS']) + offP - (lxs - lxr)
    gw, gp = PG.to_z2(w[(D, nu)], 'EE'), PG.to_z2(P[(D, nu)], 'EE')
    g = np.log(gw['reps'][:, 0, 0]/gp['reps'][:, 0, 0]) + gw['logN'] - gp['logN']
    sr = np.log(sw['repsS']/sp['repsS']) + sw['logNs'] - sp['logNs']
    pm = w[(D, nu)]['per_mode']
    print(f"# tag {t} nu={nu}  R={len(dw)}  diag: " + "; ".join(f"mode {o['mode']}: nbad {o['diag']['nbad']} wind {o['diag']['wind']:.2f} imdrift {o['diag']['imdrift']:.1e}" for o in pm))
    print("  rep   d_waves(Re)   d_PJ(Re)   log(Gw/GP)00   log(Sw/SP)   ms=diff     Im d_waves")
    for i in range(len(dw)):
        print(f"  {i:3d}  {dw[i].real:+.5f}   {dP[i].real:+.5f}   {g[i].real:+.5f}     {sr[i].real:+.5f}   {(g[i]-sr[i]).real:+.5f}   {dw[i].imag:+.5f}")
    print(f"  mean of per-rep: d_w {dw.real.mean():+.5f} +- {dw.real.std(ddof=1)/np.sqrt(len(dw)):.5f};  ms {(g - sr).real.mean():+.5f} +- {(g - sr).real.std(ddof=1)/np.sqrt(len(dw)):.5f}")
