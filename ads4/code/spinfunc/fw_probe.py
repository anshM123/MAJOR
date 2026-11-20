"""
fw_probe (FOLD_WINDOW_LARGEJ.md §2.5): precision of a single direct fold-window node (no Gamma fit): relative jackknife error of the
state/reference helicity-diagonal ratio (generalised basis, fw_fix), estimators dir and CV (exact scalar from sm_exact when available),
and the implied Gamma error  dGamma ~ m^2 |D_h(x)| delta_rel  (c2 units; the a_h normalisation is O(1)).
Usage: python3 fw_probe.py J:Delta nu
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import fw_gamma as FG
import z3_nodes as Z3

spec, nu = sys.argv[1], float(sys.argv[2])
J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
m = Z3.mphys(D, J); x = nu/m
FG.CUR[J] = dict(jk=None, D=D)
R = FG.Q._strip(FG.fwk2('EE', J), J)
Bs, Br = FG.gbasis(J, D, R), FG.gbasis(J, 2e4, R)
print(f"J={J} Delta={D:g} nu={nu} x={x:.4f} (2J+1)x={(2*J+1)*x:.2f}; bases: {FG.GB_LOG}")
for kind in ('EE', 'C3H'):
    P = FG.fwpts(kind, J)
    st = FG.FX.fix_direct(P[(D, round(nu, 9))], J, D, nu); rf = FG.FX.fix_direct(P[(2e4, round(nu, 9))], J, 2e4, nu)
    Rn = len(st['reps'])
    bst, brf = np.linalg.inv(FG.FX.bcur(D, J)[0]), np.linalg.inv(FG.FX.bcur(2e4, J)[0])
    def diag(rep, b, B):
        return np.diag(B['Ti'] @ (rep @ b.T) @ B['TiR'].T)
    cvfac = None
    try:
        import fw_scalar as FS
        FS.install(J, D)
        o = Z3.scalar_check(J, nu, {D: st, 2e4: rf}, D)
        cvfac = o
    except Exception as ex:
        print(f"   (no exact scalar for CV: {str(ex)[:60]})")
    for est in ('dir', 'cv'):
        if est == 'cv' and cvfac is None:
            continue
        def ratio(sel, jk=None):
            a = diag(st['reps'][sel].mean(0), bst, Bs); b = diag(rf['reps'][sel].mean(0), brf, Br)
            f = 1.0 if est == 'dir' else (cvfac['cv'] if jk is None else cvfac['cv_jk'][jk])
            return a/b*f
        v0 = ratio(np.ones(Rn, bool))
        jk = np.array([ratio(np.arange(Rn) != i, i) for i in range(Rn)])
        err = np.sqrt((Rn - 1)/Rn*np.sum(np.abs(jk - jk.mean(0))**2, axis=0))
        rel = err/np.abs(v0)
        Dh = np.array([abs(FG.E.Dc_fixed(h, J, x)) for h in ('++', '+-')])
        env_amp = np.sqrt(2/(np.pi*max((2*J + 1)*x, 1.0)))
        print(f"   {kind:3s} {est}: rel err of the helicity ratio (++, +-) = {np.round(rel, 4).tolist()};  "
              f"dGamma ~ m^2 |D_h| rel = {np.round(m*m*Dh*rel).tolist()} c2  (with the envelope {env_amp:.2f} instead of |D_h|: "
              f"{np.round(m*m*env_amp*rel).tolist()})")
