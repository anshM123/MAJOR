"""
pw_probe (HWIN_THIMBLE.md): precision of pw_points graviton nodes (as fw_probe): relative jackknife error of the state/reference
helicity-diagonal ratio in the generalised basis (fw_gamma.gbasis on the fw_k2 rows), estimators dir and CV (exact scalar), and the
implied dGamma ~ m^2 |D_h(x)| delta_rel (c2).  Also, where fw_points has the same node, the ratio of the two methods' central values.
Usage: python3 pw_probe.py J:Delta [PW_TAG via env]
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import fw_gamma as FG
import z3_nodes as Z3
import pw_points as PP
import pw_scal as PS


def gmat(rec, kind, sel=None):
    """state graviton matrix sum_modes reps exp(logN_mode - logN_0) (replica subset sel), and logN_0."""
    gs = rec['g'][kind]
    ref = gs[0]['logN']
    tot = sum(np.asarray(g['reps'])*np.exp(g['logN'] - ref) for g in gs)
    if sel is not None:
        tot = tot[sel]
    return tot.mean(0), ref


def main(J, D, pr=print):
    import fw_scalar as FS
    FS.install(J, D)
    FG.CUR[J] = dict(jk=None, D=D)
    Rk = FG.Q._strip(FG.fwk2('EE', J), J)
    Bs, Br = FG.gbasis(J, D, Rk), FG.gbasis(J, 2e4, Rk)
    bst, brf = np.linalg.inv(FG.FX.bcur(D, J)[0]), np.linalg.inv(FG.FX.bcur(2e4, J)[0])
    db = PP.load(J, D)
    out = {}
    for (DD, nu) in sorted(k for k in db if k[0] == D):
        if (PP.DREF, nu) not in db:
            continue
        st, rf = db[(D, nu)], db[(PP.DREF, nu)]
        m = Z3.mphys(D, J); x = nu/m
        o = PS.check(J, D, nu, st, rf)
        Rn = st['R']
        for kind in st['g']:
            def diag(M, b, B):
                return np.diag(B['Ti'] @ (M @ b.T) @ B['TiR'].T)
            def ratio(sel, est):
                a, la = gmat(st, kind, sel); b, lb = gmat(rf, kind, sel)
                r = diag(a, bst, Bs)/diag(b, brf, Br)*np.exp(la - lb)
                if est == 'cv':
                    # MC scalar ratio on the same replicas replaced by the exact one
                    lsr = PS.logY_reps(st, sel) - PS.logY_reps(rf, sel)
                    import z3_nodes as Z
                    lxs, _ = Z.lx_exact(D, J, nu); lxr, _ = Z.lx_exact(PP.DREF, J, nu)
                    r = r*np.exp((lxs - lxr) - lsr)
                return r
            for est in ('dir', 'cv'):
                v0 = ratio(np.ones(Rn, bool), est)
                jk = np.array([ratio(np.arange(Rn) != i, est) for i in range(Rn)])
                err = np.sqrt((Rn - 1)/Rn*np.sum(np.abs(jk - jk.mean(0))**2, axis=0))
                rel = err/np.abs(v0)
                Dh = np.array([abs(FG.E.Dc_fixed(h, J, x)) for h in ('++', '+-')])
                out[(nu, kind, est)] = (v0, rel)
                pr(f"  nu={nu:6.2f} {kind:3s} {est}: helicity ratio rel err (++, +-) {np.round(rel, 5).tolist()};  dGamma ~ "
                   f"{np.round(m*m*Dh*rel).tolist()} c2")
        pr(f"  nu={nu:6.2f} scalar check: Delta c2 = {o['dc2']:+.1f} +- {o['ec2']:.1f}")
    return out


if __name__ == '__main__':
    spec = sys.argv[1]
    J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
    main(J, D)
