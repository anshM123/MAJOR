"""
fw_scalar (FOLD_WINDOW_LARGEJ.md): validation of the direct fold-window method against the EXACT scalar skeleton at the same (Delta, J):
the scalar computed on the same common-random-number samples as the graviton (fw_points) must reproduce the exact state/reference ratio,
d_D - d_ref = [log Y_D - log(q_D X_D)] - [same at Delta_ref] = 0  (z3_nodes.scalar_check), with X the exact CLS-D object
X(nu) = g C(nu + i) + g~ C(nu - i) (n3_exact) built from the Arb-exact rho_2 of sm_exact (SMOOTH_BOUND §2: cache/sm_exact_<Dg>.pkl; at
Delta_g = 300 the 'n0t' n-sum configuration) for the state, and the heavy-block hybrid (rho_2 := 0) for Delta_ref = 2e4 (z3 convention).
Delta c2 = m^2 Re(d - d_ref) ex2/E2 is the error in c2 units the method would make on a known object.
Usage: python3 fw_scalar.py J:Delta[,..]  -> stdout, cache/fw_scalar.txt (appended)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import z3_nodes as Z3
import fw_gamma as FG

C = Z3.C


def sm_rec(D, J):
    Dg = int(round(D - J))
    db = pickle.load(open(os.path.join(C, 'sm_exact_%d.pkl' % Dg), 'rb'))
    keys = [k for k in db if k[0] == J and abs(k[1]) < 1e-9]
    k = sorted(keys, key=lambda k: (len(k) < 4, k))[0] if Dg == 300 else sorted(keys, key=len)[0]
    if Dg == 300:
        kk = [k_ for k_ in keys if len(k_) > 3 and k_[3] == 'n0t']
        k = kk[0] if kk else k
    r = db[k]
    st = r['st']
    rec = dict(J=J, m=st['m'], eta=st['eta'], Delta=st['Delta'], nus=np.asarray(r['nus'], float), x=np.asarray(r['x'], float),
               rho={2: np.asarray(r['rho'][2], float)}, src=k)
    E2 = (2 - rec['x']**2)/(1 - rec['x']**2)**2
    import n3_exact as NX
    rec['ex'] = {2: np.array([NX.hyb(v + 0j, rec).real for v in rec['nus']]) + E2*rec['rho'][2]}
    return rec


def install(J, D):
    db = pickle.load(open(Z3.EXPK, 'rb')) if os.path.exists(Z3.EXPK) else {}
    Z3._exact_override = {(round(float(D), 9), int(J)): sm_rec(D, J)}
    orig = Z3.exact_rec

    def exact_rec(DD, JJ):
        key = (round(float(DD), 9), int(JJ))
        if key in Z3._exact_override:
            return Z3._exact_override[key]
        if DD >= 1e4:                       # heavy reference: hybrid with rho_2 = 0 (z3 fallback), computed here (no ExactState)
            from scipy.special import j0, eval_legendre
            import mpmath as mp
            al = DD - JJ - 2.0; ab = al + 2*JJ + 1
            m = np.sqrt(al*ab); eta = (al**2 + ab**2)/(2*al*ab)
            nus = np.linspace(0.5, 60.0, 120); x = nus/m
            E2 = (2 - x*x)/(1 - x*x)**2
            P = np.array([float(mp.legendre(-0.5 + 1j*v, eta).real) for v in nus])
            hyb = E2*eval_legendre(JJ, 1 - 2*x*x) + 2*(P - j0((2*JJ + 1)*x))
            r = dict(J=JJ, m=m, eta=eta, Delta=float(DD), nus=nus, x=x, ex={2: hyb}, rho={2: 0*hyb}, heavy_rho0=True)
            Z3._exact_override[key] = r
            return r
        return orig(DD, JJ)
    Z3.exact_rec = exact_rec


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    for spec in sys.argv[1].split(','):
        J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
        install(J, D)
        pr(f"# exact scalar validation J={J} Delta={D:g} (state exact: sm_exact {Z3._exact_override[(round(D, 9), J)]['src']})")
        P = FG.fwpts('EE', J, D)
        by = {}
        for (DD, nu), r in P.items():
            by.setdefault(nu, {})[DD] = r
        for nu in sorted(by):
            dd = by[nu]
            if D not in dd or Z3.DREF not in dd:
                continue
            o = Z3.scalar_check(J, nu, dd, D)
            pr(f"  nu={nu:7.3f} x={nu/o['m']:.4f}: Re(d - d_ref) {o['d'].real:+.5f} +- {o['err'].real:.5f}, Im {o['d'].imag:+.5f} +- "
               f"{o['err'].imag:.5f};  Delta c2 = {o['dc2']:+9.1f} +- {o['ec2']:7.1f}  ({abs(o['dc2'])/max(o['ec2'], 1e-9):.1f} sigma)")
    with open(os.path.join(C, 'fw_scalar.txt'), 'a') as f:
        f.write("\n".join(lines) + "\n")
