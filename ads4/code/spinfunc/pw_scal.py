"""
pw_scal (HWIN_THIMBLE.md): exact-scalar validation of pw_points nodes.
d = log Y - log(q X) (z3_nodes convention; X exact: sm_exact rho_2 via fw_scalar.install for the state, heavy hybrid for 2e4);
Delta c2 = m^2 Re(d_state - d_ref) ex2/E2 must vanish; jackknife over the common-seed replicas.  Y_state = sum over the stored modes.
Usage: python3 pw_scal.py J:Delta[,J:Delta...]  [env PW_TAG]  -> stdout, cache/pw_scal.txt (appended)
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import numpy as np
import l1_crystal as Cr
import pw_points as PP


def logY_reps(rec, sel=None):
    """log Y from the per-mode replica sums (common normalisation), on the replica subset sel."""
    pm = rec['per_mode']
    ref = pm[0]['logN']
    tot = sum(np.asarray(o['reps'])*np.exp(o['logN'] - ref) for o in pm)
    if sel is not None:
        tot = tot[sel]
    return np.log(tot.mean()) + ref - Cr.log_beta(rec['Delta'], rec['J'])


def check(J, D, nu, st, rf):
    import z3_nodes as Z3
    lxs, rs = Z3.lx_exact(D, J, nu); lxr, _ = Z3.lx_exact(PP.DREF, J, nu)
    R = len(st['per_mode'][0]['reps'])
    def dd(sel):
        return (logY_reps(st, sel) - logY_reps(rf, sel)) - (lxs - lxr)
    v0 = dd(np.ones(R, bool))
    jk = np.array([dd(np.arange(R) != i) for i in range(R)])
    wrap = lambda z: complex(z.real, (z.imag + np.pi) % (2*np.pi) - np.pi)
    v0 = wrap(v0); jk = np.array([wrap(z - v0) for z in jk]) + v0
    e = np.sqrt((R - 1)/R*np.sum(np.abs(jk.real - jk.real.mean())**2)), np.sqrt((R - 1)/R*np.sum(np.abs(jk.imag - jk.imag.mean())**2))
    m = Z3.mphys(D, J); x = nu/m; E2 = (2 - x*x)/(1 - x*x)**2
    ex2 = float(np.interp(nu, rs['nus'], rs['ex'][2]))
    fac = m*m*ex2/E2
    return dict(d=v0, err=e, dc2=fac*v0.real, ec2=abs(fac)*e[0], m=m, cv=np.exp((lxs - lxr) - (logY_reps(st) - logY_reps(rf))))


if __name__ == '__main__':
    import fw_scalar as FS
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    for spec in sys.argv[1].split(','):
        J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
        FS.install(J, D)
        db = PP.load(J, D)
        pr(f"# pw scalar validation J={J} Delta={D:g} tag '{PP.TAG}'")
        for (DD, nu) in sorted(k for k in db if k[0] == D):
            if (PP.DREF, nu) not in db:
                continue
            st, rf = db[(D, nu)], db[(PP.DREF, nu)]
            o = check(J, D, nu, st, rf)
            rel = [abs(np.std(p['reps'])/np.mean(p['reps']))/np.sqrt(len(p['reps'])) for p in st['per_mode']]
            pr(f"  nu={nu:7.3f} x={nu/o['m']:.4f} modes {st['modes']} nu_f {st['nuf']:g} R={st['R']} 2^{st['mq']}x2: "
               f"Re d {o['d'].real:+.6f} +- {o['err'][0]:.6f}, Im {o['d'].imag:+.6f} +- {o['err'][1]:.6f};  "
               f"Delta c2 = {o['dc2']:+8.1f} +- {o['ec2']:6.1f} ({abs(o['dc2'])/max(o['ec2'], 1e-9):.1f} sigma);  "
               f"per-mode rel err {np.round(rel, 4).tolist()}")
    with open(os.path.join(PP.C, 'pw_scal.txt'), 'a') as f:
        f.write("\n".join(lines) + "\n")
