"""
n3_scalar_validate (doc 25): validation of the small-nu contour method (n3_hybrid: wall-free hybrid contour, nu_f = max(nu, 15),
adaptive RQMC) on the EXACT scalar skeleton at Delta_g = 100, at the e1 nodes nu in [NUMIN, 20] (x = nu/m < 0.2), states
m/Lam in {1.0, 1.05, 1.1, 1.2, 1.3}, J in {0, 2, ..., 10} (those present in e1).
Comparison (as n1_scalar_validate / l2): log|Y/(q X)| = a_nu + b_J, with the EXACT magic-block object X(nu) = g C(nu+i) + g~ C(nu-i) of
n3_exact (no local Taylor shift; rho_2 continued to nu +- i by interpolation), residual in c2 units Delta c2 = m^2 eps ex_2/E_2.
Also reported with l2's local-shift X (shows the size of the comparison-object error at small nu).
Usage: OMP_NUM_THREADS=1 nohup python3 n3_scalar_validate.py run       (resumable; cache/n3_scal_100.pkl)
       python3 n3_scalar_validate.py analyse                          -> cache/n3_scalar_validate.txt
"""
import os, sys, pickle, time, warnings
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
import numpy as np

PK = os.path.join(HERE, 'cache', 'n3_scal_100.pkl')
MQ = int(os.environ.get('N3_M', '10'))
NUMIN = float(os.environ.get('N3_NUMIN', '6.0'))
NUMAX = float(os.environ.get('N3_NUMAX', '20.0'))
NUF_LOW = float(os.environ.get('N3_NUF_LOW', '13'))
MFL = (1.0, 1.05, 1.1, 1.2, 1.3)
JS = (0, 2, 4, 6, 8, 10)


def pts():
    import n3_exact as X
    db = X.db(100)
    out = []
    for key in sorted(db):
        mfl, J = key
        if round(mfl, 3) not in MFL or J not in JS:
            continue
        r = db[key]
        for i, nu in enumerate(r['nus']):
            if NUMIN <= nu <= NUMAX:
                out.append((key, i, r))
    out.sort(key=lambda t: (-t[2]['nus'][t[1]], t[0][1], t[0][0]))      # nodes from the top down
    return out


def run():
    import n3_hybrid as Hy
    import l1_crystal as Cr
    sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
    from exactmom import log_q
    for key, i, r in pts():
        res = pickle.load(open(PK, 'rb')) if os.path.exists(PK) else {}
        if (key, i) in res:
            continue
        t = time.time()
        nu, Dl, J = float(r['nus'][i]), float(r['Delta']), key[1]
        nuf = max(nu, 15.0) if nu >= NUF_LOW else NUF_LOW          # same contour rule as n3_run
        o = Hy.integrate(nu, Dl, J, 'S', rho1=3.5, m=MQ, R=8, branch_from='ref', adaptive=True, nuf=nuf)
        rec = dict(mfl=key[0], J=J, m=r['m'], Delta=Dl, x=r['x'][i], nu=nu, ex2=r['ex'][2][i], rho2=r['rho'][2][i],
                   logY=o['logN'] - Cr.log_beta(Dl, J) + np.log(o['Es']), relerr=abs(o['errs']/o['Es']), reps=o['repsS'],
                   logN=o['logN'], diag=o['diag'], nuf=o['nuf'], logq=log_q(Dl, J, 3.0), mq=MQ)
        res = pickle.load(open(PK, 'rb')) if os.path.exists(PK) else {}
        res[(key, i)] = rec
        pickle.dump(res, open(PK + '.tmp%d' % os.getpid(), 'wb')); os.replace(PK + '.tmp%d' % os.getpid(), PK)
        print(f"J={J} m/Lam={key[0]} x={rec['x']:.3f} nu={nu:.2f}: rel err {rec['relerr']:.1e}, nbad {o['diag']['nbad']}, "
              f"{time.time() - t:.0f}s", flush=True)


def analyse(pr, relmax=0.05):
    import l2_scalar_validate as V
    import n3_exact as NX
    res = pickle.load(open(PK, 'rb'))
    db = NX.db(100)
    keys_all = sorted(res)
    keys = [k for k in keys_all if res[k]['relerr'] <= relmax]
    pr(f"# n3 scalar validation (Delta_g = 100): {len(keys_all)} points computed, {len(keys)} with RQMC rel. error <= {relmax:.0%}")
    for k in keys:
        q = res[k]
        q['Xn3'], _ = NX.X(q['nu'], db[k[0]])
        q['Xl2'], _ = V.Xexact(q['nu'], q['m'], q['J'], q['ex2'])
    for lab, xk in (('exact shifted X (n3_exact)', 'Xn3'), ("l2's local-shift X", 'Xl2')):
        nus = sorted(set(round(res[k]['nu'], 9) for k in keys)); iv = {v: j for j, v in enumerate(nus)}
        Js = sorted(set(res[k]['J'] for k in keys)); Jl = [J for J in Js if J != 0]
        rows = [(round(res[k]['nu'], 9), res[k]['J'], res[k]['logY'].real - res[k]['logq'] - np.log(abs(res[k][xk])), res[k]) for k in keys]
        A = np.zeros((len(rows), len(nus) + len(Jl))); b = np.array([r_[2] for r_ in rows])
        for i, (v, J, lr, q) in enumerate(rows):
            A[i, iv[v]] = 1
            if J != 0:
                A[i, len(nus) + Jl.index(J)] = 1
        w = np.array([q['m']**2*abs(q['ex2'])/V.E2(q['x']) for (v, J, lr, q) in rows])
        # weights: c2 units, divided by the RQMC error in c2 units (floor 1) so that noisy points do not dominate
        mc = np.array([q['m']**2*q['relerr']*abs(q['ex2'])/V.E2(q['x']) for (v, J, lr, q) in rows])
        sol, *_ = np.linalg.lstsq(A*(w/np.maximum(mc, 1.0))[:, None], b*w/np.maximum(mc, 1.0), rcond=None)
        eps = b - A @ sol
        dc2 = np.array([q['m']**2*e*q['ex2']/V.E2(q['x']) for e, (v, J, lr, q) in zip(eps, rows)])
        chi = dc2/np.maximum(mc, 1e-9)
        pr(f"## comparison object: {lab}")
        pr(f"   {len(rows)} points, {len(nus)} nodes nu in [{min(nus):.2f}, {max(nus):.2f}], J in {Js}; dof {len(rows) - A.shape[1]}")
        pr(f"   residual |Delta c2| (c2 units): rms {np.sqrt(np.mean(dc2**2)):.1f}, max {np.abs(dc2).max():.1f};  RQMC error (c2 units): "
           f"median {np.median(mc):.1f}, max {mc.max():.1f};  chi2/dof {np.sum(chi**2)/max(len(rows) - A.shape[1], 1):.2f}")
        for v in nus:
            s = np.array([r_[0] == v for r_ in rows])
            pr(f"     nu={v:6.2f}: {s.sum():3d} pts, |Dc2| rms {np.sqrt(np.mean(dc2[s]**2)):7.1f} max {np.abs(dc2[s]).max():7.1f}; "
               f"RQMC err median {np.median(mc[s]):6.1f}; rel err median {np.median([r_[3]['relerr'] for r_, ss in zip(rows, s) if ss]):.1e}; "
               f"exact |c2| median {np.median([abs(r_[3]['m']**2*r_[3]['rho2']) for r_, ss in zip(rows, s) if ss]):.0f}")
        for J in Js:
            s = np.array([r_[1] == J for r_ in rows])
            pr(f"     J={J:2d}: {s.sum():3d} pts, |Dc2| rms {np.sqrt(np.mean(dc2[s]**2)):7.1f}, max {np.abs(dc2[s]).max():7.1f}")
        if xk == 'Xn3':
            pr("   n(J) = " + ", ".join(f"{J}:{np.exp(sol[len(nus) + j]):.4f}" for j, J in enumerate(Jl)))


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run()
    else:
        lines = []
        def pr(s=""):
            print(s, flush=True); lines.append(s)
        analyse(pr)
        open(os.path.join(HERE, 'cache', 'n3_scalar_validate.txt'), 'w').write("\n".join(lines) + "\n")
