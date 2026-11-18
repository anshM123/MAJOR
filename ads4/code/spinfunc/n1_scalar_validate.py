"""
n1_scalar_validate (doc 24.3): validation of the non-perturbative crystal integration (n1_direct: re-centred log-coordinate contour,
path-continued logs, masks, RQMC) on the EXACT scalar skeleton, at the near-gap states that matter (Delta_g = 100, m/Lam = 1.0-1.2,
J = 0 ... 10) -- the same (state, node) points and the same comparison as l9_scalar_k3 (CHMRS magic-block X(nu) from exact e1 Mellin
data, fit log|Y/(qX)| = a_nu + b_J, residuals in c2 units Delta c2 = m^2 eps ex/E_2), with Y = the direct integral instead of the
truncated Laplace series.  Extra points outside x in [0.3, 0.7] (x in [0.2, 0.3) and (0.7, 0.85]) are added (N1_XWIDE=1).
Usage: OMP_NUM_THREADS=1 python3 n1_scalar_validate.py run J [J ...]   (resumable, cache/n1_scal_100.pkl)
       python3 n1_scalar_validate.py analyse                            -> cache/n1_scalar_validate.txt
"""
import os, sys, pickle, time, warnings
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
import numpy as np

PK = os.path.join(HERE, 'cache', 'n1_scal_100.pkl')
MQ = int(os.environ.get('N1_M', '13'))
XWIDE = (0.2, 0.85)


def pts(Js):
    import l2_scalar_validate as V
    db = pickle.load(open(os.path.join(V.EC, 'e1_100_Df3.pkl'), 'rb'))
    out = []
    for J in Js:
        for mfl in (1.0, 1.05, 1.1, 1.2):
            if (mfl, J) not in db:
                continue
            r = db[(mfl, J)]
            for i, x in enumerate(r['x']):
                if XWIDE[0] <= x <= XWIDE[1] and r['nus'][i] >= 20.0:
                    out.append(((mfl, J), i, r))
    return out


def run(Js):
    import l2_scalar_validate  # noqa: F401  (adds the exactskel / uniform paths)
    import n1_direct as D
    import l1_crystal as Cr
    from exactmom import log_q
    res = pickle.load(open(PK, 'rb')) if os.path.exists(PK) else {}
    for key, i, r in pts(Js):
        if (key, i) in res:
            continue
        t = time.time()
        nu, Dl, J = r['nus'][i], r['Delta'], key[1]
        c, wc, L = D.full_saddle(nu, Dl, J)
        E, err, nm = D.integrate_qmc(None, c, wc, L, nu, Dl, J, MQ, 8, 1.2)
        lN, _ = D.norm_w(lambda zl: D.PW.Cr.build(zl, nu, Dl, J, pref=True), c, wc, L)
        rec = dict(mfl=key[0], J=J, m=r['m'], Delta=Dl, x=r['x'][i], nu=nu, ex2=r['ex'][2][i], rho2=r['rho'][2][i],
                   logY=lN - Cr.log_beta(Dl, J) + np.log(E), relerr=abs(err/E), nmask=nm, logq=log_q(Dl, J, 3.0))
        res = pickle.load(open(PK, 'rb')) if os.path.exists(PK) else {}
        res[(key, i)] = rec
        pickle.dump(res, open(PK + '.tmp%d' % os.getpid(), 'wb')); os.replace(PK + '.tmp%d' % os.getpid(), PK)
        print(f"J={J} m/Lam={key[0]} x={rec['x']:.3f} nu={nu:.2f}: rel err {rec['relerr']:.1e}, nmask {nm}, {time.time() - t:.0f}s", flush=True)


def analyse(pr):
    import l2_scalar_validate as V
    res = pickle.load(open(PK, 'rb'))
    k3 = pickle.load(open(os.path.join(HERE, 'cache', 'l9_scal_100.pkl'), 'rb'))
    for lab, win in (('x in [0.3, 0.7] (the l9_scalar_k3 points)', (0.3, 0.7)), ('all x in [0.2, 0.85], nu >= 20', XWIDE)):
        allk = [k for k in sorted(res) if win[0] <= res[k]['x'] <= win[1]]
        keys = [k for k in allk if res[k]['relerr'] <= 0.02]
        pr(f"   ({len(allk) - len(keys)} of {len(allk)} points dropped: RQMC relative error > 2 %, all at nu <= "
           f"{max([res[k]['nu'] for k in allk if res[k]['relerr'] > 0.02] or [0]):.1f})")
        Js = sorted(set(res[k]['J'] for k in keys))
        for k in keys:
            q = res[k]; q['X'], _ = V.Xexact(q['nu'], q['m'], q['J'], q['ex2'])
        nus = sorted(set(round(res[k]['nu'], 9) for k in keys)); iv = {v: j for j, v in enumerate(nus)}
        Jl = [J for J in Js if J != 0]
        def fitres(getly, kk):
            rows = [(round(res[k]['nu'], 9), res[k]['J'], getly(k) - res[k]['logq'] - np.log(abs(res[k]['X'])), res[k]) for k in kk]
            A = np.zeros((len(rows), len(nus) + len(Jl))); b = np.array([r_[2] for r_ in rows])
            for i, (v, J, lr, q) in enumerate(rows):
                A[i, iv[v]] = 1
                if J != 0:
                    A[i, len(nus) + Jl.index(J)] = 1
            w = np.array([q['m']**2*abs(q['ex2'])/V.E2(q['x']) for (v, J, lr, q) in rows])
            sol, *_ = np.linalg.lstsq(A*w[:, None], b*w, rcond=None)
            eps = b - A @ sol
            dc2 = np.array([q['m']**2*e*q['ex2']/V.E2(q['x']) for e, (v, J, lr, q) in zip(eps, rows)])
            mc = np.array([q['m']**2*q.get('relerr', 0.0)*abs(q['ex2'])/V.E2(q['x']) for (v, J, lr, q) in rows])
            return rows, dc2, mc
        rows, dc2, mc = fitres(lambda k: res[k]['logY'].real, keys)
        pr(f"## direct vs exact skeleton, {lab}: {len(keys)} points, J in {Js}")
        pr(f"   residual |Delta c2| (c2 units): rms {np.sqrt(np.mean(dc2**2)):.2f}, max {np.abs(dc2).max():.2f};  RQMC error (c2 units): "
           f"median {np.median(mc):.2f}, max {mc.max():.2f}")
        for J in Js:
            s = np.array([r_[1] == J for r_ in rows])
            pr(f"     J={J:2d}: {s.sum():3d} pts, rms {np.sqrt(np.mean(dc2[s]**2)):6.2f}, max {np.abs(dc2[s]).max():6.2f}, "
               f"median RQMC err {np.median(mc[s]):5.2f};  exact |c2| median {np.median([abs(q['m']**2*q['rho2']) for (v, JJ, lr, q) in rows if JJ == J]):.0f}")
        # same points with K3 (linear form, l9_scalar_k3) for comparison
        kk3 = [k for k in keys if k in k3]
        if kk3 and win == (0.3, 0.7):
            for k in kk3:
                k3[k]['X'] = res[k]['X']
            rows3, d3, _ = fitres(lambda k: (k3[k]['logY0'].real + np.log(abs(np.sum(k3[k]['S'])))), kk3)
            rowsd, dd, _ = fitres(lambda k: res[k]['logY'].real, kk3)
            pr(f"   on the {len(kk3)} common points: direct rms {np.sqrt(np.mean(dd**2)):.2f} / max {np.abs(dd).max():.2f};  K3 (linear) rms "
               f"{np.sqrt(np.mean(d3**2)):.2f} / max {np.abs(d3).max():.2f}")
            for J in Js:
                s = np.array([r_[1] == J for r_ in rowsd])
                if s.any():
                    pr(f"     J={J:2d}: direct rms/max {np.sqrt(np.mean(dd[s]**2)):6.2f}/{np.abs(dd[s]).max():6.2f};  K3 {np.sqrt(np.mean(d3[s]**2)):6.2f}/{np.abs(d3[s]).max():6.2f}")


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run([int(a) for a in sys.argv[2:]])
    else:
        lines = []
        def pr(s=""):
            print(s, flush=True); lines.append(s)
        analyse(pr)
        open(os.path.join(HERE, 'cache', 'n1_scalar_validate.txt'), 'w').write("\n".join(lines) + "\n")
