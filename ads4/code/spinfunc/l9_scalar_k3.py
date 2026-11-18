"""
l9_scalar_k3 (M6 round 5, doc 23.1): accuracy of the CLS-D crystal Laplace series at J^2/m = O(1), measured on the EXACT scalar skeleton.
Same comparison as l2_scalar_validate (exact e1 Mellin data, CHMRS magic-block X(nu), fit log|Y/(qX)| = a_nu + b_J, residuals in c2
units Delta c2 = m^2 eps ex/E_2), but
  * prefactor form (Q_J polynomial prefactor, l1_crystal.build(pref=True)), orders S0, S2, S4, S6 kept separately (S6 = K3, l9_tps);
  * J = 0 ... 10 (incl. J = 8, 10, absent from l2), the near-gap states m/Lam = 1.0 ... 1.2 at Delta_g = 100 (m ~ 100 - 125).
Reported per J: rms / max |Delta c2| for truncation K1, K2, K3 in linear and cumulant (log) form; the ratio of successive orders.
Usage: L1_DEG=8 python3 l9_scalar_k3.py Dg [step]   -> cache/l9_scal_<Dg>.pkl (resumable), cache/l9_scalar_k3_<Dg>.txt
"""
import os, sys, pickle, time
os.environ.setdefault('L1_DEG', '8')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
import l1_crystal as C
import l9_tps as K
import l2_scalar_validate as V
from exactmom import log_q

JS = (0, 2, 4, 6, 8, 10)
MFL = (1.0, 1.05, 1.1, 1.2)
XW = (0.3, 0.7)


def laplace6(nu, Delta, J):
    zs, a = C.saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF, G, Q = C.build(zv, nu, Delta, J, pref=True)
    _, HF = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2, G2, Q2 = C.build(T.linsub(zs, L), nu, Delta, J, pref=True)
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    lamF2, G2, Q2 = C.build(T.linsub(zs, L), nu, Delta, J, pref=True)
    S = K.expand_pref6(lamF2, G2, Q2)
    logY0 = lamF[0] + G[0] + 3.5*np.log(2*np.pi) + np.log(np.linalg.det(L)) - C.log_beta(Delta, J)
    return logY0, np.array(S)


def run(Dg, step=1):
    db = pickle.load(open(os.path.join(V.EC, 'e1_%d_Df3.pkl' % Dg), 'rb'))
    pk = os.path.join(HERE, 'cache', 'l9_scal_%d.pkl' % Dg)
    res = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    t0 = time.time(); n = 0
    for J in JS:
        for mfl in MFL:
            key = (mfl, J)
            if key not in db:
                continue
            r = db[key]
            idx = [i for i, x in enumerate(r['x']) if XW[0] <= x <= XW[1]][::step]
            for i in idx:
                if (key, i) in res:
                    continue
                nu = r['nus'][i]
                lY, S = laplace6(nu, r['Delta'], J)
                res[(key, i)] = dict(mfl=mfl, J=J, m=r['m'], Delta=r['Delta'], x=r['x'][i], nu=nu, ex2=r['ex'][2][i],
                                     rho2=r['rho'][2][i], logY0=lY, S=S, logq=log_q(r['Delta'], J, 3.0))
                n += 1
                if n % 10 == 0:
                    pickle.dump(res, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
                    print(f"  Dg={Dg}: {n} new, {time.time() - t0:.0f}s (J={J}, m/Lam={mfl})", flush=True)
    pickle.dump(res, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    return res


def logform(S, order):
    a2, a4, a6 = S[1]/S[0], S[2]/S[0], S[3]/S[0]
    v = np.log(S[0])
    if order >= 1:
        v = v + a2
    if order >= 2:
        v = v + a4 - a2*a2/2
    if order >= 3:
        v = v + a6 - a2*a4 + a2**3/3
    return v


def analyse(Dg, res, pr, Js=JS):
    keys = [k for k in sorted(res) if res[k]['J'] in Js]
    for k in keys:
        q = res[k]; q['X'], _ = V.Xexact(q['nu'], q['m'], q['J'], q['ex2'])
    nus = sorted(set(round(res[k]['nu'], 9) for k in keys)); iv = {v: j for j, v in enumerate(nus)}
    Jl = [J for J in Js if J != 0]
    pr(f"## Dg = {Dg}: {len(keys)} (state, node) points, J in {Js}, m/Lam in {MFL}, x in {XW}")
    out = {}
    for form in ('lin', 'log'):
        for order in (1, 2, 3):
            rows = []
            for k in keys:
                q = res[k]
                if form == 'lin':
                    ly = q['logY0'].real + np.log(abs(sum(q['S'][:order + 1])))
                else:
                    ly = (q['logY0'] + logform(q['S'], order)).real
                rows.append((round(q['nu'], 9), q['J'], ly - q['logq'] - np.log(abs(q['X'])), q))
            A = np.zeros((len(rows), len(nus) + len(Jl))); b = np.array([r_[2] for r_ in rows])
            for i, (v, J, lr, q) in enumerate(rows):
                A[i, iv[v]] = 1
                if J != 0:
                    A[i, len(nus) + Jl.index(J)] = 1
            w = np.array([q['m']**2*abs(q['ex2'])/V.E2(q['x']) for (v, J, lr, q) in rows])
            sol, *_ = np.linalg.lstsq(A*w[:, None], b*w, rcond=None)
            eps = b - A @ sol
            dc2 = np.array([q['m']**2*e*q['ex2']/V.E2(q['x']) for e, (v, J, lr, q) in zip(eps, rows)])
            c2 = np.array([q['m']**2*q['rho2'] for (v, J, lr, q) in rows])
            s = f"  {form} K{order}: rms |Dc2| = {np.sqrt(np.mean(dc2**2)):7.2f};  per J (rms / max):"
            for J in Js:
                sel = np.array([r_[1] == J for r_ in rows])
                s += f"  J={J}: {np.sqrt(np.mean(dc2[sel]**2)):.1f}/{np.abs(dc2[sel]).max():.1f}"
                out[(form, order, J)] = (np.sqrt(np.mean(dc2[sel]**2)), np.abs(dc2[sel]).max(), np.median(np.abs(c2[sel])))
            pr(s)
    pr("  exact |c2| median per J: " + "  ".join(f"J={J}: {out[('lin', 3, J)][2]:.0f}" for J in Js))
    # successive-order sizes (state-dependent part is what matters; report raw |S_k/S_0| medians)
    for J in Js:
        sel = [res[k] for k in keys if res[k]['J'] == J]
        r = np.array([[abs(q['S'][o]/q['S'][0]) for o in (1, 2, 3)] for q in sel])
        pr(f"  J={J}: median |S2/S0|, |S4/S0|, |S6/S0| = " + ", ".join(f"{v:.3g}" for v in np.median(r, 0)) +
           f";  m^3 |S6/S0| median {np.median(r[:, 2]*np.array([q['m']**3 for q in sel])):.3g}")
    return out


if __name__ == '__main__':
    Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    step = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    res = run(Dg, step)
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    analyse(Dg, res, pr)
    open(os.path.join(HERE, 'cache', 'l9_scalar_k3_%d.txt' % Dg), 'w').write("\n".join(lines) + "\n")
