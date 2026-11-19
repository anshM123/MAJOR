"""
n4_validate (NEARGAP100_SMALLNU.md): scalar validation of the A+ + A- thimble sums (n4_run 'scal' points, cache/n4_scal.pkl; and the
scalar on the common samples of the graviton production points, cache/n4_raw.pkl) against the EXACT Mellin object:
   d_D = log Y_D - log(q_D X_D)      (X = g C(nu+i) + g~ C(nu-i), n3_exact; records from e1 or n4_exact)
must be a (nu, J)-constant (CLS-D normalisation N(nu) and the J-only constant of Q_J).  Reported per (J, nu) group, relative to a
reference block (Delta_ref = 2e4 if present, else the first): Re and Im of d_D - d_ref with replica errors (common random numbers ->
replica-wise differences), and Delta c2 = m_D^2 Re(d_D - d_ref) ex2/E2 (c2 units, as n3_scalar_validate).
Also: the same with A+ alone / A- alone (tests the intersection numbers).
Usage: python3 n4_validate.py  -> cache/n4_validate.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import n3_exact as NX
import n4_exact as E4


def exact_rec(D_, J):
    db = NX.db(100)
    for k, r in db.items():
        if k[1] == J and abs(r['Delta'] - D_) < 1e-6:
            return r, r.get('logq')
    r = E4.get(D_, J)
    return r, r['logq']


def logq(D_, J, r):
    sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
    from exactmom import log_q
    return log_q(D_, J, 3.0)


def groups():
    out = {}
    f = os.path.join(HERE, 'cache', 'n4_scal.pkl')
    if os.path.exists(f):
        for (D_, J, nu), raw in pickle.load(open(f, 'rb')).items():
            out.setdefault(('scal', J, nu), {})[D_] = raw
    f = os.path.join(HERE, 'cache', 'n4_raw.pkl')
    if os.path.exists(f):
        for (mode, D_, J, nu), raw in pickle.load(open(f, 'rb')).items():
            if J % 2 == 1:
                continue        # odd J: the scalar n0 = +1 integral is not the identical-scalar functional (S^0 sum) -> no exact object
            out.setdefault((mode, J, nu), {})[D_] = raw
    return out


def analyse(pr):
    G = groups()
    allres = []
    for (mode, J, nu), dd in sorted(G.items(), key=lambda t: (t[0][1], t[0][2], t[0][0])):
        Ds = sorted(dd)
        ref = 20000.0 if 20000.0 in dd else Ds[0]
        rows = {}
        for D_ in Ds:
            raw = dd[D_]
            r, _ = exact_rec(D_, J)
            lq = logq(D_, J, r)
            Xv = NX.X(nu, r)[0]
            lx = lq + np.log(complex(Xv))
            th = raw['th']
            parts = {}
            for wh in ('A+', 'A-'):
                S = th[wh]['S']
                parts[wh] = (S['logN'] - raw['logbeta'], S['repsS'])
            lp, rp = parts['A+']; lm, rm = parts['A-']
            reps = {'sum': lp + np.log(rp + rm*np.exp(lm - lp)), 'A+': lp + np.log(rp), 'A-': lm + np.log(rm)}
            x = nu/r['m']; E2 = (2 - x*x)/(1 - x*x)**2
            ex2 = float(np.interp(nu, r['nus'], r['ex'][2])) if 'ex' in r else E2
            rows[D_] = dict(reps=reps, lx=lx, m=r['m'], fac=r['m']**2*ex2/E2,
                            relerr={wh: abs(np.std(th[wh]['S']['repsS'], ddof=1)/np.sqrt(len(th[wh]['S']['repsS']))/th[wh]['S']['Es'])
                                    for wh in ('A+', 'A-')})
        pr(f"## {mode} J={J} nu={nu:.4f}: reference block Delta={ref:g}; thimble rel. errors (A+/A-): " +
           ", ".join(f"{D_:g}: {rows[D_]['relerr']['A+']:.3f}/{rows[D_]['relerr']['A-']:.3f}" for D_ in Ds))
        for lab in ('sum', 'A+', 'A-'):
            for D_ in Ds:
                if D_ == ref:
                    continue
                # replica values Y_i (complex, common seeds): ratio of replica MEANS, jackknife errors (Re and Im of the log)
                ys = np.exp(rows[D_]['reps'][lab] - rows[D_]['reps'][lab][0]); yr = np.exp(rows[ref]['reps'][lab] - rows[ref]['reps'][lab][0])
                off = rows[D_]['reps'][lab][0] - rows[ref]['reps'][lab][0] - (rows[D_]['lx'] - rows[ref]['lx'])
                def lr(sel):
                    return np.log(ys[sel].mean()/yr[sel].mean()) + off
                nR = len(ys); allm = np.ones(nR, bool)
                dm = lr(allm)
                jk = np.array([lr(np.arange(nR) != i) for i in range(nR)])
                jk = jk.real + 1j*((jk.imag - dm.imag + np.pi) % (2*np.pi) - np.pi + dm.imag)
                err = complex(np.sqrt((nR - 1)/nR*np.sum((jk.real - jk.real.mean())**2)),
                              np.sqrt((nR - 1)/nR*np.sum((jk.imag - jk.imag.mean())**2)))
                dm = complex(dm.real, (dm.imag + np.pi) % (2*np.pi) - np.pi)
                dc2 = rows[D_]['fac']*dm.real
                ec2 = abs(rows[D_]['fac'])*err.real
                pr(f"   [{lab:3s}] Delta={D_:9.4f} (m={rows[D_]['m']:7.2f}): d - d_ref = {dm.real:+.4f} {dm.imag:+.4f}i  (replica err "
                   f"{err.real:.4f}, {err.imag:.4f});  Delta c2 = {dc2:+8.1f} +- {ec2:6.1f};  m^2 (d - d_ref) = {rows[D_]['m']**2*dm.real:+.0f} {rows[D_]['m']**2*dm.imag:+.0f}i")
                allres.append((mode, J, nu, D_, lab, dm, err, dc2, ec2))
    return allres


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    analyse(pr)
    open(os.path.join(HERE, 'cache', 'n4_validate.txt'), 'w').write("\n".join(lines) + "\n")
