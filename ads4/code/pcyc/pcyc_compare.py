"""
pcyc_compare (PCYC_CLOSE.md section 3.3): the direct test of P_cyc for the EE form,
    T(nu, Delta) = <Q>_crys / <Q>_D' - 1     (P_cyc => T = 0 exactly),
from the crystal-cycle runs (pcyc_crys, cache/cr_*.pkl) and the D' quadrature (pcyc_dquad, cache/dq_*.npz; finest grid available, and
the quadrature error estimated from the difference to the next-coarser grid).  Also |T| against the failure threshold scaled from
Delta_g = 300 / 1000:  c2-constant scaling eps(Delta) = 2e-4 (300/Delta)^2 (the RD_CLOSE failure scenario), and the Delta-independent
floor 2e-5.
Usage: python3 pcyc_compare.py -> cache/pcyc_compare.txt
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np


def dq(Delta, kind):
    fs = glob.glob(os.path.join(HERE, 'cache', f'dq_D{Delta:g}_{kind}_*.npz'))
    out = []
    for f in fs:
        d = np.load(f)
        g = os.path.basename(f)[:-4].split('_')[3]
        out.append((int(np.prod([int(v) for v in g.split('-')])), g, d['nus'], d['I']*np.exp(float(d['shift']))))
    out.sort(key=lambda t: t[0])
    return [t[1:] for t in out]


def qdp(Delta, nu):
    S = dq(Delta, 'S'); E = dq(Delta, 'E')
    if not S or not E:
        return None
    def has(lst):
        res = []
        for g, nus, I in lst:
            j = np.where(np.abs(nus - nu) < 1e-9)[0]
            if j.size:
                res.append((g, I[j[0]]))
        return res
    S = has(S); E = has(E)
    if not S or not E:
        return None
    gS, IS = S[-1]; gE, IE = E[-1]
    q = IE/IS
    err = 0.0
    if len(S) > 1 and len(E) > 1:
        err = max(abs(S[-2][1]/IS - 1)*abs(q), abs(E[-2][1]/IE - 1)*abs(q))
    return q, err, gS, gE


def cvvals():
    out = {}
    fn = os.path.join(HERE, 'cache', 'pcyc_cv.txt')
    if os.path.exists(fn):
        for ln in open(fn):
            if ' deg=2:' not in ln:
                continue
            f = ln.split(' deg=2:')[0]
            cv = ln.split('CV ')[1]
            q, e = cv.split(' +- ')
            out[f] = (complex(q.replace(' ', '')), float(e))
    return out


def main():
    CV = cvvals()
    L = ["# pcyc_compare: T = <Q>_crys/<Q>_D' - 1 (crystal: n4 A+ + A- thimbles below the fold / n3 hybrid contour above; D': product rule,",
         "#   exact in Delta; sigma = crystal jackknife (x1) (+) D' grid error; 'CV' = control-variate estimate (pcyc_cv, degree 2))",
         "#  Delta    nu  contour  <Q>_D'            D' err    estimator  Re T       Im T       sigma_T    |T|/sigma  thr(c2-scaled)"]
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'cache', 'cr_D*_nu*.pkl'))):
        if 'trial' in f:
            continue
        if 'cvtest' in f and os.path.exists(f.replace('_cvtest', '')):
            continue                                   # identical seed and size as the production run: not independent
        r = pickle.load(open(f, 'rb'))
        Dl, nu = r['Delta'], r['nu']
        d = qdp(Dl, nu)
        if d is None:
            continue
        q, qerr, gS, gE = d
        cont = 'n3' if 'S' in r['th'] else 'n4'
        thr = 2.1e-4*(300.0/Dl)**2
        ests = [('plain', r['Qcrys'], r['Qerr'])]
        b = os.path.basename(f)
        if b in CV:
            ests.append(('CV', CV[b][0], CV[b][1]))
        for name, qc, qe in ests:
            T = qc/q - 1
            sT = float(np.hypot(qe/abs(q), qerr/abs(q)))
            rows.append(dict(Delta=Dl, nu=nu, cont=cont, est=name, T=T, sT=sT, f=b, qerr=qerr/abs(q)))
            L.append(f"  {Dl:6g} {nu:5g}  {cont:3s}  {q.real:+.6e}  {qerr/abs(q):7.1e}   {name:5s}  {T.real:+.2e}  {T.imag:+.2e}  {sT:8.1e}  "
                     f"{abs(T)/sT:6.2f}     {thr:.1e}   [{b}]")
    # combined per (Delta, nu, contour): CV estimates where available (else plain), over independent seeds;
    # error = max(jackknife/sqrt(n), seed scatter/sqrt(n)) (+) D' error;  x3 column = HRS_CHECK's jackknife-inflation rule
    L.append('')
    L.append('# combined over runs/seeds (n3: CV estimate; n4: plain, the CVs do not help there): Delta nu contour nruns  Re T  Im T  sigma  |T|/sigma  thr(c2-scaled)  |T|+2sigma vs thr')
    groups = {}
    for r in rows:
        k = (r['Delta'], r['nu'], r['cont'])
        groups.setdefault(k, {}).setdefault(r['f'], {})[r['est']] = r
    for k in sorted(groups):
        runs = [(d.get('CV', d.get('plain')) if k[2] == 'n3' else d['plain']) for d in groups[k].values()]   # CVs help on n3, not on n4
        Ts = np.array([r['T'] for r in runs]); sj = np.array([r['sT'] for r in runs])
        n = len(runs)
        Tm = Ts.mean()
        sc = float(np.std(np.abs(Ts - Tm))*np.sqrt(n/(n - 1))) if n > 1 else 0.0
        sig = max(float(np.sqrt(np.sum(sj**2))/n), sc/np.sqrt(n) if n > 1 else 0.0, runs[0]['qerr'])
        thr = 2.1e-4*(300.0/k[0])**2
        L.append(f'  {k[0]:6g} {k[1]:5g} {k[2]:3s} {n:2d}   {Tm.real:+.2e} {Tm.imag:+.2e}  {sig:.1e}  {abs(Tm)/sig:5.2f}   {thr:.1e}   '
                 f'{abs(Tm) + 2*sig:.1e} ({(abs(Tm) + 2*sig)/thr:.2g} x thr)')
    open(os.path.join(HERE, 'cache', 'pcyc_compare.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    return rows


if __name__ == '__main__':
    main()
