r"""
rh_t1_collinear (R4ARC_HEFT.md sec. 2): the collinear (a = b) behaviour of borderline/top-Regge light functions, and the Mellin
mechanism behind it, on the product family T' = sigma x T (LIGHTSIDE_DIST), per Wick class of the spin-2 factor.

(A) scan: symmetrised term-wise light function V(dphi) at fixed mean angle (P1's mean 0.1 and P1's (2,4) polarisations), dphi -> 0.
    Prediction (og gauge: y = 2 sin^2(dphi/2) enters only through X13 = 2y t1 t3, the (13) links and f_k):
      class 1243 (one (13) free-field link, (13)-channel twist 1):  V = y^{-3/2} (A + B log 8y) + ...   (p = 3)
      classes s, t (no (13) link, (13)-channel content double-trace like): V = y^{-1} (A + B log 8y) + ...   (p = 2)
(B) premise: the d13 integrand Gamma(d13)^2 sum_n Res_n(d13) has no singularity in -1 < Re d13 <= 2/3 except the double pole of
    Gamma(d13)^2 at 0: compare the line integrals at Re d13 = 2/3, 0.35, 0.15 and at Re d13 = -0.35 plus the residue at 0.
Usage: python3 rh_t1_collinear.py scan <k> ; python3 rh_t1_collinear.py shift <k> <dphi> ; python3 rh_t1_collinear.py report
"""
import os, sys, time, pickle
import numpy as np
from scipy.special import loggamma
import rh_common as R

PHIBAR = 0.1
DPHIS = (1.2, 0.6, 0.3, 0.15, 0.075, 0.0375)
AMPS = {2: {'1243': ('D4',), '1324': ('D4',), 's': ('D8',), 't': ('D8',)}, 4: {'s': ('D8',), 't': ('D8',)}}
PAIRS = {'1243+1324': ('1243', '1324'), 's+t': ('s', 't')}
NC = {'1243': 160, '1234': 160}


def cfg(dphi):
    _, _, p2, p4 = R.P.config('P1')
    return PHIBAR + dphi/2, PHIBAR - dphi/2, p2, p4


def line_value(Sab, Sba, k, amps, c2, Y=8.0, h=0.2):
    """symmetrised term-wise values per (class, amp) on the line Re d13 = c2 (residues computed once per amp, all classes)."""
    d13s = c2 + 1j*np.arange(-Y, Y + h/2, h)
    need = {}
    for cl, names in amps.items():
        for nm in names:
            need.setdefault(nm, []).append(cl)
    out = {}
    for d13 in d13s:
        w = np.exp(2*loggamma(d13))*h/(2*np.pi)
        for o, S in (('ab', Sab), ('ba', Sba)):
            sg = 1.0 if o == 'ab' else (-1.0)**k
            for nm, cls in need.items():
                N = max(NC.get(cl, 12) for cl in cls)
                r = R.residues(S, R.poly_D(int(nm[1:])), d13, N)
                for cl in cls:
                    out[(cl, nm)] = out.get((cl, nm), 0) + sg*w*r[cl][:NC.get(cl, 12) + 1]
    res = {}
    for key, r in out.items():
        cl = key[0]
        if cl in NC:
            pp, tl = R.tail(r, 80)
            res[key] = np.sum(r) + (tl if np.isfinite(tl) else 0)
        else:
            res[key] = np.sum(r)
    return res


def residue0(Sab, Sba, k, amps, rad=0.2, npts=64):
    """(1/2 pi i) oint_{|d13| = rad} Gamma(d13)^2 sum_n Res_n  (the double pole at d13 = 0), per (class, amp)."""
    th = 2*np.pi*np.arange(npts)/npts
    need = {}
    for cl, names in amps.items():
        for nm in names:
            need.setdefault(nm, []).append(cl)
    out = {}
    for t in th:
        d13 = rad*np.exp(1j*t)
        wgt = np.exp(2*loggamma(d13))*d13/npts            # (1/2 pi i) dd13 = d13 dtheta/(2 pi)
        for o, S in (('ab', Sab), ('ba', Sba)):
            sg = 1.0 if o == 'ab' else (-1.0)**k
            for nm, cls in need.items():
                N = max(NC.get(cl, 12) for cl in cls)
                r = R.residues(S, R.poly_D(int(nm[1:])), d13, N)
                for cl in cls:
                    out[(cl, nm)] = out.get((cl, nm), 0) + sg*wgt*np.sum(r[cl][:NC.get(cl, 12) + 1])
    return out


def scan(k):
    out = {}
    for dphi in DPHIS:
        t0 = time.time()
        a, b, p2, p4 = cfg(dphi)
        cls = tuple(AMPS[k])
        Sab = R.setup(k, a, b, p2, p4, classes=cls); Sba = R.setup(k, b, a, p2, p4, classes=cls)
        out[dphi] = line_value(Sab, Sba, k, AMPS[k], 2/3)
        print(f'k {k} dphi {dphi}: ' + ', '.join(f'{key}: {v.real:+.6e}' for key, v in out[dphi].items()) + f' [{time.time()-t0:.0f}s]', flush=True)
    pickle.dump(out, open(f'{R.CACHE}/t1_scan_k{k}.pkl', 'wb'))


def shift(k, dphi):
    a, b, p2, p4 = cfg(dphi)
    cls = tuple(AMPS[k])
    Sab = R.setup(k, a, b, p2, p4, classes=cls); Sba = R.setup(k, b, a, p2, p4, classes=cls)
    out = {c2: line_value(Sab, Sba, k, AMPS[k], c2, h=0.1) for c2 in (2/3, 0.4, -0.4)}
    out['res0'] = residue0(Sab, Sba, k, AMPS[k])
    for key in out[2/3]:
        print(f'k {k} dphi {dphi} {key}: lines 2/3 {out[2/3][key].real:+.10e}, 0.4 {out[0.4][key].real:+.10e};'
              f' -0.4 + Res0 = {(out[-0.4][key] + out["res0"][key]).real:+.10e} (Res0 {out["res0"][key].real:+.4e})', flush=True)
    pickle.dump(out, open(f'{R.CACHE}/t1_shift_k{k}_{dphi}.pkl', 'wb'))


def report():
    L = []
    pr = lambda s: (print(s), L.append(s))
    for k in (2, 4):
        try:
            D = pickle.load(open(f'{R.CACHE}/t1_scan_k{k}.pkl', 'rb'))
        except FileNotFoundError:
            continue
        keys = list(D[DPHIS[0]])
        units = [(p, cls) for p, cls in PAIRS.items() if all(any(kk[0] == c for kk in keys) for c in cls)] + [(c, (c,)) for c in sorted(set(kk[0] for kk in keys))]
        for p, cls in units:
            nm = [kk[1] for kk in keys if kk[0] == cls[0]][0]
            m = -1.5 if p in ('1243+1324', '1243', '1324') else -1.0
            pr(f'k = {k}, {"pair" if "+" in p else "class"} {p}, top-Regge term {nm}: predicted V ~ y^{m} (A + B log 8y) + O(y^({m}+1) log y)')
            ys = np.array([2*np.sin(d/2)**2 for d in DPHIS]); V = np.array([sum(D[d][(c, nm)] for c in cls).real for d in DPHIS])
            loc = [np.log(abs(V[i+1]/V[i]))/np.log(DPHIS[i+1]/DPHIS[i]) for i in range(len(DPHIS) - 1)]
            X = np.vstack([np.ones(4), np.log(8*ys[-4:]), ys[-4:], ys[-4:]*np.log(8*ys[-4:])]).T
            coef, *_ = np.linalg.lstsq(X, V[-4:]*ys[-4:]**(-m), rcond=None)
            for i, d in enumerate(DPHIS):
                pr(f'   dphi {d:7.4f}  y {ys[i]:.3e}  V {V[i]:+.6e}   V*y^{-m:+.1f} = {V[i]*ys[i]**(-m):+.6e}' +
                   (f'   local slope -d log|V|/d log dphi -> next: {-loc[i]:+.3f}' if i < len(loc) else ''))
            pr(f'   fit V y^{-m:+.1f} = A + B log 8y + C y + D y log 8y on the 4 smallest dphi: A = {coef[0]:+.5e}, B = {coef[1]:+.5e}')
    for f in sorted(os.listdir(R.CACHE)):
        if f.startswith('t1_shift_') and f.endswith('.pkl'):
            D = pickle.load(open(f'{R.CACHE}/{f}', 'rb'))
            for key in D[2/3]:
                ref = D[2/3][key]
                pr(f'{f[:-4]} {key}: I(2/3) = {ref.real:+.8e}; |I(0.4)/I(2/3) - 1| = {abs(D[0.4][key]/ref - 1):.1e};'
                   f' |(I(-0.4) + Res0)/I(2/3) - 1| = {abs((D[-0.4][key] + D["res0"][key])/ref - 1):.1e};  Res0/I(2/3) = {(D["res0"][key]/ref).real:+.5f},'
                   f' I(-0.4)/I(2/3) = {(D[-0.4][key]/ref).real:+.5f}')
    open(f'{R.CACHE}/t1_report.txt', 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'scan':
        scan(int(sys.argv[2]))
    elif sys.argv[1] == 'shift':
        shift(int(sys.argv[2]), float(sys.argv[3]))
    else:
        report()
