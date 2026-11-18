"""
l7_hel_relative: the helicity-RELATIVE graviton correction of C^_2 (clean of nu-only factors and of the scalar reference):
   R(x; m) = [d_+-(x)/d_++(x)] / [d44(theta)/((1-x^2)^4 P_J(cos theta))],  normalised at x_ref,  log R = sum_k e_k(x)/m^k  (fit over Delta),
   dc2_hel(x) := Re(e2 + e1^2/2) = dc2_+-(x) - dc2_++(x) - (same at x_ref)   [the 'hel' pattern quantity of doc 19.3],
plus stability of the fit (number of terms, Delta subset) and the helicity mixing as a second-order eigenvalue shift
   s(x) = m^2 |M~_{+,-}|^2 / |M~_++ M~_+-|   (c2-equivalent of the mixing, in units of the smaller density; order 2, basis from order 0).
Usage: python3 l7_hel_relative.py [xref] -> cache/l7_hel_relative.txt
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l7_dc2_analysis as An

def main(xref=0.05):
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# l7_hel_relative: dc2_+- - dc2_++ (x_ref = {xref}); fits log R = sum_k e_k/m^k over Delta = 100...3200")
    store = {}
    for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'l7_J*.pkl')), key=lambda f: int(f.split('_J')[1][:-4])):
        J = int(fn.split('_J')[1][:-4])
        if J < 4:
            continue
        xs, out = An.analyse(J, pickle.load(open(fn, 'rb')), xref)
        if len(out) < 6:
            pr(f"  J={J}: {len(out)} Delta rows, skipped"); continue
        Dls = sorted(out); ms = np.array([out[D]['m'] for D in Dls])
        L = np.array([np.log(out[D]['+-'][:, 2]/out[D]['++'][:, 2]) for D in Dls])
        L0 = np.array([np.log(out[D]['+-'][:, 0]/out[D]['++'][:, 0]) for D in Dls])
        res = {}
        for lab, sel, nt in (('all/4', slice(0, 6), 4), ('all/5', slice(0, 6), 5), ('D>=200/4', slice(1, 6), 4), ('all/3', slice(0, 6), 3)):
            A = np.array([ms[sel]**(-k) for k in range(1, nt + 1)]).T
            e = np.linalg.lstsq(A, L[sel], rcond=None)[0]
            res[lab] = (e[0], (e[1] + 0.5*e[0]**2))
        pr(f"  J={J}: LO helicity-ratio shape |R_0(x) - 1| max over x vs Delta: " + " ".join(f"{v:.1e}" for v in np.abs(np.exp(L0) - 1).max(axis=1)))
        for lab, (e1, c2) in res.items():
            pr(f"    fit {lab:9s}: e1(x) = " + " ".join(f"{v.real:+.2f}{v.imag:+.2f}i" for v in e1))
            pr(f"                dc2_hel(x) = " + " ".join(f"{x:.2f}:{v.real:+.0f}{v.imag:+.0f}i" for x, v in zip(xs, c2)))
        # mixing as a c2-equivalent second-order shift
        sv, vs = An.rank1([res_[:, :, 0]/np.abs(res_[:, :, 0]).max() for res_ in
                           [pickle.load(open(fn, 'rb'))[(Dls[-1], x)]['M'] for x in xs]])
        raw = pickle.load(open(fn, 'rb'))
        smax = []
        for D in Dls:
            Ms0 = [raw[(D, x)]['M'][:, :, 0] for x in xs]
            sv, vs = An.rank1([M/np.abs(M).max() for M in Ms0]); Ti = np.linalg.inv(np.array(vs).T)
            s = []
            for x in xs:
                A_ = Ti @ raw[(D, x)]['M'][:, :, 2] @ Ti.T
                s.append(raw[(D, x)]['m']**2*abs(A_[0, 1])**2/abs(A_[0, 0]*A_[1, 1]))
            smax.append(s)
        smax = np.array(smax)
        pr("    mixing c2-equivalent m^2|M~+-|^2/|M~++ M~+-| at Delta=3200 per x: " + " ".join(f"{x:.2f}:{v:.0f}" for x, v in zip(xs, smax[-1]))
           + f"   (Delta=100: max {smax[0][(xs >= 0.3) & (xs <= 0.7)].max():.0f} on [0.3,0.7])")
        store[J] = dict(xs=xs, res=res, mix=smax, ms=ms)
    pickle.dump(store, open(os.path.join(HERE, 'cache', 'l7_hel_relative.pkl'), 'wb'))
    open(os.path.join(HERE, 'cache', 'l7_hel_relative.txt'), 'w').write("\n".join(lines) + "\n")

if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.05)
