"""
nf_reprcut (NEARGAP_FILL.md §3): where does the representation probe's deviation sit?  For the odd-J direct-bulk twist edges, the
shape-fit probes (degree K = 4..8, half-samples; OLS as in hs_odd, or WLS) are split by x-window: the extrapolated window x < x50 =
50/m (nu < 50: below the first bulk node at x = 0.05, the fit is a polynomial CONTINUATION; in the final budget this window is REMOVED
from the central value and covered by sn_odd's measured envelope) and the fitted bulk [x50, 0.8).  E_repr evaluated on lam_c (window
removed) = the quantity the budget actually bounds.  Jackknife over the 8 replicas for every difference.
Usage: python3 nf_reprcut.py J [mfl,...]  -> cache/nf_reprcut_J<J>.txt / .json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nf_guard  # noqa: F401
sys.path.insert(0, os.path.join(HERE, '..', 'extfold'))
import ef_core as C
import numpy as np
from numpy.polynomial.legendre import leggauss
H = C.H
import hs_odd as HO
Dg = 1000.0


def win(J, mfl, M, a, b, nb=400):
    """bulk-model contribution to lam/N_s on x in [a, b]."""
    hs, base, Ns, st = H.base_of(Dg, J, mfl)
    t, w = leggauss(nb); xb = a + 0.5*(b - a)*(t + 1); wb = 0.5*(b - a)*w
    w2, w3 = H.state_weights(Dg, J, mfl, '+-', xb)
    return float(((wb*w2*M.k['C2']['bulk'](xb)).sum() + (wb*w3*M.k['C3H']['bulk'](xb)).sum())/Ns[0])


if __name__ == '__main__':
    J = int(sys.argv[1]); D = 1000.0 + J
    C.SET[0] = 'all'
    rowm = H.HN.mphys(D, J); Lm = H.E.lam_of(Dg); mtw = rowm/Lm; HO.XROW[0] = 1/mtw
    mfls = [mtw] + ([float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [])
    rK, rI, isdir = HO.rows(J, D); bad = H.refbad(J, D, HO.RTOL)
    reps = [HO.rows(J, D, jk=j)[1] for j in range(8)]
    dset = set(np.round(H.direct_x(J, D), 6)) - H.failed_x(J, D)
    isd = lambda rr: np.array([round(x, 6) in dset for x in rr['xs']])
    probes = {'K4': dict(kb=4), 'K5': dict(kb=5), 'K6': dict(kb=6), 'K7': dict(kb=7), 'K8': dict(kb=8),
              'half0': dict(kb=4, sub=0), 'half1': dict(kb=4, sub=1)}
    Mf = {k: HO.Model(J, rK, rI, isdir, bad, bulk='dir', **v) for k, v in probes.items()}
    Mj = {k: [HO.Model(J, rK, rr, isd(rr), bad, bulk='dir', **v) for rr in reps] for k, v in probes.items()}
    lines, res = [], {}

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# nf_reprcut J={J}: lam = full value (hs_odd 'gk'); lam_c = lam minus the bulk model on x < x50 = 50/m (the window that sn_odd removes); "
       "low = the model's own contribution on x < x50; per probe: value, difference to K6 +- jackknife sigma of the difference")
    for mfl in mfls:
        x50 = 50.0/(mfl*Lm)
        o = {}
        for k in probes:
            lam = HO.value(J, mfl, Mf[k], 'gk')[0]; low = win(J, mfl, Mf[k], 0.0, x50)
            lj = np.array([HO.value(J, mfl, m_, 'gk')[0] for m_ in Mj[k]]); wj = np.array([win(J, mfl, m_, 0.0, x50) for m_ in Mj[k]])
            o[k] = dict(lam=lam, low=low, lam_c=lam - low, lj=lj, lcj=lj - wj)
        r = {}
        for k in probes:
            d = o[k]['lam'] - o['K6']['lam']; dc = o[k]['lam_c'] - o['K6']['lam_c']
            sd = HO.jk_sigma(o[k]['lj'] - o['K6']['lj']); sdc = HO.jk_sigma(o[k]['lcj'] - o['K6']['lcj'])
            r[k] = dict(lam=o[k]['lam'], lam_c=o[k]['lam_c'], low=o[k]['low'], d=d, sd=sd, dc=dc, sdc=sdc,
                        sig=HO.jk_sigma(o[k]['lj']), sigc=HO.jk_sigma(o[k]['lcj']))
        res[f'{mfl:.7f}'] = r
        pr(f"  m/Lam={mfl:.7f} (x50 = {x50:.4f}):")
        for k, v in r.items():
            pr(f"    {k:6s}: lam {v['lam']:+.4f} (jk {v['sig']:.4f}) | low window {v['low']:+.4f} | lam_c {v['lam_c']:+.4f} (jk {v['sigc']:.4f}) | "
               f"lam-K6 {v['d']:+.4f} +- {v['sd']:.4f} | lam_c-K6 {v['dc']:+.4f} +- {v['sdc']:.4f}")
        pp = ['K4', 'K5', 'K7', 'half0', 'half1']
        pr(f"    E_repr (protocol, 3 max|lam - lam_K6|) = {3*max(abs(r[k]['d']) for k in pp):.3f};  on lam_c: {3*max(abs(r[k]['dc']) for k in pp):.3f}")
    open(os.path.join(HERE, 'cache', 'nf_reprcut_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(HERE, 'cache', 'nf_reprcut_J%d.json' % J), 'w'), indent=1, default=float)
