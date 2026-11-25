"""
b4_odd (BATCH4_FOLD.md §3; a copy of b5fold/b5_odd.py with Delta_g = 2000): odd-J near-gap states at Delta_g = 2000 (parity-odd 1x1 +-)
on their own rows Delta = 2000 + J (J = 21, 27, 29; seam J = 23), with the copies b4_hsodd (hs_odd, Delta_g and the bulk's lower end
x = 0.025 = nu 50 parametrised) and b4_snodd (sn_odd at Delta_g = 2000).  b5_odd's description, exactly as
EXT_FOLD (ef_oddrun + ef_window) and NEARGAP_FILL §3 (nf_reprcut + nf_reprsum) did for J = 23, 27, 31:
  1. highspin/hs_odd.py via its copy b4_hsodd.py (Delta_g = 2000, bulk from x = 0.025) on the merged point set
     (b4_merge: direct bulk row + blo + rz flanks, GK15 edge, small-nu nodes);
  2. smallnu/sn_odd.main unchanged (measured nu < 50 window, FULL normalisation, x3 envelope / signed);
  3. the representation term on lam_c (NEARGAP_FILL §3): every probe (K4, K5, K7, two half-samples at K4) split at x50 = 50/m, the
     probes' continuation below the first bulk node removed, jackknife on every difference; E_repr_c = 3 max |lam_c(probe) - lam_c(K6)|;
  4. L x3 = lam_c / (E_mc + E_quad + E_miss + E_repr_c + window x3), L signed = lam_signed / (... + E_signed + E_lo), both also with
     the protocol E_repr (on lam) for comparison.
Env: HS_XCUT (default 0.8), B4_XB (bulk lower end, default 0.025), B4_RZODD, B4_SET (fix | all; odd J identical).
Usage: python3 b4_odd.py J tag [mfl,...]   -> cache/odd_<tag>/{hs_odd,sn_odd,b4_odd}_J<J>.{txt,json}
"""
import os, sys, json, runpy
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_core as C
import numpy as np
from numpy.polynomial.legendre import leggauss
H = C.H
Dg = C.DG
PROBES = ('K4', 'K5', 'K7', 'half0', 'half1')


def win(HO, J, mfl, M, a, b, nb=400):
    hs, base, Ns, st = H.base_of(Dg, J, mfl)
    t, w = leggauss(nb); xb = a + 0.5*(b - a)*(t + 1); wb = 0.5*(b - a)*w
    w2, w3 = H.state_weights(Dg, J, mfl, '+-', xb)
    return float(((wb*w2*M.k['C2']['bulk'](xb)).sum() + (wb*w3*M.k['C3H']['bulk'](xb)).sum())/Ns[0])


if __name__ == '__main__':
    J, tag = int(sys.argv[1]), sys.argv[2]
    mfls = sys.argv[3] if len(sys.argv) > 3 else None
    D = Dg + J
    C.SET[0] = os.environ.get('B4_SET', 'fix')
    out = os.path.join(C.B4C, 'odd_' + tag); os.makedirs(out, exist_ok=True)
    H.HSC = out; C.S.SNC = out; C.S.HSC = out
    sys.argv = ['hs_odd.py', str(J)] + ([mfls] if mfls else [])
    runpy.run_path(os.path.join(HERE, 'b4_hsodd.py'), run_name='__main__')
    H.HSC = out
    import b4_snodd as SO        # sn_odd at Delta_g = 2000; reads out/hs_odd_J<J>.json (FULL normalisation)
    SO.main(J)
    import b4_hsodd as HO        # functions only
    hs = json.load(open(os.path.join(out, 'hs_odd_J%d.json' % J)))['states']
    sn = json.load(open(os.path.join(out, 'sn_odd_J%d.json' % J)))['states']
    rowm = H.HN.mphys(D, J); Lm = H.E.lam_of(Dg); mtw = rowm/Lm; HO.XROW[0] = 1/mtw
    rK, rI, isdir = HO.rows(J, D); bad = H.refbad(J, D, HO.RTOL)
    reps = [HO.rows(J, D, jk=j)[1] for j in range(8)]
    dset = set(np.round(H.direct_x(J, D), 6)) - H.failed_x(J, D)
    isd = lambda rr: np.array([round(x, 6) in dset for x in rr['xs']])
    pdef = {'K4': dict(kb=4), 'K5': dict(kb=5), 'K6': dict(kb=6), 'K7': dict(kb=7), 'K8': dict(kb=8),
            'half0': dict(kb=4, sub=0), 'half1': dict(kb=4, sub=1)}
    Mf = {k: HO.Model(J, rK, rI, isdir, bad, bulk='dir', **v) for k, v in pdef.items()}
    Mj = {k: [HO.Model(J, rK, rr, isd(rr), bad, bulk='dir', **v) for rr in reps] for k, v in pdef.items()}
    lines, res = [], {}

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# b4_odd J={J} tag={tag} set={C.SET[0]} XCUT={HO.XCUT} GKMIN={HO.GKMIN}: bulk nodes {Mf['K6'].nb}, edge direct nodes {len(Mf['K6'].ex)}")
    for key in sorted(hs, key=float):
        h = hs[key]; s = sn[[k for k in sn if abs(float(k) - float(key)) < 1e-6][0]]
        mfl = float(key); edge = h['main_key'].split('/')[1]
        x50 = 50.0/(mfl*Lm)
        o = {}
        for k in pdef:
            lam = HO.value(J, mfl, Mf[k], edge)[0]; low = win(HO, J, mfl, Mf[k], 0.0, x50)
            lj = np.array([HO.value(J, mfl, m_, edge)[0] for m_ in Mj[k]]); wj = np.array([win(HO, J, mfl, m_, 0.0, x50) for m_ in Mj[k]])
            o[k] = dict(lam=lam, lam_c=lam - low, lj=lj, lcj=lj - wj)
        r = {}
        for k in pdef:
            r[k] = dict(lam=o[k]['lam'], lam_c=o[k]['lam_c'], d=o[k]['lam'] - o['K6']['lam'], dc=o[k]['lam_c'] - o['K6']['lam_c'],
                        sd=HO.jk_sigma(o[k]['lj'] - o['K6']['lj']), sdc=HO.jk_sigma(o[k]['lcj'] - o['K6']['lcj']))
        Ec = 3*max(abs(r[p]['dc']) for p in PROBES); pk = max(PROBES, key=lambda p: abs(r[p]['dc']))
        Ep = 3*max(abs(r[p]['d']) for p in PROBES)
        base_c = h['E_tot'] - h['E_repr'] + Ec
        L3, Ls = s['lamc']/(h['E_tot'] + s['E_env3']), s['lam_signed']/(h['E_tot'] + s['E_signed'] + s['E_lo'])
        L3c, Lsc = s['lamc']/(base_c + s['E_env3']), s['lam_signed']/(base_c + s['E_signed'] + s['E_lo'])
        L10c = s['lamc']/(base_c + s['E_env10'])
        res[key] = dict(J=J, mfl=mfl, main=h['main'], main_key=h['main_key'], jk=h['jk'], pub=h['pub'], lamc=s['lamc'],
                        lam_signed=s['lam_signed'], E_mc=h['E_mc'], E_quad=h['E_quad'], E_miss=h['E_miss'], E_trunc=h['E_trunc'],
                        E_repr=h['E_repr'], E_repr_chk=Ep, E_repr_c=Ec, E_env3=s['E_env3'], E_env10=s['E_env10'], E_signed=s['E_signed'] + s['E_lo'],
                        E_lo=s['E_lo'], L3=L3, L3c=L3c, L10c=L10c, Ls=Ls, Lsc=Lsc, probe=pk, d=r[pk]['dc'], sd=r[pk]['sdc'],
                        budget3c=base_c + s['E_env3'], lower3c=s['lamc'] - base_c - s['E_env3'],
                        probes={k: {kk: float(vv) for kk, vv in v.items()} for k, v in r.items()})
        pr(f"  m/Lam={mfl:.7f}: main {h['main_key']} {h['main']:+.4f} +- {h['jk']:.4f} | lam_c {s['lamc']:+.4f} | signed {s['lam_signed']:+.4f} | "
           f"mc {h['E_mc']:.3f} quad {h['E_quad']:.3f} miss {h['E_miss']:.3f} trunc {h['E_trunc']:.3f} | repr on lam {h['E_repr']:.3f} "
           f"(recomputed {Ep:.3f}) -> on lam_c {Ec:.3f} [{pk} {r[pk]['dc']:+.4f} +- {r[pk]['sdc']:.4f}] | window x3 {s['E_env3']:.3f}, "
           f"x10 {s['E_env10']:.3f}, signed {s['E_signed'] + s['E_lo']:.3f}")
        pr(f"        L x3: protocol {L3:.2f} -> lam_c {L3c:.2f} | signed: {Ls:.2f} -> {Lsc:.2f} | x10: {L10c:.2f} | budget (x3, lam_c) "
           f"{base_c + s['E_env3']:.3f}, lower {s['lamc'] - base_c - s['E_env3']:+.4f}")
        for k, v in r.items():
            pr(f"        {k:6s}: lam {v['lam']:+.4f} lam_c {v['lam_c']:+.4f} | lam-K6 {v['d']:+.4f} +- {v['sd']:.4f} | lam_c-K6 {v['dc']:+.4f} +- {v['sdc']:.4f}")
    open(os.path.join(out, 'b4_odd_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(out, 'b4_odd_J%d.json' % J), 'w'), indent=1, default=float)
