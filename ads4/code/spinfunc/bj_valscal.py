"""
bj_valscal (BELOW_FOLD_SPINJ.md §2.2): the J/Delta-resummed P_J series of bj_ee (gam, Dhat from the EE vectors V = V_{5,12}, W = V_{5,34},
pref_sum with triangular truncation) applied to the SCALAR integrand (pref := 1), against bg_scalarJ (N = 12, R = 4; itself validated
against the exact Mellin skeleton to <= 0.04 c2 on nu <= 10, BELOW_FOLD2 §2.2).  Checks t(V, W) = spin_t (sign aside) and the P^(k)
machinery (k = 0 here) at J = 24-32 and large J.
Usage: python3 bj_valscal.py  -> cache/bj_valscal.txt
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_ser as S, bg_geom as GG, bg_ee as BE, bj_ee as BJ, bg_scalarJ as SJ

C = os.path.join(HERE, 'cache')


def scalar_images(Dl, J, Kr, N, ntri, nus):
    depth = Kr
    G2 = GG.geometry(Kr + 2)
    V, W = BJ.vw(BE.EEGeom(Kr + 2, G2))
    t = BJ.tser(V, W)
    if np.real(complex(t.lead('t')[0])) < 0:
        t = BJ.tser(V, [w.scale(-1.0) for w in W])
    dl = (t - 1.0).norm()
    Dh = S.ASer(dl.d - 2, dl.c)
    one = S.ASer(0, [S.Poly.const(1.0)] + [S.Poly() for _ in range(Kr)])
    P = dict(X=[[[one]]] + [[[None]] for _ in range(8)], Dh=Dh, J=J, Delta=Dl)
    pr = BJ.pref_sum(P, N, Kr, ntri)[0][0]
    logs, polys, Gs, _ = GG.scalar_pieces(Kr + 3)
    A = GG.assemble(logs, polys + [('PJ', pr)], Kr)
    Bf = GG.bfac(Kr)
    tot = np.zeros(len(nus), complex)
    for j, cj in A['corr'].items():
        Pk = S.lmul(cj, Bf.c, Kr)
        for k in range(Kr + 1):
            if len(Pk[k]):
                tot += Dl**(-k)*(0.5j*np.asarray(nus))**j*S.integrate_fast(Pk[k].E, Pk[k].C, A['E0'], A['E1'], nus)
    # compare t with spin_t
    t0s, dls = GG.spin_t(Gs, Kr + 2)
    pts = S.atom_points(5, 4)
    dmax = max(np.max(np.abs(dls.c[k].ev(pts) - (t0s.real/abs(t0s.real))*dl.c[k - dl.d].ev(pts) if k - dl.d < len(dl.c) and k >= dl.d else
                             dls.c[k].ev(pts))) for k in range(2, Kr + 3))
    return tot*np.exp(-pr.d*np.log(Dl)), dmax


if __name__ == '__main__':
    lines = []
    def pr_(s=''):
        print(s, flush=True); lines.append(s)
    r = pickle.load(open(os.path.join(C, 'bg_scalarJ_N12_R4.pkl'), 'rb'))
    nus = r['nus']
    pr_('# bj_valscal: m^2 [log(I_bj(Delta)/I_bj(ref)) - log(ratio_bgJ)] - const(nu<1.4); bg_scalarJ N=12 R=4 is the trusted resummation '
        '(<= 0.04 c2 vs the exact skeleton on nu <= 10)')
    for (Dl, J) in ((336.0, 24), (336.0, 28), (336.0, 32), (364.0, 64), (400.0, 100), (1150.0, 150)):
        for (Kr, N, ntri) in ((4, 12, 100), (3, 12, 100), (2, 12, 100), (4, 8, 100)):
            t0 = time.time()
            st, dm = scalar_images(Dl, J, Kr, N, ntri, nus)
            rf, _ = scalar_images(2e4, J, Kr, N, ntri, nus)
            bgJ = SJ.ratio(r, Dl, J)
            m = Dl
            d = m*m*(np.log(st/rf) - np.log(bgJ)).real
            d -= np.mean(d[nus < 1.4])
            pr_(f'  ({Dl:g},{J}) Kr={Kr} N={N} ntri={ntri}: max|.| nu<=4 {np.max(np.abs(d[nus <= 4])):.3f}, 4-10 {np.max(np.abs(d[(nus > 4) & (nus <= 10)])):.3f}, '
                f'10-18 {np.max(np.abs(d[nus > 10])):.3f}  (dDhat vs spin_t {dm:.1e}; {time.time() - t0:.0f}s)')
    open(os.path.join(C, 'bj_valscal.txt'), 'w').write('\n'.join(lines) + '\n')
