"""
s3_structures (OFFDIAG_KERNEL.md §3): leading-Regge diagonal and helicity off-diagonal images of the other k = 2 (j' = 0) (2,4) dual
structures (s1_offdiag_kernel.STRUCTS: the symmetric span P^2, PH, HH (q^2 = 4PP - 4PH + HH), the bulk-point-null flip-free
Nq = Hh q = 2PH - HH, and 2<->4-asymmetric ones), and the search for a combination with vanishing off-diagonal first moment.
Modes:
  scan   : scaling limit (nu = 400, 2000; t = nu beta in (0, 60]) per structure: o_i(nu, t/nu) (FLIP24), n_i(nu) (forward PRES), line-13 flip
           check; forward-normalisation ratios n_i/n_q2 vs nu (bulk-point nullity of Nq in the Regge overlap); h_-+ per structure (3-term
           fit); kappa-scan of S = q^2 + kappa Nq: G1(kappa) = int t Re Y_kappa/beta^2 dt.   -> cache/s3_scan.txt, cache/s3_scan.pkl
  table  : exact (nu, beta) tables of o_Nq/n_q2 and n_Nq/n_q2 for the LP (same nodes as s1)   -> cache/s1_table_<tag>.pkl (dict with 'oN','nN')
Usage: python3 s3_structures.py scan | table <tag> <beta indices>
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np
from scipy.integrate import cumulative_trapezoid
import s1_offdiag_kernel as S
import v3_regge_kernel as V

CACHE = os.path.join(HERE, 'cache')
NAMES = ['q2', 'PP', 'PH', 'HH', 'Nq', 'A1', 'A1m', 'A2', 'A3', 'A3m', 'VV']


def img(name, nu, betas, pols):
    return S.img_vec(nu, np.asarray(betas, float), pols, Amonos=S.A_monos_struct(name), nphi=256 if np.max(betas)*nu < 40 else 512)


def scan():
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    out = {}
    pr("# s3 scan: k = 2 (2,4) structures at leading Regge order")
    # line-13 flips and h0 per structure (nu = 6, 20; beta = 0.4)
    pr("\n(a) line-13 flips |L(1,1;-1,1)|/|L_pres| and 3-term structure of the preserving image (fit P_{l-1}, P_l, P_{l+1}):")
    for nm in NAMES:
        for nu in (6.0, 20.0):
            bs = np.linspace(0.12, min(3.0, max(0.6, 30/nu)), 14)
            vals = img(nm, nu, bs, S.PRES)
            (hm, h0, hp), res = V.fit3(nu, bs, vals)
            f13 = abs(img(nm, nu, [0.4], (1, 1, -1, 1))[0])/abs(img(nm, nu, [0.4], S.PRES)[0])
            out.setdefault('h', {})[(nm, nu)] = (hm, h0, hp)
            pr(f"   {nm:4s} nu={nu:4.0f}: flip13 {f13:.1e};  fit resid {res:.1e}, |h0|/(|h-|+|h+|) {abs(h0)/(abs(hm) + abs(hp)):.1e}, "
               f"h+/h- {hp/hm:.4f}")
    # forward normalisations relative to q2 vs nu
    pr("\n(b) forward (beta = 0) normalisation n_i(nu)/n_q2(nu)  [bulk-point null structures must -> 0 at large nu]:")
    nus = np.array([1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0])
    nf = {nm: np.array([img(nm, nu, [0.0], S.PRES)[0] for nu in nus]) for nm in NAMES}
    out['nf'] = (nus, nf)
    for nm in NAMES:
        r = nf[nm]/nf['q2']
        pr(f"   {nm:4s}: " + " ".join(f"{abs(v):.2e}" for v in r) + f"   (nu |r| at 3000: {3000*abs(r[-1]):.3f}, nu^2 |r|: {9e6*abs(r[-1]):.1f})")
    # scaling-limit off-diagonals
    pr("\n(c) scaling limit: o_i(nu, t/nu)/n_q2(nu) / beta^2 = R_i(t) + i nu F_i(t)-type; moments int t Re(o_i/n_q2)/beta^2 dt  (t <= 60)")
    ts = np.arange(0.2, 60.01, 0.2)
    sc = {}
    for nu in (400.0, 2000.0):
        nq = img('q2', nu, [0.0], S.PRES)[0]
        for nm in NAMES:
            o = np.array([img(nm, nu, [t/nu], S.FLIP24)[0] for t in ts])/nq
            be = ts/nu
            sc[(nm, nu)] = o/be**2
        for nm in NAMES:
            Z = sc[(nm, nu)]
            cR = cumulative_trapezoid(ts*Z.real, ts, initial=0); cI = cumulative_trapezoid(ts*Z.imag/nu, ts, initial=0)
            k = ts > 30
            pr(f"   nu={nu:5.0f} {nm:4s}: int t Re dt = {cR[k].mean():+9.4f} (+-{cR[k].std():.3f});  int t Im/nu dt = {cI[k].mean():+9.5f};"
               f"  Re at t=1,2,3.7,6: " + " ".join(f"{Z.real[np.argmin(abs(ts - t))]:+.3f}" for t in (1, 2, 3.7, 6)))
        pr(f"   [{time.time() - t0:.0f}s]")
    out['sc'] = (ts, sc)
    # kappa scan for S = q2 + kappa Nq (flip-free at the bulk point, same bulk-point density as q2)
    pr("\n(d) S(kappa) = q^2 + kappa Nq:  G1(kappa) := int t Re Y_kappa/beta^2 dt (scaling, nu = 2000; Y_kappa = (o_q + k o_N)/(n_q + k n_N))")
    nu = 2000.0
    nq = img('q2', nu, [0.0], S.PRES)[0]; nN = img('Nq', nu, [0.0], S.PRES)[0]
    for kap in (-2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0):
        Z = (sc[('q2', nu)] + kap*sc[('Nq', nu)])/(1 + kap*nN/nq)
        c = cumulative_trapezoid(ts*Z.real, ts, initial=0)
        pr(f"   kappa={kap:+5.2f}: G1 = {c[ts > 30].mean():+8.3f}")
    pickle.dump(out, open(os.path.join(CACHE, 's3_scan.pkl'), 'wb'))
    open(os.path.join(CACHE, 's3_scan.txt'), 'w').write(buf.getvalue())


def table(tag, idx, names=('Nq',)):
    out = {}
    t0 = time.time()
    for i in idx:
        be = S.BETAS[i]
        nus = S.nu_grid(be)
        rec = {nm: (np.zeros(len(nus), complex), np.zeros(len(nus), complex)) for nm in names}
        for k, nu in enumerate(nus):
            nq = S.img_vec(nu, np.array([0.0]), S.PRES, nphi=192)[0]
            for nm in names:
                Am = S.A_monos_struct(nm)
                rec[nm][0][k] = S.img_vec(nu, np.array([be]), S.FLIP24, nphi=192, Amonos=Am)[0]/nq
                rec[nm][1][k] = S.img_vec(nu, np.array([0.0]), S.PRES, nphi=192, Amonos=Am)[0]/nq
        out[be] = (nus, rec)
        print(f"  [{tag}] beta={be}: {len(nus)} nodes [{time.time() - t0:.0f}s]", flush=True)
        pickle.dump(out, open(os.path.join(CACHE, f's3_table_{tag}.pkl'), 'wb'))


if __name__ == '__main__':
    if sys.argv[1] == 'scan':
        scan()
    elif sys.argv[1] == 'table':
        table(sys.argv[2], [int(v) for v in sys.argv[3].split(',')])


def table2(tag, idx, name='Nq'):
    """o_name(nu, beta)/n_name(nu) on s1's nodes (2 kernel calls per node); combine with ratio_table() for n_name/n_q2."""
    out = {}
    t0 = time.time(); Am = S.A_monos_struct(name)
    for i in idx:
        be = S.BETAS[i]; nus = S.nu_grid(be)
        Z = np.zeros(len(nus), complex)
        for k, nu in enumerate(nus):
            Z[k] = S.img_vec(nu, np.array([be]), S.FLIP24, nphi=192, Amonos=Am)[0]/S.img_vec(nu, np.array([0.0]), S.PRES, nphi=192, Amonos=Am)[0]
        out[be] = (nus, Z)
        print(f"  [{tag}] beta={be}: {len(nus)} nodes [{time.time() - t0:.0f}s]", flush=True)
        pickle.dump(out, open(os.path.join(CACHE, f's3_table_{tag}.pkl'), 'wb'))


def ratio_table(name='Nq'):
    nus = np.concatenate([np.arange(0.005, 10, 0.005), np.arange(10, 100, 0.05), np.geomspace(100, 3200, 800)])
    Am = S.A_monos_struct(name)
    r = np.array([S.img_vec(nu, np.array([0.0]), S.PRES, nphi=128, Amonos=Am)[0]/S.img_vec(nu, np.array([0.0]), S.PRES, nphi=128)[0]
                  for nu in nus])
    pickle.dump(dict(nu=nus, r=r), open(os.path.join(CACHE, f's3_ratio_{name}.pkl'), 'wb'))
    print(f"ratio n_{name}/n_q2: max |r| = {np.abs(r).max():.4f} at nu = {nus[np.argmax(np.abs(r))]:.3f}")


if __name__ == '__main__' and sys.argv[1] == 'table2':
    table2(sys.argv[2], [int(v) for v in sys.argv[3].split(',')])
if __name__ == '__main__' and sys.argv[1] == 'ratio':
    ratio_table()
