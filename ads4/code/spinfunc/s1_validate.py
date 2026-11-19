"""
s1_validate (OFFDIAG_KERNEL.md §1): validation of the exact leading-Regge off-diagonal kernel (s1_offdiag_kernel) against the CFT crystal
chain (q3/l9 K2 rows: Re AND Im of Gamma_off), the direct (RQMC) ladder points of r6 (LARGEJ §2.3), and r4's universal g.
Also: (a) which helicity configurations are non-zero; (b) scaling limit nu -> oo at fixed t = nu beta: R(t) = Re Y/beta^2,
F(t) = Im Y/(beta t), the ratio R/F, and the moments int t R dt, int t F dt.
Model fitted to every CFT row (sign of the row fixed by the Im overlap; c, s per row):
    Re Gamma_off^CFT/m = Re Y + (c/nu) Im Y + s Re D,    D = forward-normalised PRESERVING image at the state (= P_l + O(1/nu) shift terms)
(s = a real change of OPE basis ++ -> ++ + s(+-)-type: a per-state shear; c = a nu-only phase).
Usage: python3 s1_validate.py  -> cache/s1_validate.txt, cache/s1_cftrows.pkl, cache/s1_scaling.pkl
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import j0
import s1_offdiag_kernel as S

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


def cft_rows():
    fn = os.path.join(CACHE, 's1_cftrows.pkl')
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    import r2_ladder                                         # noqa (r1 caches merged into q3's loader)
    import q3_analysis as Q
    out = {}
    for J in (6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32):
        g = Q.even_grid(J)
        for D, r in sorted(g.items()):
            if D < 300:
                continue
            out[(J, D)] = dict(m=r['m'], xs=r['xs'], re2=r[('C2', 'off')], im2=r[('C2', 'off', 'im')], re3=r[('C3H', 'off')],
                               im3=r[('C3H', 'off', 'im')], a=r['a'])
    pickle.dump(out, open(fn, 'wb'))
    return out


def kernel_at(nu, be):
    v = S.img_vec(nu, np.array([0.0, be]), S.PRES)
    o = S.img_vec(nu, np.array([be]), S.FLIP24)[0]
    return o/v[0], v[1]/v[0]                                  # Y (off), D (preserving, forward-normalised)


def main():
    t0 = time.time()
    pr("# s1_validate: exact leading-Regge off-diagonal vs the CFT chain")
    # (a) configurations
    pr("\n(a) helicity configurations (a1, a3; a2, a4) [z1, z3 at pbh; z2', z4' at ph], |L|/|L_pres| at (nu, beta) = (20, 0.3), (80, 0.05):")
    for nu, be in ((20.0, 0.3), (80.0, 0.05)):
        pres2 = S.img2(nu, be, S.PRES); pres3 = S.img3(nu, be, S.PRES)
        for lab, pol in (('flip 24 (a2 = a4)', (1, -1, -1, -1)), ('flip 13 (a1 = a3)', (1, 1, -1, 1)), ('double flip', (1, 1, -1, -1))):
            pr(f"   nu={nu}, beta={be}: {lab:18s}  C2: {abs(S.img2(nu, be, pol))/abs(pres2):.2e}   C3: {abs(S.img3(nu, be, pol))/abs(pres3):.2e}")
    # (b) scaling limit
    pr("\n(b) scaling limit (nu = 500, 2000; t = nu beta): R = Re Y/beta^2, F = Im Y/(beta t)")
    sc = {}
    for nu in (500.0, 2000.0):
        ts = np.arange(0.2, 80.01, 0.2)
        Y = np.array([S.img_vec(nu, np.array([t/nu]), S.FLIP24, nphi=256 if t < 40 else 512)[0] for t in ts])
        n0 = S.img_vec(nu, np.array([0.0]), S.PRES)[0]
        Y = Y/n0; be = ts/nu
        R, F = Y.real/be**2, Y.imag/(be*ts)
        sc[nu] = (ts, R, F)
        cR = cumulative_trapezoid(ts*R, ts, initial=0); cF = cumulative_trapezoid(ts*F, ts, initial=0)
        k = ts > 40
        pr(f"   nu={nu:6.0f}: R/F in [{np.min(R/F):.4f}, {np.max(R/F):.4f}];  int_0^T t R dt (T>40: mean {cR[k].mean():.4f} +- {cR[k].std():.4f});"
           f"  int t F dt: {cF[k].mean():.5f} +- {cF[k].std():.5f}")
        pr("      t : " + " ".join(f"{t:6.1f}" for t in (0.5, 1, 1.5, 2, 2.5, 3, 3.7, 4.5, 5, 6, 7, 8, 10, 12)))
        pr("      R : " + " ".join(f"{R[np.argmin(abs(ts - t))]:+6.3f}" for t in (0.5, 1, 1.5, 2, 2.5, 3, 3.7, 4.5, 5, 6, 7, 8, 10, 12)))
    pickle.dump(sc, open(os.path.join(CACHE, 's1_scaling.pkl'), 'wb'))
    # (c) CFT rows
    pr("\n(c) CFT K2 rows (q3 loader; r4's 25 rows).  Im: median and rms of Im_CFT/Im_ker on t <= 9 (sign of the row fixed by the Im overlap).")
    pr("    Re fit: Re_CFT = Re Y + (c/nu) Im Y + s Re D (t <= 9).  rms residual and rms Re_CFT, in units of beta^2.")
    d = cft_rows()
    fits = {}
    for key in sorted(d, key=lambda k: (k[1], k[0])):
        J, D = key; r = d[key]; m = r['m']; be = np.arccosh(1 + (2*J + 1)**2/(2*m*m))
        sel = r['xs']*m*be <= 9.0
        xs = r['xs'][sel]; cre = r['re2'][sel]/m; cim = r['im2'][sel]/m; nus = xs*m
        KY, KD = np.array([kernel_at(nu, be) for nu in nus]).T
        sg = np.sign(np.sum(cim*KY.imag)); cre, cim = sg*cre, sg*cim
        A = np.array([KY.imag/nus, KD.real]).T
        sol = np.linalg.lstsq(A, cre - KY.real, rcond=None)[0]
        res = cre - KY.real - A @ sol
        A1 = np.array([KD.real]).T                              # shear only (c = 0)
        sol1 = np.linalg.lstsq(A1, cre - KY.real, rcond=None)[0]; res1 = cre - KY.real - A1 @ sol1
        rat = cim/KY.imag
        fits[key] = dict(beta=be, c=sol[0], s=sol[1], s1=sol1[0], sg=sg)
        pr(f"   (J={J:2d}, Delta={D:5.0f}) beta={be:.4f}: Im ratio median {np.median(rat):.3f}, rms|Im_CFT - Im_ker|/rms Im_ker "
           f"{np.sqrt(np.mean((cim - KY.imag)**2)/np.mean(KY.imag**2)):.3f} | c = {sol[0]:+.2f}, s/beta^2 = {sol[1]/be**2:+.3f}, resid "
           f"{np.sqrt(np.mean(res**2))/be**2:.3f}; shear only: s/beta^2 = {sol1[0]/be**2:+.3f}, resid {np.sqrt(np.mean(res1**2))/be**2:.3f}"
           f" (rms Re_CFT {np.sqrt(np.mean(cre**2))/be**2:.3f}, rms Re Y {np.sqrt(np.mean(KY.real**2))/be**2:.3f})")
    sv = np.array([f['s1']/f['beta']**2 for k, f in fits.items() if k[1] >= 800])
    pr(f"   shear-only s/beta^2 over the Delta >= 800 rows: mean {sv.mean():+.3f}, sd {sv.std():.3f}, range [{sv.min():+.3f}, {sv.max():+.3f}]")
    # C3H: the kernel off-diagonal is identically 0; is the CFT one a basis effect (s3 x sigma_h P_l) or a 1/m effect?
    pr("\n    C3H rows: the leading-Regge C^_3 off-diagonal vanishes identically (a).  CFT Re Gamma_3,off/m vs m at fixed (nu, beta):")
    for (J1, D1), (J2, D2) in (((8, 400), (16, 800)), ((16, 800), (32, 1600)), ((6, 400), (12, 800)), ((12, 800), (24, 1600))):
        r1, r2 = d[(J1, D1)], d[(J2, D2)]
        t1 = r1['xs']*r1['m']*np.arccosh(1 + (2*J1 + 1)**2/(2*r1['m']**2)); t2 = r2['xs']*r2['m']*np.arccosh(1 + (2*J2 + 1)**2/(2*r2['m']**2))
        ok = (t1 > 1.5) & (t1 < 8)
        v1 = np.interp(t1[ok], t1, r1['re3']); v2 = np.interp(t1[ok], t2, r2['re3'])
        w1 = np.interp(t1[ok], t1, r1['re2']); w2 = np.interp(t1[ok], t2, r2['re2'])
        pr(f"     ({J1},{D1}) -> ({J2},{D2}) (m x {r2['m']/r1['m']:.2f}): rms Gamma_3,off ratio {np.sqrt(np.mean(v2**2)/np.mean(v1**2)):.2f}"
           f" (C2: {np.sqrt(np.mean(w2**2)/np.mean(w1**2)):.2f}; leading Regge => x {r2['m']/r1['m']:.2f}, 1/m-suppressed => x 1)")
    # (d) direct ladder points (r6 table, LARGEJ §2.3), C2, Delta_g = 300
    pr("\n(d) direct (RQMC) ladder points (LARGEJ §2.3 / cache/r6_ladder_direct.txt), C2 Gamma_off/(m beta^2):")
    direct = {(18, 336.0): [-0.75, -1.24, -1.06, -0.48, 0.03, 0.20], (36, 672.0): [-0.65, -1.18, -1.12, -0.62, -0.05, 0.27]}
    nusd = np.array([25.0, 33.4, 41.7, 50.1, 58.4, 66.8])
    for (J, D), vals in direct.items():
        al = D - J - 2.0; m = np.sqrt(al*(al + 2*J + 1)); be = np.arccosh(1 + (2*J + 1)**2/(2*m*m))
        KY, KD = np.array([kernel_at(nu, be) for nu in nusd]).T
        cre = np.array(vals)*be**2
        best = None
        for sg in (1, -1):
            A1 = np.array([KD.real]).T; s1 = np.linalg.lstsq(A1, sg*cre - KY.real, rcond=None)[0]
            rr = sg*cre - KY.real - A1 @ s1
            if best is None or np.sum(rr**2) < best[0]:
                best = (np.sum(rr**2), sg, s1[0], rr)
        _, sg, s1, rr = best
        pr(f"   ({J},{D:.0f}) beta={be:.4f}: direct {np.round(sg*np.array(vals), 2).tolist()}; kernel Re Y/b^2 {np.round(KY.real/be**2, 2).tolist()};"
           f" kernel + shear (s/beta^2 = {s1/be**2:+.3f}) {np.round((KY.real + s1*KD.real)/be**2, 2).tolist()}, rms resid {np.sqrt(np.mean(rr**2))/be**2:.3f}")
    # (e) r4's universal g vs  -(R - sigma J0)  (scaling limit, r4 sign convention)
    g4 = pickle.load(open(os.path.join(CACHE, 'r4_g.pkl'), 'rb'))
    ts, R, F = sc[2000.0]
    Rt = np.interp(g4['tg'], ts, R)
    sig = np.linalg.lstsq(np.array([j0(g4['tg'])]).T, g4['gmed'] + Rt, rcond=None)[0][0]
    pr(f"\n(e) r4 median g(t) vs -R(t) + sigma J0(t) (one constant, t = 1 ... 12): sigma = {sig:.3f}")
    pr("      t     : " + " ".join(f"{t:6.1f}" for t in g4['tg'][::2]))
    pr("      g_r4  : " + " ".join(f"{v:+6.2f}" for v in g4['gmed'][::2]))
    pr("      -R+sJ0: " + " ".join(f"{v:+6.2f}" for v in (-Rt + sig*j0(g4['tg']))[::2]))
    pr("      -R    : " + " ".join(f"{v:+6.2f}" for v in (-Rt)[::2]))
    pr(f"      rms(g_r4 - model) = {np.sqrt(np.mean((g4['gmed'] + Rt - sig*j0(g4['tg']))**2)):.3f}; r4 spread (median) {np.median(g4['gsd']):.3f}")
    pickle.dump(fits, open(os.path.join(CACHE, 's1_validate_fits.pkl'), 'wb'))
    pr(f"\n[{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, 's1_validate.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
