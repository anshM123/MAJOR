"""Result Q, Stage 6 (RESULT_Q_STAGE6.md §2.3): heavy minima per gauge under NON-COVARIANT beyond-leading-Regge stress terms.

In every model of stages 3-6 (LR, HYB, EXA) the raw image is covariant (sum_j R_j(nu) F(nu - ij)); all gauges then give identical images and
identical minima (stage6_gauge covariant test, <= 1e-10).  Gauge dependence needs a non-covariant part.  Leading Regge is exactly covariant at
every order in beta (stages 4/5, [P]), so the admissible non-covariant part is beyond LR: O(1/m^2) at fixed beta, i.e. the size of the
finite-J constants (scalar skeleton: 5; graviton flat normal forms: c4 = 1 (++), 41 (+-)).  Two structures, strength K (per-harmonic image
F ~ F0 = 2):
   (a) shift-blind  dF_j = F0 K [nu^2 - (nu - ij)^2]/m^2      -> image error F0 K s_a/m^2,  s_a = mu_1 - M_1
   (b) j-weighted   dF_j = F0 K (j/jmax)^2 nu^2/m^2           -> image error F0 K s_b/m^2,  s_b = Re sum_j int c h_j (j/jmax)^2 nu^2
Heavy-action shifts: dB4 = 2K u^3 s^(4v)/Lam^2, d(dB4) = 2K u^3 s^(4d) (the derivative functional reads Lam^2 x the image), dB5 = s_h K u^4 s^(5)/Lam^2;
dFF = -c dB4 - a (d(dB4) - dB5).  Only the derivative piece matters at O(1); value and C5 pieces are 1/Lam^2-suppressed.
Per gauge (every valid gauge of stage6_gauge at dnu = 0.04 + the SB gauge of stage6_canon): coefficients, K_tol(+/-) = largest |K| keeping all
states >= 0, and min FF/N at K = +-5.  States: exact (Stage 3 + Stage 6) states, gap states J <= 60, region-A states m <= 30 Lam, J <= 60, bhat < 30.
Output: cache/stage6_stress.txt
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import common_q as CQ
import stage6_model as M
import stage4_hj as HJ4
import stage5_hj as HJ5

CACHE = M.CACHE
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def coefficients():
    G = pickle.load(open(os.path.join(CACHE, 'stage6_gauge.pkl'), 'rb'))
    C = pickle.load(open(os.path.join(CACHE, 'stage6_canon.pkl'), 'rb'))
    F4 = HJ4.numfuncs('EE3'); F5 = HJ5.numfuncs('C5_pp')

    def s4(nu, chi, deriv):
        d = nu[1] - nu[0]; R = {j: F4[j](nu.astype(complex)) for j in F4}; Nh = sum(R.values()).real
        mu1 = d*np.sum(chi*Nh*nu**2); M1 = d*np.sum(chi*sum(R[j]*(nu - 1j*j)**2 for j in R)).real
        sb = d*np.sum((chi*sum((j*j/9.)*R[j] for j in R)).real*nu**2)
        return mu1 - M1, sb

    def s5(nu, c):
        d = nu[1] - nu[0]; R = {j: F5[j](nu.astype(complex)) for j in F5}; Ns = sum(R.values()); h = {j: R[j]/Ns for j in R}
        mu1 = d*np.sum(c.real*nu**2); M1 = d*np.sum(sum(c*h[j]*(nu - 1j*j)**2 for j in h)).real
        sb = d*np.sum((c*sum((j*j/4.)*h[j] for j in h)).real*nu**2)
        return mu1 - M1, sb
    out4 = {}
    nu4 = G[('C4', 0.04, 'nu')]
    for k in G[('C4v', 0.04)]:
        ov, od = G[('C4v', 0.04)][k], G[('C4d', 0.04)][k]
        if ov['res'] < 1e-8 and od['res'] < 1e-8:
            out4[k] = (s4(nu4, ov['chi'], False), s4(nu4, od['chi'], True))
    out4['SB (shift-blind)'] = (s4(C['C4v']['nu'], C['C4v']['chi'], False), s4(C['C4d']['nu'], C['C4d']['chi'], True))
    out5 = {}
    nu5 = G[('C5', 0.04, 'nu')]
    for k, o in G[('C5', 0.04)].items():
        if o['res'] < 1e-8:
            out5[k] = s5(nu5, o['c'])
    out5['SB (shift-blind)'] = s5(C['C5']['nu'], C['C5']['c'])
    return out4, out5


def states(name, Dg, flat, Lam):
    rows = []
    for h in ('++', '+-'):
        for (u, Js, B2f, B3f) in flat[h]:
            sel = (Js <= M.jmax_gap(u, Dg)) & ((2*Js + 1)*np.sqrt(u) < 30) & (Js <= 60) & (u >= 1/30.**2)
            for J, b2, b3 in zip(Js[sel], B2f[sel], B3f[sel]):
                rows.append((h, u, int(J), b2, b3, None))
        for (u, J, b2, b3) in M.gap_states(name, Dg, 60)[h]:
            rows.append((h, u, J, b2, b3, None))
    F, c2, c3 = M.load(name)
    for (Delta, J), rec in M.load_exact(Dg).items():
        u = 1/rec['mfl']**2
        for h in ('++', '+-'):
            if (h == '++' and J % 2) or (h == '+-' and J < 4):
                continue
            Js, B2, B3 = CQ.parts(h, u, max(J, 4), [c2], [c3]); k = list(Js).index(J)
            rows.append((h, u, J, B2[k, 0], B3[k, 0], (rec[4]['v0'], Lam**2*rec[4]['d0'], rec[5]['v0'])))
    return rows


if __name__ == '__main__':
    t0 = time.time()
    c4, c5 = coefficients()
    pr("# stage6_stress: non-covariant beyond-LR stress terms of strength K (structures a: shift-blind, b: j-weighted), per gauge")
    pr("  coefficients s = (s_a, s_b):")
    for k, (v, d) in c4.items():
        pr(f"    C4 {k:18s}: value ({v[0]:+.4f}, {v[1]:+.4f})   derivative ({d[0]:+.4f}, {d[1]:+.4f})")
    for k, v in c5.items():
        pr(f"    C5 {k:18s}: ({v[0]:+.4f}, {v[1]:+.4f})")
    from stage6_scan import XS, JCAP, BMAX
    ref4, ref5 = 'H0 min-norm', 'G0 min-norm'
    combos = [(g4, ref5) for g4 in c4] + [(ref4, g5) for g5 in c5 if g5 != ref5] + [('SB (shift-blind)', 'SB (shift-blind)')]
    summary = {}
    for delta in (0.05, 0.2):
        name = M.FUNCS[delta]; F, c2, c3 = M.load(name); a, c = F['a'], F['c']
        flat = M.flat_table(name, XS, JCAP)
        for Dg in (100, 300, 1000):
            Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
            rows = states(name, Dg, flat, Lam)
            FF = np.zeros(len(rows)); NN = np.zeros(len(rows)); U = np.zeros(len(rows)); S = np.zeros(len(rows)); lab = []
            for i, (h, u, J, b2, b3, fw) in enumerate(rows):
                f, n, _, P = M.hyb(F, h, u, float(J), b2, b3, tab, Lam, fw=fw)
                FF[i], NN[i], U[i], S[i] = f, n, u, M.SIG[h]
                lab.append((h, J, 1/np.sqrt(u), float(P['bh'])))
            pr(f"\n== delta = {delta}, Dg = {Dg}: {len(rows)} (state, channel) cases; unstressed min FF/N = {np.min(FF/NN):+.4f} [{time.time() - t0:.0f}s]")
            for st in ('a', 'b'):
                ix = 0 if st == 'a' else 1
                res = []
                for (g4, g5) in combos:
                    sv, sd, s5_ = c4[g4][0][ix], c4[g4][1][ix], c5[g5][ix]
                    dK = -c*2*U**3*sv/Lam**2 - a*(2*U**3*sd - S*U**4*s5_/Lam**2)     # dFF per unit K
                    pos = dK > 0; neg = dK < 0
                    Kp = np.min(FF[neg]/(-dK[neg])) if neg.any() else np.inf      # K > 0 lowers FF where dK < 0
                    Km = np.min(FF[pos]/dK[pos]) if pos.any() else np.inf
                    r5p = (FF + 5*dK)/NN; r5m = (FF - 5*dK)/NN
                    ip, im = int(np.argmin(r5p)), int(np.argmin(r5m))
                    res.append((g4, g5, Kp, Km, r5p[ip], lab[ip], r5m[im], lab[im]))
                    pr(f"   ({st}) C4:{g4:18s} C5:{g5:18s}: K_tol +{Kp:8.2f} / -{Km:8.2f};  min FF/N at K=+5: {r5p[ip]:+.4f} ({lab[ip][0]} J={lab[ip][1]} m/Lam={lab[ip][2]:.3f});"
                       f"  at K=-5: {r5m[im]:+.4f} ({lab[im][0]} J={lab[im][1]} m/Lam={lab[im][2]:.3f})")
                v5 = [r[4] for r in res] + [r[6] for r in res]; kt = [min(r[2], r[3]) for r in res]
                pr(f"   ({st}) over gauges: min FF/N at |K| = 5 ranges {min(v5):+.4f} .. {max(min(r[4], r[6]) for r in res):+.4f};"
                   f" K_tol (both signs) ranges {min(kt):.2f} .. {max(kt):.2f}")
                summary[(delta, Dg, st)] = res
    pickle.dump(dict(c4=c4, c5=c5, summary=summary), open(os.path.join(CACHE, 'stage6_stress.pkl'), 'wb'))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, 'stage6_stress.txt'), 'w').write(buf.getvalue())
