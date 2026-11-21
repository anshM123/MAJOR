"""Result Q, Stage 7 (RESULT_Q_STAGE7.md §1-§2): what a measurement of K has to measure.

The de-shifted forward pieces are F_chi[F] := Re sum_j int chi(nu) R_j(nu) F_j(nu) dnu (C4: 4 harmonics j = +-1, +-3, real chi).
Every valid gauge solves the same Abel equation W[chi] = Phi; two gauges differ by Dchi with W[Dchi] = 0 identically.  This script
 (1) computes the null space of the discrete Abel operator (Stage 4/6: 2 exact + 1 near-null directions) on two grids, its nu-profile,
     and its relation to the real zeros of Nhat4;
 (2) decomposes every Stage-6 gauge difference (H1..H4, SB) on that null space (residual outside the span);
 (3) defines the response of the heavy action to an ARBITRARY beyond-LR correction of the raw image, in c2-like units
        F_j(nu) = 2 [P(nu - ij) + Gamma_j(nu)/m^2]      (per-harmonic forward image; F0 = 2 as stage6_stress)
        d(dB4) = 2 u^3 Re sum_j int chi_d R_j Gamma_j,   dB4 = 2 u^3/Lam^2 Re sum_j int chi_v R_j Gamma_j
     and its amplification densities by nu-window (where would Gamma have to be known, and how well);
 (4) per binding state (++ J = 0, 2, 4; +- J = 4..8; m from the gap to 1.5 Lam; Dg = 100, 300, 1000; delta = 0.05, 0.2; HYB model of
     stage6_model) the tolerance on Gamma: uniform |Gamma| <= G (worst sign pattern), Gamma = g nu^2-type growth, and the gauge spread.
Output: cache/stage7_null.txt, cache/stage7_null.pkl
"""
import os, sys, io, time, pickle, json
os.environ.setdefault('OMP_NUM_THREADS', '2'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import stage4_deshift as DS4
import stage6_model as M
import common_q as CQ

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()
ZEROS4 = np.array([1.068482, 4.545932, 8.169433, 13.094015])     # real zeros of Nhat4 (Stage 4 §1.4)
WIN = [(0, 2.5), (2.5, 8.3), (8.3, 18), (18, 1e9)]             # nu windows: below/at thimble range/n3 range/above null support


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def quant(nu, w, qs=(0.1, 0.5, 0.9)):
    c = np.cumsum(w)/np.sum(w)
    return [float(nu[min(np.searchsorted(c, q), len(nu) - 1)]) for q in qs]


def windows(nu, dens):
    tot = np.sum(dens)
    return [float(np.sum(dens[(nu >= a) & (nu < b)])/tot) for a, b in WIN]


def null_space(S, k=3):
    U, s, Vh = S.svd()
    return s/s[0], Vh[-k:][::-1]          # ordered: least singular first


def main():
    t0 = time.time()
    G = pickle.load(open(os.path.join(CACHE, 'stage6_gauge.pkl'), 'rb'))
    C = pickle.load(open(os.path.join(CACHE, 'stage6_canon.pkl'), 'rb'))
    out = dict()
    pr("# stage7_null: null space of the C4 Abel operator, gauge decomposition, response to beyond-LR corrections")
    for dnu, N in ((0.04, 1024), (0.02, 2048)):
        S = DS4.Deshift4('EE3', dnu, N)
        sv, V = null_space(S, 3)
        pr(f"\n## grid dnu = {dnu}, N = {N} (nu_max = {N*dnu:.0f}); smallest s.v./s_max: " + " ".join(f"{x:.1e}" for x in sv[::-1][:5]))
        nu = S.nu
        for i, v in enumerate(V):
            chi = v/S.n
            dY = v**2; dN = np.abs(chi*S.Nhat); dR = np.abs(chi)*S.n
            # W of the null vector (should vanish): direct quadrature at a few t
            Wt = S.W(chi, np.array([0.0, 1.0, 3.0, 6.0]))
            pr(f"  null-{i} (s.v. {sv[::-1][i] if False else sorted(sv)[i]:.1e}): |Y|^2 quantiles nu = {quant(nu, dY)}; |chi Nhat| quantiles "
               f"{quant(nu, dN)}; |chi| sum|R_j| windows {WIN[0][1]}/{WIN[1][1]}/{WIN[2][1]}/inf: "
               + " ".join(f"{x:.3f}" for x in windows(nu, dR)) + f"; max|W(t<=6)|/|Y|_1 = {np.max(np.abs(Wt))/np.sum(np.abs(v))/dnu:.1e}")
            # local structure at the real zeros of Nhat4
            pk = [float(nu[j]) for j in np.argsort(-np.abs(v))[:1]]
            pr(f"           peak |Y| at nu = {pk[0]:.2f}; nearest Nhat4 zero {ZEROS4[np.argmin(np.abs(ZEROS4 - pk[0]))]:.3f}")
        out[('null', dnu)] = dict(sv=sv, V=V, nu=nu, n=S.n)
        # gauge differences on the null span
        for key, lab in (('C4v', 'value'), ('C4d', 'deriv')):
            if (key, dnu) not in G:
                continue
            fam = G[(key, dnu)]
            Y0 = fam['H0 min-norm']['Y']
            cand = {k: o['Y'] for k, o in fam.items() if k != 'H0 min-norm'}
            sb = C[key]['chi'] if dnu == 0.02 else C[key]['chi_04']
            cand['SB'] = sb*S.n
            pr(f"  {lab}: gauge differences Y_g - Y_H0 projected on the 3 (2 exact + near-null) null directions")
            for k, Yg in cand.items():
                d = Yg - Y0
                cf = V @ d; res = np.linalg.norm(d - V.T @ cf)/max(np.linalg.norm(d), 1e-300)
                cf2 = V[:2] @ d; res2 = np.linalg.norm(d - V[:2].T @ cf2)/max(np.linalg.norm(d), 1e-300)
                pr(f"     {k:16s}: |dY|_2 = {np.linalg.norm(d):.3e}; coefficients on null-0,1,2 = " + " ".join(f"{x:+.3e}" for x in cf) +
                   f"; residual outside span(3) {res:.1e}, outside span(2 exact) {res2:.1e}")
                out[('gdiff', key, dnu, k)] = dict(cf=cf, res=res, res2=res2)
    # ---------------------------------------------------------------- (3) response densities, dnu = 0.02 family
    dnu, N = 0.02, 2048
    S = DS4.Deshift4('EE3', dnu, N); nu = S.nu
    fams = {}
    for key in ('C4v', 'C4d'):
        fam = {k: o['chi'] for k, o in G[(key, dnu)].items() if k != 'H1 tol 1e-6'}     # H1 is not grid-stable (Stage 6)
        fam['SB'] = C[key]['chi']
        fams[key] = fam
    stored = json.load(open(os.path.join(HERE, 'functionals_q4_deshift.json')))['EE3']
    pr("\n## response of the de-shifted pieces to a beyond-LR per-harmonic correction Gamma_j(nu) (c2 units: m^2 x relative)")
    pr("   amplification A[w] = int |chi| sum_j |R_j| w(nu) dnu; windows nu < 2.5 | 2.5-8.3 | 8.3-18 | > 18 (fractions of A[1])")
    resp = {}
    for key in ('C4v', 'C4d'):
        fam = fams[key]
        chi_st = np.array(stored['chi_value' if key == 'C4v' else 'chi_deriv'])
        pr(f"  {key}: stored Stage-4 weight vs H0 (dnu = 0.02): max|chi_st - chi_H0|/max|chi_H0| = "
           f"{np.max(np.abs(chi_st - fam['H0 min-norm']))/np.max(np.abs(fam['H0 min-norm'])):.1e}")
        for k, chi in fam.items():
            ad = np.abs(chi)*S.n
            A0 = dnu*np.sum(ad); A2 = dnu*np.sum(np.abs(chi)*sum(np.abs(S.Rj[j])*np.abs(nu - 1j*j)**2 for j in S.Rj))
            spread = dnu*np.sum(np.abs(chi - fam['H0 min-norm'])*S.n)
            spread2 = dnu*np.sum(np.abs(chi - fam['H0 min-norm'])*sum(np.abs(S.Rj[j])*np.abs(nu - 1j*j)**2 for j in S.Rj))
            resp[(key, k)] = dict(A0=A0, A2=A2, spread=spread, spread2=spread2, win=windows(nu, ad),
                                  winsp=windows(nu, np.abs(chi - fam['H0 min-norm'])*S.n + 1e-300))
            pr(f"     {k:16s}: A[1] = {A0:.4f} (windows " + " ".join(f"{x:.3f}" for x in resp[(key, k)]['win']) +
               f"), A[|nu-ij|^2] = {A2:.3f}; vs H0: int|dchi| sum|R| = {spread:.4f} (windows " +
               " ".join(f"{x:.3f}" for x in resp[(key, k)]['winsp']) + f"), x|nu-ij|^2: {spread2:.3f}")
    out['resp'] = resp
    # ---------------------------------------------------------------- (4) per binding state tolerances
    from stage6_scan import BMAX
    pr("\n## tolerances at the binding near-gap states (HYB model, Stage 6).  d(dB4) = 2u^3 Re sum_j int chi_d R_j Gamma_j enters FF as -a d(dB4);")
    pr("   G_U = largest uniform |Gamma| (worst sign pattern) keeping FF >= 0 in gauge g = FF/(2 a u^3 A_d[1] + 2 c u^3 A_v[1]/Lam^2);")
    pr("   G_Q = same for |Gamma_j(nu)| <= g |nu - ij|^2 (natural O(nu^2/m^2) growth): FF/(2 a u^3 A_d[|nu-ij|^2] + ...);")
    pr("   G_S = uniform |Gamma| that can move the state by FF through the GAUGE SPREAD alone (vs H0): FF/(2 a u^3 spread_d + ...)")
    states_out = []
    gk = [k for k in fams['C4d']]
    for delta in (0.05, 0.2):
        name = M.FUNCS[delta]; F, c2, c3 = M.load(name); a, c = F['a'], F['c']
        for Dg in (100, 300, 1000):
            Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
            rows = []
            for h, Js in (('++', (0, 2, 4)), ('+-', (4, 5, 6, 7, 8))):
                for J in Js:
                    ug = M.u_gap(J, Dg); mg = 1/np.sqrt(ug)
                    ms = [mg] + [x for x in np.arange(1.0, 1.5001, 0.05) if x > mg + 1e-9]
                    for mr in ms:
                        u = 1/mr**2
                        Js_, B2, B3 = CQ.parts(h, u, max(J, 4), [c2], [c3]); kk = list(Js_).index(J)
                        f, n, _, P = M.hyb(F, h, u, float(J), B2[kk, 0], B3[kk, 0], tab, Lam)
                        rec = dict(delta=delta, Dg=Dg, h=h, J=J, m=mr, u=u, FF=float(f), N=float(n))
                        for g in gk:
                            dd, vv = resp[('C4d', g)], resp[('C4v', g)]
                            den_U = 2*a*u**3*dd['A0'] + 2*c*u**3*vv['A0']/Lam**2
                            den_Q = 2*a*u**3*dd['A2'] + 2*c*u**3*vv['A2']/Lam**2
                            den_S = 2*a*u**3*dd['spread'] + 2*c*u**3*vv['spread']/Lam**2
                            rec[('GU', g)] = f/den_U; rec[('GQ', g)] = f/den_Q
                            rec[('GS', g)] = f/den_S if den_S > 0 else np.inf
                        rows.append(rec)
            states_out += rows
            pr(f"\n  delta = {delta}, Dg = {Dg} (Lam = {Lam:.1f}, a = {a:.4f}, c = {c}):")
            pr("    state                FF/N    |  G_U (H0, H2, H3, H4, SB)              |  G_Q (H0 .. SB)                 | G_S (H2, H3, H4, SB)")
            # per (h, J): worst m
            for h, Js in (('++', (0, 2, 4)), ('+-', (4, 5, 6, 7, 8))):
                for J in Js:
                    rr = [r for r in rows if r['h'] == h and r['J'] == J]
                    r = min(rr, key=lambda r: min(r[('GU', g)] for g in gk))
                    pr(f"    {h} J={J:2d} m={r['m']:.3f}Lam  {r['FF']/r['N']:+.4f} | " + " ".join(f"{r[('GU', g)]:6.2f}" for g in gk) +
                       " | " + " ".join(f"{r[('GQ', g)]:6.3f}" for g in gk) + " | " +
                       " ".join(f"{r[('GS', g)]:7.2f}" for g in gk if g != 'H0 min-norm'))
            wr = min(rows, key=lambda r: min(r[('GU', g)] for g in gk))
            pr(f"    binding (smallest G_U over gauges): {wr['h']} J={wr['J']} m={wr['m']:.3f}Lam: G_U = "
               f"{min(wr[('GU', g)] for g in gk):.2f} c2 (best gauge {max(wr[('GU', g)] for g in gk):.2f}); "
               f"G_Q = {min(wr[('GQ', g)] for g in gk):.3f} (|Gamma| <= G_Q |nu-ij|^2: at nu = 5 that is {25*min(wr[('GQ', g)] for g in gk):.1f} c2)")
    out['states'] = states_out
    pickle.dump(out, open(os.path.join(CACHE, 'stage7_null.pkl'), 'wb'))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, 'stage7_null.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
