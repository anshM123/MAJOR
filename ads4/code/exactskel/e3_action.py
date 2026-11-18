"""
e3: AdS heavy action of Omega^U (M02; M01, F_g3 for comparison) on the e1 states with the EXACT scalar skeleton.

Same model and code as OMEGA_U (omegau/adsaction.py, read-only), with the skeleton remainder replaced:
   eps_k m^{2k} C^h_k = e_k D_h + 2 eps_k sigma^k [P - J_0] + sigma^k env_k rho_k^exact      (all 64 nodes; no control cut)
   forward: R4, R5 exact.
So  F_exact = Fh + S_exact  (normal form), margin = F_exact/Ns.  What remains [I]: B' dictionary / reading B (I1, H2), transfer of the
scalar remainder to gravitons and to +- (H4), odd-k two-pair moments and odd B' envelope (H3) -- unchanged from OMEGA_U.
Also, from the truncated series on the same states: OMEGA_U's E, U, kappa = (Fh - E)/U, and the signed exact remainder on the
series-uncontrolled nodes per unit U ("kappa_eff": the relative error the tolerance convention had to cover).
Thresholds over the scanned set (m/Lam <= MFLMAX, J <= JMAX):
   n_exact(delta) = smallest grid m/Lam such that F_exact/Ns > delta for all scanned states with m >= n Lam;
   n_series(kappa_req, delta) as in OMEGA_U u7 (F_hyb - E - kappa_req U - delta Ns > 0), restricted to the same states.
Usage: python3 e3_action.py Dg Df [Df ...]   -> cache/e3_<Dg>.txt, cache/e3_<Dg>_recs.pkl
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
os.environ.setdefault('OMEGAU_DF', '3')
import numpy as np
from common import load, FUNS, SIGMA, KS, Transforms, flat_unstripped_parts, flat_parts, gl01, peval
from adsaction import action, envk

Dg = int(sys.argv[1]); DFS = [float(a) for a in sys.argv[2:]] or [3.0, 2.0]
JCUT = int(os.environ.get('E3_JCUT', '1000'))
NP, NPF = 64, 400
pn, pw = gl01(NP)
pf, wf = gl01(NPF)
Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
Fs = {n: load(n) for n in FUNS}
TRs = {n: Transforms(PH) for n, (F, PH) in Fs.items()}
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def dodd_E(PH, h, J, rho_ex, mfl):
    """odd-k reading 'E' (I3 alternative): odd moments = eps_k x (exact k-continued scalar) instead of the B' odd envelope:
    shift = sum_{k odd} u^{k-2} int phi_k sigma (eps_k E_k - (1-y)^{-(k+1)/2}) (D_h + rho_k)."""
    from common import Dh_one, EPS
    s = SIGMA[h]; u = 1/mfl**2; x = pn*np.sqrt(u); y = x*x
    D = Dh_one(h, J, y)
    tot = 0.0
    for k in KS:
        if k % 2 == 0:
            continue
        Ek = (2 - y)/(1 - y)**(1 + k/2)
        tot += u**(k - 2)*np.sum(pw*peval(PH[k], pn)*s*(EPS[k]*Ek - 1/(1 - y)**((k + 1)/2))*(D + rho_ex[k]))
    return tot

def unc_exact(PH, h, st, sk_ser, rho_ex, mfl):
    """signed exact remainder on the series-uncontrolled nodes (per k summed), normal form."""
    s = SIGMA[h]; u = 1/mfl**2; x = pn*np.sqrt(u)
    tot = 0.0
    for k in KS:
        ctrl = sk_ser[k][2]
        wk = u**(k - 2)*peval(PH[k], pn)*pw*envk(k, x)
        tot += np.sum((s**k)*wk*rho_ex[k]*(~ctrl))
    return tot

allrecs = {}
for DF in DFS:
    db = pickle.load(open(os.path.join(HERE, 'cache', 'e1_%d_Df%g.pkl' % (Dg, DF)), 'rb'))
    fh = os.path.join(HERE, 'cache', 'e1_%d_Df%g_hiJ.pkl' % (Dg, DF))
    if os.path.exists(fh):                     # high-J sample (J = 11..28 on a coarser m grid)
        db.update(pickle.load(open(fh, 'rb')))
    recs = {(n, h): [] for n in FUNS for h in SIGMA}
    mfls = sorted(set(k[0] for k in db))
    for mfl in mfls:
        u = 1/mfl**2
        for h in SIGMA:
            Js = sorted(J for (mm, J) in db if mm == mfl and ((J % 2 == 0) if h == '++' else J >= 4))
            Js = [J for J in Js if J <= JCUT]
            if not Js:
                continue
            fl = {}
            for n, (F, PH) in Fs.items():
                Jsu, parts, FWu = flat_unstripped_parts(F, PH, h, u, max(Js), pf, wf)
                Jss, B2, B3, FWs = flat_parts(F, h, u, max(Js))
                fl[n] = (Jsu, parts, FWu, Jss, B2, B3, FWs)
            for J in Js:
                r0 = db[(mfl, J)]
                st = dict(m=r0['m'], m2=r0['m2'], eta=r0['eta'], beta=r0['beta'], Delta=r0['Delta'])
                ske = {k: (r0['rho'][k], np.zeros(NP), np.ones(NP, bool)) for k in KS}
                R4, R5 = r0['fw_exact']
                ske['fw'] = (R4, 0.0, R5, 0.0, True)
                sks = r0['ser']
                for n, (F, PH) in Fs.items():
                    Jsu, parts, FWu, Jss, B2, B3, FWs = fl[n]
                    i = int(np.searchsorted(Jsu, J)); i2 = int(np.searchsorted(Jss, J))
                    args = (F, PH, TRs[n], h, st, int(J), mfl, Lam, pn, pw)
                    flat = (parts[i], FWu[i], (B2[i2], B3[i2], FWs[i2]))
                    re = action(*args, ske, *flat)
                    rs = action(*args, sks, *flat)
                    ue = unc_exact(PH, h, st, sks, r0['rho'], mfl)
                    dE = dodd_E(PH, h, J, r0['rho'], mfl)
                    recs[(n, h)].append(dict(mfl=mfl, J=J, Ns=re['Ns'], Fh=re['Fh'], Sex=re['S'], Fex=re['Fh'] + re['S'],
                                             E=rs['E'], U=rs['U'], Sser=rs['S'], Sall=rs['Sall'], Uex=ue, NU=re['NU'], FexE=re['Fh'] + re['S'] + dE))
    allrecs[DF] = recs

pickle.dump(allrecs, open(os.path.join(HERE, 'cache', 'e3_%d_recs.pkl' % Dg), 'wb'))

GRID = sorted(set(r['mfl'] for v in allrecs.values() for vv in v.values() for r in vv))

def n_thr(states, ok):
    for n0 in GRID:
        if all(ok(r) for r in states if r['mfl'] >= n0 - 1e-9):
            return n0
    return None

pr(f"# e3: Omega^U heavy action with the exact scalar skeleton, Delta_g = {Dg}, Lam = {Lam:.4f}; states m/Lam in "
   f"[{GRID[0]:.2f}, {GRID[-1]:.2f}], J <= {max(r['J'] for v in allrecs.values() for vv in v.values() for r in vv)}")
for DF, recs in allrecs.items():
    pr(f"\n## skeleton Df = {DF:g}")
    for n in FUNS:
        for h in SIGMA:
            v = recs[(n, h)]
            if not v:
                continue
            w = min(v, key=lambda r: r['Fex']/r['Ns'])
            pr(f"  {n} {h}: {len(v)} states; min F_exact/Ns = {w['Fex']/w['Ns']:+.4f} at m/Lam={w['mfl']:.2f}, J={w['J']}"
               f"  (Fh/Ns {w['Fh']/w['Ns']:+.4f}, S_exact/Ns {w['Sex']/w['Ns']:+.4f})")
    pr("  M02 per band (both helicities): min F_exact/Ns [state] | min Fh/Ns | OMEGA_U kappa_abs min | min kappa_eff = S_exact(unc)/U")
    for lo, hi in ((1.0, 1.05), (1.05, 1.1), (1.1, 1.2), (1.2, 1.3), (1.3, 1.4), (1.4, 1.5), (1.5, 1.75), (1.75, 2.01)):
        for h in SIGMA:
            vb = [r for r in recs[('M02', h)] if lo <= r['mfl'] < hi]
            if not vb:
                continue
            w = min(vb, key=lambda r: r['Fex']/r['Ns'])
            kap = [(r['Fh'] - r['E'])/r['U'] for r in vb if r['U'] > 0]
            ke = [r['Uex']/r['U'] for r in vb if r['U'] > 0]
            wE = min(vb, key=lambda r: r['FexE']/r['Ns'])
            pr(f"   [{lo:.2f},{hi:.2f}) {h}: {w['Fex']/w['Ns']:+.4f} [m/Lam={w['mfl']:.2f}, J={w['J']}] (reading E: {wE['FexE']/wE['Ns']:+.4f} [{wE['mfl']:.2f}, {wE['J']}]) | {min(r['Fh']/r['Ns'] for r in vb):+.4f} | "
               f"{min(kap) if kap else float('inf'):.3g} | {min(ke) if ke else float('nan'):+.3g} (max |.| {max(map(abs, ke)) if ke else float('nan'):.3g})")

pr("\n## thresholds over the scanned states (worst over both helicities and the skeletons listed)")
for n in FUNS:
    for delta in (0.0, 0.03):
        st = [r for DF in allrecs for h in SIGMA for r in allrecs[DF][(n, h)]]
        ne = n_thr(st, lambda r: r['Fex'] - delta*r['Ns'] > 0)
        bad = [r for r in st if r['Fex'] - delta*r['Ns'] <= 0]
        pr(f"  {n}: n_exact(delta={delta}) = {ne}   ({len(bad)} non-positive states" +
           (f"; highest at m/Lam={max(r['mfl'] for r in bad):.2f})" if bad else ")"))
        neE = n_thr(st, lambda r: r['FexE'] - delta*r['Ns'] > 0)
        neB = n_thr(st, lambda r: min(r['FexE'], r['Fex']) - delta*r['Ns'] > 0)
        pr(f"        odd-k reading E (eps_k x continued scalar): n = {neE};  worst of both readings: n = {neB}")
    for kr in (0.1, 0.3, 1.0):
        for delta in (0.0, 0.03):
            st = [r for DF in allrecs for h in SIGMA for r in allrecs[DF][(n, h)]]
            ns = n_thr(st, lambda r: (r['Fh'] - r['E'] - kr*r['U'] - delta*r['Ns']) > 0)
            pr(f"     series criterion (OMEGA_U) on the same states: n(kappa_req={kr}, delta={delta}) = {ns}")
pr("\n## per-state M02 margins (F_exact/Ns), rows m/Lam, columns J; ++ then +-")
for DF, recs in allrecs.items():
    for h in SIGMA:
        v = recs[('M02', h)]
        Js = sorted(set(r['J'] for r in v))
        pr(f"  Df={DF:g} {h}   J: " + " ".join(f"{J:7d}" for J in Js))
        for mfl in sorted(set(r['mfl'] for r in v)):
            row = {r['J']: r['Fex']/r['Ns'] for r in v if r['mfl'] == mfl}
            pr(f"   {mfl:5.2f}  " + " ".join(f"{row[J]:+7.4f}" if J in row else "      -" for J in Js))
open(os.path.join(HERE, 'cache', 'e3_%d.txt' % Dg), 'w').write("\n".join(out) + "\n")
