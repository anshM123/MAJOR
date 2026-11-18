"""
e10: the even-kernel choice (S_ff or K*) in the near-gap thresholds, D-strong (Omega^U, M02; as e3_action.py) and D-safe (NF02 / Omega^{NF,SC};
as e7_nf.py), on the same exact-skeleton states (cache/e1_<Dg>_Df*.pkl, including *_hiJ.pkl).
Kernel option: E_KERNEL = Sff | Kstar (evenkernel.py).
 * Leading heavy action: the even moments enter normalised (C+_{4,6} -> Psi_k/N_k^K; C+_2 divided by M_K(nu;7)^2).  For a flip-free K the
   normalised heavy moment is the kernel-independent preserved density (both helicities), and the normalised graviton residue of C+_2 is
   unchanged (evenkernel.residue_factor_C2 = 1).  The code applies the kernel's normalised response ratio on every node (= 1) and recomputes;
   the thresholds must coincide with e3 / e7.
 * Finite-m TTTT corrections (SPINNING_TTTT_EXT.md sects. 6-7): the even-k heavy moment of helicity h is multiplied by 1 + c2^K(k,h)/m^2,
   with c2 the largest |m^2(R-1)| measured for kernel K at beta = 1, nu <= 3, m = 39-101 (signed with the nu <= 1 sign); both a signed
   shift and an absolute bound sum_k |c2|/m^2 |P_k| are reported (P_k = per-k even part: B' flat + curvature + skeleton remainder; forward
   C4 piece with c2(4)).  S_ff at k = 2, opposite helicity, is NOT resolved by the data (isotropic-mode noise); K*'s value is used and flagged.
 * Kernel difference: |R^K* - R^Sff| <= 3e-3 (normalised, where both are resolved) -> bound 3e-3 sum_k |P_k| on F^K* - F^Sff.
Usage: E_KERNEL=Kstar python3 e10_kernel.py Dg  -> cache/e10_<Dg>_<kernel>.txt
"""
import os, sys, pickle, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
os.environ.setdefault('OMEGAU_DF', '3')
import numpy as np
from common import load, SIGMA, KS, Transforms, flat_unstripped_parts, flat_parts, gl01, peval
from adsaction import action, envk
from nf_common import NF, TransformsNF, flat_parts_nf, action_nf
import evenkernel as EK

KERNEL = os.environ.get('E_KERNEL', 'Sff'); assert KERNEL in EK.KERNELS
# measured |c2| (largest |m^2(R-1)| over nu = 0, 1, 3 and m = 39-101 [k = 2: m = 70, 101 at nu = 3], beta = 1), sign at nu <= 1
C2 = {'Kstar': {(2, '++'): +268, (2, '+-'): +177, (4, '++'): +194, (4, '+-'): +122, (6, '++'): +144, (6, '+-'): -223},
      'Sff':   {(2, '++'): +252, (2, '+-'): +177, (4, '++'): +194, (4, '+-'): +122, (6, '++'): +144, (6, '+-'): -225}}
KDIFF = 3e-3
Dg = int(sys.argv[1])
NP, NPF = 64, 400
pn, pw = gl01(NP); pf, wf = gl01(NPF)
Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
FM, PHM = load('M02'); TRM = Transforms(PHM)
FN = NF(); TRN = TransformsNF(FN)
EVEN = [i for i, k in enumerate(KS) if k % 2 == 0]
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

# kernel's normalised response ratio on the p-nodes for each even k (normalised moment: M_K/M_K); asserted = 1
for k in (2, 4, 6):
    r = EK.Mpres(Lam*pn, k + 5, KERNEL)/EK.Mpres(Lam*pn, k + 5, KERNEL)
    assert np.allclose(r, 1.0) and np.allclose(EK.residue_factor_C2(Lam*pn, KERNEL), 1.0)

def Sk_per_k(PHk, h, mfl, ske, nf=False):
    s = SIGMA[h]; u = 1/mfl**2; x = pn*np.sqrt(u); res = []
    for k in KS:
        rho, err, ctrl = ske[k]
        phi = FN.phi(k, pn) if nf else peval(PHk[k], pn)
        res.append(np.sum((s**k)*u**(k - 2)*phi*pw*envk(k, x)*rho*ctrl))
    return np.array(res)

recs = []
for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'e1_%d_Df*.pkl' % Dg))):
    DF = os.path.basename(fn).split('_Df')[1].replace('.pkl', '').replace('_hiJ', '')
    if DF not in ('3', '2'):                    # as e3/e7 in run_all.sh (Df = 2.5 has only 16 states)
        continue
    db = pickle.load(open(fn, 'rb'))
    mfls = sorted(set(kk[0] for kk in db))
    for mfl in mfls:
        u = 1/mfl**2
        for h in SIGMA:
            Js = sorted(J for (mm, J) in db if mm == mfl and ((J % 2 == 0) if h == '++' else J >= 4))
            if not Js:
                continue
            Jsu, partsU, FWu = flat_unstripped_parts(FM, PHM, h, u, max(Js), pf, wf)
            Jss, B2, B3, FWs = flat_parts(FM, h, u, max(Js))
            pfn, wfn = gl01(max(400, 2*max(Js) + 100))
            JsN, B2N, B3N, partsN = flat_parts_nf(FN, h, u, max(Js), pfn, wfn)
            for J in Js:
                r0 = db[(mfl, J)]
                st = dict(m=r0['m'], m2=r0['m2'], eta=r0['eta'], beta=r0['beta'], Delta=r0['Delta'])
                ske = {k: (r0['rho'][k], np.zeros(NP), np.ones(NP, bool)) for k in KS}
                R4, R5 = r0['fw_exact']; ske['fw'] = (R4, 0.0, R5, 0.0, True)
                i = int(np.searchsorted(Jsu, J)); i2 = int(np.searchsorted(Jss, J))
                re = action(FM, PHM, TRM, h, st, int(J), mfl, Lam, pn, pw, ske, partsU[i], FWu[i], (B2[i2], B3[i2], FWs[i2]))
                Pk = partsU[i] + re['curv'] + Sk_per_k(PHM, h, mfl, ske)
                FW = re['Fh'] + re['S'] - (B2[i2] + B3[i2] + FWs[i2] + re['curv'].sum()) - np.sum(Sk_per_k(PHM, h, mfl, ske))  # forward part (hybrid + exact)
                iN = int(np.searchsorted(JsN, J))
                rn = action_nf(FN, TRN, h, st, int(J), mfl, Lam, pn, pw, ske, partsN[iN], B2N[iN], B3N[iN])
                PkN = partsN[iN] + rn['curv'] + Sk_per_k(None, h, mfl, ske, nf=True)
                recs.append(dict(Df=DF, h=h, mfl=mfl, J=J, m=st['m'], U=dict(F=re['Fh'] + re['S'], Ns=re['Ns'], P=Pk, FW=FW),
                                 N=dict(F=rn['Fh'] + rn['S'], Ns=rn['Ns'], P=PkN, FW=0.0), bh=(2*J + 1)/mfl))
pickle.dump(recs, open(os.path.join(HERE, 'cache', 'e10_%d_%s_recs.pkl' % (Dg, KERNEL)), 'wb'))

def shifts(r, fx):
    m2 = r['m']**2; d = r[fx]; c = C2[KERNEL]
    sgn = sum(c[(k, r['h'])]/m2*d['P'][i] for i, k in zip(EVEN, (2, 4, 6))) + c[(4, r['h'])]/m2*d['FW']
    ab = sum(abs(c[(k, r['h'])])/m2*abs(d['P'][i]) for i, k in zip(EVEN, (2, 4, 6))) + abs(c[(4, r['h'])])/m2*abs(d['FW'])
    kd = KDIFF*(sum(abs(d['P'][i]) for i in EVEN) + abs(d['FW']))
    return sgn, ab, kd

GRID = sorted(set(r['mfl'] for r in recs))
def n_thr(ok, sel=lambda r: True):
    for n0 in GRID:
        if all(ok(r) for r in recs if sel(r) and r['mfl'] >= n0 - 1e-9):
            return n0
    return None

pr(f"# e10: even kernel {KERNEL} in the near-gap thresholds, Delta_g = {Dg}, Lam = {Lam:.4f}; {len(recs)} (state, helicity, skeleton) records;"
   f" m/Lam in [{GRID[0]:.2f}, {GRID[-1]:.2f}]")
pr(f"  normalised response ratio and C+_2 residue factor on all p-nodes: 1 (asserted).  c2 table ({KERNEL}): " + ", ".join(f"k={k}{h}:{v:+d}" for (k, h), v in C2[KERNEL].items())
   + ("  [k = 2 +-: K* value; S_ff unresolved]" if KERNEL == 'Sff' else ""))
for fx, name in (('U', 'D-strong Omega^U (M02)'), ('N', 'D-safe NF02 / Omega^{NF,SC}')):
    pr(f"\n## {name}")
    for delta in (0.0, 0.03):
        n0 = n_thr(lambda r: r[fx]['F'] - delta*r[fx]['Ns'] > 0)
        ns = n_thr(lambda r: r[fx]['F'] + shifts(r, fx)[0] - delta*r[fx]['Ns'] > 0)
        na = n_thr(lambda r: r[fx]['F'] - shifts(r, fx)[1] - delta*r[fx]['Ns'] > 0)
        nk = n_thr(lambda r: r[fx]['F'] - shifts(r, fx)[2] - delta*r[fx]['Ns'] > 0)
        pr(f"  delta = {delta}: n (leading, = e3/e7) = {n0};  with signed c2/m^2 shift: {ns};  with |c2|/m^2 bound: {na};  with the K*-S_ff difference bound (3e-3): {nk}")
    sh = [(shifts(r, fx), r) for r in recs]
    w = max(sh, key=lambda t: t[0][1]/t[1][fx]['Ns'])
    pr(f"  largest |c2| bound / N_s: {w[0][1]/w[1][fx]['Ns']:.3g} at m/Lam = {w[1]['mfl']:.2f}, J = {w[1]['J']}, {w[1]['h']} (F/N_s there {w[1][fx]['F']/w[1][fx]['Ns']:+.4f});"
       f" signed shift / N_s range [{min(t[0][0]/t[1][fx]['Ns'] for t in sh):+.2e}, {max(t[0][0]/t[1][fx]['Ns'] for t in sh):+.2e}];"
       f" largest kernel-difference bound / N_s {max(t[0][2]/t[1][fx]['Ns'] for t in sh):.2e}")
    for lo, hi in ((1.0, 1.05), (1.05, 1.1), (1.1, 1.2), (1.2, 1.5), (1.5, 3.01)):
        v = [r for r in recs if lo <= r['mfl'] < hi]
        if not v:
            continue
        wm = min(v, key=lambda r: r[fx]['F']/r[fx]['Ns'])
        pr(f"   m/Lam [{lo:.2f},{hi:.2f}): min F/N_s {wm[fx]['F']/wm[fx]['Ns']:+.4f} [{wm['h']} J={wm['J']} m/Lam={wm['mfl']:.2f}];"
           f" min (F - |c2| bound)/N_s {min((r[fx]['F'] - shifts(r, fx)[1])/r[fx]['Ns'] for r in v):+.4f};  max |c2| bound/N_s {max(shifts(r, fx)[1]/r[fx]['Ns'] for r in v):.2e}")
open(os.path.join(HERE, 'cache', 'e10_%d_%s.txt' % (Dg, KERNEL)), 'w').write("\n".join(out) + "\n")
