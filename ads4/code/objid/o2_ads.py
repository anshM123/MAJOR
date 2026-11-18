"""
o2: object identity (H14) in the near-gap AdS action (Delta_g = 100), as a bulk-point-factor model on top of the exact-skeleton action.
For object X the moment is modelled as R_X(k, x) x (exactskel model moment), i.e. the exact scalar C-skeleton moment multiplied by the
ratio of bulk-point densities (o1_flat.R) [M: leading order in 1/m at fixed x; R_X = 1 + O(x^4) in the Regge regime].  Then
   F^X = F_ex(e3/e7) + dF_X,   dF_X = sum_k sum_i (R_X(k,x_i) - 1) u^{k-2} phi_k(p_i) w_i [ek D_h + 2 eps_k sigma^k (P - J0) + sigma^k env_k rho_k],
on the 64 nodes of e1 (exact rho_k, Df = 3 and 2).  Forward terms (nu = 0, R = 1) unchanged.  Validation: the same node sum with R = 1 against
e7/e3's F_ex (which use 400-node flat parts and exact transforms).
Functionals: NF02 (D-safe, Omega^{NF,SC}) with X in {C, mix (k=2 C, k>=3 Psi), Psi, PsiPure, Pi}; M02 (D-strong, Omega^U) with {C, Pi, Psi, mix}.
Thresholds n_X(delta) over all e1 states (J <= 10 and the high-J sample), and the bounds at n_X (formulas of exactskel e6/e9, recomputed here).
Usage: python3 o2_ads.py  -> cache/o2_ads.txt, cache/o2_recs.pkl
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ('omegau', 'exactskel'):
    sys.path.insert(0, os.path.join(HERE, '..', p))
sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
from common import gl01, load, peval, Dh_one, ek, EPS, SIGMA, KS
from adsaction import envk
from nf_common import NF
import remk
from o1_flat import R

Dg = 100; Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
pn, pw = gl01(64)
NFF = NF(); M02 = load('M02')[1]
OBJ = {'NF02': ('C', 'mix', 'Psi', 'PsiPure', 'Pi'), 'M02': ('C', 'Pi', 'Psi', 'mix')}
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def integrand(fun, h, r0):
    mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]; x = pn*np.sqrt(u); y = x*x
    P = remk.Pders_vec(Lam*pn, r0['eta'], nmax=0)[0]; Dh = Dh_one(h, J, y); jz = j0((2*J + 1)*x)
    I = {}
    for k in KS:
        ph = NFF.phi(k, pn) if fun == 'NF02' else peval(M02[k], pn)
        I[k] = u**(k - 2)*ph*pw*(ek(k, y, s)*Dh + 2*EPS[k]*s**k*(P - jz) + s**k*envk(k, x)*r0['rho'][k])
    return x, I

def main():
    e7 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e7_%d_recs.pkl' % Dg), 'rb'))
    e3 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e3_%d_recs.pkl' % Dg), 'rb'))
    recs = []; val = []
    for fn in sorted(glob.glob(os.path.join(HERE, '..', 'exactskel', 'cache', 'e1_%d_Df*.pkl' % Dg))):
        DF = float(os.path.basename(fn).split('_Df')[1].replace('.pkl', '').replace('_hiJ', ''))
        db = pickle.load(open(fn, 'rb'))
        for (mfl, J), r0 in db.items():
            for h in ('++', '+-'):
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                for fun in ('NF02', 'M02'):
                    if fun == 'NF02':
                        ref = [r for r in e7 if r['h'] == h and abs(r['mfl'] - mfl) < 1e-9 and r['J'] == J and float(r['Df']) == DF]
                    else:
                        ref = [r for r in e3.get(DF, {}).get(('M02', h), []) if abs(r['mfl'] - mfl) < 1e-9 and r['J'] == J]
                    if not ref:
                        continue
                    ref = ref[0]
                    x, I = integrand(fun, h, r0)
                    base = sum(I[k].sum() for k in KS)
                    rec = dict(fun=fun, h=h, mfl=mfl, J=J, Df=DF, Ns=ref['Ns'], Fex=ref['Fex'], base64=base)
                    for obj in OBJ[fun]:
                        rec['dF_' + obj] = sum(((R(obj, k, x) - 1)*I[k]).sum() for k in KS)
                    recs.append(rec)
    pickle.dump(recs, open(os.path.join(HERE, 'cache', 'o2_recs.pkl'), 'wb'))
    return recs

def report(recs):
    pr("# o2: H14 object-identity model on the exact-skeleton near-gap action, Delta_g = 100 (Df = 3, 2; e1 states incl. high-J sample)")
    # validation of the 64-node smeared sum: NF02 has no forward terms -> base64 ~ Fex (flat part 400 vs 64 nodes)
    v = [r for r in recs if r['fun'] == 'NF02']
    d = np.array([abs(r['base64'] - r['Fex'])/r['Ns'] for r in v])
    pr(f"  validation NF02: |64-node sum - F_ex|/N_s: median {np.median(d):.1e}, 95% {np.quantile(d, 0.95):.1e}, max {d.max():.1e} "
       f"(max at m/Lam={v[int(d.argmax())]['mfl']}, J={v[int(d.argmax())]['J']})")
    pr("  binding-state table (Df = 3): F/N_s for each object")
    for fun in ('NF02', 'M02'):
        for h in ('++', '+-'):
            for mfl in (1.0, 1.04, 1.1, 1.3):
                for J in (0, 2, 4, 6):
                    rr = [r for r in recs if r['fun'] == fun and r['h'] == h and abs(r['mfl'] - mfl) < 1e-9 and r['J'] == J and r['Df'] == 3.0]
                    if not rr:
                        continue
                    r = rr[0]
                    pr(f"   {fun} {h} m/Lam={mfl:.2f} J={J}: " + " ".join(f"{o} {(r['Fex'] + r['dF_' + o])/r['Ns']:+.4f}" for o in OBJ[fun]))
    grid = sorted(set(r['mfl'] for r in recs))
    pr("  thresholds n_X(delta) over all scanned states (both Df), and the worst state")
    res = {}
    for fun in ('NF02', 'M02'):
        for obj in OBJ[fun]:
            for delta in (0.0, 0.03):
                v = [r for r in recs if r['fun'] == fun]
                F = lambda r: r['Fex'] + r['dF_' + obj]
                n = next((n0 for n0 in grid if all(F(r) - delta*r['Ns'] > 0 for r in v if r['mfl'] >= n0 - 1e-9)), None)
                w = min(v, key=lambda r: F(r)/r['Ns'])
                bad = [r for r in v if F(r) - delta*r['Ns'] <= 0]
                res[(fun, obj, delta)] = n
                pr(f"   {fun} {obj:7s} delta={delta}: n = {n}; min F/N_s = {F(w)/w['Ns']:+.4f} ({w['h']} m/Lam={w['mfl']} J={w['J']} Df={w['Df']}); "
                   f"{len(bad)} non-positive" + (f", highest m/Lam {max(r['mfl'] for r in bad)}" if bad else ""))
    return res

def bounds(res):
    from scipy.integrate import quad
    import json
    from fractions import Fraction as Fr
    import evenkernel as EK
    pr("  bounds at the thresholds (formulas of exactskel e9 / e6; S_ff kernel factor)")
    Jd = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_certified.json')))
    c = [float(Fr(x)) for x in Jd['M02']['psi2_exact']]; k0 = 10884217/396000
    for (fun, obj, delta), n in sorted(res.items()):
        if n is None or delta != 0.0:
            continue
        L = Lam/n
        if fun == 'NF02':
            LA = quad(lambda t: float(EK.residue_factor_C2(Lam*t, 'Sff'))*float(NFF.psi2(t))/(t*t + 2.25/L**2), 0, 1, limit=400, points=[1.5/L])[0]
            X = n**8*16*LA
        else:
            LA = quad(lambda t: float(EK.residue_factor_C2(Lam*t, 'Sff'))*sum(ci*t**i for i, ci in enumerate(c))/(t*t + 2.25/L**2), 0, 1, limit=400, points=[1.5/L])[0]
            X = n**8*16*LA/k0
        pr(f"   {fun} {obj:7s}: n = {n}: sqrt(t4^2+aT^2) <= {240*np.sqrt(X)/Lam**4:.3g}")

if __name__ == "__main__":
    recs = main()
    res = report(recs)
    bounds(res)
    open(os.path.join(HERE, 'cache', 'o2_ads.txt'), 'w').write("\n".join(out) + "\n")
