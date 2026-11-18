"""
a2: heavy-action shift of the functionals when the SKELETON-MODEL moments are replaced by EXACT <TTTT> moments on the computed low nodes.

Model of exactskel/e3 (Omega^U, M02) and e7 (Omega^{NF,SC}, NF02), per node p_i (nu = Lam p_i, x = nu/m, y = x^2) and k, in units of
eps_k m^{2k} C^h_k / sigma^k:
   M^h_k = ek(k,y) D_h(1-2y) + 2 eps_k [P_{-1/2+i nu}(cosh beta) - J_0((2J+1)x)] + env_k rho^skel_k        (D_++ = P_J, D_+- = d^J_44)
with rho^skel_k = exact scalar-skeleton remainder (exactskel e1, Df = 3 / 2).  Exact graviton replacement, with the functional's own
normalisation (normalised replacement = leading heavy formula exactly 2P/m^{2k} |c|^2 by construction):
   G^h_k = eps_k 2 P R_grav(k, nu; wave h, kernel K_k),   K_k = FF (even k, NF) or E (even k, M02);  O (odd k, both).
Shift of the action on the computed nodes:  dF = sum_k sum_{i computed} sigma^k u^{k-2} phi_k(p_i) w_i (G - M).
Rescale-invariant version: R_grav -> R_grav / c*, c* = median over computed (k, nu) of R_grav/R_skel (an overall factor common to all k, nu
multiplies F and N_s alike and cannot change a sign).
Coverage: share of sum_i |u^{k-2} phi_k w_i env_k| (the remainder norm of the model) carried by the computed nodes.
Usage: python3 a2_action.py cache/a1_<tag>.pkl [...]   -> cache/a2_summary.txt (appends per call: printed)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU); sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
TT = os.path.join(HERE, '..', 'fourpt', 'tttt')
for p in (TT, os.path.join(TT, '..'), os.path.join(TT, '..', 'oddk'), os.path.join(TT, '..', '..', 'e2e')):
    sys.path.insert(0, p)
import numpy as np
from scipy.special import j0
from common import gl01, load, peval, Dh_one, ek, EPS, SIGMA, KS
from adsaction import envk
from nf_common import NF
from t16_kernels_ext import rel

pn, pw = gl01(64)
NFF = NF(); M02 = load('M02')[1]
WAVE = {'++': ('same', 'same'), '+-': ('opp', 'opp')}


def weights(fun, k, u):
    x = pn*np.sqrt(u)
    ph = NFF.phi(k, pn) if fun == 'NF02' else peval(M02[k], pn)
    return u**(k - 2)*ph*pw, envk(k, x)


def analyse(a1, fun, h, Df, e1, e7, keven=None):
    Dg, mfl, J = a1['Dg'], a1['mfl'], a1['J']
    u = 1/mfl**2; s = SIGMA[h]; Lam = a1['Lam']; m = a1['m']
    r0 = e1[(round(mfl, 6), J)]
    assert abs(r0['Delta'] - a1['Delta']) < 1e-6, (r0['Delta'], a1['Delta'])
    inodes = a1['inodes']; nus = a1['NUS']
    P = a1['P']
    wave = WAVE[h]
    dF = 0.0; Fpart = 0.0; cov_num = cov_den = 0.0; ratios = []; rows = []
    for k in KS:
        wk, env = weights(fun, k, u)
        cov_den += np.sum(np.abs(wk*env)); kern = 'O' if k % 2 else (keven or ('FF' if fun == 'NF02' else 'E'))
        for j, i in enumerate(inodes):
            nu = nus[j]; x = nu/m; y = x*x
            Rg = rel(*a1['grav'][(k, nu)][(wave, kern)]).real
            Rs = a1[('skel', Df)][(k, nu)]['R_ex']
            Dh = Dh_one(h, J, y)
            M = ek(k, y)*Dh + 2*EPS[k]*(P[j] - j0((2*J + 1)*x)) + env[i]*r0['rho'][k][i]
            G = EPS[k]*2*P[j]*Rg
            dF += s**k*wk[i]*(G - M); Fpart += s**k*wk[i]*M
            cov_num += abs(wk[i]*env[i]); ratios.append(Rg/Rs)
            rows.append((k, nu, Rg, Rs, (G - M)/env[i]))
    c = float(np.median(ratios))
    # uniform pointwise remainder |delta rho| <= d on ALL 64 nodes that would just cancel the skeleton margin: d_crit = F_ex / sum |w env|
    full_norm = sum(np.sum(np.abs(weights(fun, k, u)[0]*weights(fun, k, u)[1])) for k in KS)
    dF_shape = dF - (c - 1)*Fpart
    rec = [r for r in e7 if r['h'] == h and abs(r['mfl'] - mfl) < 1e-6 and r['J'] == J and float(r['Df']) == Df] if e7 is not None else []
    return dict(full_norm=full_norm, dF=dF, dF_shape=dF_shape, c=c, Fpart=Fpart, cov=cov_num/cov_den, rows=rows, e7=rec[0] if rec else None)


if __name__ == "__main__":
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    for fn in sys.argv[1:]:
        a1 = pickle.load(open(fn, 'rb'))
        Dg = a1['Dg']
        pr(f"## {os.path.basename(fn)}: Dg={Dg} m/Lam={a1['mfl']} J={a1['J']} m={a1['m']:.3f} beta={a1['beta']:.4f}; computed nodes nu <= {max(a1['NUS'][:len(a1['inodes'])]):.2f} ({len(a1['inodes'])} of 64)")
        for Df in (3.0, 2.0):
            e1 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e1_%d_Df%g.pkl' % (Dg, Df)), 'rb'))
            e7 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e7_%d_recs.pkl' % Dg), 'rb'))
            e3 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e3_%d_recs.pkl' % Dg), 'rb'))
            for fun in ('NF02', 'M02'):
                for h in ('++', '+-'):
                    if (h == '++' and a1['J'] % 2) or (h == '+-' and a1['J'] < 4):
                        continue
                    if fun == 'NF02':
                        recs = [dict(r, Df=float(r['Df'])) for r in e7]
                    else:
                        recs = [dict(r, Df=Df, h=h) for r in e3.get(Df, {}).get(('M02', h), [])]
                    for keven in ((None, 'KS') if fun == 'NF02' else (None,)):
                        res = analyse(a1, fun, h, Df, e1, recs, keven)
                        ref = res['e7']
                        lab = fun + ('[K*]' if keven else '')
                        mar = f"F/Ns(skel) = {ref['Fex']/ref['Ns']:+.4f}, Ns = {ref['Ns']:.4g}" if ref else "(margin: see e3)"
                        Ns = ref['Ns'] if ref else np.nan
                        pr(f"  Df={Df:g} {lab} {h}: dF = {res['dF']:+.3e} ({res['dF']/Ns:+.2e} Ns)  c* = {res['c']:.5f}  dF_shape = {res['dF_shape']:+.3e} ({res['dF_shape']/Ns:+.2e} Ns)"
                           f"  coverage {res['cov']:.3f};  {mar}" + (f";  d_crit = {ref['Fex']/res['full_norm']:.2e}" if ref else ""))
                        dr = {}
                        for (k, nu, Rg, Rs, d) in res['rows']:
                            dr.setdefault(k, []).append((nu, d))
                        pr("      delta rho_k (graviton - model, per node; nu: value): " + " | ".join(
                            f"k{k}: " + " ".join(f"{d:+.1e}" for nu, d in sorted(v)[::3]) + f" [max|.| {max(abs(d) for nu, d in v):.1e}]" for k, v in sorted(dr.items())))
    open(os.path.join(HERE, 'cache', 'a2_summary.txt'), 'a').write("\n".join(out) + "\n")
