"""
Shared helpers for spinfunc/: bulk-point density ratios R_X(k, x) of the graviton functionals (derived in f1_bulkpoint.py), and the
near-gap integrand of the exact-skeleton model (generalisation of objid/o2_ads.integrand to Delta_g = 100, 300, 1000).
Imports omegau/, exactskel/, objid/ read-only.

Objects (density relative to the B' density = CLS-D's flat sum rule, leading order in 1/m at fixed x = nu/m):
  'A'      route A: CLS-D weights f^{+-}_{k;2} + eta-inversion (H0_SPINNING_FUNCTIONAL.md); R = 1 for all k        [Q CLS-D + D-sketch]
  'PsiG'   the functional as currently defined (Omega^{NF,SC}): C+_2 kept (R = 1, assumed), Psi_3, Psi_7 = CLS-D weights (R = 1),
           Psi_4, Psi_6 (weight r^{k-3}) and Psi_5 (weight (u'+v')/(32u'v')): R_k of f1                                [D given CLS-D]
  'PsiG2'  as PsiG but C+_2 taken as the bare Regge moment Pi-hat_2 (sensitivity only; R_2 = (1-r^2)^2/(1+r^2))      [M]
  'PsiS'   OBJECT_IDENTITY's scalar-analogue 'mix' (k = 2 -> C, k >= 3 -> scalar Psi^A), for comparison                  [superseded]
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ('omegau', 'exactskel', 'objid'):
    sys.path.insert(0, os.path.join(HERE, '..', p))
import numpy as np

def rstar(x):
    s = np.sqrt(np.maximum(1 - x*x, 0.0)); return (1 - s)/(1 + s)

def R_grav(k, x):
    """bulk-point density of the project's graviton Psi_k relative to B' (f1_bulkpoint.py)."""
    r2 = rstar(x)**2
    if k == 4:
        return (1 - r2)**2/(1 + r2)
    if k == 6:
        return (1 - r2)**4/(1 + r2)
    if k == 5:
        return (1 + 6*r2 + r2*r2)/(1 - r2)**2
    return np.ones_like(x)

def R(obj, k, x):
    if obj == 'A':
        return np.ones_like(x)
    if obj == 'PsiG':
        return R_grav(k, x)
    if obj == 'PsiG2':
        if k == 2:
            r2 = rstar(x)**2; return (1 - r2)**2/(1 + r2)
        return R_grav(k, x)
    if obj == 'PsiS':
        from o1_flat import R as Rs
        return Rs('mix', k, x)
    raise ValueError(obj)

OBJS = ('A', 'PsiG', 'PsiG2', 'PsiS')


class NearGap:
    """per-k, per-node pieces of the NF02 (Omega^{NF,SC}) near-gap action in the exact-skeleton model (exactskel e7 / objid o2):
    I_k(p_i) = u^{k-2} phi_k(p_i) w_i [ e_k D_h  +  2 eps_k sigma^k (P_{-1/2+i nu}(cosh beta) - J0((2J+1)x))  +  sigma^k env_k rho_k ]."""
    def __init__(self, Dg):
        from common import gl01
        from nf_common import NF
        self.Dg = Dg
        self.Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
        self.pn, self.pw = gl01(64)
        self.F = NF()

    def pieces(self, h, r0):
        from scipy.special import j0
        from common import Dh_one, ek, EPS, SIGMA, KS
        from adsaction import envk
        import remk
        mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]; pn, pw = self.pn, self.pw
        x = pn*np.sqrt(u); y = x*x
        P = remk.Pders_vec(self.Lam*pn, r0['eta'], nmax=0)[0]; Dh = Dh_one(h, J, y); jz = j0((2*J + 1)*x)
        out = {}
        for k in KS:
            pre = u**(k - 2)*self.F.phi(k, pn)*pw
            out[k] = dict(dens=pre*ek(k, y, s)*Dh, curv=pre*2*EPS[k]*s**k*(P - jz), rem=pre*s**k*envk(k, x)*r0['rho'][k])
        return x, out


def load_states(Dg):
    """all exact-skeleton states (e1_<Dg>_Df*.pkl incl. hiJ) joined with the e7 records (F_ex, N_s)."""
    import pickle, glob
    C = os.path.join(HERE, '..', 'exactskel', 'cache')
    e7 = pickle.load(open(os.path.join(C, 'e7_%d_recs.pkl' % Dg), 'rb'))
    idx = {(r['h'], round(r['mfl'], 9), r['J'], float(r['Df'])): r for r in e7}
    res = []
    for fn in sorted(glob.glob(os.path.join(C, 'e1_%d_Df*.pkl' % Dg))):
        DF = float(os.path.basename(fn).split('_Df')[1].replace('.pkl', '').replace('_hiJ', ''))
        for (mfl, J), r0 in pickle.load(open(fn, 'rb')).items():
            for h in ('++', '+-'):
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                ref = idx.get((h, round(mfl, 9), J, DF))
                if ref is None:
                    continue
                res.append((h, DF, r0, ref))
    return res
