"""
edger/er_core.py (EDGE_R.md): shared helpers for the graviton (EDGE-R) quantity of AG_PROOF Lemma P4W-chi,

    E_chi[G_e] = int_{p_w}^1 psi_k chi G_e dp,   chi(p) = erfc((p_e - p)/sigma)/2,   cost_E,e = 2|E_chi,2| + u|E_chi,3|,

evaluated from REAL-AXIS graviton data in the EE pipeline's convention (HIGHSPIN / ee1000): Gamma_{k,h}(x) in c2 units (= m^2 x the
correction relative to the leading envelope env_{h,k}(x), kappa-curvature image subtracted = exact - hybrid), so that the k-th cost of the
entry is exactly the EE correction restricted by chi:

    cost_k = int chi(p(x)) w_k(x) Gamma_k(x) dx,   w_k = psi_k(p) env_{h,k}(x) [u] / m^2 * dp/dx     (hs_core.state_weights)

(2 int psi_2 G = int psi_2 env_2 Gamma_2/m^2 dp: G_e = env_2 Gamma_2/(2 m^2) in P4W units; the k = 3 cost u int psi_3 G_3 with
G_3 = env_3 Gamma_3/m^2).  Units: P4W's H (= f++ f+-/(f++ + f+-) even J, f+- odd J; f = N_flat/5 + min(0, curv), lsw_common.floors)
and the EE pipeline's N_s share psi_k and the leading-image normalisation, so every cost below is directly comparable with
AG_PROOF's tolerance eps_E tol = (1 - d_chi) H and with S_E (plain) / S_E,t (+- transferred) of ag_chi.scale_SE.
All project modules are imported read-only; er_guard (imported first) confines writes to edger/cache.
"""
import er_guard  # noqa: F401  (first)
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
CACHE = os.path.join(HERE, 'cache')
for sub in ('highspin', 'agproof', 'p4win', 'lsw', 'spinfunc'):
    p = os.path.join(CODE, sub)
    if p not in sys.path:
        sys.path.append(p)
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.special import erfc

Dg = 1000
_CHI = None


def chi_rows():
    """AG_PROOF's 47 Theorem W rows (ag_chi cache): cutoff (p_e, sigma, T), d_chi/H, H, tolerances, S_E plain / transferred."""
    global _CHI
    if _CHI is None:
        _CHI = json.load(open(os.path.join(CODE, 'agproof', 'cache', 'chi.json')))
    return _CHI


def w_state(J, mfl):
    for r in chi_rows():
        if r['J'] == J and abs(r['mfl'] - mfl) < 2e-3:
            return r
    return None


def nearest_cut(J, mfl):
    """the cutoff of the nearest Theorem W state (same parity class first, then |dJ|, then |d mfl|)."""
    rs = sorted(chi_rows(), key=lambda r: (abs(r['J'] - J), abs(r['mfl'] - mfl)))
    return rs[0]


def chi(p, pe, sg):
    return 0.5*erfc((pe - np.asarray(p, float))/sg)


def lam():
    import ee1000_common as E
    return E.lam_of(Dg)


def beta_of(J, mfl):
    return 2*np.arcsinh((J + 0.5)/(mfl*lam()))


def H_of(J, mfl):
    """P4W's H for an arbitrary state (lsw_common.floors; identical to jstar_full's H at the sampled states)."""
    import lsw_common as LW
    f = {h: v['floor'] for h, v in LW.floors(J, mfl).items()}
    if J % 2:
        return f['+-'], f
    return f['++']*f['+-']/(f['++'] + f['+-']), f


def SE_units(J, mfl):
    """(S_E plain, S_E,t transferred) = (2*2 V_Fa2 + u*2 V_Fa3) e_sc(1) [x (1 - x1^2)^-2]  (ag_chi.scale_SE, Lemma F profile)."""
    import ag_flat as AF
    import p4w_contour as K
    o = AF.edge_compare(J, mfl, ny=0) if 'ny' in AF.edge_compare.__code__.co_varnames else AF.edge_compare(J, mfl)
    beta = beta_of(J, mfl); u = 1/mfl**2
    e1 = float(K.e_sc(beta, mfl)(np.array([1.0]))[0])
    SE = (2*2*o['VF2'] + u*2*o['VF3'])*e1
    return SE, SE*(1 - 1/mfl**2)**-2, o


def weights(J, mfl, h, x):
    """hs_core.state_weights: per unit Gamma (c2) at x, the k = 2 and k = 3 cost densities in x."""
    import hs_core as HC
    return HC.state_weights(Dg, J, mfl, h, x)


def edge_cost(J, mfl, h, x, wx, g2, g3, pe, sg):
    """(cost_2, cost_3) = sum_i wx_i chi(p_i) w_k(x_i) Gamma_k(x_i), restricted to p >= p_w = 3/(beta Lam) and p <= 1."""
    L_ = lam(); m = mfl*L_; p = np.asarray(x)*m/L_
    pw = 3.0/(beta_of(J, mfl)*L_)
    sel = (p >= pw) & (p <= 1.0)
    w2, w3 = weights(J, mfl, h, x)
    c = chi(p, pe, sg)*sel*wx
    return float((c*w2*g2).sum()), float((c*w3*g3).sum())


def gl(a, b, n):
    t, w = np.polynomial.legendre.leggauss(n)
    return a + 0.5*(b - a)*(t + 1), 0.5*(b - a)*w
