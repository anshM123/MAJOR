"""
sn_window (SMALLNU_WINDOW.md): sign-free cost of the small-nu window nu < NU_W (= 50) of a near-gap hybrid state at Delta_g = 1000
for a nu-dependent envelope |Gamma(nu)| <= G(nu), in EXACTLY NEARGAP2's E_low convention (neargap2/ng_bind.errmat / absdens /
gauge_factor, imported read-only): per node |weight| per unit |Gamma| (c2) of the k = 2 and k = 3 parts, times the min-norm gauge
factor max(1, |Re X_min/psi2|) (<= 6.8 at nu < 5, <= 1.24 on 5-12, 1 above 25), normalised by N_s; for even J the 2x2 (++, +-) matrix
with the off-diagonal (|c_off| = m|Gamma_off| <= G), Frobenius norm.  1024 hybrid nodes (NEARGAP2 §1.1).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sn_core as S
import numpy as np
import ng_bind as NB
E = S.E
DG = 1000.0
NU_W = 50.0
_DC = {}


def _d(J, mfl, npn=1024):
    k = (J, round(mfl, 9), npn)
    if k not in _DC:
        Lm = E.lam_of(DG)
        hs = ['+-'] if J % 2 else ['++', '+-']
        st = {h: E.hyb(DG, J, mfl, h, npn) for h in hs}
        nu = st[hs[0]]['pn']*Lm
        _DC[k] = dict(hs=hs, Ns=np.array([st[h]['Ns'] for h in hs]), nu=nu, abs=NB.absdens(DG, J, mfl, npn),
                      gf=NB.gauge_factor(DG, nu), m=mfl*Lm)
    return _DC[k]


def cost(J, mfl, G, numin=0.0, numax=NU_W, gauge=True):
    """sign-free cost (lam/N_s units) of |Gamma| <= G(nu) on numin <= nu < numax.  G: scalar or callable of nu."""
    d = _d(J, mfl)
    nu = d['nu']
    Gv = G(nu) if callable(G) else np.full_like(nu, float(G))
    Gv = Gv*(d['gf'] if gauge else 1.0)
    sel = (nu >= numin) & (nu < numax)
    return NB.errmat(d, sel, Gv)


def weight(J, mfl, numin=0.0, numax=NU_W, gauge=True):
    """W = cost per unit uniform |Gamma| (c2)."""
    return cost(J, mfl, 1.0, numin, numax, gauge)


def envelope(xn, gn, k=1.0, g_lo=1e4):
    """piecewise-constant MEASURED envelope in x: on [x_i, x_{i+1}] the larger of the two node bounds g_i (= |Gamma_i| + 3 sigma_i,
    max over kinds/entries), times k; below the lowest node g_lo (not scaled); beyond the last node the last bound (times k).
    Returns a function of x."""
    o = np.argsort(xn); xn = np.asarray(xn, float)[o]; gn = np.asarray(gn, float)[o]
    seg = np.maximum(gn[:-1], gn[1:]) if len(gn) > 1 else gn

    def f(x):
        x = np.asarray(x, float)
        i = np.clip(np.searchsorted(xn, x, side='right') - 1, 0, max(len(seg) - 1, 0))
        g = k*np.where(x >= xn[-1], gn[-1], seg[i] if len(seg) else gn[0])
        return np.where(x < xn[0], g_lo, g)
    return f


def cost_x(J, mfl, fx, numin=0.0, numax=NU_W, gauge=True):
    """cost of an envelope given as a function of x (= nu/m_state; the chain's rows are used at fixed x)."""
    m = mfl*E.lam_of(DG)
    return cost(J, mfl, lambda nu: fx(nu/m), numin, numax, gauge)
