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
_HC = {}


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


def cost_entries(J, mfl, gfun, numin=0.0, numax=NU_W, gauge=True):
    """as cost(), but with a separate envelope per matrix entry: gfun(entry, x) for entry in ('+-',) (odd) or ('++', '+-', 'off')
    (even; 'off' bounds |c_off| = m|Gamma_off|), x = nu/m_state.  Same errmat structure (Frobenius of the normalised 2x2)."""
    d = _d(J, mfl)
    nu = d['nu']; x = nu/d['m']
    sel = (nu >= numin) & (nu < numax)
    g = d['gf'] if gauge else 1.0
    hs = d['hs']; Ns = d['Ns']; n = len(hs)
    M = np.zeros((n, n))
    for i, h in enumerate(hs):
        a2, a3 = d['abs'][h]
        M[i, i] = ((a2 + a3)*gfun(h, x)*g)[sel].sum()/Ns[i]
    if n == 2:
        a2, a3 = d['abs']['off']
        M[0, 1] = M[1, 0] = ((a2 + a3)*gfun('off', x)*g)[sel].sum()/np.sqrt(Ns[0]*Ns[1])
    return float(np.sqrt((M*M).sum()))


def signed_dF(J, mfl, gam, numin=0.0, numax=NU_W):
    """signed window contribution: dF matrix (normalised by sqrt(Ns_i Ns_j)) of Gamma given as gam(key, x) for key in
    ('C2', h) / ('C3H', h) (h = '++', '+-', 'off'; 'off' = c_off) or 'C2' / 'C3H' (odd), on numin <= nu < numax, with the chain's
    own density weights (ng_chain.dens: psi_k env_k / m^2; off: eo2 c_off / m^2).  Returns (dF raw, base, Ns)."""
    Q = S.Q
    Lm = E.lam_of(DG); m = mfl*Lm
    hs = ['+-'] if J % 2 else ['++', '+-']
    kk = (J, round(mfl, 9))
    if kk not in _HC:
        _HC[kk] = {h: E.hyb(DG, J, mfl, h, 1024) for h in hs}
    st = _HC[kk]
    s0 = st[hs[0]]; x, u, y = s0['x'], s0['u'], s0['y']
    nu = s0['pn']*Lm
    sel = (nu >= numin) & (nu < numax)
    env = Q.P4.env if J % 2 else Q.A.env
    n = len(hs); dF = np.zeros((n, n))
    for i, h in enumerate(hs):
        k2 = 'C2' if J % 2 else ('C2', h); k3 = 'C3H' if J % 2 else ('C3H', h)
        v = st[h]['p2']*env(h, x, 2)*gam(k2, x)/m**2 + st[h]['p3']*u*env(h, x, 3)*gam(k3, x)/m**2
        dF[i, i] = v[sel].sum()
    if n == 2:
        eo2 = np.sqrt(Q.A.env('++', x, 2)*Q.A.env('+-', x, 2))
        v = s0['p2']*eo2*gam(('C2', 'off'), x)/m**2 + s0['p3']*u*eo2/(2 - y)*gam(('C3H', 'off'), x)/m**2
        dF[0, 1] = dF[1, 0] = v[sel].sum()
    return dF, np.array([st[h]['base'] for h in hs]), np.array([st[h]['Ns'] for h in hs])
