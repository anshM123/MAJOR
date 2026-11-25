"""
qb_bind (QBIND_COROLLARY.md §1): error budget of the BINDING states of Result Q's upper functional Omega^Q_up (delta = 0.2,
'Q3-full-delta0.2', the published functional) at Delta_g = 1000, with the protocol NEARGAP2 used for Omega^EE (neargap2/ng_bind.py),
extended to the pieces Omega^Q has and Omega^EE does not (the forward k = 4, 5 pieces).

Per-state value (RESULT_Q_STATEMENT §1, RESULT_Q_REVISION §1.1, BELOW_FOLD2 §2, DESHIFT_CLOSE §1.3):
    FF_h = HYB_h(B2, B3, FW, B4)                    flat with exact Wigner functions + leading-Regge AdS curvature (stage6_model.hyb)
           - a 2u^3 s_d(h, J)                       measured k = 4 derivative correction (analytic-radial engine; FWD_CLOSE, PARITY_ODD_CLOSE)
           - c (2u^3/Lam^2) s_v                     measured k = 4 value correction (BELOW_FOLD2 §2: s_v ~ -13.4)
           + dB2_h + dB3_h                          smeared k = 2, 3 graviton corrections from the K2 rows (Result D's objects, Q's psi2, psi3)
    lambda = min eig of  N^{-1/2} [diag(FF_h) + dF_offdiag] N^{-1/2}   (1x1 for J < 4 and odd J, 2x2 (++, +-) for even J >= 4)
The smeared corrections are computed on the OWN-MASS K2 rows (qb_rows.py, plus NEARGAP2's own rows J = 4, 6, read-only) through the
UNCHANGED EE_NEARGAP chain (qb_chain = ng_chain), exactly as ng_bind, with Q's weights in place of NF02's.

Error budget (sign-free; x3 on truncation / interpolation / measured forward terms as in ng_bind; G conventions stated):
  E_trunc  3 x sum_nodes |w_Q| |Gamma_K2 - Gamma_K1| on the own row's x-range                         (as ng_bind)
  E_interp 3 x max(|lam - lam_linear-x|, |lam - lam_every-other-node|)                                (as ng_bind)
  E_low    G_low = 1e4 c2 x W_Q(nu < max(50, x0 m)) x max(1, |Re X^Q_minnorm / psi2^Q|) (Q's own min-norm Abel weight, rv_k2eng_Y)
  E_edge   G_edge = max(1e4, 10 max|Gamma| last 4 own-row nodes) x W_Q(x > x_max)                        (as ng_bind)
  E_sd     3 x (|truncation of s_d| + m-transfer spread) x a 2u^3/N                                     (measured, engine; HRS_CLOSE §6)
  E_sv     3 x c (2u^3/Lam^2) |s_v| / N
  E_C5     a u^3 A5 G5/(m^2 N), G5 = 1e4 c2 (no k = 5 engine; DESHIFT_CLOSE §1.3 bound), A5 = 1.80 (XP = published gauge; 7.66 worst
           Stage-6 gauge also printed)
  E_OO     (OO states) 0.12 c2 x W_Q(all nodes) + a 2u^3 x 5e-6 / N  (PARITY_ODD_CLOSE §1 (F), RESULT_Q_FWD_CLOSE §2.3)
  E_off    (even J >= 4) the forward off-diagonal |s_off| (measured bound) enters the 2x2 as a 2u^3 |s_off|; included in the matrix
  Base     flat value (polynomial quadrature) and LR curvature (tables); Lemma C bound on the curvature printed separately.
Usage: python3 qb_bind.py tag J:h:mfl[,mfl...]:sd:dsd[:soff] ...    (h = ++, +-, OO)   -> cache/qb_bind_<tag>.{txt,json}
"""
import os, sys, json, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qb_chain as NC                                   # chdir(spinfunc), guard installed
RQ = os.path.join(HERE, '..', 'resultq')
sys.path.append(RQ)
import numpy as np
E, Q, W = NC.E, NC.Q, NC.W
NGC = NC.NGC
DG = 1000.0
G_LOW = 1.0e4
NU_LOW = 50.0
G5 = 1.0e4
A5 = {'XP': 1.80, 'new': 1.91, 'worst': 7.66}
S_V = 13.4                                              # |s_v| (BELOW_FOLD2 §2)
QNAME = 'Q3-full-delta0.2'

import stage6_model as M                               # read-only (tables cached in resultq/cache)
import common_q as CQ
FQ, C2Q, C3Q = M.load(QNAME)
AQ, CQc = FQ['a'], FQ['c']
LAM = M.lam(DG)
TAB = M.Tables(QNAME, LAM, 2300.0)


def psi2(p):
    return np.polynomial.polynomial.polyval(p, C2Q)


def psi3(p):
    return np.polynomial.polynomial.polyval(p, C3Q)


# ---------------------------------------------------------------------------------------------------------------- Q base
def qbase(J, mfl, h):
    """HYB pieces of Omega^Q at (J, m/Lam, h): FF, N, parts, flat FF/N and the curvature part; h in ('++', '+-')."""
    u = 1/mfl**2
    Js, B2, B3 = CQ.parts(h, u, max(J, 4), [C2Q], [C3Q]); k = list(Js).index(J)
    FF, N, _, parts = M.hyb(FQ, h, u, float(J), B2[k, 0], B3[k, 0], TAB, LAM)
    # flat (no curvature): the same normal form with flat forward pieces
    c4 = M.C4[h]; s = M.SIG[h]; JJ = J*(J + 1.0)
    B4f = 2*u*u; dB4f = (c4 - 2*JJ)*u**3; B5f = s*u**3
    FWf = -AQ*(dB4f - B5f)
    FFf = B2[k, 0] + B3[k, 0] + FWf - CQc*B4f
    Nf = abs(B2[k, 0]) + abs(B3[k, 0]) + abs(FWf) + CQc*abs(B4f)
    return dict(FF=float(FF), N=float(N), B2=float(parts['B2']), B3=float(parts['B3']), FW=float(parts['FW']), B4=float(parts['B4']),
                FFf=float(FFf), Nf=float(Nf), B2f=float(B2[k, 0]), B3f=float(B3[k, 0]), u=u)


def lemmaC(J, mfl, h, npn=256):
    """sign-free Lemma C bound (POS_ERRORS §2.1) on the LR curvature of the smeared parts with Q's weights, and the forward P0 bound."""
    pn, pw = Q.gl01(npn); u = 1/mfl**2; m = mfl*LAM
    beta = 2*np.arcsinh((J + .5)/m)
    dP = (1 - 1/np.cosh(beta/2)) + np.minimum(2.0, pn*LAM*(2*np.sinh(beta/2) - beta))
    sm = (np.abs(psi2(pn))*pw*(2 + u*pn*pn)*dP).sum() + (np.abs(psi3(pn))*pw*u*dP).sum()
    dP0 = 1 - 1/np.cosh(beta/2)
    fw = (CQc*2*u*u + AQ*(M.C4[h] + 1)*u**3)*dP0
    return float(sm), float(fw)


# ---------------------------------------------------------------------------------------------------------------- weights
def qnodes(J, mfl, npn=256):
    pn, pw = Q.gl01(npn)
    u = 1/mfl**2; m = mfl*LAM
    x = pn*np.sqrt(u); y = x*x
    return dict(pn=pn, pw=pw, u=u, m=m, x=x, y=y, nu=pn*LAM, p2=psi2(pn)*pw, p3=psi3(pn)*pw)


def env_of(J):
    return Q.P4.env if J % 2 else Q.A.env


def qdens(J, mfl, f, npn=256):
    """per-node smeared correction densities of a state with Q's weights (as ng_chain.dens with NF02 -> Q)."""
    s = qnodes(J, mfl, npn); x, m, u, y = s['x'], s['m'], s['u'], s['y']
    hs = ['+-'] if J % 2 else (['++'] if J < 4 else ['++', '+-'])
    G = dict(f(x, m)); env = env_of(J)
    out = dict(hs=hs, nu=s['nu'], x=x, s=s)
    if J % 2:
        out['+-'] = (s['p2']*env('+-', x, 2)*G['C2']/m**2, s['p3']*u*env('+-', x, 3)*G['C3H']/m**2)
    else:
        if len(hs) == 2:
            G[('C2', 'off')] = G[('C2', 'off')]/m; G[('C3H', 'off')] = G[('C3H', 'off')]/m
            eo2 = np.sqrt(Q.A.env('++', x, 2)*Q.A.env('+-', x, 2))
            out['off'] = (s['p2']*eo2*G[('C2', 'off')]/m, s['p3']*u*eo2/(2 - y)*G[('C3H', 'off')]/m)
        for h in hs:
            out[h] = (s['p2']*Q.A.env(h, x, 2)*G[('C2', h)]/m**2, s['p3']*u*Q.A.env(h, x, 3)*G[('C3H', h)]/m**2)
    return out


def absdens(J, mfl, npn=256):
    s = qnodes(J, mfl, npn); x, m, u, y = s['x'], s['m'], s['u'], s['y']
    hs = ['+-'] if J % 2 else (['++'] if J < 4 else ['++', '+-'])
    env = env_of(J)
    out = {h: (np.abs(s['p2'])*env(h, x, 2)/m**2, u*np.abs(s['p3'])*env(h, x, 3)/m**2) for h in hs}
    if len(hs) == 2:
        eo2 = np.sqrt(env('++', x, 2)*env('+-', x, 2))
        out['off'] = (np.abs(s['p2'])*eo2/m**2, u*np.abs(s['p3'])*eo2/(2 - y)/m**2)
    return out


_GF = {}


def gauge_factor(nu):
    """max(1, |Re X^Q_minnorm / psi2^Q|) in a +-0.5 window, nu <= 25 (Q's own min-norm Abel weight, resultq/cache/rv_k2eng_Y_1000.pkl)."""
    if 'g' not in _GF:
        y = pickle.load(open(os.path.join(RQ, 'cache', 'rv_k2eng_Y_1000.pkl'), 'rb'))
        with np.errstate(divide='ignore', invalid='ignore'):
            r = np.abs((y['a']*y['Y']).real/y['phi'].real)
        _GF['g'] = (y['nu'], np.where(np.isfinite(r), r, 1.0))
    g, r = _GF['g']
    out = np.ones_like(nu)
    for i, v in enumerate(nu):
        if v <= 25.0:
            sel = np.abs(g - v) <= 0.5
            out[i] = max(1.0, float(r[sel].max()) if sel.any() else 1.0)
    return out


# ---------------------------------------------------------------------------------------------------------------- chain helpers (ng_bind copies)
def grid_for(J, K, noown):
    NC.NOOWN = noown
    if J % 2:
        NC.KK, NC.PK = K, None
    else:
        NC.KK, NC.PK = 2, K
    W.DESPIKE_LOG.clear()
    g = Q.odd_grid(J) if J % 2 else Q.even_grid(J)
    return g, list(W.DESPIKE_LOG)


def model_of(J, g):
    grids = {'even': {}, 'odd': {}}
    grids['odd' if J % 2 else 'even'][J] = g
    em, om = NC.models(grids)
    return (om if J % 2 else em)[J][0]


def lin_model(J, g):
    from scipy.interpolate import CubicSpline as CS

    class _Lin:
        def __init__(self, x, y, bc_type=None):
            self.x, self.y = np.asarray(x), np.asarray(y)

        def __call__(self, t):
            return np.interp(np.asarray(t, float), self.x, self.y)
    Q.NG.CubicSpline = _Lin; Q.CubicSpline = _Lin
    try:
        return model_of(J, g)
    finally:
        Q.NG.CubicSpline = CS; Q.CubicSpline = CS


def sub(r, keep):
    n = len(r['xs']); sel = np.zeros(n, bool); sel[keep::2] = True; sel[0] = True; sel[-1] = True
    return {k: (v[sel] if isinstance(v, np.ndarray) and v.shape == (n,) else v) for k, v in r.items()}


def interp_nodes(g, m, x, key_fn):
    Ds = sorted(g); ms = np.array([g[D]['m'] for D in Ds])
    vals = [np.interp(np.minimum(x, g[D]['xs'][-1]), g[D]['xs'], key_fn(g[D])) for D in Ds]
    if m <= ms[0]:
        return vals[0]
    if m >= ms[-1]:
        return vals[-1]
    i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
    t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
    return (1 - t)*vals[i0] + t*vals[i1]


# ---------------------------------------------------------------------------------------------------------------- the state
class State:
    def __init__(self, J, hq, mfl, sd, dsd, soff=0.0, sd_pp=None):
        """hq: '++', '+-' or 'OO' (OO = the parity-odd even-J channel, image = ++ (PARITY_ODD_CLOSE)); sd/dsd: measured s_d of the
        channel of interest and its error (truncation + m-transfer); soff: |s_off| bound (even J >= 4); sd_pp: s_d of the other
        diagonal entry for even J >= 4 (the ++ channel when hq = '+-')."""
        self.J, self.hq, self.mfl, self.sd, self.dsd, self.soff, self.sd_pp = J, hq, mfl, sd, dsd, soff, sd_pp


def lam_matrix(J, mfl, hs, sds, dF, soff_signed=0.0):
    """normalised Q matrix and its min eigenvalue; returns (lam, FFh, Nh)."""
    u = 1/mfl**2
    FFh, Nh = [], []
    for h, sd in zip(hs, sds):
        b = qbase(J, mfl, h)
        FFh.append(b['FF'] - AQ*2*u**3*sd)          # s_v (Lam^-2-suppressed, a_h-convention dependent) is not in the central value: E_sv
        Nh.append(b['N'])
    FFh, Nh = np.array(FFh), np.array(Nh)
    Fm = np.diag(FFh) + dF
    if len(hs) == 2:
        Fm[0, 1] += -AQ*2*u**3*soff_signed; Fm[1, 0] = Fm[0, 1]
    n = 1/np.sqrt(Nh)
    return float(np.linalg.eigvalsh(n[:, None]*Fm*n[None, :])[0]), FFh, Nh


def dF_of(d, sel=None):
    hs = d['hs']; n = len(hs)
    sel = np.ones_like(d['nu'], bool) if sel is None else sel
    dF = np.zeros((n, n))
    for i, h in enumerate(hs):
        dF[i, i] = d[h][0][sel].sum() + d[h][1][sel].sum()
    if n == 2:
        dF[0, 1] = dF[1, 0] = d['off'][0][sel].sum() + d['off'][1][sel].sum()
    return dF


def errmat(ab, hs, Nh, sel, G):
    n = len(hs)
    Mx = np.zeros((n, n))
    Gv = np.broadcast_to(G, sel.shape)
    for i, h in enumerate(hs):
        a2, a3 = ab[h]
        Mx[i, i] = ((a2 + a3)*Gv)[sel].sum()/Nh[i]
    if n == 2:
        a2, a3 = ab['off']
        Mx[0, 1] = Mx[1, 0] = ((a2 + a3)*Gv)[sel].sum()/np.sqrt(Nh[0]*Nh[1])
    return float(np.sqrt((Mx*Mx).sum()))


def run_state(st, pr, npn=256):
    J, mfl = st.J, st.mfl
    m = mfl*LAM; u = 1/mfl**2
    oo = st.hq == 'OO'
    hs = ['+-'] if J % 2 else (['++'] if J < 4 else ['++', '+-'])
    # s_d per diagonal entry
    if J % 2 or J < 4:
        sds = [st.sd]
    else:
        sds = [st.sd_pp, st.sd] if st.hq == '+-' else [st.sd, st.sd_pp]
    ii = hs.index('++' if oo else st.hq)                                   # the entry of interest
    g2, dsp2 = grid_for(J, 2, False)
    g1, _ = grid_for(J, 1, False)
    g0m, _ = grid_for(J, 2, True)
    f2, f1, f0m = model_of(J, g2), model_of(J, g1), model_of(J, g0m)
    d2, d1, dm = qdens(J, mfl, f2, npn), qdens(J, mfl, f1, npn), qdens(J, mfl, f0m, npn)
    nu, x = d2['nu'], d2['x']
    # own row = the row nearest the state's mass (it is AT the mass for every state run here)
    Ds = sorted(g2); ms = np.array([g2[D]['m'] for D in Ds]); Down = Ds[int(np.argmin(np.abs(ms - m)))]
    own = g2[Down]; x0, xmax = own['xs'][0], own['xs'][-1]
    zero = np.zeros((len(hs), len(hs)))
    lam_hyb, FFh, Nh = lam_matrix(J, mfl, hs, [0.0]*len(hs), zero)               # K = 0 (HYB) value
    lam_sd, _, _ = lam_matrix(J, mfl, hs, sds, zero, st.soff)                    # + measured forward s_d (published level)
    lam2 = lam_matrix(J, mfl, hs, sds, dF_of(d2), st.soff)[0]
    lamm = lam_matrix(J, mfl, hs, sds, dF_of(dm), st.soff)[0]
    lamK1 = lam_matrix(J, mfl, hs, sds, dF_of(d1), st.soff)[0]
    # the off-diagonal forward term enters with unknown sign: take the worse sign
    if len(hs) == 2:
        lam2 = min(lam2, lam_matrix(J, mfl, hs, sds, dF_of(d2), -st.soff)[0])
    # truncation (node-wise, sign-free)
    tr = {}
    for kd in ('C2', 'C3H'):
        for h in (hs + (['off'] if len(hs) == 2 else [])):
            key = kd if J % 2 else (kd, h)
            gd = {}
            for D in g2:
                if D not in g1:
                    continue
                x2r, x1r = np.round(g2[D]['xs'], 6), np.round(g1[D]['xs'], 6)
                xc = np.intersect1d(x2r, x1r)
                if len(xc) < 2:
                    continue
                v2 = np.array([g2[D][key][int(np.argmin(np.abs(x2r - v)))] for v in xc])
                v1 = np.array([g1[D][key][int(np.argmin(np.abs(x1r - v)))] for v in xc])
                gd[D] = dict(m=g2[D]['m'], xs=xc, v=np.abs(v2 - v1))
            if isinstance(key, tuple) and key[1] == 'off':
                gd = {D: dict(v, v=v['v']*g2[D]['m']) for D, v in gd.items()}
            tr[(kd, h)] = interp_nodes(gd, m, x, lambda r: r['v'])
    ab = absdens(J, mfl, npn)
    win = (x >= x0) & (x <= xmax)
    n = len(hs); Mt = np.zeros((n, n))
    for i, h in enumerate(hs):
        a2, a3 = ab[h]
        Mt[i, i] = (a2*tr[('C2', h)] + a3*tr[('C3H', h)])[win].sum()/Nh[i]
    if n == 2:
        a2, a3 = ab['off']
        Mt[0, 1] = Mt[1, 0] = (a2*tr[('C2', 'off')] + a3*tr[('C3H', 'off')])[win].sum()/np.sqrt(Nh.prod())
    E_trunc = 3*float(np.sqrt((Mt*Mt).sum()))
    E_tl = 3*abs(lam2 - lamK1)
    # interpolation
    laml = lam_matrix(J, mfl, hs, sds, dF_of(qdens(J, mfl, lin_model(J, g2), npn)), st.soff)[0]
    lsub = []
    for keep in (0, 1):
        gs = dict(g2); gs[Down] = sub(own, keep)
        lsub.append(lam_matrix(J, mfl, hs, sds, dF_of(qdens(J, mfl, model_of(J, gs), npn)), st.soff)[0])
    E_interp = 3*max(abs(lam2 - laml), *[abs(lam2 - v) for v in lsub])
    # windows
    lowsel = nu < max(NU_LOW, x0*m)
    gf = gauge_factor(nu)
    E_low = errmat(ab, hs, Nh, lowsel, G_LOW*gf)
    gedge = 1e4
    for kd in ('C2', 'C3H'):
        for h in hs:
            key = kd if J % 2 else (kd, h)
            gedge = max(gedge, 10*np.max(np.abs(np.asarray(own[key])[-4:])))
    edsel = x > xmax
    E_edge = errmat(ab, hs, Nh, edsel, gedge)
    lamc = lam_matrix(J, mfl, hs, sds, dF_of(d2, ~(lowsel | edsel)), st.soff)[0]
    if len(hs) == 2:
        lamc = min(lamc, lam_matrix(J, mfl, hs, sds, dF_of(d2, ~(lowsel | edsel)), -st.soff)[0])
    # forward pieces (diagonal entry of interest; for the 2x2 the bound is applied to every diagonal entry: Frobenius)
    fw = np.array([AQ*2*u**3/Nh[i] for i in range(n)])
    E_sd = 3*st.dsd*float(np.sqrt((fw**2).sum())) if n == 2 else 3*st.dsd*float(fw[0])
    E_sv = 3*CQc*(2*u**3/LAM**2)*S_V*float(np.sqrt(((1/Nh)**2).sum()))
    E_C5 = {k: AQ*u**3*v*G5/(m**2)*float(np.sqrt(((1/Nh)**2).sum())) for k, v in A5.items()}
    E_OO = 0.0
    if oo:
        allsel = np.ones_like(nu, bool)
        E_OO = errmat(ab, hs, Nh, allsel, 0.12) + AQ*2*u**3*5e-6/Nh[ii]
    tot = E_trunc + E_interp + E_low + E_edge + E_sd + E_sv + E_C5['XP'] + E_OO
    tot_w = tot - E_C5['XP'] + E_C5['worst']
    # base
    lc = {h: lemmaC(J, mfl, h, npn) for h in hs}
    base = {}
    for i, h in enumerate(hs):
        b = qbase(J, mfl, h)
        base[h] = dict(flat=b['FFf']/b['Nf'], hyb=b['FF']/b['N'], curv=(b['FF'] - b['FFf'])/b['N'],
                       lemmaC=(lc[h][0] + lc[h][1])/b['N'], B2=b['B2'], B3=b['B3'], FW=b['FW'], B4=b['B4'], N=b['N'])
    out = dict(J=J, h=st.hq, mfl=mfl, Down=float(Down), rows=[float(D) for D in Ds], rows_model=sorted(float(D) for D in g0m),
               sd=st.sd, dsd=st.dsd, soff=st.soff, lam_hyb=lam_hyb, lam_sd=lam_sd, lam_model=lamm, lam_K2=lam2, lam_K1=lamK1,
               lam_central=lamc, lam_lin=laml, lam_sub=lsub, E_trunc=E_trunc, E_trunc_lam=E_tl, E_interp=E_interp, E_low=E_low,
               E_edge=E_edge, G_edge=float(gedge), E_sd=E_sd, E_sv=E_sv, E_C5=E_C5, E_OO=E_OO, E_tot=tot, lower=lamc - tot,
               L=lamc/tot, E_tot_worstC5=tot_w, L_worstC5=lamc/tot_w, base=base, despike=dsp2, x0=float(x0), xmax=float(xmax),
               dF_K2=dF_of(d2).tolist(), dF_model=dF_of(dm).tolist(), N=Nh.tolist(), FF=FFh.tolist())
    pr(f"  {st.hq} J={J} m/Lam={mfl:.5f} own row {Down:g} (rows {out['rows']}; published {out['rows_model']}): "
       f"HYB {lam_hyb:+.5f} -> +s_d {lam_sd:+.5f} -> +dB (published rows) {lamm:+.5f} | own row K2 {lam2:+.5f} "
       f"(K1 {lamK1:+.5f}; lin {laml:+.5f}; sub {lsub[0]:+.5f}/{lsub[1]:+.5f}) | central (windows cut) {lamc:+.5f}")
    pr(f"        smeared dF (K2, own row) {np.array(out['dF_K2']).round(8).tolist()}  N {np.round(Nh, 4).tolist()}  s_d {sds} +- {st.dsd}"
       + (f"  |s_off| {st.soff}" if n == 2 else ''))
    pr(f"        errors: trunc {E_trunc:.2e}  interp {E_interp:.2e}  low(nu<{max(NU_LOW, x0*m):.0f}, G=1e4 x gauge<= {gf[lowsel].max():.2f}) "
       f"{E_low:.2e}  edge(x>{xmax:.3f}, G={gedge:.2e}) {E_edge:.2e}  s_d {E_sd:.2e}  s_v {E_sv:.2e}  C5(G=1e4, A5=1.80) "
       f"{E_C5['XP']:.2e}" + (f"  OO {E_OO:.2e}" if oo else '') +
       f"  -> total {tot:.2e}; lower bound {lamc - tot:+.5f}; L = {lamc/tot:.1f}  [worst-gauge C5 (A5=7.66): total {tot_w:.2e}, L = {lamc/tot_w:.1f}]")
    pr(f"        lambda-level K2 truncation 3|K2-K1| = {E_tl:.2e} (node-wise used above: {E_trunc:.2e})")
    for h, b in base.items():
        pr(f"        base {h}: flat {b['flat']:+.6f}, HYB {b['hyb']:+.6f} (LR curvature {b['curv']:+.2e} N; Lemma-C bound {b['lemmaC']:.2e} N) "
           f"B2 {b['B2']:+.4f} B3 {b['B3']:+.4f} FW {b['FW']:+.4f} B4 {b['B4']:+.4f} N {b['N']:.4f}")
    if dsp2:
        pr(f"        despiked: {dsp2}")
    return out


if __name__ == '__main__':
    tag = sys.argv[1]
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    pr(f"# qb_bind Dg={DG:g} functional {QNAME} (a = {AQ:.4f}, c = {CQc}): per state HYB -> +s_d -> +smeared K2 (own row); error budget")
    res = []
    for spec in sys.argv[2:]:
        p = spec.split(':')
        J, hq, mfs, sd, dsd = int(p[0]), p[1], p[2], float(p[3]), float(p[4])
        soff = float(p[5]) if len(p) > 5 else 0.0
        sdpp = float(p[6]) if len(p) > 6 else None
        for v in mfs.split(','):
            res.append(run_state(State(J, hq, float(v), sd, dsd, soff, sdpp), pr))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(NGC, 'qb_bind_%s.txt' % tag), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(NGC, 'qb_bind_%s.json' % tag), 'w'))
