"""
q3 (S1 + S4 fixes, S1_S4_FIXES.md): near-gap graviton corrections for even J <= 20 (2x2) and odd J = 5 ... 13 (1x1, parity-odd +-) with
C^_2' and the FIXED C^_3, and n at Delta_g = 300, 1000 -- the l9/p4 K2 linear chain (H0 §23.2, PARITY_ODD_ROUTEA §6), read-only imports.

Changes relative to the chains of v3_heavy_impact / v4_odd_impact (the C^_2' baseline of the review):
  * k = 2 local de-shift: spinning (1 +- u_T)/2 with the exact rational u_T (v4_chi.u_spin)                    [as v3/v4: C^_2']
  * k = 3 (C^_3): Q3_C3SHIFT = 'none' (DEFAULT, the S4 decision: the leading-Regge heavy image of Psi_3[o_H q] is l-isotypic,
    sigma_h c3(nu) P_l -- q1/q1b), or 'scalar' (the old chains: scalar-g de-shift applied to k = 3; for comparison only)
  * J = 11 ... 20: crystal points from q2_run (cache/q2_*), merged with the l8/l7 C2 references the l9 engine reuses
  * every record: logK -> logK - i pi (Delta + J) (PARITY_ODD_ROUTEA §12 prescription; asserted on Im logK mod 2 pi; a no-op for the
    state/reference ratios on even-integer rows, which all rows here are)
Near gap (one combined pass over ALL exact-skeleton states, e1/e7 of spf_common.load_states): even J -> l9_neargap 2x2 formula with the
grid model (cubic spline in x, linear in 1/m between rows, held outside); odd J -> p6 1x1 formula.  n, min lambda/N_s, t* (common scale
on all computed corrections), and per J: min lambda, max |dF_hh|/N_h, max |dF_+-|/sqrt(N N), and t*_J (scale at which a state of that J
first fails).  Row variants (env Q3_ROWS): 'all' (every computed row <= 3200), 'no200' (J >= 11 without the Delta = 200 row),
'no100_200' (no Delta = 100, 200 rows for any J).
Usage: Q3_C3SHIFT=none Q3_ROWS=all python3 q3_analysis.py [Dg ...]  -> cache/q3_analysis_<shift>_<rows>.txt, cache/q3_grids_<shift>.pkl
"""
import os, sys, io, glob, pickle, time
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np
from scipy.special import j0
from scipy.interpolate import CubicSpline
import l9_analysis as A
import l9_neargap as NG
import p4_odd_analysis as P4
from spf_common import NearGap, load_states
from common import gl01, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk
from v4_chi import u_spin

C3SHIFT = os.environ.get('Q3_C3SHIFT', 'none')
ROWS = os.environ.get('Q3_ROWS', 'all')
JEVEN = [int(v) for v in os.environ.get('Q3_JEVEN', '0,2,4,6,8,10,12,14,16,18,20').split(',')]
JODD = [int(v) for v in os.environ.get('Q3_JODD', '5,7,9,11,13').split(',')]
CACHE = os.path.join(HERE, 'cache')
PHASE_LOG = []
NONFINITE = []


# ------------------------------------------------------------------------------------------------ records (q2 caches + phase strip)
def _strip(R, J):
    out = {}
    for k, r in R.items():
        D = k[0]
        if not all(np.all(np.isfinite(r[o])) for o in ('M0', 'M2', 'M4', 'logC', 'logK') if o in r):
            PHASE_LOG.append(0.0); NONFINITE.append((J, k)); continue          # dropped (overflow in the l1 engine)
        dev = (r['logK'].imag - np.pi*(D + J)) % (2*np.pi)
        dev = min(dev, 2*np.pi - dev)
        PHASE_LOG.append(dev)
        assert dev < 1e-6, (k, J, r["logK"])
        out[k] = dict(r, logK=r['logK'] - 1j*np.pi*(D + J))
    return out


def _q2_even(kind, J):
    R = {}
    for fn in sorted(glob.glob(os.path.join(CACHE, 'q2_pts_%s_J%d_k2*.pkl' % (kind, J)))):
        R.update(pickle.load(open(fn, 'rb')))
    return R


_orig_load = A._load
_orig_records = A.L.records


def _load(kind, J, sfx):
    if J >= 11:
        return _strip(_q2_even(kind, J), J) if sfx == '_k2' else {}
    return _strip(_orig_load(kind, J, sfx), J)


def _records(J):
    return _strip(_orig_records(J), J)


A._load = _load
A.L.records = _records


def _q2_odd(kind, J, anti=False, k3=False):
    if J < 11:
        R = _orig_p4load_raw(kind, J, k3)
    else:
        R = {}
        for fn in sorted(glob.glob(os.path.join(CACHE, 'q2_p3_%s_J%d_*.pkl' % (kind, J)))):
            R.update(pickle.load(open(fn, 'rb')))
    R = {(k[0], k[1]): v for k, v in R.items() if k[2] == 1.0}
    return _strip(R, J)


def _orig_p4load_raw(kind, J, k3):
    fn = os.path.join(CACHE, 'p3_pts_%s_J%d%s.pkl' % (kind, J, '_k3' if k3 else ''))
    return pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}


P4.load = _q2_odd


# ------------------------------------------------------------------------------------------------ de-shifts
def _shift_generic(env, Dc, SIG):
    import l2_scalar_validate as L2
    def shift(J, h, k, xs, m):
        sg = SIG[h] if k == 3 else 1.0
        f = lambda t: sg*env(h, t, k)*Dc(h, J, t)
        out = []
        for x in xs:
            if k == 2:
                u = complex(u_spin(x*m)); gp, gm = 0.5*(1 + u), 0.5*(1 - u)            # C^_2' (spinning, exact)
            elif C3SHIFT == 'scalar':
                gp, gm = L2.gfac(x*m); s = gp + gm; gp, gm = gp/s, gm/s                  # old chains (scalar g on k = 3)
            else:
                out.append(0.0); continue                                                # C^_3: no de-shift (S4 decision)
            out.append((gp*f(x + 1j/m) + gm*f(x - 1j/m) - f(x)).real)
        return np.array(out)
    return shift


A.shift = _shift_generic(A.env, A.Dc, A.SIG)
P4.shift = _shift_generic(P4.env, P4.Dc, P4.SIG)


# ------------------------------------------------------------------------------------------------ grids
def row_filter(J, Ds):
    if ROWS == 'no200' and J >= 11:
        return [D for D in Ds if D != 200.0]
    if ROWS == 'no100_200':
        return [D for D in Ds if D not in (100.0, 200.0)]
    return Ds


ALIGN = os.environ.get('Q3_ALIGN', '1') == '1'
OFFSCALE = os.environ.get('Q3_OFFSCALE', '1') == '1'
XCUT = float(os.environ.get('Q3_XCUT', '0'))   # 1: interpolate/hold c_off = m Gamma_off (validated m-independent); 0: Gamma_off (l9_neargap)
ALIGN_LOG = {}


def align_rows(J, g):
    """The sign of each row's helicity basis vectors (l8_analysis.basis -> realvec: phase fixed by the largest component) is arbitrary,
    so the off-diagonal Gamma_off of different Delta rows can come with OPPOSITE signs; l9_neargap interpolates them in 1/m as they are
    (a partial cancellation between rows).  Fix: per row, a sign s_D such that the robust, m-independent IMAGINARY C2 mixing
    m Im Gamma_off (H0 §22.2; asymptotically Delta-independent) overlaps positively with that of the most asymptotic row; s_D multiplies
    the off-diagonals of both kinds (common basis).  Diagnostic: the normalised overlap (should be close to +-1)."""
    Ds = sorted(g)
    if J < 4 or len(Ds) < 2:
        return g
    ref = g[Ds[-1]]
    xr = ref['xs']; ir = ref[('C2', 'off', 'im')]
    for D in Ds:
        r = g[D]
        xi = np.intersect1d(np.round(r['xs'], 6), np.round(xr, 6))
        a = np.array([r[('C2', 'off', 'im')][np.argmin(np.abs(r['xs'] - x))] for x in xi])
        b = np.array([ir[np.argmin(np.abs(xr - x))] for x in xi])
        ov = float(a @ b/np.sqrt((a @ a)*(b @ b)))
        sD = 1.0 if (ov >= 0 or not ALIGN) else -1.0
        ALIGN_LOG[(J, D)] = (ov, sD)
        if sD < 0:
            for kd in ('C2', 'C3H'):
                r[(kd, 'off')] = -r[(kd, 'off')]
                r[(kd, 'off', 'im')] = -r[(kd, 'off', 'im')]
    return g


DIRECT400 = os.environ.get('Q3_DIRECT', '0') == '1'      # replace the Delta = 400 row by q4's direct (RQMC) row where computed


def even_grid(J):
    g = A.grid(J, 2, ('C2', 'C3H'))
    if DIRECT400:
        fn = os.path.join(CACHE, 'q4_direct_grids.pkl')
        dg = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
        if J in dg:
            g[400.0] = dg[J]['direct']
        fr = os.path.join(CACHE, 'q4_direct_grids_rows.pkl')
        for (Jr, Dr), v in (pickle.load(open(fr, 'rb')) if os.path.exists(fr) else {}).items():
            if Jr == J:
                g[Dr] = v['direct']
    Ds = row_filter(J, sorted(D for D in g if D <= 3200.0))
    return align_rows(J, {D: g[D] for D in Ds})


def odd_grid(J):
    g = P4.grid(J, K=2, kinds=('C2', 'C3H'))
    Ds = row_filter(J, sorted(D for D in g if D <= 3200.0))
    return {D: g[D] for D in Ds}


def odd_model(g):
    Ds = sorted(g); ms = np.array([g[D]['m'] for D in Ds])
    spl = {D: {kd: CubicSpline(g[D]['xs'], g[D][kd], bc_type='natural') for kd in ('C2', 'C3H')} for D in Ds}
    xmax = {D: g[D]['xs'][-1] for D in Ds}
    def f(x, m):
        if m <= ms[0]:
            i0 = i1 = 0; t = 0.0
        elif m >= ms[-1]:
            i0 = i1 = len(Ds) - 1; t = 0.0
        else:
            i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
            t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
        return {kd: (1 - t)*spl[Ds[i0]][kd](np.minimum(x, xmax[Ds[i0]])) + t*spl[Ds[i1]][kd](np.minimum(x, xmax[Ds[i1]]))
                for kd in ('C2', 'C3H')}
    return f, Ds


# ------------------------------------------------------------------------------------------------ near gap (combined)
def neargap(Dg, emods, omods):
    NGp = NearGap(Dg); F = NGp.F; pn, pw, Lm = NGp.pn, NGp.pw, NGp.Lam
    blocks = {}
    states = list(load_states(Dg))
    # optional extension (Q3_EXTRA_M, e.g. "1.4,1.5,1.75,2"): hybrid states (flat + leading-Regge curvature, scalar remainder := 0, as
    # ROUTEA_LARGE_B's F_hyb) at m/Lam beyond the exact-skeleton set, for the even/odd J with crystal grids; DF tag 'hyb' [M]
    for mfl in [float(v) for v in os.environ.get('Q3_EXTRA_M', '').split(',') if v]:
        for J in sorted(set(emods) | set(omods)):
            if J < 11 or (mfl <= 1.3 + 1e-9 and J <= 20):          # J <= 20: exact-skeleton states exist up to 1.3 Lam
                continue
            mm = mfl*Lm
            alpha = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*mm*mm))/2
            if alpha < Dg - 2 - 1e-9:                             # twist gap (omegau/common.state convention: tau >= Delta_g)
                continue
            r0 = dict(mfl=mfl, J=J, m=mfl*Lm, eta=1.0 + ((2*J + 1)/mfl)**2/(2*Lm*Lm), rho={k: np.zeros_like(pn) for k in range(2, 8)})
            for h in ('++', '+-'):
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                states.append((h, 'hyb', r0, None))
    for h, DF, r0, ref in states:
        mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]
        x = pn*np.sqrt(u); y = x*x
        P = remk.Pders_vec(Lm*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
        pf2, wf2 = gl01(max(400, 2*J + 100))
        Js, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
        p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
        cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
        fac = (1 - y)**2 if h == '++' else (1 - y)**-2
        rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
        Ns = abs(B2 + cv2) + abs(B3 + cv3)
        blocks.setdefault((round(mfl, 9), J, DF), {})[h] = dict(base=B2 + B3 + cv2 + cv3 + rem, Ns=Ns, r0=r0, x=x, y=y, p2=p2, p3=p3,
                                                                 u=u, m=r0['m'])
    recs = []
    for (mfl, J, DF), b in blocks.items():
        hs = [h for h in ('++', '+-') if h in b]
        base = np.array([b[h]['base'] for h in hs]); Ns = np.array([b[h]['Ns'] for h in hs])
        dF = np.zeros((len(hs), len(hs)))
        e = b[hs[0]]; x, y, p2, p3, u, m = e['x'], e['y'], e['p2'], e['p3'], e['u'], e['m']
        rho2, rho3 = e['r0']['rho'][2], e['r0']['rho'][3]
        r20 = np.polyfit(x[:3], rho2[:3], 2)[-1] if np.any(rho2) else 0.0
        cov = False
        if J % 2 == 0 and J in emods:
            cov = True
            G = emods[J][0](x, m)
            if OFFSCALE and len(hs) == 2:
                G = dict(G); G[('C2', 'off')] = G[('C2', 'off')]/m; G[('C3H', 'off')] = G[('C3H', 'off')]/m
            if XCUT > 0:                       # diagnostic: graviton correction switched off below x = XCUT (model remainder kept there)
                G = {k: np.where(x < XCUT, (A.SIG[k[1]] if (k[0] == 'C3H' and k[1] in A.SIG) else 1.0)*0*v + (
                    m*m*(rho2 - A.Dsh(k[1], J, x)*r20) if k[1] in ('++', '+-') and k[0] == 'C2' else
                    m*m*A.SIG[k[1]]*(rho3 - A.Dsh(k[1], J, x)*r20) if k[1] in ('++', '+-') else 0*v), v) for k, v in G.items()}
            for kk, h in enumerate(hs):
                Dh = A.Dsh(h, J, x)
                dF[kk, kk] = (p2*A.env(h, x, 2)*(G[('C2', h)]/m**2 - (rho2 - Dh*r20))).sum()
                dF[kk, kk] += (p3*u*A.env(h, x, 3)*(G[('C3H', h)]/m**2 - A.SIG[h]*(rho3 - Dh*r20))).sum()
            if len(hs) == 2:
                eo2 = np.sqrt(A.env('++', x, 2)*A.env('+-', x, 2))
                o = (p2*eo2*G[('C2', 'off')]).sum()/m + (p3*u*eo2/(2 - y)*G[('C3H', 'off')]).sum()/m
                dF[0, 1] = dF[1, 0] = o
        elif J % 2 == 1 and J in omods:
            cov = True
            G = omods[J][0](x, m)
            Dh = P4.Dc('+-', J, x)
            dF[0, 0] = (p2*P4.env('+-', x, 2)*(G['C2']/m**2 - (rho2 - Dh*r20))).sum()
            dF[0, 0] += (p3*u*P4.env('+-', x, 3)*(G['C3H']/m**2 - P4.SIG['+-']*(rho3 - Dh*r20))).sum()
        recs.append(dict(mfl=mfl, J=J, DF=DF, hs=hs, base=base, Ns=Ns, dF=dF, cov=cov))
    return recs


def lam(r, t):
    Fm = np.diag(r['base']) + t*r['dF']; n = 1/np.sqrt(r['Ns'])
    return np.linalg.eigvalsh(n[:, None]*Fm*n[None, :])[0]


def tstar(recs, cond):
    out = []
    for sgn in (1, -1):
        lo, hi = 0.0, 1.0
        while cond(sgn*hi) and hi < 1e4:
            lo, hi = hi, 2*hi
        if hi >= 1e4:
            out.append(sgn*np.inf); continue
        for _ in range(40):
            mid = 0.5*(lo + hi)
            lo, hi = (mid, hi) if cond(sgn*mid) else (lo, mid)
        out.append(sgn*lo)
    return out[1], out[0]


def summarize(Dg, recs, pr):
    grid = sorted(set(r['mfl'] for r in recs))
    def nval(t):
        return next((n0 for n0 in grid if all(lam(r, t) > 0 for r in recs if r['mfl'] >= n0 - 1e-9)), None)
    w1 = min(recs, key=lambda r: lam(r, 1.0)); w0 = min(recs, key=lambda r: lam(r, 0.0))
    ts = tstar(recs, lambda t: nval(t) == grid[0])
    pr(f"  Delta_g = {Dg}: n = {nval(1.0)} (uncorrected {nval(0.0)}); min lam/N_s = {lam(w1, 1.0):+.4f} [m/Lam={w1['mfl']} J={w1['J']} "
       f"Df={w1['DF']} {w1['hs']}] (uncorr. {lam(w0, 0.0):+.4f} [J={w0['J']}]); t* = ({ts[0]:+.1f}, {ts[1]:+.1f})")
    res = dict(n=nval(1.0), min=lam(w1, 1.0), w=(w1['mfl'], w1['J'], w1['DF']), tstar=ts, perJ={})
    for J in sorted(set(r['J'] for r in recs)):
        rj = [r for r in recs if r['J'] == J]
        if not any(r['cov'] for r in rj) and J > 3:
            pr(f"      J={J:2d} (uncovered): min lam/N_s {min(lam(r, 0.0) for r in rj):+.4f}")
            continue
        wj = min(rj, key=lambda r: lam(r, 1.0))
        dg = max(abs(r['dF'][k, k])/r['Ns'][k] for r in rj for k in range(len(r['hs'])))
        of = max((abs(r['dF'][0, 1])/np.sqrt(r['Ns'].prod()) if len(r['hs']) == 2 else 0) for r in rj)
        tj = tstar(rj, lambda t: all(lam(r, t) > 0 for r in rj))
        res['perJ'][J] = dict(min=lam(wj, 1.0), unc=min(lam(r, 0.0) for r in rj), at=wj['mfl'], diag=dg, off=of, tJ=tj)
        pr(f"      J={J:2d}: min lam/N_s {lam(wj, 1.0):+.4f} (uncorr. {min(lam(r, 0.0) for r in rj):+.4f}) at m/Lam={wj['mfl']}; "
           f"max |dF_hh|/N_h {dg:.4f}; max |dF_+-|/sqrt(NN) {of:.4f}; t*_J = ({tj[0]:+.1f}, {tj[1]:+.1f})")
    return res


def run(Dgs, pr):
    emods, omods, grids = {}, {}, {'even': {}, 'odd': {}}
    for J in JEVEN:
        try:
            g = even_grid(J)
        except Exception as ex:
            pr(f"  even J={J}: no grid ({ex})"); continue
        if len(g) >= 1:
            grids['even'][J] = g
            if OFFSCALE and J >= 4:                           # model c_off := m Gamma_off (m-independent, §2.1) instead of Gamma_off
                g = {D: {**r, **{(kd, 'off'): r[(kd, 'off')]*r['m'] for kd in ('C2', 'C3H')}} for D, r in g.items()}
            if len(g) == 1:                                   # single row: duplicate (held everywhere)
                D = list(g)[0]; g = {D: g[D], D + 1e-6: dict(g[D], m=g[D]['m']*(1 + 1e-9))}
            emods[J] = NG.model(g, J, ('C2', 'C3H'))
    for J in JODD:
        try:
            g = odd_grid(J)
        except Exception as ex:
            pr(f"  odd J={J}: no grid ({ex})"); continue
        if len(g) >= 1:
            grids['odd'][J] = g
            if len(g) == 1:
                D = list(g)[0]; g = {D: g[D], D + 1e-6: dict(g[D], m=g[D]['m']*(1 + 1e-9))}
            omods[J] = odd_model(g)
    pr("# off-diagonal row alignment (J, Delta: overlap of m Im Gamma_off,C2 with the largest row, applied sign): " +
       "; ".join(f"{J},{D:.0f}: {ov:+.2f}{'' if sD > 0 else ' FLIP'}" for (J, D), (ov, sD) in sorted(ALIGN_LOG.items())))
    pr(f"# q3 [C3 shift: {C3SHIFT}; rows: {ROWS}; align: {ALIGN}; off-diagonal model: {'c_off = m Gamma_off' if OFFSCALE else 'Gamma_off (as l9)'}; direct Delta=400 rows: {DIRECT400}]  even rows: " +
       "; ".join(f"J={J}: " + ",".join(f"{D:.0f}" for D in sorted(grids['even'][J])) for J in sorted(grids['even'])) +
       "  | odd rows: " + "; ".join(f"J={J}: " + ",".join(f"{D:.0f}" for D in sorted(grids['odd'][J])) for J in sorted(grids['odd'])))
    out = {}
    for Dg in Dgs:
        recs = neargap(Dg, emods, omods)
        out[Dg] = summarize(Dg, recs, pr)
        out[Dg]['recs'] = recs
    return grids, out


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    Dgs = [int(a) for a in sys.argv[1:]] or [300, 1000]
    grids, out = run(Dgs, pr)
    pr(f"# non-finite records dropped: {NONFINITE}")
    pr(f"# phase strip: {len(PHASE_LOG)} records, max |Im logK - pi(Delta + J)| mod 2 pi = {max(PHASE_LOG):.1e}")
    pr(f"[{time.time() - t0:.0f}s]")
    tag = f"{C3SHIFT}_{ROWS}" + ('' if ALIGN else '_noalign') + ('' if OFFSCALE else '_gammaoff') + ('_direct400' if DIRECT400 else '') + os.environ.get('Q3_TAG', '')
    open(os.path.join(CACHE, 'q3_analysis_%s.txt' % tag), 'w').write(buf.getvalue())
    pickle.dump(dict(grids=grids, out={Dg: {k: v for k, v in o.items() if k != 'recs'} for Dg, o in out.items()}),
                open(os.path.join(CACHE, 'q3_grids_%s.pkl' % tag), 'wb'))
